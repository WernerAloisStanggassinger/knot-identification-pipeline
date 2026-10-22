# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""Pipeline-wide subprocess helper that kills the whole process group on timeout.

Background (Disc #222, 2026-05-21):
  `subprocess.run([..., 'conda', 'run', ...], timeout=N)` sends SIGTERM on
  timeout only to the direct child (conda's run wrapper). Conda's
  subprocess handling does NOT pass the signal on to the inner python
  process that does the actual work (e.g. knotpy.kauffman).

  Result: on timeout the conda front end ends, the inner python becomes
  an orphan and keeps running, sometimes for hours at 90% CPU. On
  2026-05-21 we observed 17 such orphans after ~2h, blocking 10 cores.

Fix idea:
  start_new_session=True puts the subprocess tree into its own process
  group. On timeout os.killpg() kills ALL processes in the group:
  conda, the inner python, all children. Clean cleanup, no orphans.

Usage (drop-in for subprocess.run calls that invoke conda):

  # Before (orphan risk):
  res = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

  # After (process-group killing):
  from pipeline_subprocess import run_capture_with_group_kill
  res = run_capture_with_group_kill(cmd, timeout=600)

The function has the same signature as subprocess.run with capture_output=True
and text=True. On timeout it raises subprocess.TimeoutExpired like the
original (the caller's try/except stays the same).
"""
from __future__ import annotations
import os
import signal
import subprocess
import threading

# ─────────────────────────────────────────────────────────────────────────
# Active registry + external signal killer (OOM 2026-07-01, "hole B").
#
# The timeout group kill below fires ONLY in the except-TimeoutExpired branch.
# Two cases otherwise stay open and leave orphans (conda/sage grandchildren
# + their multiprocessing.Pool workers, which thanks to start_new_session=True
# live in their OWN session and do NOT see an external SIGTERM sent to the
# parent):
#   1. OOM: the kernel kills single workers long before the timeout.
#   2. mem_watchdog.sh / manual kill: SIGTERM hits the Python parent,
#      but never the detached grandchild session.
# Fix: every running subprocess registers its pgid; install_signal_killer()
# installs a handler that, on an external signal, FIRST kills all active groups
# via killpg and then ends itself through the default handler.
# ─────────────────────────────────────────────────────────────────────────
_active_lock = threading.Lock()
_active_pgids: set[int] = set()


def _register_pgid(pid):
    try:
        pgid = os.getpgid(pid)
    except (ProcessLookupError, OSError):
        return None
    with _active_lock:
        _active_pgids.add(pgid)
    return pgid


def _unregister_pgid(pgid):
    if pgid is None:
        return
    with _active_lock:
        _active_pgids.discard(pgid)


def kill_all_active(sig=signal.SIGKILL):
    """Kill ALL currently registered subprocess groups (conda/sage/Pool grandchildren).

    Meant for signal handlers: closes the orphan path that the plain
    timeout group kill does not cover (OOM / external SIGTERM)."""
    with _active_lock:
        pgids = list(_active_pgids)
    for pgid in pgids:
        try:
            os.killpg(pgid, sig)
        except (ProcessLookupError, OSError):
            pass
    return pgids


def install_signal_killer(signals=(signal.SIGTERM, signal.SIGINT, signal.SIGHUP)):
    """Install a handler that, on an external kill, first kills all subprocess groups,
    then restores the default behaviour and ends the process.

    Idempotent enough to be called at the start of main(). mem_watchdog.sh sends SIGTERM,
    waits 3s, then SIGKILL; within that window this handler kills the grandchildren."""
    def _handler(signum, frame):
        kill_all_active(signal.SIGKILL)
        # Restore the default handler and send the same signal to ourselves ->
        # correct exit status, no double handling.
        try:
            signal.signal(signum, signal.SIG_DFL)
            os.kill(os.getpid(), signum)
        except (ValueError, OSError):
            os._exit(143)
    for s in signals:
        try:
            signal.signal(s, _handler)
        except (ValueError, OSError):
            # e.g. in a non-main thread -> skip silently
            pass


def run_capture_with_group_kill(cmd, timeout, input=None):
    """Equivalent of subprocess.run that kills the whole process group on timeout.

    Args:
        cmd: command + args as a list, as for subprocess.run
        timeout: timeout in seconds
        input: optional stdin text (added 2026-07-28 for aggressive_simplify.py,
               which sends its items via stdin). Without this parameter exactly that
               one call site stayed on raw subprocess.run, and a bug class that is
               99 % repaired is not repaired (Disc #586).
               Default None => stdin is inherited from the parent as before; behaviour for all
               existing callers unchanged.

    Returns:
        subprocess.CompletedProcess with stdout/stderr as str

    Raises:
        subprocess.TimeoutExpired: on timeout (after process group SIGTERM/SIGKILL).
    """
    p = subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE if input is not None else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    # Hole B: register the pgid so that an external signal handler
    # (install_signal_killer) also kills the grandchild session on OOM/kill.
    reg_pgid = _register_pgid(p.pid)
    try:
        stdout, stderr = p.communicate(input=input, timeout=timeout)
        return subprocess.CompletedProcess(p.args, p.returncode, stdout, stderr)
    except subprocess.TimeoutExpired:
        # Phase 1: graceful SIGTERM to the whole process group
        pgid = reg_pgid
        try:
            if pgid is None:
                pgid = os.getpgid(p.pid)
            os.killpg(pgid, signal.SIGTERM)
        except (ProcessLookupError, OSError):
            pass
        # Phase 2: wait 2s to see whether a graceful exit works
        try:
            p.communicate(timeout=2)
        except subprocess.TimeoutExpired:
            # Phase 3: SIGKILL to the group
            if pgid is not None:
                try:
                    os.killpg(pgid, signal.SIGKILL)
                except (ProcessLookupError, OSError):
                    pass
            # Phase 4: wait 2s + final reap
            try:
                p.communicate(timeout=2)
            except subprocess.TimeoutExpired:
                pass
        raise subprocess.TimeoutExpired(cmd, timeout)
    finally:
        _unregister_pgid(reg_pgid)
