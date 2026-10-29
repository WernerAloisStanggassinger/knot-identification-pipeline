# SPDX-License-Identifier: GPL-3.0-or-later
# Knot-Identification Pipeline — Werner Alois Stanggassinger & Claude Code (Anthropic).
# Counter-case for run_jobs() from calibrate_full_corpus_robust.py (Disc #704).
# Four cases. The predecessor (mp.Pool.imap_unordered) fails at B and C:
# it hangs silently. Exactly that must NOT happen any more here.
import os, signal, sys, time
import site_config
sys.path.insert(0, site_config.SIXD_DIR)
from calibrate_full_corpus_robust import run_jobs

MARKER = "/tmp/countercase_runjobs_marker_"

def well_behaved(job):
    return job[0], 0.1, False, f"CONVERGED {job[0]}"

def dies_once(job):
    if job[0] == "lethal":
        m = MARKER + "once"
        if not os.path.exists(m):
            open(m, "w").close()
            time.sleep(0.3)
            os.kill(os.getpid(), signal.SIGTERM)   # earlyoom equivalent
            time.sleep(5)
        return job[0], 0.1, False, f"CONVERGED {job[0]} (zweiter Versuch)"
    return well_behaved(job)

def dies_always(job):
    if job[0] == "unsalvageable":
        time.sleep(0.2)
        os.kill(os.getpid(), signal.SIGTERM)
        time.sleep(5)
    return well_behaved(job)

def raises(job):
    if job[0] == "raiser":
        raise ValueError("broken curve")
    return well_behaved(job)

def jobs(names):
    return [(n, "/no/path", 0) for n in names]

if __name__ == "__main__":
    failures = 0
    def check(name, actual, expected):
        global failures
        ok = actual == expected
        failures += (not ok)
        print(f"   {'✅' if ok else '🔴'} {name}: {actual}  (expected {expected})")

    for f in ("once",):
        try: os.remove(MARKER + f)
        except OSError: pass

    print("\n=== A  all well-behaved — no restart, no losses ===")
    r = run_jobs(jobs(["a", "b", "c", "d", "e"]), 1, fn=well_behaved)
    check("(ok, fail, given_up, restarts)", r, (5, 0, 0, 0))

    print("\n=== B  a worker dies ONCE — task must come back ===")
    t0 = time.time()
    r = run_jobs(jobs(["v", "lethal", "w", "x"]), 1, fn=dies_once)
    print(f"   duration {time.time()-t0:.1f}s")
    check("(ok, fail, given_up, restarts)", r, (4, 0, 0, 1))

    print("\n=== C  task kills the worker ALWAYS — must give up, not hang ===")
    t0 = time.time()
    r = run_jobs(jobs(["p", "unsalvageable", "q"]), 1, max_attempts=3, fn=dies_always)
    print(f"   duration {time.time()-t0:.1f}s")
    check("ok",           r[0], 2)
    check("given_up",   r[2], 1)
    check("restarts>=2", r[3] >= 2, True)

    print("\n=== D  exception IN the worker — count as fail, not as death ===")
    r = run_jobs(jobs(["y", "raiser", "z"]), 1, fn=raises)
    check("(ok, fail, given_up, restarts)", r, (2, 1, 0, 0))

    print(f"\n{'✅ ALL FOUR GREEN' if not failures else f'🔴 {failures} FAILED'}")
    sys.exit(1 if failures else 0)
