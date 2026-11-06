# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
# build_external_theta_db.sage — pre-compute Theta_numeric for the
# external_tools/bar_natan_theta/knots.sage knots and write
# external_tools/theta_db/knots_sage_theta_db.json (Drors-QR-Phase-2, Discovery #173).
#
# Usage:
#   conda run -n sage_stable sage build_external_theta_db.sage \
#       [--c-min 3] [--c-max 14] [--max-knots N] [--timeout 60] [--save-every 200]
#
# Resume support: existing JSON is read and skips already computed
# knots. Crash recovery via periodic save.
#
# Convention: every DB entry = {
#   'name': '12n_123',          # knots.sage name
#   'c': 12,
#   'theta_at_22_7__21_13': '...',     # str (Frac or Int)
#   'delta1_at_22_7': '...',
#   'duration_s': 0.42,
#   'error': None or message
# }

import os, re, sys, json, time
from pathlib import Path

import site_config
THETA_DIR = site_config.THETA_DIR
OUT_PATH  = site_config.ROOT + "/external_tools/theta_db/knots_sage_theta_db.json"

# ----- ARG-PARSE -------------------------------------------------------------
# (Sage's argparse handling is peculiar — we use sys.argv directly)
def _arg(name, default):
    for i, a in enumerate(sys.argv):
        if a == name and i+1 < len(sys.argv):
            return sys.argv[i+1]
    return default

C_MIN     = int(_arg("--c-min", 3))
C_MAX     = int(_arg("--c-max", 14))
MAX_KNOTS = int(_arg("--max-knots", 0))   # 0 = unlimited
TIMEOUT   = int(_arg("--timeout", 60))
SAVE_EVERY = int(_arg("--save-every", 200))
OUT_OVR   = _arg("--out", None)
if OUT_OVR:
    OUT_PATH = OUT_OVR

print(f"[build] c_range=[{C_MIN},{C_MAX}], max_knots={MAX_KNOTS or 'inf'}, timeout={TIMEOUT}s, save_every={SAVE_EVERY}")
print(f"[build] output: {OUT_PATH}")

# ----- LOAD KNOTS + THETA ----------------------------------------------------
print(f"[build] loading knots.sage ...")
t_load = time.time()
load(f"{THETA_DIR}/knots.sage")
print(f"[build] knots.sage loaded in {time.time()-t_load:.1f}s, {len(Knot)} entries")

print(f"[build] loading Theta_numeric.sage ...")
load(f"{THETA_DIR}/Theta_numeric.sage")
print(f"[build] Theta_numeric.sage loaded")

# ----- LOAD EXISTING DB (resume) --------------------------------------------
db = {}
if Path(OUT_PATH).exists():
    try:
        with open(OUT_PATH) as fh:
            db = json.load(fh)
        print(f"[build] resume: loaded {len(db)} existing entries")
    except Exception as e:
        print(f"[build] could not load existing DB ({e}); starting fresh")
        db = {}

# ----- FILTER KNOT NAMES -----------------------------------------------------
def name_to_c(n):
    m = re.match(r"^(\d+)", n)
    return int(m.group(1)) if m else None

names_todo = []
for n in Knot.keys():
    c = name_to_c(n)
    if c is None: continue
    if c < C_MIN or c > C_MAX: continue
    if n in db and db[n].get("theta_at_22_7__21_13") is not None and not db[n].get("error"):
        continue
    names_todo.append((c, n))
names_todo.sort()
if MAX_KNOTS > 0:
    names_todo = names_todo[:MAX_KNOTS]
print(f"[build] {len(names_todo)} knots to compute (already done: {len(db)})")

# ----- T1, T2 default -----
T1 = QQ(22)/7
T2 = QQ(21)/13

# ----- ITERATE + COMPUTE -----------------------------------------------------
t_global = time.time()
n_ok = 0; n_err = 0
last_save = 0

import signal
class _Timeout(Exception): pass
def _alarm_handler(signum, frame):
    raise _Timeout()
signal.signal(signal.SIGALRM, _alarm_handler)

def _save_db():
    global last_save
    Path(OUT_PATH).parent.mkdir(parents=True, exist_ok=True)
    tmp = OUT_PATH + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(db, fh, indent=None)
    os.replace(tmp, OUT_PATH)
    last_save = len(db)

for k_idx, (c, name) in enumerate(names_todo):
    Cs, phi = Knot[name]
    t0 = time.time()
    signal.alarm(TIMEOUT)
    try:
        D1, theta, theta_total, expo = Theta_numeric([Cs, phi], t1=T1, t2=T2)
        signal.alarm(0)
        entry = {
            "name": name, "c": c,
            "theta_at_22_7__21_13": str(theta),
            "delta1_at_22_7": str(D1),
            "duration_s": float(time.time() - t0),
            "error": None,
        }
        n_ok += 1
    except _Timeout:
        signal.alarm(0)
        entry = {"name": name, "c": c,
                 "theta_at_22_7__21_13": None, "delta1_at_22_7": None,
                 "duration_s": float(time.time() - t0),
                 "error": f"timeout after {TIMEOUT}s"}
        n_err += 1
    except Exception as e:
        signal.alarm(0)
        entry = {"name": name, "c": c,
                 "theta_at_22_7__21_13": None, "delta1_at_22_7": None,
                 "duration_s": float(time.time() - t0),
                 "error": repr(e)[:300]}
        n_err += 1
    db[name] = entry

    if (k_idx + 1) % 50 == 0 or k_idx == len(names_todo) - 1:
        elapsed = time.time() - t_global
        rate = (k_idx + 1) / max(elapsed, 0.001)
        eta = (len(names_todo) - k_idx - 1) / max(rate, 0.001)
        print(f"[build] {k_idx+1}/{len(names_todo)}  ok={n_ok} err={n_err}  "
              f"rate={rate:.1f}/s  ETA={eta/60:.1f}min  last={name}")

    if len(db) - last_save >= SAVE_EVERY:
        _save_db()
        print(f"[build]   ... saved DB ({len(db)} entries)")

_save_db()
elapsed = time.time() - t_global
print(f"[build] DONE. Total elapsed: {elapsed/60:.2f}min. Final DB size: {len(db)} entries (ok={n_ok}, err={n_err}).")
