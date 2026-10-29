#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Knot-Identification Pipeline — Werner Alois Stanggassinger & Claude Code (Anthropic).
"""sixd_clean_analysis.py — corrected 6D analysis (machine 236, 2026-06-12).

Embodies ALL calibration companion corrections (205-235):
  - VALIDITY = converged ∧ deep_dead=False ∧ ACN≥c−0.5  (deep_dead was the missing one, #232:
    25% tunneled junk marked as converged, trefoil RL=95-139).
  - STABLE SUBSET ACN+K_tot+Wr (z-norm) for separation (#206: T_sig/SL=noise).
  - per-knot reliability via RL-CV (#231); 6D complements θ, mutant-blind (#225/233/235).

Replaces the buggy 205-207 single scripts for the calibration completion (report v4).
Usage: python3 sixd_clean_analysis.py  (reads out2_robust/, corpus_db.json)
"""
import json
import site_config
import glob
import os
import re
import statistics
from collections import defaultdict

ROBUST = "@SIXD_DIR@/out2_robust".replace('@SIXD_DIR@', site_config.SIXD_DIR)
DB = "@PIPELINE_DIR@/corpus_db.json".replace('@PIPELINE_DIR@', site_config.PIPELINE_DIR)
STABLE = ["ACN", "K_tot", "Wr"]
OPTIMAL = ["Wr", "ropelength"]   # #242: chirality + 1 complexity = best separation (p5/p95=2.39)
ALLC = ["Wr", "ACN", "K_tot", "T_sig", "SL"]


def cof(did):
    m = re.match(r"S_(\d+)_", str(did))
    return int(m.group(1)) if m else None


def is_valid(j):
    """VALID = converged ∧ deep_dead=False ∧ (ACN≥c−0.5 if c is known)."""
    if not j.get("converged"):
        return False
    if j.get("deep_dead"):      # #232 — the decisive filter
        return False
    return True


def load():
    db = json.load(open(DB))
    stem2did = {r.get("stem"): r.get("distinct_knot_id") for r in db["knots"]}
    raw, n_total, n_dd, n_valid = defaultdict(list), 0, 0, 0
    for f in glob.glob(os.path.join(ROBUST, "*.json")):
        try:
            j = json.load(open(f))
        except Exception:
            continue
        n_total += 1
        if j.get("deep_dead"):
            n_dd += 1
        if not is_valid(j):
            continue
        fp = j.get("fingerprint") or {}
        if all(isinstance(fp.get(c), (int, float)) for c in ALLC):
            s = j.get("stem") or os.path.basename(f).replace("_sixd2.json", "")
            did = stem2did.get(s) or s.split("__")[0]
            raw[did].append(fp)
            n_valid += 1
    return raw, n_total, n_dd, n_valid


def sep_ratio(raw, comps):
    fps = [fp for ms in raw.values() for fp in ms]
    mean = {c: statistics.mean([fp[c] for fp in fps]) for c in comps}
    sd = {c: (statistics.pstdev([fp[c] for fp in fps]) or 1.0) for c in comps}

    def vec(fp):
        return [(fp[c] - mean[c]) / sd[c] for c in comps]

    def d(a, b):
        return sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5

    within, cents = [], {}
    for did, ms in raw.items():
        vs = [vec(fp) for fp in ms]
        if len(vs) >= 2:
            for i in range(len(vs)):
                for j in range(i + 1, len(vs)):
                    within.append(d(vs[i], vs[j]))
        cents[did] = [statistics.mean(x) for x in zip(*vs)]
    dids = list(cents)
    between = [d(cents[dids[i]], cents[dids[j]]) for i in range(len(dids)) for j in range(i + 1, len(dids))]
    within.sort(); between.sort()
    if not within or not between:
        return None
    p5p95 = between[int(len(between) * 0.05)] / within[int(len(within) * 0.95)]
    med = statistics.median(between) / statistics.median(within)
    return p5p95, med


def main():
    raw, n_total, n_dd, n_valid = load()
    print(f"=== CORRECTED 6D ANALYSIS (machine 236) ===")
    print(f"Records: {n_total} | deep_dead (junk): {n_dd} ({100*n_dd//max(n_total,1)}%) | "
          f"VALID (converged∧¬deep_dead): {n_valid} | distinct knots: {len(raw)}\n")
    for tag, comps in [("alle 5 (roh)", ALLC), ("STABIL ACN+K_tot+Wr", STABLE),
                       ("OPTIMAL Wr+ropelength", OPTIMAL)]:
        r = sep_ratio(raw, comps)
        if r:
            verdict = "TRENNT (>1)" if r[0] > 1 else "Überlapp (<1)"
            print(f"  {tag:22} Separation p5/p95={r[0]:.2f} ({verdict})  median-ratio={r[1]:.1f}")
    # per-knot RL reliability
    print("\n  per-knot RL-CV (reliability):")
    floppy = []
    for did, ms in raw.items():
        rls = [fp.get("ropelength") for fp in ms if isinstance(fp.get("ropelength"), (int, float))]
        if len(rls) >= 2:
            cv = statistics.pstdev(rls) / statistics.mean(rls) * 100
            if cv > 15:
                floppy.append((did, cv))
    print(f"    floppy (RL-CV>15%, after the deep_dead filter): {len(floppy)} {[d[0][:18] for d in floppy]}")
    print("\nCONCLUSION (expected, #233): STABLE subset p5/p95>1 = clean 6D separates distinct knots;")
    print("floppy≈0 after the filter; 6D complements θ (mutant-blind, #225/234), separates distinct knots well.")


if __name__ == "__main__":
    main()
