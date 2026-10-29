#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Knot-Identification Pipeline — Werner Alois Stanggassinger & Claude Code (Anthropic).
"""sixd_crosscheck.py — THE STANDARD DOOR for the 6D geometry counter-check (Werner 2026-07-02).

ROLE (Werner's decision 2026-07-02): 6D is STANDARD from now on, but as a fast
geometric **pre-filter / 4th counter-check** — NOT as a distinctness referee.
The algebra (det/vol/θ/HOMFLY, corpus_db.distinct_knot_id) remains the truth; 6D
CONFIRMS it geometrically where it separates, and FLAGS the pairs where it is 6D-degenerate
(mutant-like) — there the algebra decides. No outside claim until mail-ready.

TRUTH ABOUT COSTS (verified 2026-07-02): the robust 6D scan (SONO-canonical_flow,
`eich_vollkorpus_robust.py`) takes 6–13 min/knot — NOT faster than the algebra. What is
fast is the COMPARISON of already computed fingerprints (ms). This tool therefore reads
the existing store (`out2_robust/`, 269 fingerprints, full calibration 2026-06-13); it does not
SCAN. Scanning new knots = a separate optional, cost-flagged stage.

BASIS: imports the corrected, calibrated logic from `sixd_machine/sixd_clean_analysis.py`
(machine 236 — validity converged∧¬deep_dead, stable subset, all corrections 205-235).
Report v4 / Disc #276 FINAL: Wr+ropelength p5/p95=1.41 separates, floppy=0, mutant-blind.

Modes:
  (default)          Corpus counter-check: separation report + disagreement report (agreement 6D↔algebra).
  --disagreements    ONLY the 6D-degenerate distinct pairs (mutant-like → algebra decides).
  --knot <stem|did>  Per knot: fingerprint + nearest 6D neighbours + whether the algebra agrees.
  --subset {optimal,stable,all}  Distance subset (default optimal = Wr+ropelength).
"""
import argparse
import site_config
import os
import statistics
import sys

SCANNER_DIR = site_config.SIXD_DIR
sys.path.insert(0, SCANNER_DIR)
import sixd_clean_analysis as sca  # noqa: E402  — source of the calibrated corrections (236)

SUBSETS = {"optimal": sca.OPTIMAL, "stable": sca.STABLE, "all": sca.ALLC}


def _znorm(raw, comps):
    """z-norm vectoriser + distance, identical to sca.sep_ratio (same metric)."""
    fps = [fp for ms in raw.values() for fp in ms]
    mean = {c: statistics.mean([fp[c] for fp in fps]) for c in comps}
    sd = {c: (statistics.pstdev([fp[c] for fp in fps]) or 1.0) for c in comps}

    def vec(fp):
        return [(fp[c] - mean[c]) / sd[c] for c in comps]

    def dist(a, b):
        return sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5
    return vec, dist


def _centroids(raw, comps):
    vec, dist = _znorm(raw, comps)
    cents, within = {}, []
    for did, ms in raw.items():
        vs = [vec(fp) for fp in ms]
        cents[did] = [statistics.mean(x) for x in zip(*vs)]
        for i in range(len(vs)):
            for j in range(i + 1, len(vs)):
                within.append(dist(vs[i], vs[j]))
    within.sort()
    # 95th percentile of the within-spread = collision threshold (two distinct knots
    # whose centroids lie CLOSER than the typical realisation spread → 6D-indistinguishable).
    thr = within[int(len(within) * 0.95)] if within else 0.0
    return cents, dist, thr


def cmd_disagreements(raw, comps, verbose=True):
    cents, dist, thr = _centroids(raw, comps)
    dids = list(cents)
    pairs = []
    for i in range(len(dids)):
        for j in range(i + 1, len(dids)):
            dd = dist(cents[dids[i]], cents[dids[j]])
            if dd < thr:
                pairs.append((dd, dids[i], dids[j]))
    pairs.sort()
    if verbose:
        print(f"  Collision threshold (95% within-spread): {thr:.3f} [z units, subset {comps}]")
        print(f"  6D-DEGENERATE distinct pairs (algebra is referee): {len(pairs)}")
        for dd, a, b in pairs:
            print(f"    ⚠ {dd:.3f}  {a}  ≈6D≈  {b}")
        if not pairs:
            print("    (none — 6D separates ALL scanned distinct knots geometrically)")
    return pairs, thr, len(dids)


def cmd_knot(raw, comps, target):
    # target = stem or distinct_knot_id; the store is grouped by distinct_knot_id.
    cents, dist, thr = _centroids(raw, comps)
    if target not in cents:
        # Fuzzy: prefix/substring over distinct_knot_ids
        hits = [d for d in cents if target in d or d in target or d.startswith(target)]
        if len(hits) == 1:
            target = hits[0]
        elif hits:
            print(f"ambiguous — candidates: {hits[:8]}")
            return 1
        else:
            print(f"'{target}' has NO valid 6D fingerprint in the store "
                  f"(not scanned or filtered out as deep_dead/invalid).")
            return 1
    ms = raw[target]
    fp0 = {c: round(statistics.mean([fp[c] for fp in ms]), 4) for c in ("Wr", "ropelength", "ACN", "K_tot")}
    print(f"=== 6D counter-check: {target} ({len(ms)} realisation(s)) ===")
    print(f"  Fingerprint (centroid): {fp0}")
    nn = sorted(((dist(cents[target], cents[d]), d) for d in cents if d != target))[:5]
    print(f"  Nearest 6D neighbours (subset {comps}, threshold {thr:.3f}):")
    for dd, d in nn:
        flag = "  ⚠ 6D COLLISION → algebra decides" if dd < thr else ""
        print(f"    {dd:.3f}  {d}{flag}")
    print("  Algebra (distinct_knot_id) is and remains the referee of this identity.")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--disagreements", action="store_true")
    ap.add_argument("--knot", type=str, default=None)
    ap.add_argument("--subset", choices=list(SUBSETS), default="optimal")
    args = ap.parse_args()
    comps = SUBSETS[args.subset]

    raw, n_total, n_dd, n_valid = sca.load()
    if args.knot:
        sys.exit(cmd_knot(raw, comps, args.knot))

    if args.disagreements:
        print(f"=== 6D DISAGREEMENTS (subset {args.subset}={comps}, {len(raw)} distinct scanned) ===")
        cmd_disagreements(raw, comps, verbose=True)
        sys.exit(0)

    print("=== 6D GEOMETRY COUNTER-CHECK (standard door, Werner 2026-07-02) ===")
    print(f"Store: {n_total} fingerprints | deep_dead junk: {n_dd} | VALID: {n_valid} "
          f"| distinct knots (scanned): {len(raw)}")
    print(f"Subset for distance: {args.subset} = {comps}\n")
    print("Separation (p5/p95 > 1 = 6D confirms algebra geometrically):")
    for tag, cs in [("all 5 (raw)", sca.ALLC), ("STABLE ACN+K_tot+Wr", sca.STABLE),
                    ("OPTIMAL Wr+ropelength", sca.OPTIMAL)]:
        r = sca.sep_ratio(raw, cs)
        if r:
            v = "SEPARATES (>1)" if r[0] > 1 else "Overlap (<1)"
            print(f"  {tag:22} p5/p95={r[0]:.2f} ({v})  median-ratio={r[1]:.1f}")
    print("\n6D↔algebra agreement:")
    pairs, thr, ndid = cmd_disagreements(raw, comps, verbose=True)
    agree = ndid - len(set([p[1] for p in pairs] + [p[2] for p in pairs]))
    print(f"\nCONCLUSION: {agree}/{ndid} distinct knots are unambiguously separated in 6D "
          f"(algebra confirmed geometrically); {len(pairs)} degenerate pairs remain "
          f"a matter for the algebra (mutant-blind, Disc #276). 6D = counter-check, not referee.")


if __name__ == "__main__":
    main()
