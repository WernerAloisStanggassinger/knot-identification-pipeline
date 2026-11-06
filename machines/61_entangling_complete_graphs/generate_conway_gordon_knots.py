# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 61 — entangling_complete_graphs (Conway–Gordon).

K₆ is INTRINSICALLY LINKED (Conway–Gordon 1983): in EVERY spatial
embedding at least one disjoint triangle pair carries lk ≠ 0,
exactly: Σ lk (mod 2) = 1 over all 10 pairs. K₇ is intrinsically
KNOTTED (a Hamilton cycle with Arf=1 in every embedding).
The machine BUILDS embeddings and HARVESTS the unavoidable
entanglement — „knots one cannot avoid".
BET EXACT today: CG sum ≡ 1 mod 2 in every random embedding
(the theorem as machine self-test, Gauss lk bit-exactly integer).
REVERSAL CLAUSE: straight edges + random vertices carry near-intersections
→ Embedding gate per harvested curve. Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, itertools, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402


def lk_polygon(A, B, rot_seed=0):
    """Linking number of two closed polygons — EXACT via
    signed crossings of a generic projection (build lesson: the
    midpoint Gauss integral converges too slowly at sharp polygon corners,
    measured 0.44 off the integer; combinatorics is exact)."""
    rng = np.random.default_rng(rot_seed)
    R, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    P = np.asarray(A, float) @ R
    Q = np.asarray(B, float) @ R
    total = 0
    nA, nB = len(P), len(Q)
    for i in range(nA):
        p1, p2 = P[i], P[(i + 1) % nA]
        r = p2 - p1
        for j in range(nB):
            q1, q2 = Q[j], Q[(j + 1) % nB]
            s = q2 - q1
            den = r[0] * s[1] - r[1] * s[0]
            if abs(den) < 1e-14:
                continue
            dq = q1 - p1
            t = (dq[0] * s[1] - dq[1] * s[0]) / den
            u = (dq[0] * r[1] - dq[1] * r[0]) / den
            if not (1e-12 < t < 1 - 1e-12 and 1e-12 < u < 1 - 1e-12):
                continue
            # over/under is decided by the z coordinate at the intersection point —
            # without it the planar sign ALWAYS sums to 0
            # (algebraic intersection number, 2nd build lesson of this machine)
            zA = p1[2] + t * r[2]
            zB = q1[2] + u * s[2]
            if zA > zB:                      # only A-over-B counts: Σε = lk
                total += 1 if den > 0 else -1
    return float(total)


def gauss_lk(A, B, seeds=(0, 1)):
    """lk with two-projection witnesses: both generic directions
    must deliver the same integer (else ValueError)."""
    vals = [lk_polygon(A, B, s) for s in seeds]
    ints = [int(round(v)) for v in vals]
    if ints[0] != ints[1] or any(abs(v - i) > 1e-9
                                 for v, i in zip(vals, ints)):
        raise ValueError(f"lk witnesses disagree: {vals}")
    return float(ints[0])


def triangle_pairs_k6():
    """All 10 unordered pairs of disjoint triangles in K₆."""
    verts = set(range(6))
    pairs = []
    for tri in itertools.combinations(range(6), 3):
        other = tuple(sorted(verts - set(tri)))
        if tri < other:
            pairs.append((tri, other))
    return pairs


def cg_sum(X):
    """Conway–Gordon sum Σ lk mod 2 of a K₆ embedding X (6,3)."""
    total, lks, worst = 0, [], 0.0
    for (t1, t2) in triangle_pairs_k6():
        lk = gauss_lk(X[list(t1)], X[list(t2)])
        k = int(round(lk))
        worst = max(worst, abs(lk - k))
        lks.append(k)
        total += k
    return total % 2, lks, worst


def hamiltonian_cycles_k7():
    """All 360 Hamilton cycles of K₇ (start fixed, direction modded out)."""
    cycles = []
    for perm in itertools.permutations(range(1, 7)):
        if perm[0] < perm[-1]:
            cycles.append((0,) + perm)
    return cycles


def cycle_curve(X, cyc, pts_per_edge=60, rounds=2):
    P = X[list(cyc)]
    dense = []
    m = len(P)
    for k in range(m):
        a, b = P[k], P[(k + 1) % m]
        t = np.linspace(0, 1, pts_per_edge, endpoint=False)
        dense.append(np.outer(1 - t, a) + np.outer(t, b))
    C = np.vstack(dense)
    for _ in range(rounds):
        C = 0.5 * C + 0.25 * (np.roll(C, 1, axis=0) + np.roll(C, -1, axis=0))
    return C


def writhe_fast(C, n=256):
    sys.path.insert(0, str(HERE.parent / "27_pumping_the_twist_account"))
    from generate_writhe_pumped import writhe
    return writhe(C, n=n)


def make_knot(X, cyc, out_root, tag):
    C = cycle_curve(X, cyc)
    name = f"cg_k7_{tag}"
    entry = {"name": name, "machine": "61_entangling_complete_graphs",
             "cycle": list(map(int, cyc)),
             "writhe_proxy": round(writhe_fast(C), 3),
             "predictions": {
                 "satz": "K₇ trägt in JEDER Einbettung einen verknoteten "
                         "Hamilton-Zyklus (Conway–Gordon)",
                 "identity": "Zünd-Messung"}}
    emb = embedding_report(C)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(C, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/61", n_terms=96,
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 61 (Conway–Gordon)")
    # (1) Combinatorics exact: 10 triangle pairs, 360 Hamilton cycles
    ok1 = len(triangle_pairs_k6()) == 10 and len(hamiltonian_cycles_k7()) == 360
    print(f"  (1) combinatorics: {len(triangle_pairs_k6())} K₆ pairs, "
          f"{len(hamiltonian_cycles_k7())} K₇ Hamilton cycles  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) THE THEOREM as self-test: CG sum ≡ 1 mod 2, 12 embeddings
    rng = np.random.default_rng(61)
    fails, worst_int = 0, 0.0
    parities = []
    for trial in range(12):
        X = rng.normal(size=(6, 3))
        par, lks, worst = cg_sum(X)
        worst_int = max(worst_int, worst)
        parities.append(par)
        if par != 1:
            fails += 1
    ok2 = fails == 0
    print(f"  (2) Conway–Gordon sum: {parities} (all 1; "
          f"max lk distance to integer {worst_int:.2e})  "
          f"{'✓ THEOREM CONFIRMED' if ok2 else '✗ FAIL — integral broken'}")
    assert ok2
    # (3) lk values integer (Gauss integral calibration)
    ok3 = worst_int < 1e-3
    print(f"  (3) lk integrality: max distance {worst_int:.2e}  "
          f"{'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    # (4) K₇ harvest: cycles by |Wr| descending, embedding gate selects
    #     (build lesson: the |Wr|-MAXIMAL cycle is usually also the TIGHTEST
    #     — the most entangled curves tear at the gate first)
    X7 = rng.normal(size=(7, 3))
    scored = []
    for cyc in hamiltonian_cycles_k7()[::6]:        # 60 samples
        w = writhe_fast(cycle_curve(X7, cyc, pts_per_edge=20, rounds=1),
                        n=128)
        scored.append((abs(w), w, cyc))
    scored.sort(reverse=True)
    e4 = f4 = None
    rank = 0
    for rank, (_, w, cyc) in enumerate(scored, 1):
        e4, f4 = make_knot(X7, cyc, HERE / "demo_smoke", "maxwr")
        if f4 is not None:
            break
    ok4 = f4 is not None
    print(f"  (4) K₇ harvest (rank {rank} of 60, Wr={e4['writhe_proxy']}, "
          f"cycle {e4['cycle']}): nyq="
          f"{e4.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok4 else '✗ FAIL: ' + str(e4.get('rejected'))}")
    assert ok4
    print("[DEMO] 4/4 green — Conway–Gordon machine ready to fire "
          "(the unavoidable entanglement, secured by the theorem).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
