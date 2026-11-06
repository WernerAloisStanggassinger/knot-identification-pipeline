#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""
257_widmark_polystix / HEXASTIX — knots from Hart's 72-pencil packing.

Sister of polystix_gen.py (there: 6 rods, tetrastix/Borromean, 3,840 closures
fully enumerated). Here: the FOUR-direction packing, from which Hart "72 Pencils"
and Widmark builds his "Hexastix ... Knot" glass works.

QUELLEN
  G. Hart, "72 Pencils", georgehart.com/sculpture/pencils.html — verbatim:
    "an arrangement of four intersecting hexagonal tubes that penetrate each other"
    "The four hexagonal tubes are hollow ... the cavity is the rhombic dodecahedron"
  A. Widmark, Bridges 2021, pp. 293-296 — connection rule: "twists" between
    neighbors on DIFFERENT axes, "half loops" between neighbors on the SAME one.
  ⚠️ Neither Hart nor Widmark identifies one of these knots. The "72" is the
     rod count, NOT a crossing number (receipt: BELEG_was_72_bedeutet_20260726.md).

THE GEOMETRY — derived, not guessed (26.07.2026)
  Four tubes along the space diagonals (1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1).
  Per tube a HOLLOW hexagon ring in the triangular lattice of the perpendicular plane:
  hex radius R -> 6R rods.  R=3 -> 18 per tube -> 4x18 = 72 = Hart's number.

  Built concentrically, ALL tube pairs intersect. The distance of two axes is
  |(c2-c1)*n|/|n| with n = d1 x d2; for the four space diagonals the six n
  are exactly the six FACE diagonals, and the lattice part is always INTEGER
  (recomputed) — so it hits zero. An offset per tube must shift it to 1/2.

  🔑 These six conditions are linear and the matrix has rank 5 with 6 rows; its
  left null space is EXACTLY the all-ones vector. Hence necessarily:
      the six pair offsets must add up to ZERO
      -> not all six can be +1/2; three +1/2, three -1/2.
  Solution (each 1/4 of a face diagonal):
      u0=(0,0,0)  u1=(0,-1/4,+1/4)  u2=(+1/4,0,-1/4)  u3=(-1/4,+1/4,0)
  Measured on the built object: smallest axis distance 1/(2*sqrt(2)) = 0.353553,
  parallel rods sqrt(2/3) = 0.816497.  Rod radius therefore < 0.1768.

THE COMBINATORICS — why NOT enumerated here
  Closures to one cycle = (n-1)!*2^(n-1):
      n= 6 ->    3,840          (polystix_gen.py enumerates them completely)
      n=24 -> 2.2e29
      n=72 -> 2.0e123
  So: neighbourhood bound (Widmark's rule) + directed search, no counting.
  Sought is a Hamilton cycle in the rod graph: each rod exactly once, entered
  at one end, left at the other, connectors only between neighbor ends.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "the_knotfamily_machine"))
from fseries_io import polyline_to_fourier, write_fseries  # noqa: E402

from polystix_gen import (_bezier, check_rods_disjoint, min_self_distance)  # noqa: E402

RAW = [np.array(v, float) for v in [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]]
DIRS = [m / np.linalg.norm(m) for m in RAW]
# derived in hexsolve2 (see header); each 1/4 of a face diagonal
OFFS = [np.array(v, float) for v in
        [(0, 0, 0), (0, -0.25, 0.25), (0.25, 0, -0.25), (-0.25, 0.25, 0)]]


def basis_perp(d):
    """Projection of the cubic lattice onto the plane perpendicular to d = triangular lattice."""
    E = np.eye(3)
    P = [e - np.dot(e, d) * d for e in E]
    return P[0], P[1]


def hex_ring(R: int):
    """Points of the triangular lattice with hex distance EXACTLY R  ->  6R pieces (hollow tube)."""
    return [(i, j) for i in range(-R, R + 1) for j in range(-R, R + 1)
            if max(abs(i), abs(j), abs(i + j)) == R]


def rods_hexastix(ring: int = 3, half_len: float = 1.25, lam: float = 1.0):
    """Hart's four hollow hexagonal tubes. ring=3 -> 72 rods (4x18).

    ⚠️ half_len is NOT free. Measured on the finished object (26.07.):
        half_len 6.00 -> ends with a neighbor in a DIFFERENT tube:   0 of 144
        half_len 1.25 -> ends with a neighbor in a different tube: 144 of 144
    With Hart's long pencils the ends do not see each other; a circle across
    all four tubes then cannot exist at all, only half loops in the own tube.
    That is the factual difference between Hart's PACKING and Widmark's KNOTS:
    to knot, the rods must be short.
    """
    rods, tube = [], []
    for k, d in enumerate(DIRS):
        a, b = basis_perp(d)
        for i, j in hex_ring(ring):
            c = lam * (i * a + j * b) + OFFS[k]
            rods.append((c - half_len * d, c + half_len * d))
            tube.append(k)
    return rods, tube


# ----------------------------------------------------------------------------
# neighbourhood and search
# ----------------------------------------------------------------------------

def end_data(rods):
    ends, tang, owner = [], [], []
    for i, (p0, p1) in enumerate(rods):
        u = (p1 - p0) / np.linalg.norm(p1 - p0)
        ends += [p0, p1]
        tang += [-u, +u]
        owner += [i, i]
    return np.array(ends), np.array(tang), np.array(owner)


def neighbour_edges(ends, owner, max_conn: float):
    """Allowed connectors: ends of different rods closer than max_conn."""
    D = np.linalg.norm(ends[:, None, :] - ends[None, :, :], axis=-1)
    adj = {e: [] for e in range(len(ends))}
    for a, b in itertools.combinations(range(len(ends)), 2):
        if owner[a] != owner[b] and D[a, b] <= max_conn:
            adj[a].append(b)
            adj[b].append(a)
    return adj, D


def find_cycle(adj, owner, n_rods, rng, budget=400_000):
    """Randomized depth-first search with backtracking: Hamilton cycle over all rods.

    State: current end, used rods, chosen connectors.
    From one end it ALWAYS goes first through the rod (end ^ 1), then via a
    free connector to an unused rod. At the very end the connector must
    back to the start end must exist.
    """
    start = 0
    steps = [0]

    def rec(cur_end, used, chosen):
        steps[0] += 1
        if steps[0] > budget:
            return None
        other = cur_end ^ 1                       # through the rod
        if len(used) == n_rods:                   # all rods used -> close
            return chosen + [(other, start)] if start in adj[other] else None
        cands = [e for e in adj[other] if owner[e] not in used]
        rng.shuffle(cands)
        for e in cands:
            r = rec(e, used | {owner[e]}, chosen + [(other, e)])
            if r is not None:
                return r
        return None

    return rec(start, {owner[start]}, [])


def closure_polyline(rods, matching, per_rod=26, per_arc=30, bulge=0.42):
    ends, tang, _ = end_data(rods)
    partner = {}
    for a, b in matching:
        partner[a] = b
        partner[b] = a
    pts, cur = [], 0
    for _ in range(len(rods)):
        other = cur ^ 1
        seg = np.linspace(0.0, 1.0, per_rod, endpoint=False)[:, None]
        pts.append(ends[cur] + seg * (ends[other] - ends[cur]))
        nxt = partner[other]
        pts.append(_bezier(ends[other], tang[other], ends[nxt], tang[nxt], bulge, per_arc))
        cur = nxt
    return np.vstack(pts)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ring", type=int, default=3, help="hex radius: 1->24, 2->48, 3->72 rods")
    ap.add_argument("--radius", type=float, default=0.15, help="rod radius (bound 0.1768)")
    ap.add_argument("--half-len", type=float, default=1.25,
                    help="half rod length; >1.5 separates the tubes from each other (see rods_hexastix)")
    ap.add_argument("--max-conn", type=float, default=1.6, help="neighbour bound of the connectors")
    ap.add_argument("--emit", type=int, default=0)
    ap.add_argument("--terms", type=int, default=96)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--budget", type=int, default=400_000, help="search steps per attempt")
    ap.add_argument("--outdir", default="out_hex")
    ap.add_argument("--stats-only", action="store_true")
    args = ap.parse_args()

    rods, tube = rods_hexastix(args.ring, half_len=args.half_len)
    gap = check_rods_disjoint(rods, args.radius)
    ends, tang, owner = end_data(rods)
    print(f"tubes: 4 · rods per tube: {len(rods)//4} · total {len(rods)} · half rod length {args.half_len}")
    print(f"smallest axis distance {gap:.6f} (rod radius {args.radius} -> intersection-free: "
          f"{'YES' if gap > 2*args.radius else 'NO'})")

    adj, D = neighbour_edges(ends, owner, args.max_conn)
    grad = [len(adj[e]) for e in range(len(ends))]
    print(f"connector bound {args.max_conn}: degree per end min {min(grad)} / "
          f"median {int(np.median(grad))} / max {max(grad)}")
    tub = np.array([tube[o] for o in owner])
    n_tw = sum(1 for e in adj for f in adj[e] if tub[e] != tub[f]) // 2
    n_hl = sum(1 for e in adj for f in adj[e] if tub[e] == tub[f]) // 2
    print(f"allowed connectors: {n_tw} twist (cross-tube) · {n_hl} half_loop (within tube)")
    if n_tw == 0:
        print("⚠️  NO cross-tube connector — a circle across all tubes is impossible.")
    if min(grad) == 0:
        print("⚠️  Ends without any neighbor — bound too small, no circle possible.")
    if args.stats_only:
        return

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    manifest, written, versuche = [], 0, 0
    while written < args.emit and versuche < 200:
        versuche += 1
        m = find_cycle(adj, owner, len(rods), rng, args.budget)
        if m is None:
            continue
        P = closure_polyline(rods, m)
        dmin = min_self_distance(P)
        if dmin <= 1e-3:
            print(f"  attempt {versuche}: self-intersection ({dmin:.5f}) — discarded")
            continue
        kinds = ["half_loop" if tube[owner[a]] == tube[owner[b]] else "twist" for a, b in m]
        name = (f"widmark_hex{len(rods)}_{args.seed}_{written:03d}_"
                f"tw{kinds.count('twist')}_hl{kinds.count('half_loop')}")
        a_, b_ = polyline_to_fourier(P, n_terms=args.terms)
        write_fseries(outdir / f"{name}.fseries", a_, b_, knot_name=name,
                      source="257_widmark_polystix/hexastix")
        manifest.append(dict(name=name, ring=args.ring, n_rods=len(rods),
                             matching=[list(p) for p in m], kinds=kinds,
                             min_self_dist=round(float(dmin), 5)))
        written += 1
        print(f"  ✓ {name}  (min. self-distance {dmin:.4f})")

    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=1))
    print(f"written: {written} .fseries in {outdir}/ (after {versuche} search attempts)")


if __name__ == "__main__":
    main()
