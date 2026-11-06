#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""
257_widmark_polystix — knot machine after Anduriel Widmark's construction method.

SOURCE (third-party method, cleanly separated from ours):
  A. Widmark, "Sculpture Design with Hexastix and Related Non-Intersecting
  Cylinder Packings", Bridges 2021, pp. 293-296.
  archive.bridgesmathart.org/2021/bridges2021-293.pdf
  Verbatim: "I started with six rods in a tetrastix/Borromean arrangement and
  found there are 3,840 unique paths or ways to connect the ends into a
  continuous loop. [...] The knots are constructed by using either 'twists' to
  join proximal neighbors that cross on different axes or 'half loops' to join
  neighbors on the same axis."

MACHINE FORM (Stanggassinger family):
  atomic module  = non-intersecting rod packing (tetrastix / hexastix)
  symmetry       = cubic space group of the packing
  closure constraint = Widmark's neighbour connection rule (twist / half loop)

CALIBRATION: the number 3,840 must fall out of the combinatorics, otherwise we have not
understood his system.  n rods to ONE closed cycle
connecting = (n-1)! * 2^(n-1);  n=6 -> 120 * 32 = 3840.

Output: .fseries (Fourier centreline) — the format the pipeline eats.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "the_knotfamily_machine"))
from fseries_io import polyline_to_fourier, write_fseries  # noqa: E402

# ----------------------------------------------------------------------------
# 1. The rod packing
# ----------------------------------------------------------------------------

def rods_borromean(half_len: float = 1.75, offset: float = 1.0):
    """Six rods, tetrastix/Borromean: two each parallel to x, y, z.

    Cyclic offset pattern (x rods offset in y, y rods in z,
    z rods in x). Pairwise axis distance = `offset`, hence intersection-free
    for rod radius < offset/2.

    Returns: list of (P0, P1) — start/end of each rod.
    """
    L, d = half_len, offset
    rods = []
    for s in (+1, -1):
        rods.append((np.array([-L, s * d, 0.0]), np.array([+L, s * d, 0.0])))  # || x
    for s in (+1, -1):
        rods.append((np.array([0.0, -L, s * d]), np.array([0.0, +L, s * d])))  # || y
    for s in (+1, -1):
        rods.append((np.array([s * d, 0.0, -L]), np.array([s * d, 0.0, +L])))  # || z
    return rods


def check_rods_disjoint(rods, radius: float) -> float:
    """Smallest axis distance between two rods (counter-check of the packing)."""
    best = np.inf
    for (a0, a1), (b0, b1) in itertools.combinations(rods, 2):
        u, v = a1 - a0, b1 - b0
        n = np.cross(u, v)
        nn = np.linalg.norm(n)
        if nn < 1e-12:                      # parallel
            w = b0 - a0
            dist = np.linalg.norm(w - np.dot(w, u) / np.dot(u, u) * u)
        else:
            dist = abs(np.dot(b0 - a0, n / nn))
        best = min(best, dist)
    if best <= 2 * radius:
        raise ValueError(f"rods collide: axis distance {best:.3f} <= 2r {2*radius:.3f}")
    return best


# ----------------------------------------------------------------------------
# 2. The closures (Widmark's 3,840)
# ----------------------------------------------------------------------------

def enumerate_single_cycle_closures(n_rods: int):
    """All pairings of the 2n rod ends that give EXACTLY ONE closed cycle.

    End (i,0)/(i,1) -> global index 2i / 2i+1.
    Yields per closure a tuple of end-point pairs.
    """
    ends = list(range(2 * n_rods))

    def matchings(rest):
        if not rest:
            yield ()
            return
        a = rest[0]
        for j in range(1, len(rest)):
            b = rest[j]
            remain = rest[1:j] + rest[j + 1:]
            for m in matchings(remain):
                yield ((a, b),) + m

    for m in matchings(ends):
        # run the cycle: rod i connects 2i<->2i+1, then the pairing
        partner = {}
        for a, b in m:
            partner[a] = b
            partner[b] = a
        seen, cur, steps = set(), 0, 0
        while True:
            rod_other = cur ^ 1          # through the rod
            seen.add(cur); seen.add(rod_other)
            cur = partner[rod_other]     # across the connection
            steps += 1
            if cur == 0:
                break
            if steps > 2 * n_rods:
                break
        if len(seen) == 2 * n_rods and steps == n_rods:
            yield m


# ----------------------------------------------------------------------------
# 3. Geometry of a closure
# ----------------------------------------------------------------------------

def _bezier(p0, t0, p1, t1, k: float, n: int):
    """Cubic Bézier that leaves p0 in direction t0 and reaches p1 from -t1.

    Models the bent glass rod: the wire runs out of the rod axis
    and bends into the neighbour's axis — Widmark's 'twist'/'half loop'.
    """
    d = np.linalg.norm(p1 - p0)
    c0 = p0 + t0 * k * d
    c1 = p1 + t1 * k * d
    s = np.linspace(0.0, 1.0, n, endpoint=False)[:, None]
    return ((1 - s) ** 3 * p0 + 3 * (1 - s) ** 2 * s * c0
            + 3 * (1 - s) * s ** 2 * c1 + s ** 3 * p1)


def closure_polyline(rods, matching, per_rod=40, per_arc=44, bulge=0.42):
    """Build the closed 3D polyline of a closure."""
    ends, tangents = [], []
    for p0, p1 in rods:
        u = (p1 - p0) / np.linalg.norm(p1 - p0)
        ends += [p0, p1]
        tangents += [-u, +u]              # outward-pointing rod direction
    partner = {}
    for a, b in matching:
        partner[a] = b
        partner[b] = a

    pts, cur = [], 0
    for _ in range(len(rods)):
        other = cur ^ 1
        seg = np.linspace(0.0, 1.0, per_rod, endpoint=False)[:, None]
        pts.append(ends[cur] + seg * (ends[other] - ends[cur]))     # rod
        nxt = partner[other]
        pts.append(_bezier(ends[other], tangents[other],
                           ends[nxt], tangents[nxt], bulge, per_arc))  # connector
        cur = nxt
    return np.vstack(pts)


def min_self_distance(P: np.ndarray, skip: int = 12) -> float:
    """Smallest distance between non-adjacent points of the curve."""
    n = len(P)
    d = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=-1)
    idx = np.arange(n)
    ring = np.abs(idx[:, None] - idx[None, :])
    ring = np.minimum(ring, n - ring)
    d[ring <= skip] = np.inf
    return float(d.min())


# ----------------------------------------------------------------------------
# 4. Neighbourhood (Widmark's restriction)
# ----------------------------------------------------------------------------

def end_geometry(rods):
    ends, axis = [], []
    for i, (p0, p1) in enumerate(rods):
        ends += [p0, p1]
        ax = int(np.argmax(np.abs(p1 - p0)))
        axis += [ax, ax]
    return np.array(ends), np.array(axis)


def classify_pair(ends, axis, a, b):
    """'half loop' = same axis, 'twist' = different axes."""
    return "half_loop" if axis[a] == axis[b] else "twist"


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--radius", type=float, default=0.30, help="rod radius")
    ap.add_argument("--max-conn", type=float, default=None,
                    help="Neighbor bound: only connectors shorter than this distance")
    ap.add_argument("--emit", type=int, default=0, help="how many .fseries to write")
    ap.add_argument("--terms", type=int, default=48, help="Fourier terms")
    ap.add_argument("--outdir", default="out")
    ap.add_argument("--stats-only", action="store_true")
    args = ap.parse_args()

    rods = rods_borromean()
    gap = check_rods_disjoint(rods, args.radius)
    ends, axis = end_geometry(rods)
    D = np.linalg.norm(ends[:, None, :] - ends[None, :, :], axis=-1)

    print(f"rods: {len(rods)} · smallest axis distance {gap:.3f} "
          f"(rod radius {args.radius}, hence intersection-free)")

    # --- calibration: Widmark's number-----------------------------------------
    closures = list(enumerate_single_cycle_closures(len(rods)))
    expected = math.factorial(len(rods) - 1) * 2 ** (len(rods) - 1)
    ok = "✓" if len(closures) == expected else "✗"
    print(f"Closures into ONE cycle: {len(closures)}  "
          f"(expected (n-1)!·2^(n-1) = {expected})  {ok}")
    print(f"CALIBRATION against Widmark 2021: he counts 3,840 — we count {len(closures)}. "
          f"{'AGREE' if len(closures) == 3840 else 'DEVIATION!'}")

    # --- connector lengths and neighbourhood---------------------------
    lens = sorted({round(float(D[a, b]), 4)
                   for m in closures[:400] for a, b in m})
    print(f"\noccurring connector distances (selection): {lens[:8]} …")

    if args.stats_only:
        return

    # --- geometric selection--------------------------------------------
    kept, rejected = [], 0
    for m in closures:
        if args.max_conn is not None and any(D[a, b] > args.max_conn for a, b in m):
            rejected += 1
            continue
        kept.append(m)
    print(f"\nafter neighbour bound {args.max_conn}: {len(kept)} of {len(closures)} "
          f"({rejected} discarded)")

    if not args.emit:
        return

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    written, collided = 0, 0
    manifest = []
    for i, m in enumerate(kept):
        if written >= args.emit:
            break
        P = closure_polyline(rods, m)
        dmin = min_self_distance(P)
        # Two SEPARATE questions, formerly wrongly mixed:
        #  (a) is the curve embedded?  -> decides the knot type
        #  (b) does it fit in glass of radius r? -> pure build question, later
        if dmin <= 1e-3:                     # real self-intersection
            collided += 1
            continue
        kinds = [classify_pair(ends, axis, a, b) for a, b in m]
        name = (f"widmark_bor6_{i:05d}_"
                f"tw{kinds.count('twist')}_hl{kinds.count('half_loop')}")
        a, b = polyline_to_fourier(P, n_terms=args.terms)
        write_fseries(outdir / f"{name}.fseries", a, b, knot_name=name,
                          source="257_widmark_polystix")
        manifest.append({"name": name, "matching": [list(p) for p in m],
                         "kinds": kinds, "min_self_dist": round(dmin, 4)})
        written += 1

    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=1))
    print(f"written: {written} .fseries in {outdir}/  "
          f"({collided} discarded due to self-contact)")


if __name__ == "__main__":
    main()
