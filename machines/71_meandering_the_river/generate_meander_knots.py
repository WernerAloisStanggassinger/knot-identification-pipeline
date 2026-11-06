# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 71 — meandering_the_river (meander knots).

Closed meander = Jordan curve × line with 2n intersections =
pair of non-crossing matchings (top/bottom) forming ONE cycle.
Counting sequence 1, 2, 8, 42 (A005315) = self-test of the own enumeration.
KNOT PATH: meander arcs + river line as a closed arc,
spliced at ONE crossing (two components → one curve), then
lifted alternatingly (inherited from 39) → 2n−1 crossings.
Nobody has mapped the meander census onto knot types.
REVERSAL CLAUSE: semicircle geometry = vertex trap (odd
subdivision, 46er lesson). Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from itertools import combinations
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "39_scribbling_then_lifting"))
sys.path.insert(0, str(HERE.parent))
from generate_scribble_knots import crossings, lift  # noqa: E402
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402


def noncrossing_matchings(n):
    """All non-crossing perfect matchings of 2n points
    (Catalan many), divide & conquer."""
    pts = list(range(2 * n))
    def rec2(rest):
        if not rest:
            return [[]]
        a = rest[0]
        res = []
        for i in range(1, len(rest), 2):
            b = rest[i]
            for mi in rec2(rest[1:i]):
                for mo in rec2(rest[i + 1:]):
                    res.append([(a, b)] + mi + mo)
        return res
    return rec2(pts)


def is_single_cycle(upper, lower, n):
    nxt_u = {}
    for (a, b) in upper:
        nxt_u[a] = b; nxt_u[b] = a
    nxt_l = {}
    for (a, b) in lower:
        nxt_l[a] = b; nxt_l[b] = a
    seen, cur, use_u = set(), 0, True
    while cur not in seen:
        seen.add(cur)
        cur = (nxt_u if use_u else nxt_l)[cur]
        use_u = not use_u
        if cur in seen and len(seen) < 2 * n:
            return False
    return len(seen) == 2 * n


def closed_meanders(n):
    M = noncrossing_matchings(n)
    return [(u, l) for u in M for l in M if is_single_cycle(u, l, n)]


def semicircle(a, b, up, n_pts=41):
    """Semicircle FROM a TO b (direction counts — build lesson: direction-
    blind arcs produce jump phantoms, 9 ghost crossings);
    ODD subdivision against the vertex trap."""
    c, r = (a + b) / 2.0, abs(b - a) / 2.0
    if up:
        th = np.linspace(np.pi, 0, n_pts) if a < b \
            else np.linspace(0, np.pi, n_pts)
    else:
        th = np.linspace(np.pi, 2 * np.pi, n_pts) if a < b \
            else np.linspace(2 * np.pi, np.pi, n_pts)
    return np.column_stack([c + r * np.cos(th),
                            r * np.sin(th)])


def meander_curve(upper, lower, n):
    """Planar meander curve as polyline (arc sequence in the cycle)."""
    nxt_u = {}
    for (a, b) in upper:
        nxt_u[a] = b; nxt_u[b] = a
    nxt_l = {}
    for (a, b) in lower:
        nxt_l[a] = b; nxt_l[b] = a
    pts, cur, use_u = [], 0, True
    for _ in range(2 * n):
        nx = (nxt_u if use_u else nxt_l)[cur]
        pts.append(semicircle(float(cur), float(nx), use_u)[:-1])
        cur, use_u = nx, not use_u
    return np.vstack(pts)


def river_loop(n, depth=4.5, y0=0.13):
    """The ,river‘ line as a closed arc: axis segment
    (OFFSET by y0 — vertex trap, 4th sighting: the line ran exactly
    through the 2n meander axis points, Gauss parity tore) +
    return arc BELOW all meander arcs (depth > n)."""
    x = np.linspace(-1.0, 2.0 * n, 161)
    line = np.column_stack([x, np.full_like(x, y0)])
    th = np.linspace(0, np.pi, 121)[1:-1]
    c, r = (2.0 * n - 1.0) / 2, (2.0 * n + 1.0) / 2
    back = np.column_stack([c + r * np.cos(th), y0 - depth * np.sin(th)])
    return np.vstack([line, back])


def splice(meander, river):
    """ORIENTED SMOOTHING at a true transversal crossing
    (build lesson: nearest-points seam lays two strands almost on top of each other
    → Gauss parity tears; smoothing a crossing is the clean
    merging of two components into one)."""
    for i in range(len(meander) - 1):
        p, r = meander[i], meander[i + 1] - meander[i]
        for j in range(len(river) - 1):
            q, s = river[j], river[j + 1] - river[j]
            den = r[0] * s[1] - r[1] * s[0]
            if abs(den) < 1e-12:
                continue
            dq = q - p
            t = (dq[0] * s[1] - dq[1] * s[0]) / den
            u = (dq[0] * r[1] - dq[1] * r[0]) / den
            if 0.15 < t < 0.85 and 0.15 < u < 0.85:
                # crossing found → merge with orientation
                return np.vstack([meander[:i + 1], river[j + 1:],
                                  river[:j + 1], meander[i + 1:]])
    return None


def make_knot(upper, lower, n, out_root, tag):
    M = meander_curve(upper, lower, n)
    R = river_loop(n)
    C2 = splice(M, R)
    if C2 is None:
        return {"name": f"me_n{n}_{tag}",
                "rejected": "keine transversale Spleiß-Kreuzung"}, None
    cross = crossings(C2)
    name = f"me_n{n}_{tag}"
    entry = {"name": name, "machine": "71_meandering_the_river",
             "n_pairs": n, "n_crossings": len(cross),
             "predictions": {
                 "alternating": "per 39-Hebung",
                 "zensus": "Mäander→Knoten-Karte = Zünd-Programm"}}
    if len(cross) < 3:
        entry["rejected"] = f"nur {len(cross)} Kreuzungen"
        return entry, None
    P = lift(C2, cross)
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(P, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/71",
                              n_terms=max(96, 12 * len(cross)),
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 71 (meanders)")
    # (1) CENSUS SELF-TEST: 1, 2, 8, 42 (A005315)
    counts = [len(closed_meanders(n)) for n in (1, 2, 3, 4)]
    ok1 = counts == [1, 2, 8, 42]
    print(f"  (1) meander count n=1..4: {counts} (literature: 1,2,8,42)  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) Meander is SIMPLE: planar curve has 0 self-crossings
    ms = closed_meanders(3)
    M = meander_curve(*ms[4], 3)
    n_self = len(crossings(M))
    ok2 = n_self == 0
    print(f"  (2) Jordan property: {n_self} self-crossings of the "
          f"meander curve  {'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) Splice + lift: knot built, crossings counted
    e3 = f3 = None
    for idx, (u, l) in enumerate(closed_meanders(4)):
        e3, f3 = make_knot(u, l, 4, HERE / "demo_smoke", f"i{idx}")
        if f3 is not None:
            break
    ok3 = f3 is not None
    print(f"  (3) {e3['name']}: {e3['n_crossings']} crossings after "
          f"splice, nyq={e3.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok3 else '✗ FAIL: ' + str(e3.get('rejected'))}")
    assert ok3
    print("[DEMO] 3/3 green — meander machine ready to fire "
          "(census 1,2,8,42 self-confirmed; knot map = ignition).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
