# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 46 — lacing_the_boots (lacings, padding rail).

Two eyelet columns, lacing = order permutation of the eyelets
(closed: loop arc at the top). Plane curve from rounded
segments → self-crossings → lifted ALTERNATING (39 heritage incl.
Gauss parity check). IGNITION QUESTION with everyday punch: which lacing
is the first knotted one? (Standard lacings presumably tightenable =
unknots — exactly that is what Monday measures.) Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "39_scribbling_then_lifting"))
sys.path.insert(0, str(HERE.parent))
from generate_scribble_knots import crossings, lift  # noqa: E402
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402


def eyelets(n, width=2.0, dy=1.0):
    L = [(0.0, i * dy) for i in range(n)]
    R = [(width, i * dy) for i in range(n)]
    return L, R


def lacing_path(sequence, n, width=2.0, dy=1.0, pts_per_seg=41):
    # 41 ODD: crossing at fraction 1/2 would otherwise fall exactly on vertices
    # of both diagonals (vertex trap, third sighting — cf. machine 02)
    """sequence: list ('L'|'R', index) — closed path incl. top arc."""
    Lh, Rh = eyelets(n, width, dy)
    P = [ (Lh[i] if side == 'L' else Rh[i]) for (side, i) in sequence ]
    P = np.array(P, float)
    # Round corners: Catmull-Rom-like smoothing via dense interpolation
    dense = []
    m = len(P)
    for k in range(m):
        a, b = P[k], P[(k + 1) % m]
        t = np.linspace(0, 1, pts_per_seg, endpoint=False)
        seg = np.outer(1 - t, a) + np.outer(t, b)
        dense.append(seg)
    C = np.vstack(dense)
    # slight smoothing (defuse corners)
    for _ in range(3):
        C = 0.5 * C + 0.25 * (np.roll(C, 1, axis=0) + np.roll(C, -1, axis=0))
    return C


def crisscross(n):
    """True criss-cross lacing: up on one diagonal (L0,R1,L2,…),
    down on the other (…,L3,R2,L1,R0) — each eyelet exactly once,
    diagonals cross in the middle field. (Build lesson: row-by-row lacing
    = crossing-free ladder.)"""
    up = [('L' if i % 2 == 0 else 'R', i) for i in range(n)]
    down = [('R' if i % 2 == 0 else 'L', i) for i in range(n - 1, -1, -1)]
    return up + down


def zigzag_skip(n, skip=2):
    """Variant: skip lacing (skip offset) — candidate for knotting."""
    seq = []
    side = 'L'
    visited_L = [i for i in range(0, n)]
    order = [(i * skip) % n for i in range(n)]
    for i in order:
        seq.append(('L', i))
        seq.append(('R', (i + skip) % n))
    return seq


def make_knot(kind, n, out_root, skip=2):
    seq = crisscross(n) if kind == "crisscross" else zigzag_skip(n, skip)
    C = lacing_path(seq, n)
    cross = crossings(C)
    name = f"la_{kind}{skip if kind != 'crisscross' else ''}_n{n}"
    entry = {"name": name, "machine": "46_lacing_the_boots",
             "kind": kind, "n_eyelets_per_col": n,
             "n_crossings": len(cross),
             "predictions": {
                 "alternating": "per Konstruktion (39-Hebung)",
                 "zuend_frage": "erste verknotete Schnürung? (Standard "
                                "vermutlich Unknot — Montag misst)"}}
    if len(cross) < 1:
        entry["rejected"] = "kreuzungsfrei"
        return entry, None
    P = lift(C, cross)
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(P, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/46",
                              n_terms=max(96, 12 * len(cross)),
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 46 (lacings)")
    seq = crisscross(5)
    ok1 = len(seq) == 10
    print(f"  (1) Criss-Cross n=5: {len(seq)} eyelet visits  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    C = lacing_path(seq, 5)
    cr = crossings(C)
    ok2 = len(cr) >= 3
    print(f"  (2) crossings in the middle field: {len(cr)}  "
          f"{'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    e, f = make_knot("crisscross", 5, HERE / "demo_smoke")
    ok3 = f is not None
    print(f"  (3) criss-cross built: {e['n_crossings']} crossings, "
          f"embedded={e.get('embedding',{}).get('embedded')}  "
          f"{'✓' if ok3 else '✗ FAIL: ' + str(e.get('rejected'))}")
    assert ok3
    e4, f4 = make_knot("skip", 5, HERE / "demo_smoke", skip=2)
    status = "gebaut" if f4 else f"abgelehnt ({e4.get('rejected')})"
    print(f"  (4) Skip lacing skip=2: {e4['n_crossings']} crossings → "
          f"{status}  ✓ (gate decides)")
    print("[DEMO] 4/4 green — lacing machine ready to fire "
          "(ignition question: the first knotted lacing).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
