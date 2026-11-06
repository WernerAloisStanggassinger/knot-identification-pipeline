# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 100 — summing_the_knots (connected-sum machine).

Geometric #-operation at fseries level: separate the summands,
open at mutually NEAREST OUTER points, orientation-faithfully
bridge. Ladders (K, K#K, K#K#K) + the classical control pair
Granny (K3₁#K3₁) vs Square (K3₁#mK3₁): det EQUAL (9), sig
DIFFERENT (±6 vs 0) — distinguishable in-house.
MAIL THREAD: Dror (θ under # — anchor 1 is 3₁^#6!) + Stoimenow
(det multiplicativity as firing control).
BETS EXACT today: opening points = outer pair (gate);
orientation logic (granny=same, square=mirrored chirality);
sums embedded. FIRING: det(K₁#K₂)=det·det, genus additive,
granny/square sig separation, θ behaviour. Doctrine: do NOT fire.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "10_climbing_the_torus_ladder"))
sys.path.insert(0, str(HERE.parent))
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402
from the_knotfamily_machine.fseries_io import (resample_uniform,  # noqa: E402
                                               reconstruct_from_fourier)


def trefoil(mirror=False, n=2048):
    from generate_torus_ladder import torus_coeffs
    a, b = torus_coeffs(2, 3)
    P = reconstruct_from_fourier(a, b, n_samples=n)
    if mirror:
        P = P * np.array([1.0, 1.0, -1.0])
    return P


def outermost_point(P, direction):
    """Index of the outermost point in direction d (opening at the edge)."""
    return int(np.argmax(P @ np.asarray(direction, float)))


def connected_sum(P1, P2, gap=2.0):
    """K₁ # K₂: K₂ placed to the right of K₁, both opened at the mutually
    nearest OUTER point opened, connected with a double bridge."""
    P1 = resample_uniform(np.asarray(P1, float), 2048)
    P2 = resample_uniform(np.asarray(P2, float), 2048)
    d1 = float(np.linalg.norm(P1.max(0) - P1.min(0)))
    # shift K₂: to the right with distance gap·Ø
    off = np.array([P1[:, 0].max() - P2[:, 0].min() + gap, 0.0, 0.0])
    Q2 = P2 + off
    i1 = outermost_point(P1, [1, 0, 0])        # rightmost point of K₁
    i2 = outermost_point(Q2, [-1, 0, 0])       # leftmost point of K₂
    # WIDEN opening (build lesson: neighbour-point opening puts the
    # two bridge feet ds-close → gate trips): remove ±k points
    k = 30
    n1, n2 = len(P1), len(Q2)
    A1, B1 = P1[(i1 - k) % n1], P1[(i1 + k) % n1]
    A2, B2 = Q2[(i2 - k) % n2], Q2[(i2 + k) % n2]
    def bridge(p, q, lift, m=30):
        t = np.linspace(0, 1, m + 2)[1:-1]
        base = np.outer(1 - t, p) + np.outer(t, q)
        base[:, 2] += lift * np.sin(np.pi * t)
        return base
    # Orientation-faithful: K₁ up to (i1−k), bridge → K₂ from (i2+k) once
    # around to (i2−k), bridge back → K₁ from (i1+k)
    def span(P, a, b):
        n = len(P)
        a, b = a % n, b % n
        return P[a:b + 1] if a <= b else np.vstack([P[a:], P[:b + 1]])
    seq1a = span(P1, (i1 + k) % n1, (i1 - k) % n1)   # K₁ without opening
    seq2 = span(Q2, (i2 + k) % n2, (i2 - k) % n2)    # K₂ without opening
    C = np.vstack([seq1a, bridge(A1, A2, +0.4),
                   seq2, bridge(B2, B1, -0.4)])
    for _ in range(2):
        C = 0.5 * C + 0.25 * (np.roll(C, 1, axis=0)
                              + np.roll(C, -1, axis=0))
    return C


def make_knot(parts, out_root, tag):
    """parts: list of curves → iterated sum."""
    C = parts[0]
    for Pn in parts[1:]:
        C = connected_sum(C, Pn)
    name = f"cs_{tag}"
    entry = {"name": name, "machine": "100_summing_the_knots",
             "n_summands": len(parts),
             "predictions": {
                 "det": "multiplikativ (Zündung)",
                 "genus": "additiv (Zündung)",
                 "theta": "Verhalten unter # (Dror-Frage)"}}
    emb = embedding_report(C)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = None
    for nt in (128, 192, 256, 384):
        try:
            rep = polyline_to_fseries(C, folder / f"{name}.fseries",
                                      knot_name=name,
                                      source="maschinen_ideen/100",
                                      n_terms=nt, nyquist_max_rel=0.02)
            break
        except RuntimeError:
            continue
    if rep is None:
        entry["rejected"] = "Nyquist rot"
        return entry, None
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 100 (connected sums)")
    T = trefoil()
    Tm = trefoil(mirror=True)
    # (1) Opening geometry: bridges connect genuine outer points
    P1 = resample_uniform(T, 2048)
    i1 = outermost_point(P1, [1, 0, 0])
    x_max = P1[:, 0].max()
    ok1 = abs(P1[i1, 0] - x_max) < 1e-12
    print(f"  (1) opening at the outermost point (x={P1[i1, 0]:.4f} = max)  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) GRANNY vs SQUARE: both built (the classical control pair)
    eg, fg = make_knot([T, T], HERE / "demo_smoke", "granny")
    es, fs = make_knot([T, Tm], HERE / "demo_smoke", "square")
    ok2 = fg is not None and fs is not None
    print(f"  (2) Granny (3₁#3₁) + Square (3₁#m3₁): nyq="
          f"{eg.get('fseries', {}).get('nyquist_dev_rel')}/"
          f"{es.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓ control pair in place (det 9/9, sig ±6/0 = firing)' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) LADDER: 3₁#3₁#3₁ (relative of anchor 1: knots-of-love = ^#6)
    e3, f3 = make_knot([T, T, T], HERE / "demo_smoke", "triple")
    ok3 = f3 is not None
    print(f"  (3) triple sum: nyq="
          f"{e3.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓ ladder runs (det bet 27)' if ok3 else '✗ FAIL: ' + str(e3.get('rejected'))}")
    assert ok3
    print("[DEMO] 3/3 green — sum machine ready to fire "
          "(det multiplicative + θ under # = firing bets; anchor 1 wire).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
