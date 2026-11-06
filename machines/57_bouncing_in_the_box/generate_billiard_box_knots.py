# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 57 — bouncing_in_the_box (3D billiard in a cuboid).

Billiard path in the unit cube: direction (p,q,r) coprime, position
per axis = triangle wave tri(p·t+φ) — the unfolding makes Lissajous
(01) the SMOOTH shadow of this angular machine (Jones–Przytycki
rail as literature anchor, lit check before outside communication).
EXACT BET: path closes after period 2 with EXACTLY 2p/2q/2r wall
touches per axis (countable). REVERSAL CLAUSE (vertex trap, 4th
sighting): corner/edge hit = two axes at the same time at the
wall → phase gate checks minimum separation. Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from math import gcd
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402

N_SAMPLES = 8192
EDGE_GAP_MIN = 1e-3     # minimum time gap between two wall touches


def tri(u):
    """Triangle wave, period 2, values [0,1]: 0→1→0 (wall at 0 and 1)."""
    v = np.mod(u, 2.0)
    return np.where(v <= 1.0, v, 2.0 - v)


def billiard(freqs, phases, n=N_SAMPLES):
    t = np.linspace(0, 2.0, n, endpoint=False)
    return np.column_stack([tri(f * t + ph)
                            for f, ph in zip(freqs, phases)])


def wall_times(f, ph, eps=1e-12):
    """Times of the wall touches of tri(f·t+φ) in [0,2):
    f·t+φ ≡ 0 or 1 (mod 2) → t = (k − φ)/f."""
    ts = []
    for k in range(int(np.ceil(ph - eps)), int(np.floor(2 * f + ph)) + 1):
        tt = (k - ph) / f
        if -eps <= tt < 2.0 - eps:
            ts.append(tt % 2.0)
    return np.array(sorted(ts))


def edge_gap(freqs, phases):
    """Smallest time gap between wall touches of DIFFERENT axes
    (0 = edge/corner hit = vertex trap)."""
    W = [wall_times(f, p) for f, p in zip(freqs, phases)]
    g = np.inf
    for i in range(3):
        for j in range(i + 1, 3):
            d = np.abs(W[i][:, None] - W[j][None, :])
            d = np.minimum(d, 2.0 - d)
            g = min(g, float(d.min()))
    return g


def make_knot(freqs, phases, out_root):
    p, q, r = freqs
    assert gcd(gcd(p, q), r) == 1, "Richtungs-Tripel muss teilerfremd sein"
    gap = edge_gap(freqs, phases)
    name = f"bb_{p}_{q}_{r}"
    entry = {"name": name, "machine": "57_bouncing_in_the_box",
             "freqs": list(freqs), "phases": [round(x, 4) for x in phases],
             "edge_gap": round(gap, 6),
             "wall_hits_predicted": [2 * p, 2 * q, 2 * r],
             "predictions": {
                 "lit_anker": "Jones–Przytycki Billard-Knoten "
                              "(lit check before outside)",
                 "identity": "Zünd-Messung"}}
    if gap < EDGE_GAP_MIN:
        entry["rejected"] = f"Vertex-Gate (edge_gap={gap:.2e})"
        return entry, None
    P = billiard(freqs, phases)
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(P, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/57",
                              n_terms=max(96, 8 * (p + q + r)),
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 57 (billiard in the box)")
    freqs, phases = (3, 5, 7), (0.171, 0.523, 0.317)
    # (1) Closure exact: P(t) has period 2
    t0 = np.array([0.0, 0.7, 1.3])
    d = np.abs(np.column_stack([tri(f * (t0 + 2.0) + ph) for f, ph
                                in zip(freqs, phases)])
               - np.column_stack([tri(f * t0 + ph) for f, ph
                                  in zip(freqs, phases)])).max()
    ok1 = d < 1e-12
    print(f"  (1) period-2 closure: max|Δ|={d:.2e}  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) Wall touches = 2p/2q/2r counted exactly
    counts = [len(wall_times(f, p)) for f, p in zip(freqs, phases)]
    ok2 = counts == [2 * f for f in freqs]
    print(f"  (2) wall touches: {counts} = {[2*f for f in freqs]}  "
          f"{'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) Vertex gate engages on edge hit: φ=0 on two axes ⇒
    #     at t=0 BOTH touch the wall (build lesson: same phase ≠
    #     hit for coprime f — (k−φ)/p=(m−φ)/q is then unsolvable;
    #     the real hit sits at the common zero point)
    e3, f3 = make_knot((3, 5, 7), (0.0, 0.0, 0.317), HERE / "demo_smoke")
    ok3 = f3 is None and "Vertex-Gate" in str(e3.get("rejected", ""))
    print(f"  (3) vertex trap (φx=φy=0): {e3.get('rejected')}  "
          f"{'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    # (4) (3,5,7) billiard knot built + gates green
    e4, f4 = make_knot(freqs, phases, HERE / "demo_smoke")
    ok4 = f4 is not None
    print(f"  (4) bb_3_5_7: edge_gap={e4['edge_gap']}, "
          f"nyq={e4.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok4 else '✗ FAIL: ' + str(e4.get('rejected'))}")
    assert ok4
    print("[DEMO] 4/4 green — billiard machine ready to fire "
          "(Lissajous = smooth shadow, angular world = new territory).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
