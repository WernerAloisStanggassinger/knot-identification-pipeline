# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 34 — folding_panels, stage 1: BAND EDGE (Möbius path).

Closed spine + band with m half twists (m ODD) → the EDGE of the
band is ONE closed curve = the knot. Bishop frame
(rotation-minimising) + holonomy compensation (crack 2: twist = mπ − α,
α = measured frame holonomy) — the half-twist count is EXACT.

BETS (textbook-exact): circle spine ⇒ edge = T(2,m) torus knot
(m=1 unknot, m=3 trefoil det 3, m=5 det 5). Knot spine ⇒ edge =
2-cable relative (identity = ignition measurement; 06 satellite algebra
as frame). Fold-angle word = twist DISTRIBUTION (workshop link);
real planar panels = expansion stage. Doctrine: do NOT ignite before Monday.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "06_eating_own_output"))
sys.path.insert(0, str(HERE.parent))
from generate_satellite_knots import trefoil_curve  # noqa: E402
from machine_common import (embedding_report, polyline_to_fseries,  # noqa: E402
                            resample_uniform)


def bishop_frame(S):
    """Rotation-minimising frame along a closed polyline.

    Returns (T, N1, N2, holonomie_alpha): N1(end→start) rotation angle α.
    """
    n = len(S)
    T = np.roll(S, -1, axis=0) - S
    T = T / np.linalg.norm(T, axis=1)[:, None]
    N1 = np.zeros_like(T)
    v = np.cross(T[0], [0, 0, 1.0])
    if np.linalg.norm(v) < 1e-9:
        v = np.cross(T[0], [0, 1.0, 0])
    N1[0] = v / np.linalg.norm(v)
    for i in range(1, n):
        # parallel transport: rotate N1 along about the axis T[i-1]×T[i]
        n1 = N1[i - 1]
        n1 = n1 - (n1 @ T[i]) * T[i]
        N1[i] = n1 / max(np.linalg.norm(n1), 1e-12)
    # holonomy: transport from N1[-1] after step back to T[0]
    n1 = N1[-1] - (N1[-1] @ T[0]) * T[0]
    n1 /= max(np.linalg.norm(n1), 1e-12)
    N2_0 = np.cross(T[0], N1[0])
    alpha = float(np.arctan2(n1 @ N2_0, n1 @ N1[0]))
    N2 = np.cross(T, N1)
    return T, N1, N2, alpha


def band_boundary(spine, m_half, width_rel=0.08):
    """Edge of the band with EXACTLY m half twists (holonomy compensated)."""
    if m_half % 2 == 0:
        raise ValueError("m even ⇒ edge = 2 components (crack 3)")
    S = resample_uniform(spine, 1024)
    diam = float(np.linalg.norm(S.max(0) - S.min(0)))
    w = width_rel * diam
    T, N1, N2, alpha = bishop_frame(S)
    n = len(S)
    t = np.arange(n) / n
    phi = (m_half * np.pi - alpha) * t
    U = np.cos(phi)[:, None] * N1 + np.sin(phi)[:, None] * N2
    edge_plus = S + (w / 2) * U
    edge_minus = S - (w / 2) * U
    # m odd: U(end)≈−U(start) ⇒ edge+ continues into edge−
    P = np.vstack([edge_plus, edge_minus])
    gap = float(np.linalg.norm(edge_plus[-1] - edge_minus[0]))
    return P, gap, w


def make_knot(spine_kind, m_half, out_root):
    spine = (trefoil_curve() if spine_kind == "trefoil" else
             np.column_stack([np.cos(np.linspace(0, 2 * np.pi, 512,
                                                 endpoint=False)),
                              np.sin(np.linspace(0, 2 * np.pi, 512,
                                                 endpoint=False)),
                              np.zeros(512)]))
    P, gap, w = band_boundary(spine, m_half)
    name = f"bb_{spine_kind}_m{m_half}"
    bets = {"circle": {1: "Unknot", 3: "Trefoil T(2,3), det 3",
                       5: "T(2,5), det 5", 7: "T(2,7), det 7"}}
    entry = {"name": name, "machine": "34_folding_panels",
             "spine": spine_kind, "half_twists": m_half,
             "edge_gap_over_width": round(gap / w, 4),
             "predictions": {
                 "identity_bet": bets.get(spine_kind, {}).get(
                     m_half, "2-Kabel-Verwandter des Spines (Zünd-Messung)"),
                 "basis": "Rand des m-halbgetwisteten Bandes, Holonomie "
                          "kompensiert (Lehrbuch für Kreis-Spine)",
             }}
    if gap / w > 0.2:
        entry["rejected"] = "Rand schließt nicht (Holonomie-Bug?)"
        return entry, None
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate (width_rel?)"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(P, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/34_folding_panels",
                              n_terms=max(96, 24 * m_half + 48),
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 34 (band edge)")
    # 1) m even refused (crack 3)
    try:
        band_boundary(trefoil_curve(), 2)
        ok1 = False
    except ValueError:
        ok1 = True
    print(f"  (1) m=2 refused (link)  {'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # 2) circle m=1: edge closes (holonomy compensation works)
    e2, f2 = make_knot("circle", 1, HERE / "demo_smoke")
    ok2 = f2 is not None and e2["edge_gap_over_width"] < 0.2
    print(f"  (2) circle m=1 (unknot bet): gap/width="
          f"{e2['edge_gap_over_width']} embedded="
          f"{e2['embedding']['embedded']}  "
          f"{'✓' if ok2 else '✗ FAIL: ' + str(e2.get('rejected'))}")
    assert ok2
    # 3) circle m=3/5: T(2,3)/T(2,5) bets built
    built = []
    for m in (3, 5):
        e, f = make_knot("circle", m, HERE / "demo_smoke")
        if f is not None:
            built.append((m, e["predictions"]["identity_bet"]))
    ok3 = len(built) == 2
    print(f"  (3) torus ladder bets: {built}  {'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    # 4) trefoil spine m=1: 2-cable relative, gates
    e4, f4 = make_knot("trefoil", 1, HERE / "demo_smoke")
    ok4 = f4 is not None
    print(f"  (4) trefoil spine m=1: embedded={e4['embedding']['embedded']} "
          f"Nyquist={e4.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok4 else '✗ FAIL: ' + str(e4.get('rejected'))}")
    assert ok4
    print("[DEMO] 4/4 green — band edge machine ready to fire "
          "(panel stage = extension).")


def main():
    ap = argparse.ArgumentParser(description="Machine 34 — band edge")
    sub = ap.add_subparsers(dest="cmd", required=True)
    pm = sub.add_parser("make")
    pm.add_argument("spine", choices=["circle", "trefoil"])
    pm.add_argument("--m", type=int, required=True)
    sub.add_parser("demo")
    args = ap.parse_args()
    if args.cmd == "demo":
        cmd_demo(args)
    else:
        man_path = HERE / "MANIFEST.json"
        man = json.loads(man_path.read_text()) if man_path.exists() \
            else {"knots": []}
        entry, folder = make_knot(args.spine, args.m, HERE / "out_fseries")
        man["knots"] = [k for k in man["knots"]
                        if k["name"] != entry["name"]] + [entry]
        man_path.write_text(json.dumps(man, indent=1))
        print(f"[MAKE] {entry['name']} → "
              f"{folder if folder else entry.get('rejected')}")


if __name__ == "__main__":
    main()
