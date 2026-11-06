# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 44 — cabling_the_spine (cables in the Bishop tube frame).

(p,q) cable: parameter runs p times around the spine, cross-section angle
ψ = q·t − α·(p·t)/(2π) (Bishop holonomy α compensated, 34 mechanics).
BETS EXACT (06 algebra): Δ_cable = Δ_pattern(t)·Δ_spine(t^w);
(2,1) cable: pattern Δ=1, w=2 ⇒ det = |Δ_spine(1)| = 1 — a sharp bet!
CALIBRATION: circle spine (2,3) ⇒ T(2,3) (against machine 10). Do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "34_folding_panels"))
sys.path.insert(0, str(HERE.parent / "06_eating_own_output"))
sys.path.insert(0, str(HERE.parent))
from generate_band_boundary import bishop_frame  # noqa: E402
from generate_satellite_knots import (trefoil_curve, predict_satellite,  # noqa: E402
                                      TREFOIL_DELTA)
from machine_common import (embedding_report, polyline_to_fseries,  # noqa: E402
                            resample_uniform)


def cable(spine, p, q, r_rel=0.18):
    S = resample_uniform(spine, 1024)
    diam = float(np.linalg.norm(S.max(0) - S.min(0)))
    r = r_rel * diam
    T, N1, N2, alpha = bishop_frame(S)
    n_tot = len(S) * p
    idx = np.arange(n_tot) % len(S)
    wrap = np.arange(n_tot) // len(S)
    t = np.arange(n_tot) / n_tot * 2 * np.pi
    # Frame continuous across wraps: per wrap keep turning by −α; redistribute the rest of the
    # holonomy (α·p mod 2π) as a linear ramp, so that
    # ψ closes in the FIXED start frame mod 2π (build finding: knot spine
    # otherwise did not close — circle has α=0 and noticed nothing)
    # STEPPED per wrap: the stored frame N1[0] attaches to the continuously
    # transported one, lagging by +α (holonomy measurement 34) ⇒ ψ must ADD
    # +α at the seam (build finding: −α doubled the jump);
    # the mod-2π remainder is redistributed as a continuous ramp so that ψ in the
    # fixed start frame closes.
    resid = (-alpha * p) % (2 * np.pi)
    psi = q * t + alpha * wrap + resid * (np.arange(n_tot) / n_tot)
    base = S[idx]
    U = (np.cos(psi)[:, None] * N1[idx] + np.sin(psi)[:, None] * N2[idx])
    return base + r * U


def make_knot(spine_kind, p, q, out_root, r_rel=0.18):
    spine = (trefoil_curve() if spine_kind == "trefoil" else
             np.column_stack([np.cos(np.linspace(0, 2*np.pi, 512, endpoint=False)),
                              np.sin(np.linspace(0, 2*np.pi, 512, endpoint=False)),
                              np.zeros(512)]))
    P = cable(spine, p, q, r_rel)
    gap = float(np.linalg.norm(P[-1] - P[0]))
    seg = float(np.linalg.norm(P[1] - P[0]))
    name = f"cb_{spine_kind}_p{p}_q{q}"
    preds = {"closure_gap_over_seg": round(gap / max(seg, 1e-12), 2)}
    if spine_kind == "circle":
        preds["identity_bet"] = f"T({p},{q}) (Eichung gegen Maschine 10)"
    elif p == 2 and q % 2 == 1:
        sat, det = predict_satellite([1], TREFOIL_DELTA, 2)
        preds["det_bet"] = f"{det} (Δ_Spine(t²), w=2 — 06-Formel exakt)"
    entry = {"name": name, "machine": "44_cabling_the_spine",
             "spine": spine_kind, "p": p, "q": q, "predictions": preds}
    if gap > 20 * seg:
        entry["rejected"] = "Kabel schließt nicht (Holonomie-Buchführung)"
        return entry, None
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate (r_rel vs Spine-Krümmung)"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(P, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/44",
                              n_terms=max(96, 24 * (p + q)),
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 44 (cables)")
    e1, f1 = make_knot("circle", 2, 3, HERE / "demo_smoke")
    ok1 = f1 is not None
    print(f"  (1) circle (2,3) → T(2,3) calibration: gap/seg="
          f"{e1['predictions']['closure_gap_over_seg']} embedded="
          f"{e1.get('embedding',{}).get('embedded')}  "
          f"{'✓' if ok1 else '✗ FAIL: ' + str(e1.get('rejected'))}")
    assert ok1
    sat, det = predict_satellite([1], TREFOIL_DELTA, 2)
    ok2 = det == 1
    print(f"  (2) 06 formula: (2,·) cable of the trefoil, w=2 ⇒ det = {det} "
          f"(companion vanishes)  {'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    e3, f3 = make_knot("trefoil", 2, 1, HERE / "demo_smoke", r_rel=0.10)
    ok3 = f3 is not None
    print(f"  (3) trefoil (2,1) cable: det bet "
          f"{e3['predictions'].get('det_bet')} embedded="
          f"{e3.get('embedding',{}).get('embedded')}  "
          f"{'✓' if ok3 else '✗ FAIL: ' + str(e3.get('rejected'))}")
    assert ok3
    print("[DEMO] 3/3 green — cable machine ready to fire (06×34 triangle closed).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
