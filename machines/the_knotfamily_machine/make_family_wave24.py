# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave24.py — Wave 24: kf26-Foundation-Probe (outer_spiral).

Cycle 3, wave 2.

kf26_outer_spiral (c=18, det=8281=7^2*13^2, D6, sig=0, fib=T) is a singleton
with a **double prime square** pattern (very rare). Hypothesis:
S_18_d8281_D6_amphi-Foundation.

Operator: make_outer_module_spiral(scale, twist) -> build_outer (S_6, 6
module copies). Mini grid 3x3 around (scale=0.5, twist=2.0) = kf26 position,
(0.5, 2.0) left out -> 8 new knots kf276..kf283.

  scale in {0.45, 0.50, 0.55}  (Step 0.05)
  twist in {1.75, 2.00, 2.25}  (Step 0.25)
"""
from __future__ import annotations
import site_config
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from knot_machine import build_outer, s6_generator  # noqa: E402
from fseries_io import polyline_to_fourier, write_fseries  # noqa: E402

OUT_BASE = Path((site_config.GENERATED_KNOTS_DIR))


def make_outer_module_spiral(scale: float, twist: float) -> np.ndarray:
    """Identical to wave3.make_outer_module_spiral — 6-vertex S_6 module with
    spiral path from M[0] to g_1(M[0])."""
    G = s6_generator()
    M0 = np.array([scale, 0.0, 0.4])
    M5 = G @ M0
    M = np.zeros((6, 3))
    M[0] = M0
    M[5] = M5
    for i in range(1, 5):
        t = i / 5
        base = (1 - t) * M0 + t * M5
        r = scale * (0.5 + 0.3 * np.sin(t * np.pi))
        ang = t * twist * 2 * np.pi
        M[i] = np.array([r * np.cos(ang), r * np.sin(ang), base[2]])
    return M


# 3x3 mini grid around (0.5, 2.0) = kf26, centre left out
SCALES = [0.45, 0.50, 0.55]
TWISTS = [1.75, 2.00, 2.25]
SKIP = (0.50, 2.00)  # kf26 position


def _fmt(x):
    return f"{x:.2f}".rstrip("0").rstrip(".")


def build_variants():
    out = []
    kf = 276
    for s in SCALES:
        for tw in TWISTS:
            if (s, tw) == SKIP:
                continue
            name = f"kf{kf:03d}_outersp_s{_fmt(s)}_t{_fmt(tw)}"
            out.append((name, s, tw))
            kf += 1
    return out


def main():
    OUT_BASE.mkdir(parents=True, exist_ok=True)
    index = []
    variants = build_variants()
    print(f"{'name':<32}  {'scale':>6}  {'twist':>6}  {'pts':>5}  {'bbox':>6}"
          f"  {'odd':>7}  {'even':>7}")

    G = s6_generator()
    for name, s, tw in variants:
        M = make_outer_module_spiral(s, tw)
        ce = float(np.linalg.norm(G @ M[0] - M[-1]))
        # build_outer accepts closure_tol; with pure make_outer_module_spiral
        # the closure is exact (M[5] = G @ M[0] by construction), tol generous.
        P = build_outer(M, drop_seam_duplicates=True, closure_tol=ce + 1e-4)
        try:
            a, b = polyline_to_fourier(P, n_terms=20, n_samples=2048, center=True)
        except Exception as e:
            print(f"  {name}: Fourier error: {e}")
            index.append({"name": name, "error": str(e)})
            continue
        out_dir = OUT_BASE / name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{name}.fseries"
        descr = f"outer_module_spiral_then_build_outer {{'scale': {s}, 'twist': {tw}}}"
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.v1 wave24 ({descr})")
        amps = np.sqrt((a ** 2 + b ** 2).sum(axis=1))
        odd_sum = float(amps[::2].sum())
        even_sum = float(amps[1::2].sum())
        bbox = float(np.abs(P).max())
        index.append({
            "name": name, "scale": s, "twist": tw, "descr": descr,
            "closure_err_module": ce,
            "polyline_points": len(P), "bbox": bbox,
            "odd_sum": odd_sum, "even_sum": even_sum,
            "amps_top5": [float(x) for x in amps[:5]],
            "out": str(out_path),
        })
        print(f"  {name:<32}  {s:>6.2f}  {tw:>6.2f}  {len(P):>5}"
              f"  {bbox:>6.3f}  {odd_sum:>7.3f}  {even_sum:>7.3f}")

    idx_path = OUT_BASE / "wave24_index.json"
    idx_path.write_text(json.dumps(index, indent=2))
    print(f"\n-> {idx_path}")
    print(f"-> {len(variants)} fseries files in erzeugte_knoten/kf{{276..283}}_*/")


if __name__ == "__main__":
    main()
