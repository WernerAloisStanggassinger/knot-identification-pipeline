# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave23.py — Wave 23: kf54 foundation probe (global z-twist).

Cycle 3, wave 1.

kf54_global_twist_3pi (c=15, det=961=31^2, D3, sig=0, fib=F) is a singleton
post cycle 2. Pattern: prime-square determinant like the nitro family
(S_24_d961_D6_amphi). Hypothesis from project_singletons_analysis_post_zyklus1:
**S_15_d961_D3_chiral** as D3 sister of the nitro family (D6 at the same det).

Mini-grid in global_z_twist around total_twist=3pi (kf54 position), step pi/8,
8 new knots (kf268..kf275), range ±0.5pi (i.e. 2.5pi..3.5pi, kf54 itself
left out).

Operator: polyline_post — build_inner (C_3) first, then global_z_twist on
finished polyline (exactly like wave 5 for kf52/53/54).
"""
from __future__ import annotations
import site_config
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from knot_machine import build_inner  # noqa: E402
from fseries_io import polyline_to_fourier, write_fseries  # noqa: E402

OUT_BASE = Path((site_config.GENERATED_KNOTS_DIR))
MODULE_NPY = ROOT / "double3_module_extracted.npy"


def global_z_twist(P: np.ndarray, total_twist: float) -> np.ndarray:
    """Identical to wave5.global_z_twist — rotate points about the z axis with
    angle proportional to the z height (linear z mapping)."""
    P = P.copy()
    z = P[:, 2]
    z_min, z_max = z.min(), z.max()
    if z_max - z_min < 1e-12:
        return P
    angles = total_twist * (z - z_min) / (z_max - z_min)
    cos_a = np.cos(angles)
    sin_a = np.sin(angles)
    x_new = cos_a * P[:, 0] - sin_a * P[:, 1]
    y_new = sin_a * P[:, 0] + cos_a * P[:, 1]
    P[:, 0] = x_new
    P[:, 1] = y_new
    return P


# total_twist in units of pi, step 1/8 around 3pi, 8 values (3pi itself = kf54 left out)
TWIST_FACTORS = [2.500, 2.625, 2.750, 2.875, 3.125, 3.250, 3.375, 3.500]


def _fmt(x):
    return f"{x:.3f}".rstrip("0").rstrip(".")


def build_variants():
    out = []
    kf = 268
    for tf in TWIST_FACTORS:
        name = f"kf{kf:03d}_gtw_{_fmt(tf)}pi"
        out.append((name, tf))
        kf += 1
    return out


def main():
    if not MODULE_NPY.exists():
        print(f"MISSING: {MODULE_NPY}", file=sys.stderr)
        sys.exit(2)
    M = np.load(MODULE_NPY)

    # Build C_3 polyline once — all variants share this base,
    # the global twist acts on it.
    P0 = build_inner(M, drop_seam_duplicates=True)

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    index = []
    variants = build_variants()
    print(f"{'name':<28}  {'tw/pi':>7}  {'pts':>5}  {'bbox':>6}"
          f"  {'odd':>7}  {'even':>7}")

    for name, tf in variants:
        total_twist = tf * np.pi
        P = global_z_twist(P0, total_twist)
        try:
            a, b = polyline_to_fourier(P, n_terms=20, n_samples=2048, center=True)
        except Exception as e:
            print(f"  {name}: Fourier error: {e}")
            index.append({"name": name, "error": str(e)})
            continue
        out_dir = OUT_BASE / name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{name}.fseries"
        descr = f"polyline_post_global_z_twist {{'total_twist_pi': {tf}}}"
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.v1 wave23 ({descr})")
        amps = np.sqrt((a ** 2 + b ** 2).sum(axis=1))
        odd_sum = float(amps[::2].sum())
        even_sum = float(amps[1::2].sum())
        bbox = float(np.abs(P).max())
        index.append({
            "name": name,
            "total_twist_pi": tf,
            "total_twist": total_twist,
            "descr": descr,
            "polyline_points": len(P), "bbox": bbox,
            "odd_sum": odd_sum, "even_sum": even_sum,
            "amps_top5": [float(x) for x in amps[:5]],
            "out": str(out_path),
        })
        print(f"  {name:<28}  {tf:>7.3f}  {len(P):>5}"
              f"  {bbox:>6.3f}  {odd_sum:>7.3f}  {even_sum:>7.3f}")

    idx_path = OUT_BASE / "wave23_index.json"
    idx_path.write_text(json.dumps(index, indent=2))
    print(f"\n-> {idx_path}")
    print(f"-> {len(variants)} fseries files in erzeugte_knoten/kf{{268..275}}_*/")


if __name__ == "__main__":
    main()
