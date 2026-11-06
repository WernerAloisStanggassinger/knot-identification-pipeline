# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave21.py — Wave 21: sh+1.25 column gap fill.

The sh+1.25 column is already heavily covered (wave 13/15/16):
  tw = -2     → kf178 → 11n_133
  tw = -1.5   → kf192 → 11n_133
  tw = -1     → kf179 → genuine singleton
  tw = -0.75  → kf193 → 11n_133
  tw = -0.5   → kf194 → 11n_133
  tw =  0     → kf155 → 12n_837
  tw = +1     → kf180 → S_16_d289_D9_chiral
  tw = +1.5   → kf190 → 12n_837
  tw = +2     → kf181 → S_18_d243_D3_chiral

Wave 21 fills the 8 large gaps between these support points:
  tw = -1.75, -1.25, -0.25, +0.25, +0.5, +0.75, +1.25, +1.75

8 knots kf251..kf258, all at sh=+1.25.

Hypotheses:
  - The negative tw gaps (-1.75/-1.25/-0.25) should all fall into
    11n_133 (stable wide region) — if not, new basins.
  - The positive tw gaps may show basin boundaries between
    12n_837 (tw=0,1.5) and S_16_d289 (tw=1) and S_18_d243 (tw=2).
    stage-6 failures at bifurcation points to be expected.
  - kf179 (singleton at tw=-1) is an island — wave 21 checks whether the
    neighbours tw=-0.75/-1.25 really fall into 11n_133 or are new
    island members are of kf179.
"""
from __future__ import annotations
import site_config
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from knot_machine import build_inner, rot_z  # noqa: E402
from fseries_io import polyline_to_fourier, write_fseries  # noqa: E402

OUT_BASE = Path((site_config.GENERATED_KNOTS_DIR))
MODULE_NPY = ROOT / "double3_module_extracted.npy"


def perturb_zshear(M, slope):
    Mp = M.copy()
    Mp[1:-1, 2] = M[1:-1, 2] + slope * M[1:-1, 0]
    return Mp


def perturb_helix_twist(M, twist_total):
    Mp = M.copy()
    L = len(M) - 1
    for i in range(1, L):
        ang = twist_total * (i / L)
        R = rot_z(ang)
        Mp[i] = R @ M[i]
    return Mp


def perturb_shear_then_twist(M, slope, twist_total):
    return perturb_helix_twist(perturb_zshear(M, slope), twist_total)


def _fmt(x):
    if x == 0:
        return "0"
    return f"{x:+.3f}".rstrip("0").rstrip(".")


SH_FIXED = +1.25
TWISTS = [-1.75, -1.25, -0.25, +0.25, +0.5, +0.75, +1.25, +1.75]


def build_variants():
    out = []
    kf = 251
    for tw in TWISTS:
        name = f"kf{kf:03d}_sh{_fmt(SH_FIXED)}_tw{_fmt(tw)}"
        out.append((name, SH_FIXED, tw))
        kf += 1
    return out


def main():
    if not MODULE_NPY.exists():
        print(f"MISSING: {MODULE_NPY}", file=sys.stderr)
        sys.exit(2)
    M = np.load(MODULE_NPY)
    R120 = rot_z(2 * np.pi / 3)

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    index = []
    variants = build_variants()
    print(f"{'name':<32}  {'sh':>5}  {'tw':>5}  {'pts':>5}  {'bbox':>6}"
          f"  {'odd':>7}  {'even':>7}")

    for name, sh, tw in variants:
        Mp = perturb_shear_then_twist(M, sh, tw)
        ce = float(np.linalg.norm(R120 @ Mp[0] - Mp[-1]))
        P = build_inner(Mp, drop_seam_duplicates=True, closure_tol=ce + 1e-4)
        try:
            a, b = polyline_to_fourier(P, n_terms=20, n_samples=2048, center=True)
        except Exception as e:
            print(f"  {name}: Fourier error: {e}")
            index.append({"name": name, "error": str(e)})
            continue
        out_dir = OUT_BASE / name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{name}.fseries"
        descr = f"module_shear_then_twist {{'slope': {sh}, 'twist_total': {tw}}}"
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.v1 wave21 ({descr})")
        amps = np.sqrt((a ** 2 + b ** 2).sum(axis=1))
        odd_sum = float(amps[::2].sum())
        even_sum = float(amps[1::2].sum())
        bbox = float(np.abs(P).max())
        index.append({
            "name": name, "slope": sh, "twist_total": tw, "descr": descr,
            "polyline_points": len(P), "bbox": bbox,
            "odd_sum": odd_sum, "even_sum": even_sum,
            "amps_top5": [float(x) for x in amps[:5]],
            "out": str(out_path),
        })
        print(f"  {name:<32}  {sh:+.2f}  {tw:+.2f}  {len(P):>5}"
              f"  {bbox:>6.3f}  {odd_sum:>7.3f}  {even_sum:>7.3f}")

    idx_path = OUT_BASE / "wave21_index.json"
    idx_path.write_text(json.dumps(index, indent=2))
    print(f"\n→ {idx_path}")
    print(f"→ {len(variants)} fseries files in erzeugte_knoten/kf{{251..258}}_*/")


if __name__ == "__main__":
    main()
