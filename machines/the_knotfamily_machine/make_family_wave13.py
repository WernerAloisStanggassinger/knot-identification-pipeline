# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave13.py — Wave 13: pos-shear high resolution × wave-9 twist grid.

Goal: maximum KI-hit yield at minimum failure risk. Wave-10 finding
(Discovery #180 MAIN FINDING #3): pos-shear has a 92% success rate, neg-shear 23%.
Wave 13 exploits this and fills up the pos-shear 2D grid.

Axis design (gap filler for wave 9/10/11):
  shear ∈ {+0.5, +0.75, +1.25, +1.5}  (4 values; sh=+1 already covered in wave 10)
  twist ∈ {-2, -1, +1, +2}             (4 values from the wave-9 coarse grid)
  → 4×4 = 16 variants kf170..kf185

Not overlapping with:
  Wave 9:  sh∈{±1, ±1.5} × tw∈{-2..+2 step 1}  (sh=+1 and sh=+1.5 newly combined here with +0.5/+0.75/+1.25)
  Wave 10: sh∈{±1}      × tw∈{-3..+3 step 0.5} (sh=+1 already exhausted there)
  Wave 11: sh∈{±1.5,±1.25,±1,±0.75,±0.5} × tw=0 (tw=0 deliberately absent here)

Expected findings:
  - sh+0.5/+0.75 × tw{±1,±2}: never probed yet — could fall into 12n_837/839/706
    clusters or open a new topological region
  - sh+1.25/+1.5 × tw{±1,±2}: extends the S_16_d289_D9_chiral region
    (wave-9 Discovery #179: K12a_1019 at sh+1.5)

Convention: perturb_shear_then_twist as in wave 8/9/10/11/12.

Output: erzeugte_knoten/kf{170..185}_<descr>/<name>.fseries + wave13_index.json.
Sweep via separate sweep_wave13.sh.
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


def perturb_zshear(M: np.ndarray, slope: float) -> np.ndarray:
    Mp = M.copy()
    Mp[1:-1, 2] = M[1:-1, 2] + slope * M[1:-1, 0]
    return Mp


def perturb_helix_twist(M: np.ndarray, twist_total: float) -> np.ndarray:
    Mp = M.copy()
    L = len(M) - 1
    for i in range(1, L):
        ang = twist_total * (i / L)
        R = rot_z(ang)
        Mp[i] = R @ M[i]
    return Mp


def perturb_shear_then_twist(M: np.ndarray, slope: float,
                             twist_total: float) -> np.ndarray:
    Mp = perturb_zshear(M, slope)
    Mp = perturb_helix_twist(Mp, twist_total)
    return Mp


def _fmt(x: float) -> str:
    if x == 0:
        return "0"
    return f"{x:+.2f}".rstrip("0").rstrip(".")


SHEARS = [+0.5, +0.75, +1.25, +1.5]
TWISTS = [-2.0, -1.0, +1.0, +2.0]


def build_variants():
    out = []
    kf = 170
    for sh in SHEARS:
        for tw in TWISTS:
            name = f"kf{kf:03d}_sh{_fmt(sh)}_tw{_fmt(tw)}"
            out.append((name, sh, tw))
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
            a, b = polyline_to_fourier(P, n_terms=20, n_samples=2048,
                                       center=True)
        except Exception as e:
            print(f"  {name}: Fourier error: {e}")
            index.append({"name": name, "error": str(e)})
            continue

        out_dir = OUT_BASE / name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{name}.fseries"
        descr = (f"module_shear_then_twist "
                 f"{{'slope': {sh}, 'twist_total': {tw}}}")
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.v1 wave13 ({descr})")

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

    idx_path = OUT_BASE / "wave13_index.json"
    idx_path.write_text(json.dumps(index, indent=2))
    print(f"\n→ {idx_path}")
    print(f"→ {len(variants)} fseries files in erzeugte_knoten/kf{{170..185}}_*/")


if __name__ == "__main__":
    main()
