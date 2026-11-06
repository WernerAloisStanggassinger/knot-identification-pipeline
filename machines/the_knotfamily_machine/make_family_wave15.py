# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave15.py — Wave 15: det=27 completion via the kf179 neighbourhood.

Per Disc #185 TODO d + corpus inventory 2026-05-16: kf179 (sh+1.25, tw-1) is
the only remaining det=27 singleton — Z/2 sym, vol=15.91, c=15, sig=+6.
Differs from:
  S_18_d27_D3_chiral (vol=18.02, D3, sig=+6, 1 member kf17)
  S_18_d27_D3_chiral_m (vol=18.02, D3, sig=-6, 1 member kf43, mirror-pair)
  S_18_d27_D3_chiral_v14 (vol=14.85, D3, sig=+6, 7 members)

Hypothesis: kf179 could have its own basin, which becomes visible with step 0.125-0.25
around (sh+1.25, tw-1). Wave 15 tests this with 6 knots:

  kf192 (sh+1.25,  tw-1.5)
  kf193 (sh+1.25,  tw-0.75)
  kf194 (sh+1.25,  tw-0.5)
  kf195 (sh+1.125, tw-1)
  kf196 (sh+1.375, tw-1)
  kf197 (sh+1.5,   tw-0.75)

Goal: confirm basin growth (>=2 matches on kf179 theta) or as
verify it as a genuine point island (all 6 fall into other known clusters).

Pos-shear region (wave-13 finding: low OOM risk).
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
    return f"{x:+.3f}".rstrip("0").rstrip(".")


VARIANTS_SPEC = [
    ("kf192", +1.25,  -1.5),
    ("kf193", +1.25,  -0.75),
    ("kf194", +1.25,  -0.5),
    ("kf195", +1.125, -1.0),
    ("kf196", +1.375, -1.0),
    ("kf197", +1.5,   -0.75),
]


def build_variants():
    return [(f"{kf}_sh{_fmt(sh)}_tw{_fmt(tw)}", sh, tw)
            for kf, sh, tw in VARIANTS_SPEC]


def main():
    if not MODULE_NPY.exists():
        print(f"MISSING: {MODULE_NPY}", file=sys.stderr)
        sys.exit(2)
    M = np.load(MODULE_NPY)
    R120 = rot_z(2 * np.pi / 3)

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    index = []
    variants = build_variants()
    print(f"{'name':<32}  {'sh':>6}  {'tw':>6}  {'pts':>5}  {'bbox':>6}"
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
                      source=f"knotfamily_machine.v1 wave15 ({descr})")

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
        print(f"  {name:<32}  {sh:+.3f}  {tw:+.3f}  {len(P):>5}"
              f"  {bbox:>6.3f}  {odd_sum:>7.3f}  {even_sum:>7.3f}")

    idx_path = OUT_BASE / "wave15_index.json"
    idx_path.write_text(json.dumps(index, indent=2))
    print(f"\n→ {idx_path}")
    print(f"→ {len(variants)} fseries files in erzeugte_knoten/kf{{192..197}}_*/")


if __name__ == "__main__":
    main()
