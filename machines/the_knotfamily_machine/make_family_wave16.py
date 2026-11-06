# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave16.py — Wave 16: double cluster-foundation sweep
(Disc #193 TODO b) around kf191 + kf196 with step 0.0625.

Goal: for both wave-14/15 singletons with conspicuous algebra profile
find 1-3 theta-bit-identical members, to confirm the S-cluster
foundation.

Foundation A — kf191 (sh+1.375, tw+1.75, c=24/det=961/sig=0/sym=D6/amphi):
  Hypothesis: cluster `S_24_d961_D6_amphi` parallel to S_24_d81_D6_amphi
  with a different determinant (31² instead of 3⁴).
  Sweep: 3x3-Grid sh∈{1.3125, 1.375, 1.4375} x tw∈{1.6875, 1.75, 1.8125}
  centre (1.375, 1.75) = kf191 (skip) -> 8 new knots kf198..kf205.

Foundation B — kf196 (sh+1.375, tw-1, c=18/det=1/sig=0/sym=D3):
  Hypothesis: cluster `S_18_d1_D3_chiral` — det=1 is very rare,
  a Stanggassinger construction with trivial determinant.
  Sweep: 3x3-Grid sh∈{1.3125, 1.375, 1.4375} x tw∈{-1.0625, -1, -0.9375}
  centre (1.375, -1) = kf196 (skip) -> 8 new knots kf206..kf213.

Total: 16 knots kf198..kf213. Sweep -P 6, ETA ~50-60 min.
Convention: perturb_shear_then_twist as in wave 8..15.
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
    return f"{x:+.4f}".rstrip("0").rstrip(".")


VARIANTS_SPEC = [
    ("kf198", +1.3125, +1.6875),
    ("kf199", +1.3125, +1.75),
    ("kf200", +1.3125, +1.8125),
    ("kf201", +1.375,  +1.6875),
    ("kf202", +1.375,  +1.8125),
    ("kf203", +1.4375, +1.6875),
    ("kf204", +1.4375, +1.75),
    ("kf205", +1.4375, +1.8125),
    ("kf206", +1.3125, -1.0625),
    ("kf207", +1.3125, -1.0),
    ("kf208", +1.3125, -0.9375),
    ("kf209", +1.375,  -1.0625),
    ("kf210", +1.375,  -0.9375),
    ("kf211", +1.4375, -1.0625),
    ("kf212", +1.4375, -1.0),
    ("kf213", +1.4375, -0.9375),
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
    print(f"{'name':<36}  {'sh':>7}  {'tw':>7}  {'pts':>5}  {'bbox':>6}"
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
                      source=f"knotfamily_machine.v1 wave16 ({descr})")

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
        print(f"  {name:<36}  {sh:+.4f}  {tw:+.4f}  {len(P):>5}"
              f"  {bbox:>6.3f}  {odd_sum:>7.3f}  {even_sum:>7.3f}")

    idx_path = OUT_BASE / "wave16_index.json"
    idx_path.write_text(json.dumps(index, indent=2))
    print(f"\n→ {idx_path}")
    print(f"→ {len(variants)} fseries files in erzeugte_knoten/kf{{198..213}}_*/")


if __name__ == "__main__":
    main()
