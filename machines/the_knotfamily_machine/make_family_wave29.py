# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave29.py — Wave 29: poly_s08-Foundation-Probe (S_8 polygon).

Cycle 5, wave 1.

poly_s08 (c=24, det=25721=17^2*89, sym=D8, sig=0, amphi) is a singleton —
D8 symmetry so far unrepresented in the corpus between D7 (S_21_d169, S_39_d7)
and D9 (S_16_d289). Hypothesis: S_24_d25721_D8_amphi as D8 foundation.

Operator: build_polygon_symmetric(M, n=8, mirror=True) on the perturbed
module. snap_closure re-anchors M[-1] = polygon_generator(8, mirror=True).
5x4 mini grid (shear, twist) around (0, 0) = poly_s08 position, centre
left out -> 19 new knots kf321..kf339.

  shear in {-1.0, -0.5, 0.0, +0.5, +1.0}  (Step 0.5)
  twist in {-1.0, 0.0, +1.0, +2.0}        (Step 1.0)

Hypotheses:
  - If the variants fall consistently into S_24_d25721_D8_amphi → D8 foundation
    confirmed (3rd D family in a row: D7, D8, D9).
  - If they scatter → poly_s08 structurally isolated like poly_c07 before W25.
"""
from __future__ import annotations
import site_config
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from knot_machine import (  # noqa: E402
    build_polygon_symmetric, polygon_generator, rot_z,
)
from fseries_io import polyline_to_fourier, write_fseries  # noqa: E402
from make_family_polygon import snap_closure  # noqa: E402

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


SHEARS = [-1.0, -0.5, 0.0, +0.5, +1.0]
TWISTS = [-1.0, 0.0, +1.0, +2.0]
SKIP = (0.0, 0.0)  # = poly_s08


def _fmt(x):
    if x == 0:
        return "0"
    return f"{x:+.2f}".rstrip("0").rstrip(".")


def build_variants():
    out = []
    kf = 321
    for sh in SHEARS:
        for tw in TWISTS:
            if (sh, tw) == SKIP:
                continue
            name = f"kf{kf:03d}_polys08_sh{_fmt(sh)}_tw{_fmt(tw)}"
            out.append((name, sh, tw))
            kf += 1
    return out


def main():
    if not MODULE_NPY.exists():
        print(f"MISSING: {MODULE_NPY}", file=sys.stderr)
        sys.exit(2)
    M_base = np.load(MODULE_NPY)
    G8m = polygon_generator(8, True)

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    index = []
    variants = build_variants()
    print(f"{'name':<40}  {'sh':>5}  {'tw':>5}  {'pts':>5}")

    for name, sh, tw in variants:
        Mp = perturb_shear_then_twist(M_base, sh, tw)
        Ms = snap_closure(Mp, G8m)
        ce = float(np.linalg.norm(G8m @ Ms[0] - Ms[-1]))
        try:
            P = build_polygon_symmetric(
                Ms, n=8, mirror=True,
                drop_seam_duplicates=True,
                closure_tol=max(ce + 1e-4, 1e-4),
            )
        except Exception as e:
            print(f"  {name}: build error: {e}")
            index.append({"name": name, "error": str(e)})
            continue
        try:
            a, b = polyline_to_fourier(P, n_terms=20, n_samples=2048, center=True)
        except Exception as e:
            print(f"  {name}: Fourier error: {e}")
            index.append({"name": name, "error": f"fourier: {e}"})
            continue
        dest_dir = OUT_BASE / name
        dest_dir.mkdir(parents=True, exist_ok=True)
        fpath = dest_dir / f"{name}.fseries"
        write_fseries(fpath, a, b, knot_name=name, source="make_family_wave29.py")
        print(f"  {name:<40}  {sh:>+5.2f}  {tw:>+5.2f}  {len(P):>5d}  -> {fpath.name}")
        index.append({"name": name, "shear": sh, "twist": tw,
                      "n_points": len(P), "closure_residual": ce,
                      "fseries": str(fpath)})

    idx_path = OUT_BASE / "wave29_index.json"
    json.dump(index, open(idx_path, "w"), indent=2)
    print(f"\nWrote {idx_path} ({len(index)} variants)")


if __name__ == "__main__":
    main()
