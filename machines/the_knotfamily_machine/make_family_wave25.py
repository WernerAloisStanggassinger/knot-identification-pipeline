# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave25.py — Wave 25: poly_c07-Foundation-Probe (D7 polygon).

Cycle 3, wave 3.

poly_c07 (c=39, det=7, D7, sig=+14, fib=T) is a singleton — D7 symmetry
so far unrepresented in the corpus. Hypothesis: D7 foundation. If a family
exists, a new sym class is established.

Operator: build_polygon_symmetric(M, n=7, mirror=False) on a 3x3-
perturbed module. snap_closure re-anchors M[-1] = polygon_generator(7).
3x3 mini grid (shear, twist) around (0, 0) = poly_c07 position, centre
left out -> 8 new knots kf284..kf291.

  shear in {-0.5, 0.0, +0.5}  (Step 0.5)
  twist in {-1.0, 0.0, +1.0}  (Step 1.0)

Hypotheses:
  - If the variants fall into a consistent cluster S_(c)_d(det)_D7_chiral,
    the D7 foundation is confirmed.
  - If they scatter, poly_c07 is truly structurally isolated
    (sampling difficulty or real singularity).
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
    """z-shear: z_new = z + slope * x, only on interior."""
    Mp = M.copy()
    Mp[1:-1, 2] = M[1:-1, 2] + slope * M[1:-1, 0]
    return Mp


def perturb_helix_twist(M, twist_total):
    """Helix twist about the z axis, linear over interior vertex index."""
    Mp = M.copy()
    L = len(M) - 1
    for i in range(1, L):
        ang = twist_total * (i / L)
        R = rot_z(ang)
        Mp[i] = R @ M[i]
    return Mp


def perturb_shear_then_twist(M, slope, twist_total):
    return perturb_helix_twist(perturb_zshear(M, slope), twist_total)


SHEARS = [-0.5, 0.0, +0.5]
TWISTS = [-1.0, 0.0, +1.0]
SKIP = (0.0, 0.0)  # = poly_c07 (snap_closure only, no perturb)


def _fmt(x):
    if x == 0:
        return "0"
    return f"{x:+.2f}".rstrip("0").rstrip(".")


def build_variants():
    out = []
    kf = 284
    for sh in SHEARS:
        for tw in TWISTS:
            if (sh, tw) == SKIP:
                continue
            name = f"kf{kf:03d}_polyc07_sh{_fmt(sh)}_tw{_fmt(tw)}"
            out.append((name, sh, tw))
            kf += 1
    return out


def main():
    if not MODULE_NPY.exists():
        print(f"FEHLT: {MODULE_NPY}", file=sys.stderr)
        sys.exit(2)
    M_base = np.load(MODULE_NPY)
    G7 = polygon_generator(7, False)

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    index = []
    variants = build_variants()
    print(f"{'name':<36}  {'sh':>5}  {'tw':>5}  {'pts':>5}  {'bbox':>6}"
          f"  {'odd':>7}  {'even':>7}")

    for name, sh, tw in variants:
        # 1) Module-Perturbation (shear-then-twist)
        Mp = perturb_shear_then_twist(M_base, sh, tw)
        # 2) re-anchor closure for n=7 polygon-symmetric
        Ms = snap_closure(Mp, G7)
        ce = float(np.linalg.norm(G7 @ Ms[0] - Ms[-1]))
        # 3) Build C_7 polygon-symmetric (7 module copies)
        try:
            P = build_polygon_symmetric(
                Ms, n=7, mirror=False,
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
            index.append({"name": name, "error": str(e)})
            continue
        out_dir = OUT_BASE / name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{name}.fseries"
        descr = f"poly_c07_shear_then_twist {{'slope': {sh}, 'twist_total': {tw}}}"
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.v1 wave25 ({descr})")
        amps = np.sqrt((a ** 2 + b ** 2).sum(axis=1))
        odd_sum = float(amps[::2].sum())
        even_sum = float(amps[1::2].sum())
        bbox = float(np.abs(P).max())
        index.append({
            "name": name, "slope": sh, "twist_total": tw, "n": 7,
            "mirror": False, "descr": descr,
            "closure_err_after_snap": ce,
            "polyline_points": len(P), "bbox": bbox,
            "odd_sum": odd_sum, "even_sum": even_sum,
            "amps_top5": [float(x) for x in amps[:5]],
            "out": str(out_path),
        })
        print(f"  {name:<36}  {sh:>+5.2f}  {tw:>+5.2f}  {len(P):>5}"
              f"  {bbox:>6.3f}  {odd_sum:>7.3f}  {even_sum:>7.3f}")

    idx_path = OUT_BASE / "wave25_index.json"
    idx_path.write_text(json.dumps(index, indent=2))
    print(f"\n-> {idx_path}")
    print(f"-> {len(variants)} fseries-Files in erzeugte_knoten/kf{{284..291}}_*/")


if __name__ == "__main__":
    main()
