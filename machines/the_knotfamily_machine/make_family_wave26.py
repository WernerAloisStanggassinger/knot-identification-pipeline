# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave26.py — Wave 26: kf90-Foundation-Probe (shear+twist).

Cycle 4, wave 1.

kf90_sh-0.5_tw+2 (c=30, det=81=3^4, D3, sig=0, fib=T) is a singleton
post cycle 3 with prime-square pattern (3^4 after 3^3=27, 3^5=243 already
as S-cluster). Hypothesis: S_30_d81_D3_chiral as the next 3^k foundation
in the D3 stem.

Operator: perturb_shear_then_twist + build_inner (C_3), Wave-8-Konvention.
3x4 mini grid around (sh=-0.5, tw=+2.0) = kf90 position, centre left out
-> 11 new knots kf292..kf302.

  shear in {-1.0, -0.75, -0.5, -0.25}  (Step 0.25)
  twist in {+1.5, +2.0, +2.5}          (Step 0.5)
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
    """Linear-Shear z' = z + slope·x for inner vertices (wave-2 convention)."""
    Mp = M.copy()
    Mp[1:-1, 2] = M[1:-1, 2] + slope * M[1:-1, 0]
    return Mp


def perturb_helix_twist(M: np.ndarray, twist_total: float) -> np.ndarray:
    """Rotate interior vertex i in xy by angle twist_total*(i/L), L=len(M)-1
    (wave 4-7 convention)."""
    Mp = M.copy()
    L = len(M) - 1
    for i in range(1, L):
        ang = twist_total * (i / L)
        R = rot_z(ang)
        Mp[i] = R @ M[i]
    return Mp


def perturb_shear_then_twist(M: np.ndarray, slope: float,
                             twist_total: float) -> np.ndarray:
    """First shear, then twist — order is not commutative. Wave 8."""
    return perturb_helix_twist(perturb_zshear(M, slope), twist_total)


SHEARS = [-1.0, -0.75, -0.5, -0.25]
TWISTS = [+1.5, +2.0, +2.5]
SKIP = (-0.5, +2.0)  # = kf90-Position


def _fmt(x):
    if x == 0:
        return "0"
    return f"{x:+.3f}".rstrip("0").rstrip(".")


def build_variants():
    out = []
    kf = 292
    for sh in SHEARS:
        for tw in TWISTS:
            if (sh, tw) == SKIP:
                continue
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
    print(f"{'name':<32}  {'sh':>6}  {'tw':>5}  {'pts':>5}  {'bbox':>6}"
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
        descr = f"module_shear_then_twist {{'slope': {sh}, 'twist_total': {tw}}}"
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.v1 wave26 ({descr})")
        amps = np.sqrt((a ** 2 + b ** 2).sum(axis=1))
        odd_sum = float(amps[::2].sum())
        even_sum = float(amps[1::2].sum())
        bbox = float(np.abs(P).max())
        index.append({
            "name": name, "slope": sh, "twist_total": tw, "descr": descr,
            "closure_err_module": ce,
            "polyline_points": len(P), "bbox": bbox,
            "odd_sum": odd_sum, "even_sum": even_sum,
            "amps_top5": [float(x) for x in amps[:5]],
            "out": str(out_path),
        })
        print(f"  {name:<32}  {sh:>+6.3f}  {tw:>+5.2f}  {len(P):>5}"
              f"  {bbox:>6.3f}  {odd_sum:>7.3f}  {even_sum:>7.3f}")

    idx_path = OUT_BASE / "wave26_index.json"
    idx_path.write_text(json.dumps(index, indent=2))
    print(f"\n-> {idx_path}")
    print(f"-> {len(variants)} fseries files in erzeugte_knoten/kf{{292..302}}_*/")


if __name__ == "__main__":
    main()
