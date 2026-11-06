# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave8.py — Wave 8: combined shear+twist matrix.

Wave 2 had a shear series (slope ±0.5, ±1.0, ±1.5) without twist, waves 4-7
helix twist without shear. Wave 8 fills the 2D grid:

  shear ∈ {-0.5, 0.0, +0.5}  ×  twist ∈ {-2.0, -1.0, 0.0, +1.0, +2.0}

= 15 variants kf86..kf100. Order: first shear (z'=z+slope·x),
then helix twist (R_z(twist·i/L) on inner vertex i).

Goal: transitions between the shear basin (wave 2 → mostly K_24 cluster) and
helix-twist classes (waves 4-7 → 12n_839, 12n_837, 11n_133, trefoil).
Expected bifurcations on the diagonals, new classes where both
perturbations collide.

(0,0) = identical reproduction of the baseline knot (kf12) stays in the grid
as control/self-test.

Output: erzeugte_knoten/kf{86..100}_<descr>/<name>.fseries + wave8_index.json.
Sweep via separate sweep_wave8.sh.
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
    """Linear shear z' = z + slope·x for inner vertices (wave-2 convention)."""
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
    """First shear, then twist — the order is not commutative."""
    Mp = perturb_zshear(M, slope)
    Mp = perturb_helix_twist(Mp, twist_total)
    return Mp


def _fmt(x: float) -> str:
    """+0.5 / -1.0 / 0 — compact signed form for folder names."""
    if x == 0:
        return "0"
    return f"{x:+.2f}".rstrip("0").rstrip(".")


SHEARS = [-0.5, 0.0, +0.5]
TWISTS = [-2.0, -1.0, 0.0, +1.0, +2.0]


def build_variants():
    out = []
    kf = 86
    for sh in SHEARS:
        for tw in TWISTS:
            name = f"kf{kf:02d}_sh{_fmt(sh)}_tw{_fmt(tw)}"
            out.append((name, sh, tw))
            kf += 1
    return out


def main():
    if not MODULE_NPY.exists():
        print(f"FEHLT: {MODULE_NPY}", file=sys.stderr)
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
                      source=f"knotfamily_machine.v1 wave8 ({descr})")

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

    idx_path = OUT_BASE / "wave8_index.json"
    idx_path.write_text(json.dumps(index, indent=2))
    print(f"\n→ {idx_path}")
    print(f"→ {len(variants)} fseries-Files in erzeugte_knoten/kf{{86..100}}_*/")


if __name__ == "__main__":
    main()
