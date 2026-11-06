# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
Wave 2: 12 new variants with more aggressive perturbations.

Classes:
  A) Strong shears (slope ±1.0, ±1.5)
  B) Extreme pinches (factor 0.2, 2.0)
  C) Single-vertex displacement (breaks C₃ completely)
  D) Random gaussian scatter (different seeds)
  E) Combined: shear + xy-rotate

Output: $KNOT_PIPELINE_ROOT/erzeugte_knoten/<name>/<name>.fseries
"""
from __future__ import annotations
import site_config
import numpy as np
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from knot_machine import build_inner, rot_z  # noqa
from fseries_io import polyline_to_fourier, write_fseries  # noqa

OUT_BASE = Path((site_config.GENERATED_KNOTS_DIR))
MODULE_NPY = ROOT / "double3_module_extracted.npy"


def perturb_zshear(M, slope):
    Mp = M.copy()
    Mp[1:-1, 2] = M[1:-1, 2] + slope * M[1:-1, 0]
    return Mp


def perturb_pinch(M, factor):
    Mp = M.copy()
    z = M[1:-1, 2]
    z_abs_norm = np.abs(z) / max(np.abs(M[:, 2]).max(), 1e-9)
    scale = 1.0 + (factor - 1.0) * z_abs_norm
    Mp[1:-1, 0:2] = M[1:-1, 0:2] * scale[:, None]
    return Mp


def perturb_single_vertex(M, vidx, dx, dy, dz):
    """Displace a single interior vertex; breaks C₃ completely if
    the rest of the path (through C₃ chaining) stays symmetric — means
    the final knot has the displacement repeated 3× rotated."""
    Mp = M.copy()
    if 0 < vidx < len(M) - 1:
        Mp[vidx] = Mp[vidx] + np.array([dx, dy, dz])
    return Mp


def perturb_gaussian(M, sigma, seed):
    rng = np.random.default_rng(seed)
    Mp = M.copy()
    noise = rng.normal(scale=sigma, size=(len(M) - 2, 3))
    Mp[1:-1] = M[1:-1] + noise
    return Mp


def perturb_shear_plus_rot(M, slope, angle_deg):
    """First shear, then rotate."""
    Mp = perturb_zshear(M, slope)
    R = rot_z(np.deg2rad(angle_deg))
    Mp[1:-1] = Mp[1:-1] @ R.T
    return Mp


VARIANTS = [
    # class A: strong shears
    ("kf13_shear+1.0",   perturb_zshear,         dict(slope=1.0),                  "z-shear +1.0·x"),
    ("kf14_shear+1.5",   perturb_zshear,         dict(slope=1.5),                  "z-shear +1.5·x"),
    ("kf15_shear-1.5",   perturb_zshear,         dict(slope=-1.5),                 "z-shear -1.5·x"),
    # class B: extreme pinches
    ("kf16_pinch0.1",    perturb_pinch,          dict(factor=0.1),                 "radial pinch 0.1"),
    ("kf17_pinch2.5",    perturb_pinch,          dict(factor=2.5),                 "radial pinch 2.5"),
    # class C: single vertex (breaks C₃)
    ("kf18_v6_dz+0.5",   perturb_single_vertex,  dict(vidx=6, dx=0, dy=0, dz=0.5), "Vertex 6 +z 0.5"),
    ("kf19_v3_dy+0.4",   perturb_single_vertex,  dict(vidx=3, dx=0, dy=0.4, dz=0), "Vertex 3 +y 0.4"),
    # class D: random gaussian
    ("kf20_gauss_s42",   perturb_gaussian,       dict(sigma=0.15, seed=42),        "gauss σ=0.15 seed=42"),
    ("kf21_gauss_s7",    perturb_gaussian,       dict(sigma=0.15, seed=7),         "gauss σ=0.15 seed=7"),
    ("kf22_gauss_s99",   perturb_gaussian,       dict(sigma=0.25, seed=99),        "gauss σ=0.25 seed=99"),
    # class E: combined
    ("kf23_shear_rot",   perturb_shear_plus_rot, dict(slope=0.7, angle_deg=20),   "shear+rotate"),
    ("kf24_shear_rot-",  perturb_shear_plus_rot, dict(slope=-0.7, angle_deg=-20), "shear-rotate-"),
]


def main():
    M = np.load(MODULE_NPY)
    R120 = rot_z(2 * np.pi / 3)
    closure_orig = float(np.linalg.norm(R120 @ M[0] - M[-1]))
    print(f"module closure error (original): {closure_orig:.6f}\n")

    summary = []
    for name, fn, kwargs, descr in VARIANTS:
        Mp = fn(M, **kwargs)
        # Closure must be preserved (boundary unchanged)
        ce = float(np.linalg.norm(R120 @ Mp[0] - Mp[-1]))
        P = build_inner(Mp, drop_seam_duplicates=True, closure_tol=closure_orig + 1e-6)
        a, b = polyline_to_fourier(P, n_terms=20, n_samples=2048, center=True)

        out_dir = OUT_BASE / name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{name}.fseries"
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.v1 wave2 ({descr})")

        amps = np.sqrt((a ** 2 + b ** 2).sum(axis=1))
        bbox = float(np.abs(P).max())
        summary.append({
            "name": name, "descr": descr, "closure_err": ce,
            "bbox": bbox,
            "odd_sum": float(amps[::2].sum()),
            "even_sum": float(amps[1::2].sum()),
            "amps_top5": amps[:5].tolist(),
            "out": str(out_path),
        })
        print(f"  {name:<20} bbox={bbox:.3f} odd={summary[-1]['odd_sum']:.3f} "
              f"even={summary[-1]['even_sum']:.3f} closure={ce:.5f}  ({descr})")

    import json
    with open(OUT_BASE / "wave2_index.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nWave2-Index: {OUT_BASE / 'wave2_index.json'}")


if __name__ == "__main__":
    main()
