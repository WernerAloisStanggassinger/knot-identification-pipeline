# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
Shear-axis fine sampling: 5 intermediate values to localize the cluster jumps.

Observed jumps in wave1+2:
  -1.5 (kf15, Singleton Z/3)  ──??──  -0.5 (kf09, Cluster A D6)
   0.0 (kf12, Cluster A D6)   ──??──  +0.5 (kf08, Cluster B D3)
  +0.5 (kf08, Cluster B D3)   ──??──  +1.0 (kf13, Cluster C D6)
  +1.0 (kf13, Cluster C D6)   ──??──  +1.5 (kf14, Singleton D6 det=361)

New slope values (to localize the bifurcation points):
"""
from __future__ import annotations
import site_config
import numpy as np
import sys
import json
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


VARIANTS = [
    ("kf37_shear-1.2",  -1.2,  "z-shear -1.2·x  (zw. -1.5 Singleton und -0.5 Cluster A)"),
    ("kf38_shear-0.8",  -0.8,  "z-shear -0.8·x  (Bifurkation A↔Singleton suchen)"),
    ("kf39_shear+0.25", +0.25, "z-shear +0.25·x (Bifurkation A↔B)"),
    ("kf40_shear+0.75", +0.75, "z-shear +0.75·x (Bifurkation B↔C)"),
    ("kf41_shear+1.25", +1.25, "z-shear +1.25·x (Bifurkation C↔Singleton)"),
]


def main():
    M = np.load(MODULE_NPY)
    R120 = rot_z(2 * np.pi / 3)
    closure_orig = float(np.linalg.norm(R120 @ M[0] - M[-1]))
    print(f"Module closure error (original): {closure_orig:.6f}\n")

    summary = []
    for name, slope, descr in VARIANTS:
        Mp = perturb_zshear(M, slope)
        ce = float(np.linalg.norm(R120 @ Mp[0] - Mp[-1]))
        P = build_inner(Mp, drop_seam_duplicates=True, closure_tol=closure_orig + 1e-6)
        a, b = polyline_to_fourier(P, n_terms=20, n_samples=2048, center=True)

        out_dir = OUT_BASE / name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{name}.fseries"
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.v1 shearsweep ({descr})")

        amps = np.sqrt((a ** 2 + b ** 2).sum(axis=1))
        bbox = float(np.abs(P).max())
        summary.append({
            "name": name, "slope": slope, "descr": descr, "closure_err": ce,
            "bbox": bbox,
            "odd_sum": float(amps[::2].sum()),
            "even_sum": float(amps[1::2].sum()),
            "amps_top5": amps[:5].tolist(),
            "out": str(out_path),
        })
        print(f"  {name:<20} slope={slope:+.2f} bbox={bbox:.3f} closure={ce:.5f}  ({descr})")

    with open(OUT_BASE / "shearsweep_index.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nShearsweep index: {OUT_BASE / 'shearsweep_index.json'}")


if __name__ == "__main__":
    main()
