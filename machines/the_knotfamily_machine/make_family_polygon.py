# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
make_family_polygon.py — symmetry sweep along (n, mirror) on a
base module.

Instead of perturbing the module (like make_family.py / wave2 / wave3 / shearsweep),
this script runs through different SYMMETRIES of the same module shape:

    C_3, C_5, C_7, C_9          (pure cyclic, n unique modules)
    S_4, S_6, S_8, S_10         (Schoenflies with σ_h, n/2 unique + mirror)

Workflow per variant:
  1. Load base module M (default: double3_module_extracted.npy, C_3-conformal)
  2. Re-anchor closure: M_new with linear snap to M_new[-1] = g·M[0],
     where g = polygon_generator(n, mirror). Inner vertices get a
     linear drift, so that closure constraint is exactly satisfied and the
     module geometry is preserved minimally invasively (cf. close_seam snap).
  3. build_polygon_symmetric → closed 3D polyline
  4. Fourier(20) → .fseries → erzeugte_knoten/<name>/<name>.fseries

NO pipeline run here. Pipeline sweep separately via pipeline_run.sh
or the make_family_*.py analogue.
"""
from __future__ import annotations
import site_config
import json
import numpy as np
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from knot_machine import (  # noqa
    build_polygon_symmetric, polygon_generator, make_module_endpoint,
)
from fseries_io import polyline_to_fourier, write_fseries  # noqa

OUT_BASE = Path((site_config.GENERATED_KNOTS_DIR))
MODULE_NPY = ROOT / "double3_module_extracted.npy"


def snap_closure(M_orig: np.ndarray, generator: np.ndarray,
                 *, keep_first: bool = True) -> np.ndarray:
    """Re-anchor the closure condition M[-1] = g·M[0] for a given module.

    keep_first=True: M_new[0] = M_orig[0] (anchor at position 0).
    Linear drift is distributed over the polyline, so that M_new[-1] = g·M_orig[0].

    Convention parallel to knot_machine.close_seam(mode='snap'): minimally invasive,
    smooths the boundary offset over the whole module stretch.
    """
    M = np.asarray(M_orig, dtype=float).copy()
    L = M.shape[0]
    target_end = generator @ M[0]
    delta = target_end - M[-1]                            # offset that has to be closed
    weights = np.linspace(0.0, 1.0, L).reshape(-1, 1)     # 0 at the start, 1 at the end
    return M + weights * delta


# ──────────────────────────────────────────────────────────────────────
# family table: per entry (name, n, mirror, descr).
# ──────────────────────────────────────────────────────────────────────
VARIANTS = [
    # pure cyclic C_n (n unique modules, odd n)
    ("poly_c03",  3, False, "C_3 (Baseline = double-3-Anker)"),
    ("poly_c05",  5, False, "C_5 (5 unique Module)"),
    ("poly_c07",  7, False, "C_7 (7 unique Module)"),
    ("poly_c09",  9, False, "C_9 (9 unique Module)"),
    # Schoenflies S_n with σ_h (n/2 unique + mirror, even n)
    ("poly_s04",  4, True,  "S_4 (2 unique + Mirror)"),
    ("poly_s06",  6, True,  "S_6 (3 unique + Mirror, analog Outer-Regel)"),
    ("poly_s08",  8, True,  "S_8 (4 unique + Mirror)"),
    ("poly_s10", 10, True,  "S_10 (5 unique + Mirror)"),
]


def main():
    M_orig = np.load(MODULE_NPY)
    print(f"base module: {MODULE_NPY.name}")
    print(f"  shape: {M_orig.shape}")
    print(f"  M[0]:  {M_orig[0]}")
    print(f"  M[-1]: {M_orig[-1]}")
    print(f"  z-range: [{M_orig[:,2].min():.3f}, {M_orig[:,2].max():.3f}]")
    print()

    OUT_BASE.mkdir(parents=True, exist_ok=True)

    summary = []
    print(f"{'name':<12}  {'n':>2} {'mir':>4}  {'L':>3}  {'closure':>9}  "
          f"{'bbox':>6}  {'odd':>7}  {'even':>7}  descr")
    print("-" * 90)
    for name, n, mirror, descr in VARIANTS:
        G = polygon_generator(n, mirror)
        # re-anchor closure
        M_new = snap_closure(M_orig, G)
        ce_before = float(np.linalg.norm(G @ M_orig[0] - M_orig[-1]))
        ce_after  = float(np.linalg.norm(G @ M_new[0]  - M_new[-1]))
        # polyline from n rigid copies
        P = build_polygon_symmetric(
            M_new, n=n, mirror=mirror,
            drop_seam_duplicates=True,
            closure_tol=1e-6,
        )
        # Fourier
        a, b = polyline_to_fourier(P, n_terms=20, n_samples=2048, center=True)
        # output folder + .fseries
        out_dir = OUT_BASE / name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{name}.fseries"
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.polygon_symmetric n={n} "
                             f"mirror={mirror} ({descr})")
        amps = np.sqrt((a ** 2 + b ** 2).sum(axis=1))
        odd_sum  = float(amps[::2].sum())
        even_sum = float(amps[1::2].sum())
        bbox = float(np.abs(P).max())
        summary.append({
            "name": name, "n": n, "mirror": mirror, "descr": descr,
            "module_len": int(M_new.shape[0]),
            "polyline_len": int(P.shape[0]),
            "closure_err_before_snap": ce_before,
            "closure_err_after_snap":  ce_after,
            "bbox": bbox,
            "amps_top5": amps[:5].tolist(),
            "odd_sum": odd_sum, "even_sum": even_sum,
            "out": str(out_path),
        })
        print(f"{name:<12}  {n:>2} {str(mirror):>4}  {M_new.shape[0]:>3}  "
              f"{ce_after:>9.2e}  {bbox:6.3f}  "
              f"{odd_sum:7.4f}  {even_sum:7.4f}  {descr}")

    idx_path = OUT_BASE / "polygon_family_index.json"
    with open(idx_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nFamily-Index: {idx_path}")
    print(f"{len(summary)} knots generated under {OUT_BASE}")


if __name__ == "__main__":
    main()
