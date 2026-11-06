# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_polygon_wave2.py — extension of the symmetry sweep.

Wave 2 of the polygon family. While the first wave (C_3/5/7/9 + S_4/6/8/10)
runs in the pipeline, we generate in parallel CPU-only new .fseries for
further symmetries:

  Even cyclic C_n (NEW in the corpus — wave 1 had only odd C):
    C_2, C_4, C_6, C_8, C_10
  Further odd C_n:
    C_11, C_13
  Further Schoenflies S_n:
    S_2 (smallest — pure inversion symmetry)
    S_12, S_14

Output writes to erzeugte_knoten/poly_<name>/ and collects in
polygon_family_wave2_index.json (separate from the wave-1 index).

Performance: 10 variants in parallel via ProcessPoolExecutor
(polyline_to_fourier with 2048 samples is CPU-bound — parallelization
mandatory according to memory).
"""
from __future__ import annotations
import site_config
import json
import traceback
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from knot_machine import (
    build_polygon_symmetric, polygon_generator,
)
from fseries_io import polyline_to_fourier, write_fseries
from make_family_polygon import snap_closure  # Reuse

OUT_BASE = Path((site_config.GENERATED_KNOTS_DIR))
MODULE_NPY = ROOT / "double3_module_extracted.npy"

VARIANTS_WAVE2 = [
    # Even cyclic C_n (purely cyclic without sigma_h)
    ("poly_c02",  2, False, "C_2 (2 Module, minimale zyklische)"),
    ("poly_c04",  4, False, "C_4 (4 Module, gerade)"),
    ("poly_c06",  6, False, "C_6 (6 Module, gerade)"),
    ("poly_c08",  8, False, "C_8 (8 Module, gerade)"),
    ("poly_c10", 10, False, "C_10 (10 Module, gerade)"),
    # Higher odd C_n
    ("poly_c11", 11, False, "C_11 (11 unique Module)"),
    ("poly_c13", 13, False, "C_13 (13 unique Module)"),
    # Further Schoenflies S_n
    ("poly_s02",  2, True,  "S_2 (pure Inversion mit z-flip)"),
    ("poly_s12", 12, True,  "S_12 (6 unique + Mirror)"),
    ("poly_s14", 14, True,  "S_14 (7 unique + Mirror)"),
]


def _build_one(args):
    """Worker for ProcessPoolExecutor — generates a polygon variant."""
    name, n, mirror, descr, module_npy_path, out_base = args
    try:
        M_orig = np.load(module_npy_path)
        G = polygon_generator(n, mirror)
        M_new = snap_closure(M_orig, G)
        ce_after = float(np.linalg.norm(G @ M_new[0] - M_new[-1]))
        P = build_polygon_symmetric(
            M_new, n=n, mirror=mirror,
            drop_seam_duplicates=True,
            closure_tol=1e-6,
        )
        a, b = polyline_to_fourier(P, n_terms=20, n_samples=2048, center=True)
        out_dir = Path(out_base) / name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{name}.fseries"
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.polygon_symmetric n={n} "
                             f"mirror={mirror} ({descr}, wave2)")
        amps = np.sqrt((a ** 2 + b ** 2).sum(axis=1))
        bbox = float(np.abs(P).max())
        return {
            "name": name, "n": n, "mirror": mirror, "descr": descr,
            "module_len": int(M_new.shape[0]),
            "polyline_len": int(P.shape[0]),
            "closure_err_after_snap": ce_after,
            "bbox": bbox,
            "odd_sum":  float(amps[::2].sum()),
            "even_sum": float(amps[1::2].sum()),
            "out": str(out_path),
            "_status": "ok",
            "_log": (f"{name:<12}  {n:>2} {str(mirror):>4}  {M_new.shape[0]:>3}  "
                     f"{ce_after:>9.2e}  {bbox:6.3f}  {P.shape[0]:>6}  {descr}"),
        }
    except Exception as e:
        return {
            "name": name, "n": n, "mirror": mirror, "descr": descr,
            "error": str(e),
            "traceback": traceback.format_exc(),
            "_status": "error",
            "_log": f"{name:<12}  FEHLER: {e}",
        }


def main():
    M_orig = np.load(MODULE_NPY)
    print(f"Base module: {MODULE_NPY.name}")
    print(f"  shape: {M_orig.shape}")
    print(f"  parallel workers: {min(len(VARIANTS_WAVE2), 10)}")
    print()

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    print(f"{'name':<12}  {'n':>2} {'mir':>4}  {'L':>3}  {'closure':>9}  "
          f"{'bbox':>6}  {'L_poly':>6}  descr")
    print("-" * 95)

    # worker args: (name, n, mirror, descr, module_path, out_base)
    work_items = [(name, n, mirror, descr, str(MODULE_NPY), str(OUT_BASE))
                  for name, n, mirror, descr in VARIANTS_WAVE2]
    summary = []
    with ProcessPoolExecutor(max_workers=min(len(work_items), 10)) as pool:
        for res in pool.map(_build_one, work_items):
            print(res["_log"])
            res.pop("_log", None)
            summary.append(res)

    idx_path = OUT_BASE / "polygon_family_wave2_index.json"
    with open(idx_path, "w") as f:
        json.dump(summary, f, indent=2)
    n_ok = sum(1 for s in summary if s.get("_status") == "ok")
    print(f"\nFamily index wave 2: {idx_path}")
    print(f"{n_ok}/{len(summary)} knots generated under {OUT_BASE}")


if __name__ == "__main__":
    main()
