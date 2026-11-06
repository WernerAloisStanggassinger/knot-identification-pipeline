# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave6.py — Wave 6: helix-twist fine grid.

Continuation of wave 4 (kf48 t=1.0, kf49 t=2.0) and wave 5 (kf58 t=3.0,
kf59 t=4.0). Discovery #162 showed two interesting findings in the helix
twist range:

  - kf58 (t=3.0) → bekannt/K11n_133m  (KnotInfo hit)
  - kf59 (t=4.0) → bekannt/K3_1       (collapse to the trefoil)

Between t=1.0 and t=4.0 the intermediate values are missing. Wave 6 fills the
fine grid in 0.25 steps:

  range [1, 2]:  t ∈ {1.25, 1.5, 1.75}
  range [2, 3]:  t ∈ {2.25, 2.5, 2.75}
  range [3, 4]:  t ∈ {3.25, 3.5, 3.75}

= 9 variants kf65-kf73.

Hypotheses:
  - At [2.5, 3.0] presumably further KnotInfo hits (same mechanism
    as kf58 = K11n_133m).
  - At [3.5, 3.75] presumably transition from the KnotInfo range into the
    trefoil collapse region (kf59 t=4.0 is a trivial collapse).
  - At [1, 2] presumably new cluster members for K12n_837 or
    K_24 (Wave-4-Cluster).

Output: erzeugte_knoten/kf{65..73}_<descr>/<name>.fseries + wave6_index.json.
Sweep via separate sweep_wave6.sh.
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


def perturb_helix_twist(M: np.ndarray, twist_total: float) -> np.ndarray:
    """Identical to wave 4/5 perturb_helix_twist. Rotate interior vertex i
    in xy by angle twist_total*(i/L), L=len(M)-1."""
    Mp = M.copy()
    L = len(M) - 1
    for i in range(1, L):
        ang = twist_total * (i / L)
        R = rot_z(ang)
        Mp[i] = R @ M[i]
    return Mp


VARIANTS = [
    # range [1, 2]
    ("kf65_helix_twist_t1.25", 1.25),
    ("kf66_helix_twist_t1.5",  1.5),
    ("kf67_helix_twist_t1.75", 1.75),
    # range [2, 3]
    ("kf68_helix_twist_t2.25", 2.25),
    ("kf69_helix_twist_t2.5",  2.5),
    ("kf70_helix_twist_t2.75", 2.75),
    # range [3, 4]
    ("kf71_helix_twist_t3.25", 3.25),
    ("kf72_helix_twist_t3.5",  3.5),
    ("kf73_helix_twist_t3.75", 3.75),
]


def main():
    M = np.load(MODULE_NPY)
    R120 = rot_z(2 * np.pi / 3)
    closure_err = float(np.linalg.norm(R120 @ M[0] - M[-1]))
    print(f"Loaded module: {M.shape}, original closure = {closure_err:.6f}")
    print()

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    summary = []
    print(f"{'name':<28}  {'t':>5}  {'pts':>5}  {'bbox':>6}  {'odd':>7}  {'even':>7}")

    for name, t in VARIANTS:
        Mp = perturb_helix_twist(M, twist_total=t)
        ce = float(np.linalg.norm(R120 @ Mp[0] - Mp[-1]))
        P = build_inner(Mp, drop_seam_duplicates=True, closure_tol=ce + 1e-4)

        try:
            a, b = polyline_to_fourier(P, n_terms=20, n_samples=2048, center=True)
        except Exception as e:
            print(f"  {name}: Fourier error: {e}")
            summary.append({"name": name, "error": str(e)})
            continue

        out_dir = OUT_BASE / name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{name}.fseries"
        descr = f"module_perturb_then_inner {{'twist_total': {t}}}"
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.v1 wave6 ({descr})")

        amps = np.sqrt((a ** 2 + b ** 2).sum(axis=1))
        odd_sum = float(amps[::2].sum())
        even_sum = float(amps[1::2].sum())
        bbox = float(np.abs(P).max())

        summary.append({
            "name": name, "twist_total": t, "descr": descr,
            "polyline_points": int(P.shape[0]),
            "bbox": bbox, "odd_sum": odd_sum, "even_sum": even_sum,
            "amps_top5": amps[:5].tolist(),
            "out": str(out_path),
        })
        print(f"{name:<28}  {t:5.2f}  {P.shape[0]:>5}  "
              f"{bbox:6.3f}  {odd_sum:7.4f}  {even_sum:7.4f}")

    with open(OUT_BASE / "wave6_index.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nWave6-Index: {OUT_BASE / 'wave6_index.json'}")
    print(f"{len(summary)} modules generated under {OUT_BASE}")


if __name__ == "__main__":
    main()
