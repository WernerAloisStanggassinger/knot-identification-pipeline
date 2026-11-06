# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave4.py — Wave 4: so far unoccupied corners of phase space.

Continuation of waves 1-3 + shearsweep (Discoveries #124-126, #130). Waves 1-3
generated kf01-kf41; wave 4 occupies kf42-kf51 (10 variants) along
axes not explored so far:

  A) 2D phase diagram shear × pinch (6 variants):
       slope ∈ {-0.3, +0.3, +0.6} × pinch ∈ {0.4, 1.6}
     Extends the shearsweep (#130) by a second pinch dimension —
     Saturn/double-3-like D6-thin region vs Z/3 is as
     known as the shear axis; what happens when one pinches at the same time?

  B) Helix twist (2 variants):
       arclength-dependent xy rotation of the interior vertices.
     Breaks C3 symmetry continuously along the module — so far not
     in the corpus.

  C) Asymmetric interior perturbation (1 variant):
       Only one vertex pair perturbed (vidx=4, opposite vidx=8).
     Extends two-vertex (#129) to targeted-asymmetric.

  D) Double-Gauss anti-phase (1 variant):
       Two correlated Gauss perturbations with opposite sign
       on opposite module halves. Unlike single-gauss
       (#125 → trefoil) where scattering is disruptive.

Output: erzeugte_knoten/kf{42..51}_<descr>/<name>.fseries + wave4_index.json.

Knot machine + fseries_io compatible with the existing pipeline (one-command
pipeline_run.sh, Discovery #117). Sweep via the separate sweep_wave4.sh.
"""
from __future__ import annotations
import site_config
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from knot_machine import build_inner, rot_z  # noqa
from fseries_io import polyline_to_fourier, write_fseries  # noqa

OUT_BASE = Path((site_config.GENERATED_KNOTS_DIR))
MODULE_NPY = ROOT / "double3_module_extracted.npy"


# --------------------------------------------------------------- A) Shear × Pinch -

def perturb_shear_pinch(M, slope, pinch):
    """Combined z-shear (z' = z + slope*x) and radial pinch
    (xy *= 1 + (pinch-1)*|z|/zmax)."""
    Mp = M.copy()
    z = M[1:-1, 2]
    z_abs_norm = np.abs(z) / max(np.abs(M[:, 2]).max(), 1e-9)
    scale = 1.0 + (pinch - 1.0) * z_abs_norm
    Mp[1:-1, 0:2] = M[1:-1, 0:2] * scale[:, None]
    Mp[1:-1, 2] = z + slope * Mp[1:-1, 0]
    return Mp


# --------------------------------------------------------------- B) Helix-Twist -

def perturb_helix_twist(M, twist_total):
    """Rotate interior vertex i in xy by the angle twist_total*(i/L), L=len(M)-1.
    Smooth phase twist along the module; breaks C3 continuously."""
    Mp = M.copy()
    L = len(M) - 1
    for i in range(1, L):
        ang = twist_total * (i / L)
        R = rot_z(ang)
        Mp[i] = R @ M[i]
    return Mp


# --------------------------------------------------------------- C) Asymmetric -

def perturb_one_pair(M, vidx_a, vidx_b, amp):
    """Shift only two opposite vertices in +z resp. -z."""
    Mp = M.copy()
    Mp[vidx_a, 2] += amp
    Mp[vidx_b, 2] -= amp
    return Mp


# --------------------------------------------------------------- D) Double-Gauss anti-phase -

def perturb_double_gauss_anti(M, sigma, seed):
    """Two correlated Gauss perturbations with opposite sign
    on the module halves (1..mid) and (mid+1..L-1). Boundary M[0]/M[-1] fixed."""
    rng = np.random.default_rng(seed)
    Mp = M.copy()
    L = len(M)
    n_inner = L - 2  # interior count
    half = n_inner // 2
    delta_pos = rng.normal(scale=sigma, size=(half, 3))
    delta_neg = -rng.normal(scale=sigma, size=(n_inner - half, 3))
    Mp[1:1 + half] = M[1:1 + half] + delta_pos
    Mp[1 + half:L - 1] = M[1 + half:L - 1] + delta_neg
    return Mp


# --------------------------------------------------------------- VARIANTS -

VARIANTS = [
    # A) Shear × Pinch (6)
    ("kf42_shear-0.3_pinch0.4", perturb_shear_pinch, dict(slope=-0.3, pinch=0.4)),
    ("kf43_shear-0.3_pinch1.6", perturb_shear_pinch, dict(slope=-0.3, pinch=1.6)),
    ("kf44_shear+0.3_pinch0.4", perturb_shear_pinch, dict(slope=+0.3, pinch=0.4)),
    ("kf45_shear+0.3_pinch1.6", perturb_shear_pinch, dict(slope=+0.3, pinch=1.6)),
    ("kf46_shear+0.6_pinch0.4", perturb_shear_pinch, dict(slope=+0.6, pinch=0.4)),
    ("kf47_shear+0.6_pinch1.6", perturb_shear_pinch, dict(slope=+0.6, pinch=1.6)),
    # B) Helix-Twist (2)
    ("kf48_helix_twist_t1.0",   perturb_helix_twist,  dict(twist_total=1.0)),
    ("kf49_helix_twist_t2.0",   perturb_helix_twist,  dict(twist_total=2.0)),
    # C) Asymmetric pair (1)
    ("kf50_pair_v4_v8_a0.5",    perturb_one_pair,     dict(vidx_a=4, vidx_b=8, amp=0.5)),
    # D) Double-Gauss anti (1)
    ("kf51_double_gauss_anti",  perturb_double_gauss_anti, dict(sigma=0.25, seed=2026)),
]


def main():
    M = np.load(MODULE_NPY)
    R120 = rot_z(2 * np.pi / 3)
    closure_err = float(np.linalg.norm(R120 @ M[0] - M[-1]))
    print(f"Loaded module: {M.shape}, original closure = {closure_err:.6f}")
    print()

    OUT_BASE.mkdir(parents=True, exist_ok=True)

    summary = []
    print(f"{'name':<28}  {'bbox':>6}  {'odd':>7}  {'even':>7}  {'closure':>8}  descr")
    for name, fn, kwargs in VARIANTS:
        Mp = fn(M, **kwargs)
        ce = float(np.linalg.norm(R120 @ Mp[0] - Mp[-1]))
        P = build_inner(Mp, drop_seam_duplicates=True, closure_tol=ce + 1e-4)
        a, b = polyline_to_fourier(P, n_terms=20, n_samples=2048, center=True)

        out_dir = OUT_BASE / name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{name}.fseries"
        descr = f"{kwargs}"
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.v1 wave4 ({descr})")

        amps = np.sqrt((a ** 2 + b ** 2).sum(axis=1))
        odd_sum = float(amps[::2].sum())
        even_sum = float(amps[1::2].sum())
        bbox = float(np.abs(P).max())
        summary.append({
            "name": name, "descr": descr,
            "closure_err": ce, "bbox": bbox,
            "polyline_points": int(P.shape[0]),
            "odd_sum": odd_sum, "even_sum": even_sum,
            "amps_top5": amps[:5].tolist(),
            "out": str(out_path),
        })
        print(f"{name:<28}  {bbox:6.3f}  {odd_sum:7.4f}  {even_sum:7.4f}  "
              f"{ce:8.5f}  {descr}")

    with open(OUT_BASE / "wave4_index.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nWave4-Index: {OUT_BASE / 'wave4_index.json'}")
    print(f"{len(summary)} modules generated under {OUT_BASE}")


if __name__ == "__main__":
    main()
