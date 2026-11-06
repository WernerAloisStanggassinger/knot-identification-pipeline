# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
Wave 3: 12 more variants with targeted strategies:

A) OUTER S₆ generator validation: 4 synthetic 6-vertex modules,
   build_outer (S₆ Schoenflies, 6 module copies). Closure: M[-1] = g₁·M[0]
   with g₁ = R₆₀(z) ∘ σ_h. So far never tested on a construction.

B) Multi-vertex perturbations: deliberate disturbance patterns on several
   interior vertices simultaneously (alternating, spiral, counter-rotating).

C) Extreme combinations: shear + rotate + pinch in one,
   as well as module-reverse (M[::-1]) for inversed traversal.
"""
from __future__ import annotations
import site_config
import numpy as np
import sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from knot_machine import build_inner, build_outer, rot_z, s6_generator, make_module_endpoint  # noqa
from fseries_io import polyline_to_fourier, write_fseries  # noqa

OUT_BASE = Path((site_config.GENERATED_KNOTS_DIR))
MODULE_NPY = ROOT / "double3_module_extracted.npy"


# --------------------------------------------------------------- A) OUTER -

def make_outer_module_triangle(scale=0.5, height=0.3):
    """6-vertex module with M[0] and M[-1]=g₁(M[0]), interior forms
    a triangle shape with height offset."""
    G = s6_generator()
    M0 = np.array([scale, 0.0, height])
    M5 = G @ M0
    # 4 interior points: linear interpolation + triangle tip
    M = np.zeros((6, 3))
    M[0] = M0
    M[5] = M5
    t_vals = np.linspace(0.2, 0.8, 4)
    for i, t in enumerate(t_vals, start=1):
        base = (1 - t) * M0 + t * M5
        # shift triangle tip outward (radially in xy)
        r = np.sqrt(base[0]**2 + base[1]**2)
        if r > 0:
            offset = 0.3 * np.array([base[0]/r, base[1]/r, 0])
            # move the tip along in z as well
            z_offset = 0.5 * np.sin(t * np.pi)
            M[i] = base + offset + np.array([0, 0, z_offset])
        else:
            M[i] = base
    return M


def make_outer_module_spiral(scale=0.5, twist=2.0):
    """Spiral 6-vertex path from M[0] to g₁(M[0])."""
    G = s6_generator()
    M0 = np.array([scale, 0.0, 0.4])
    M5 = G @ M0
    M = np.zeros((6, 3))
    M[0] = M0
    M[5] = M5
    for i in range(1, 5):
        t = i / 5
        base = (1 - t) * M0 + t * M5
        r = scale * (0.5 + 0.3 * np.sin(t * np.pi))
        ang = t * twist * 2 * np.pi
        M[i] = np.array([r * np.cos(ang), r * np.sin(ang), base[2]])
    return M


def make_outer_module_zigzag(scale=0.5, amp=0.4):
    G = s6_generator()
    M0 = np.array([scale, 0.0, 0.3])
    M5 = G @ M0
    M = np.zeros((6, 3))
    M[0] = M0
    M[5] = M5
    for i in range(1, 5):
        t = i / 5
        base = (1 - t) * M0 + t * M5
        z_off = amp * (1 if i % 2 == 0 else -1)
        M[i] = base + np.array([0, 0, z_off])
    return M


def make_outer_module_smooth(scale=0.5):
    """Smooth Bezier-like path from M[0] to g₁(M[0]) without jumps."""
    G = s6_generator()
    M0 = np.array([scale, -scale*0.3, 0.4])
    M5 = G @ M0
    M = np.zeros((6, 3))
    M[0] = M0
    M[5] = M5
    # Quadratic Bezier with control point outside
    P_ctrl = np.array([0.0, 0.0, -0.3])
    for i in range(1, 5):
        t = i / 5
        # Quadratic Bezier B(t) = (1-t)^2 P0 + 2(1-t)t P_c + t^2 P5
        M[i] = (1-t)**2 * M0 + 2*(1-t)*t*P_ctrl + t**2 * M5
    return M


# --------------------------------------------------------------- B) MULTI -

def perturb_multi_alternating(M, amp=0.3):
    """Alternating ±amp on z of the interior vertices."""
    Mp = M.copy()
    for i in range(1, len(M) - 1):
        sign = 1 if i % 2 == 0 else -1
        Mp[i, 2] = M[i, 2] + sign * amp
    return Mp


def perturb_multi_two_vertex(M, v1=3, v2=8, amp=0.4):
    """Two interior vertices counter-rotating in y."""
    Mp = M.copy()
    if 0 < v1 < len(M) - 1:
        Mp[v1, 1] += amp
    if 0 < v2 < len(M) - 1:
        Mp[v2, 1] -= amp
    return Mp


def perturb_spiral_overlay(M, amp=0.2, freq=2):
    """Spiral disturbance on all interior vertices in xy."""
    Mp = M.copy()
    n_int = len(M) - 2
    for i in range(1, len(M) - 1):
        t = (i - 1) / max(1, n_int - 1)
        ang = t * freq * 2 * np.pi
        Mp[i, 0] += amp * np.cos(ang)
        Mp[i, 1] += amp * np.sin(ang)
    return Mp


def perturb_radial_outward(M, amp=0.3):
    """Push all interior vertices radially outward."""
    Mp = M.copy()
    for i in range(1, len(M) - 1):
        x, y, z = M[i]
        r = np.sqrt(x*x + y*y)
        if r > 1e-9:
            Mp[i, 0] += amp * x / r
            Mp[i, 1] += amp * y / r
    return Mp


# --------------------------------------------------------------- C) EXTREME -

def perturb_compound(M, slope=0.6, angle_deg=15, pinch=0.7):
    """shear + rotate + pinch in one go."""
    Mp = M.copy()
    # shear
    Mp[1:-1, 2] = M[1:-1, 2] + slope * M[1:-1, 0]
    # rotate xy
    R = rot_z(np.deg2rad(angle_deg))
    Mp[1:-1] = Mp[1:-1] @ R.T
    # radial pinch
    z_abs_norm = np.abs(Mp[1:-1, 2]) / max(np.abs(Mp[:, 2]).max(), 1e-9)
    scale = 1.0 + (pinch - 1.0) * z_abs_norm
    Mp[1:-1, 0:2] = Mp[1:-1, 0:2] * scale[:, None]
    return Mp


def perturb_huge_gauss(M, sigma=0.4, seed=2026):
    rng = np.random.default_rng(seed)
    Mp = M.copy()
    Mp[1:-1] = M[1:-1] + rng.normal(scale=sigma, size=(len(M) - 2, 3))
    return Mp


def perturb_module_reverse(M):
    """Traverse module backwards (M[::-1] is NOT closure-valid,
    but we can rebuild it consistently by swapping M[0]<->M[-1]
    and heeding the new closure condition: the reversed module has
    M_rev[-1] = R_-120·M_rev[0] (inverse rotation)."""
    Mp = M[::-1].copy()
    return Mp


# --------------------------------------------------------------- VARIANTS -

VARIANTS = []

# A) OUTER S6 — 4 synthetic modules
for name, factory, kwargs in [
    ("kf25_outer_triangle", make_outer_module_triangle, dict(scale=0.5, height=0.3)),
    ("kf26_outer_spiral",   make_outer_module_spiral,   dict(scale=0.5, twist=2.0)),
    ("kf27_outer_zigzag",   make_outer_module_zigzag,   dict(scale=0.5, amp=0.4)),
    ("kf28_outer_smooth",   make_outer_module_smooth,   dict(scale=0.5)),
]:
    VARIANTS.append((name, "outer", factory, kwargs))

# B) multi-vertex perturbations on the double-3 module
for name, fn, kwargs in [
    ("kf29_alt_z+0.3",       perturb_multi_alternating,  dict(amp=0.3)),
    ("kf30_two_vertex",      perturb_multi_two_vertex,   dict(v1=3, v2=8, amp=0.4)),
    ("kf31_spiral_overlay",  perturb_spiral_overlay,     dict(amp=0.2, freq=2)),
    ("kf32_radial_outward",  perturb_radial_outward,     dict(amp=0.3)),
]:
    VARIANTS.append((name, "inner", fn, kwargs))

# C) Extreme combinations
for name, fn, kwargs in [
    ("kf33_compound",        perturb_compound,           dict(slope=0.6, angle_deg=15, pinch=0.7)),
    ("kf34_huge_gauss_s2026", perturb_huge_gauss,        dict(sigma=0.4, seed=2026)),
    ("kf35_huge_gauss_s333", perturb_huge_gauss,         dict(sigma=0.5, seed=333)),
    ("kf36_module_reverse",  perturb_module_reverse,     dict()),
]:
    VARIANTS.append((name, "inner", fn, kwargs))


def main():
    M_inner = np.load(MODULE_NPY)
    R120 = rot_z(2 * np.pi / 3)
    G6 = s6_generator()
    closure_inner = float(np.linalg.norm(R120 @ M_inner[0] - M_inner[-1]))

    summary = []
    for name, mode, factory_or_fn, kwargs in VARIANTS:
        if mode == "outer":
            # synthetic module construction
            M = factory_or_fn(**kwargs)
            # Validate closure
            ce = float(np.linalg.norm(G6 @ M[0] - M[-1]))
            if ce > 1e-3:
                print(f"  WARN {name}: outer closure error {ce:.5f}")
            P = build_outer(M, drop_seam_duplicates=True, closure_tol=ce + 1e-4)
        else:
            Mp = factory_or_fn(M_inner, **kwargs)
            ce = float(np.linalg.norm(R120 @ Mp[0] - Mp[-1]))
            P = build_inner(Mp, drop_seam_duplicates=True, closure_tol=ce + 1e-4)

        # Fourier
        a, b = polyline_to_fourier(P, n_terms=20, n_samples=2048, center=True)
        out_dir = OUT_BASE / name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{name}.fseries"
        descr = f"{mode} | {kwargs}"
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.v1 wave3 ({descr})")

        amps = np.sqrt((a ** 2 + b ** 2).sum(axis=1))
        bbox = float(np.abs(P).max())
        summary.append({
            "name": name, "mode": mode, "descr": descr,
            "closure_err": ce, "bbox": bbox,
            "polyline_points": int(P.shape[0]),
            "odd_sum": float(amps[::2].sum()),
            "even_sum": float(amps[1::2].sum()),
            "amps_top5": amps[:5].tolist(),
        })
        print(f"  {name:<22} mode={mode}  bbox={bbox:.3f}  pts={P.shape[0]:>4}  "
              f"odd={summary[-1]['odd_sum']:.3f}  even={summary[-1]['even_sum']:.3f}  "
              f"closure={ce:.4f}")

    with open(OUT_BASE / "wave3_index.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nWave3-Index: {OUT_BASE / 'wave3_index.json'}")


if __name__ == "__main__":
    main()
