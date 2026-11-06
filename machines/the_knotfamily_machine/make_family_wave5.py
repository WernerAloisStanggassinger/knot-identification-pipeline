# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave5.py — Wave 5: four new axes.

Follow-up to wave 1-4 + polygon sweeps + shearsweep (discoveries #124-#130,
#144-#146). Waves 1-4 produced kf01-kf51; wave 5 occupies kf52-kf61
(10 variants) along so far unexplored axes:

  A) Global z twist about the polyline axis (3 variants):
       total_twist ∈ {π, 2π, 3π}
     So far (wave 4) there was only "Helix-Twist" on module interior vertices.
     Here the WHOLE polyline (after C_3 build) is rotated about z, with angle
     proportional to the z coordinate — global helixing.

  B) Connected-sum tests for the knots-of-love conjecture (3 variants):
       2x, 3x, 6x trefoil compound along the x axis with bridge-arc
       chaining. Discovery #161 had conjectured: knots-of-love (c=18 d=729)
       could be a 6x K3_1 connected sum (det = 3^6 = 729 ✓).
     These variants test the conjecture algebraically.

  C) Extreme helix twist (2 variants):
       twist_total ∈ {3.0, 4.0}  — wave 4 went up to 2.0.
     Extends the wave-4 scheme into the denser-crossings range.

  D) Bending along the xy axis (2 variants):
       bending of x or y as a function of z (S curve).
     So far only the z axis was perturbed (squash/stretch/flat). Bending breaks
     mirror symmetry and changes the topological position without
     changing the number of crossings much.

Output: erzeugte_knoten/kf{52..61}_<descr>/<name>.fseries + wave5_index.json.
Sweep via the separate sweep_wave5.sh.
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


# ════════════════════════════════════════════════════════════════════════
# A) Global z twist about the polyline axis
# ════════════════════════════════════════════════════════════════════════

def global_z_twist(P: np.ndarray, total_twist: float) -> np.ndarray:
    """Rotate every point of the polyline about the z axis with angle
    proportional to the z height: angle(z) = total_twist * (z - z_min)/(z_max - z_min).

    Acts on the already built C_3 polyline (not on the module).
    Breaks C_3 symmetry continuously — global helixing.
    """
    P = P.copy()
    z = P[:, 2]
    z_min, z_max = z.min(), z.max()
    if z_max - z_min < 1e-12:
        return P
    angles = total_twist * (z - z_min) / (z_max - z_min)
    cos_a = np.cos(angles)
    sin_a = np.sin(angles)
    x_new = cos_a * P[:, 0] - sin_a * P[:, 1]
    y_new = sin_a * P[:, 0] + cos_a * P[:, 1]
    P[:, 0] = x_new
    P[:, 1] = y_new
    return P


# ════════════════════════════════════════════════════════════════════════
# B) Connected-sum tests (knots-of-love conjecture)
# ════════════════════════════════════════════════════════════════════════

def trefoil_polyline(n_pts: int = 200, scale: float = 1.0,
                     center=(0.0, 0.0, 0.0)) -> np.ndarray:
    """Classical 3D trefoil parametrization.

    x(t) = sin(t) + 2 sin(2t)
    y(t) = cos(t) - 2 cos(2t)
    z(t) = -sin(3t)

    Returns (n_pts, 3) — closed polyline, seam open.
    """
    t = np.linspace(0, 2 * np.pi, n_pts, endpoint=False)
    P = np.stack([
        np.sin(t) + 2 * np.sin(2 * t),
        np.cos(t) - 2 * np.cos(2 * t),
        -np.sin(3 * t),
    ], axis=1)
    return P * scale + np.array(center)


def connectsum_chain(k: int, scale: float = 0.4, x_sep: float = 3.0,
                     y_dip: float = -4.0, n_per: int = 200,
                     bridge_n: int = 25, closing_n: int = 80) -> np.ndarray:
    """k trefoils along x axis, chained via short bridges, closed via
    a long arc below the trefoils (y_dip).

    Each trefoil keeps its main part (middle 140 pts of 200),
    the first/last points are interpreted as "cut faces"
    and connected by linear bridges to the neighboring trefoil.

    The closing bridge goes from the rightmost trefoil back to the leftmost
    via a parabola arc that runs below all trefoils (y=y_dip) —
    so no crossing with the trefoils themselves.
    """
    trefs = [trefoil_polyline(n_pts=n_per, scale=scale,
                              center=(i * x_sep, 0.0, 0.0))
             for i in range(k)]

    # main body of each trefoil (middle 70 %)
    cut = n_per // 7
    bodies = [T[cut:n_per - cut] for T in trefs]

    chain = [bodies[0]]
    # linear bridges between consecutive trefoils
    for i in range(1, k):
        p_start = bodies[i - 1][-1]
        p_end = bodies[i][0]
        bridge_t = np.linspace(0.0, 1.0, bridge_n + 2)[1:-1].reshape(-1, 1)
        bridge = (1 - bridge_t) * p_start + bridge_t * p_end
        chain.append(bridge)
        chain.append(bodies[i])

    # Closing arc: from the last trefoil end back to the first trefoil start,
    # via a parabola through (mid_x, y_dip, 0)
    p_start = bodies[-1][-1]
    p_end = bodies[0][0]
    mid_x = (p_start[0] + p_end[0]) / 2.0
    t = np.linspace(0.0, 1.0, closing_n + 2)[1:-1]
    # quadratic Bezier-like curve through control point (mid_x, y_dip, 0)
    ctrl = np.array([mid_x, y_dip, 0.0])
    closing = ((1 - t)[:, None] ** 2) * p_start \
              + (2 * (1 - t) * t)[:, None] * ctrl \
              + (t[:, None] ** 2) * p_end
    chain.append(closing)

    return np.vstack(chain)


# ════════════════════════════════════════════════════════════════════════
# C) Extreme helix twist (vertices, module level — like wave 4 only stronger)
# ════════════════════════════════════════════════════════════════════════

def perturb_helix_twist(M: np.ndarray, twist_total: float) -> np.ndarray:
    """Identical to wave 4 perturb_helix_twist. Rotate interior vertex i
    in xy by angle twist_total*(i/L), L=len(M)-1."""
    Mp = M.copy()
    L = len(M) - 1
    for i in range(1, L):
        ang = twist_total * (i / L)
        R = rot_z(ang)
        Mp[i] = R @ M[i]
    return Mp


# ════════════════════════════════════════════════════════════════════════
# D) Bending along the xy axis (in the module)
# ════════════════════════════════════════════════════════════════════════

def perturb_bend(M: np.ndarray, axis: str, amp: float) -> np.ndarray:
    """Bending of the module: x (or y) is shifted by a function of z.

    bend_x: x_new = x + amp * sin(pi * (z - z_min)/(z_max - z_min))
    bend_y: y_new = y + amp * sin(pi * (z - z_min)/(z_max - z_min))

    Acts only on interior vertices (module endpoints stay fixed for closure).
    """
    Mp = M.copy()
    z_inner = M[1:-1, 2]
    z_min, z_max = M[:, 2].min(), M[:, 2].max()
    if z_max - z_min < 1e-12:
        return Mp
    phase = np.pi * (z_inner - z_min) / (z_max - z_min)
    delta = amp * np.sin(phase)
    if axis == "x":
        Mp[1:-1, 0] = M[1:-1, 0] + delta
    elif axis == "y":
        Mp[1:-1, 1] = M[1:-1, 1] + delta
    else:
        raise ValueError(f"axis must be 'x' or 'y' (got: {axis})")
    return Mp


# ════════════════════════════════════════════════════════════════════════
# variants table: (name, mode, params) — mode ∈
#   "module_perturb_then_inner": perturb fn on double3 module, then build_inner
#   "polyline_post":             build_inner first, then post-perturb on polyline
#   "direct_polyline":           a polyline directly, without module build
# ════════════════════════════════════════════════════════════════════════

VARIANTS = [
    # A) Global z twist (post on the C_3 polyline)
    ("kf52_global_twist_pi",   "polyline_post", dict(fn=global_z_twist, total_twist=np.pi)),
    ("kf53_global_twist_2pi",  "polyline_post", dict(fn=global_z_twist, total_twist=2*np.pi)),
    ("kf54_global_twist_3pi",  "polyline_post", dict(fn=global_z_twist, total_twist=3*np.pi)),
    # B) Connected-sum tests
    ("kf55_connectsum_2tref",  "direct_polyline", dict(fn=connectsum_chain, k=2)),
    ("kf56_connectsum_3tref",  "direct_polyline", dict(fn=connectsum_chain, k=3)),
    ("kf57_connectsum_6tref",  "direct_polyline", dict(fn=connectsum_chain, k=6)),
    # C) Extreme helix twist (module level, like wave 4)
    ("kf58_helix_twist_t3.0",  "module_perturb_then_inner", dict(fn=perturb_helix_twist, twist_total=3.0)),
    ("kf59_helix_twist_t4.0",  "module_perturb_then_inner", dict(fn=perturb_helix_twist, twist_total=4.0)),
    # D) Bending xy axis
    ("kf60_bend_x_amp0.3",     "module_perturb_then_inner", dict(fn=perturb_bend, axis="x", amp=0.3)),
    ("kf61_bend_y_amp0.3",     "module_perturb_then_inner", dict(fn=perturb_bend, axis="y", amp=0.3)),
]


def main():
    M = np.load(MODULE_NPY)
    R120 = rot_z(2 * np.pi / 3)
    closure_err = float(np.linalg.norm(R120 @ M[0] - M[-1]))
    print(f"Loaded module: {M.shape}, Original closure = {closure_err:.6f}")
    print()

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    summary = []
    print(f"{'name':<28}  {'mode':<28}  {'pts':>5}  {'bbox':>6}  {'odd':>7}  {'even':>7}")

    for name, mode, kwargs in VARIANTS:
        fn = kwargs.pop("fn")

        if mode == "module_perturb_then_inner":
            Mp = fn(M, **kwargs)
            ce = float(np.linalg.norm(R120 @ Mp[0] - Mp[-1]))
            P = build_inner(Mp, drop_seam_duplicates=True, closure_tol=ce + 1e-4)
        elif mode == "polyline_post":
            # First C_3 build on unperturbed module, then fn on polyline
            P0 = build_inner(M, drop_seam_duplicates=True)
            P = fn(P0, **kwargs)
        elif mode == "direct_polyline":
            P = fn(**kwargs)
        else:
            print(f"  {name}: unknown mode {mode}, skip")
            continue

        try:
            a, b = polyline_to_fourier(P, n_terms=20, n_samples=2048, center=True)
        except Exception as e:
            print(f"  {name}: Fourier error: {e}")
            summary.append({"name": name, "error": str(e)})
            continue

        out_dir = OUT_BASE / name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{name}.fseries"
        descr = f"{mode} {kwargs}"
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.v1 wave5 ({descr})")

        amps = np.sqrt((a ** 2 + b ** 2).sum(axis=1))
        odd_sum = float(amps[::2].sum())
        even_sum = float(amps[1::2].sum())
        bbox = float(np.abs(P).max())
        summary.append({
            "name": name, "mode": mode, "descr": descr,
            "polyline_points": int(P.shape[0]),
            "bbox": bbox, "odd_sum": odd_sum, "even_sum": even_sum,
            "amps_top5": amps[:5].tolist(),
            "out": str(out_path),
        })
        print(f"{name:<28}  {mode:<28}  {P.shape[0]:>5}  "
              f"{bbox:6.3f}  {odd_sum:7.4f}  {even_sum:7.4f}")

    with open(OUT_BASE / "wave5_index.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nwave5 index: {OUT_BASE / 'wave5_index.json'}")
    print(f"{len(summary)} modules generated under {OUT_BASE}")


if __name__ == "__main__":
    main()
