# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_connectsum_ring.py — ring-shaped connected-sum construction.

Follow-up to Discovery #162 (wave 5 connected sums kf55/56/57 with linear
x chain + parabolic closing arc were topologically unclean, all exit=3).

New approach: k trefoils arranged on a circle of radius R, each trefoil
locally oriented tangential to the ring (seam direction = ring tangent at the
trefoil centre). Consecutive trefoils are connected by arc bridges
connected along the ring — the bridges run along a
circular-arc curve around the centre (instead of cutting through the interior).

Topologically, k trefoils connected in a ring: K3_1 # K3_1 # ... # K3_1
(k-fold connected sum). Expected algebra:
  c = 3k        (each trefoil knot = 3 crossings)
  det = 3^k     (determinant multiplies under connected sum)
  signature = -2k (or +2k depending on convention; trefoil has sigma=-2)

knots-of-love conjecture (Discovery #161): c=18 d=729 = 6x trefoil.
  c = 6*3 = 18 ✓
  det = 3^6 = 729 ✓
If kf64_csring_6tref delivers the same algebra as knots-of-love, the
conjecture is verified.

Output: erzeugte_knoten/kf{62,63,64}_csring_<k>tref/<name>.fseries.
Sweep separately via pipeline_run.sh.
"""
from __future__ import annotations
import site_config
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from knot_machine import rot_z  # noqa: E402
from fseries_io import polyline_to_fourier, write_fseries  # noqa: E402

OUT_BASE = Path((site_config.GENERATED_KNOTS_DIR))


def trefoil_polyline(n_pts: int = 200, scale: float = 1.0,
                     center=(0.0, 0.0, 0.0),
                     rot_z_angle: float = 0.0) -> np.ndarray:
    """Classical 3D trefoil parametrisation with optional z rotation.

    x(t) = sin(t) + 2 sin(2t)
    y(t) = cos(t) - 2 cos(2t)
    z(t) = -sin(3t)

    Tangent at t=0 points in the +x direction — when rotating by z-angle theta
    the tangent at t=0 points in the direction (cos theta, sin theta, 0).
    """
    t = np.linspace(0, 2 * np.pi, n_pts, endpoint=False)
    P = np.stack([
        np.sin(t) + 2 * np.sin(2 * t),
        np.cos(t) - 2 * np.cos(2 * t),
        -np.sin(3 * t),
    ], axis=1)
    P = P * scale
    if rot_z_angle != 0.0:
        R = rot_z(rot_z_angle)
        P = P @ R.T
    return P + np.array(center)


def arc_bridge(p_start: np.ndarray, p_end: np.ndarray,
               ring_center: np.ndarray, n_pts: int = 25) -> np.ndarray:
    """Arc bridge from p_start to p_end that runs around ring_center in the xy
    plane. Interpolates linearly in (polar angle, radius, z).

    So the bridge follows the ring instead of going through its interior —
    Prerequisite that the trefoils lie near the ring radius and the
    neighboring trefoils are on the same ring.
    """
    cx, cy = ring_center[0], ring_center[1]
    v_start = np.array([p_start[0] - cx, p_start[1] - cy])
    v_end = np.array([p_end[0] - cx, p_end[1] - cy])
    r_start = np.linalg.norm(v_start)
    r_end = np.linalg.norm(v_end)
    theta_start = np.arctan2(v_start[1], v_start[0])
    theta_end = np.arctan2(v_end[1], v_end[0])
    # angle difference: force the shortest arc direction
    dtheta = theta_end - theta_start
    if dtheta > np.pi:
        dtheta -= 2 * np.pi
    elif dtheta < -np.pi:
        dtheta += 2 * np.pi
    t = np.linspace(0.0, 1.0, n_pts + 2)[1:-1]
    angles = theta_start + dtheta * t
    radii = r_start + (r_end - r_start) * t
    z_vals = p_start[2] + (p_end[2] - p_start[2]) * t
    return np.stack([
        cx + radii * np.cos(angles),
        cy + radii * np.sin(angles),
        z_vals,
    ], axis=1)


def connectsum_ring(k: int, scale: float = 0.4, R: float = 2.5,
                    n_per_trefoil: int = 200, bridge_n: int = 30,
                    cut_frac: float = 1.0 / 8.0,
                    z_jitter: float = 0.0) -> np.ndarray:
    """k trefoils in a ring on a circle of radius R, connected via arc bridges.

    Each trefoil is rotated so that its seam tangent (in +x at t=0)
    points in the ring tangent direction at the trefoil centre. Seams between
    consecutive trefoils are thereby geometrically flat.

    cut_frac: fraction of the first/last points of the trefoil that is left out
              becomes (so that the bridge cleanly at the remaining places
              attaches). Default 1/8 — main body = 75% of the trefoil.

    z_jitter: optional small offset per trefoil in z (avoids possible
              same-plane crossings). Default 0.
    """
    ring_center = np.array([0.0, 0.0, 0.0])
    bodies = []
    for i in range(k):
        theta_ring = 2 * np.pi * i / k
        # trefoil tangent at t=0 points in +x; after rotation by
        # rot_z_angle = theta_ring + pi/2 it points in
        # (-sin theta_ring, cos theta_ring, 0) = tangential ring direction.
        rot_angle = theta_ring + np.pi / 2
        center = (R * np.cos(theta_ring), R * np.sin(theta_ring),
                  i * z_jitter)
        T = trefoil_polyline(n_pts=n_per_trefoil, scale=scale,
                             center=center, rot_z_angle=rot_angle)
        cut = max(1, int(n_per_trefoil * cut_frac))
        body = T[cut:n_per_trefoil - cut]
        bodies.append(body)

    chain = []
    for i in range(k):
        chain.append(bodies[i])
        # bridge to the next trefoil (cyclic)
        next_i = (i + 1) % k
        p_end = bodies[i][-1]
        p_start = bodies[next_i][0]
        bridge = arc_bridge(p_end, p_start, ring_center, n_pts=bridge_n)
        chain.append(bridge)

    return np.vstack(chain)


VARIANTS = [
    # (name, k, scale, R)
    # scale small compared to R, so that trefoils stay local + do not
    # overlap with neighbors. At k=6, 2*scale < 2*pi*R/6 = pi*R/3 must hold
    # → R > 6*scale/pi ≈ 0.95 for scale=0.5.
    ("kf62_csring_2tref", 2, 0.5, 2.0),
    ("kf63_csring_3tref", 3, 0.5, 2.2),
    ("kf64_csring_6tref", 6, 0.4, 3.0),
]


def main():
    OUT_BASE.mkdir(parents=True, exist_ok=True)
    summary = []
    print(f"{'name':<24}  {'k':>3}  {'scale':>6}  {'R':>5}  {'pts':>5}  "
          f"{'bbox':>6}  {'odd':>7}  {'even':>7}")

    for name, k, scale, R in VARIANTS:
        P = connectsum_ring(k=k, scale=scale, R=R)
        try:
            a, b = polyline_to_fourier(P, n_terms=20, n_samples=2048,
                                       center=True)
        except Exception as e:
            print(f"  {name}: Fourier error: {e}")
            summary.append({"name": name, "error": str(e)})
            continue

        out_dir = OUT_BASE / name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{name}.fseries"
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.v1 csring (k={k} scale={scale} R={R})")

        amps = np.sqrt((a ** 2 + b ** 2).sum(axis=1))
        odd_sum = float(amps[::2].sum())
        even_sum = float(amps[1::2].sum())
        bbox = float(np.abs(P).max())
        summary.append({
            "name": name, "k_trefoils": k, "scale": scale, "R": R,
            "polyline_points": int(P.shape[0]),
            "bbox": bbox, "odd_sum": odd_sum, "even_sum": even_sum,
            "amps_top5": amps[:5].tolist(),
            "out": str(out_path),
            "expected_algebra": {
                "c": 3 * k, "det": 3 ** k, "amphi": False,
                "vol_hyperbolic": "non-hyperbolic (connected sum)",
            },
        })
        print(f"{name:<24}  {k:>3}  {scale:>6.2f}  {R:>5.2f}  "
              f"{P.shape[0]:>5}  {bbox:6.3f}  {odd_sum:7.4f}  {even_sum:7.4f}")

    with open(OUT_BASE / "csring_index.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nIndex: {OUT_BASE / 'csring_index.json'}")
    print(f"{len(summary)} modules generated.")


if __name__ == "__main__":
    main()
