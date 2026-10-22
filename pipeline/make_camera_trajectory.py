#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
make_camera_trajectory.py
=========================
Builds a generic camera flight along a Fibonacci spiral from pole to pole.
N camera points, numbered sequentially t_step = 0..N-1; the time stamp is the
index along the spiral.

A generic, knot-independent tool, generated once per N. The scanner loads it as
input instead of computing the trajectory inline, so the camera flight is an
artifact of its own.

Convention:
  theta(i) grows monotonically from ~0 (north pole) to ~pi (south pole).
  phi(i) advances by the golden angle: 2*pi*i/golden mod 2*pi.
  -> one continuous spiral; every index has a unique predecessor and successor.

Outputs (two files with a common stem):
  <stem>.csv   human-readable:
                 t_step, theta, phi, cam_x, cam_y, cam_z
  <stem>.npz   compact, for the scanner:
                 t_step, theta, phi, cam (N x 3), rotmat (N x 3 x 3)

The rotation matrix R(i) maps Z = (0,0,1) to cam(i). The scanner applies
p -> R(i) p to the curve and projects onto the xy-plane, i.e. frame i views the
fixed knot along R(i)^T e_z.
"""

import sys
import os
import math
import csv
import argparse
import numpy as np


# ──────────────────────────────────────────────
# 1. FIBONACCI SPIRAL pole -> pole
# ──────────────────────────────────────────────

def fibonacci_trajectory(n_frames):
    """
    Returns the sequential pole-to-pole spiral.

    Returns
    -------
    theta : (N,) float64
    phi   : (N,) float64
    cam   : (N, 3) float32   unit vectors of the camera directions
    """
    golden = (1 + math.sqrt(5)) / 2
    indices = np.arange(n_frames, dtype=np.float64)
    theta = np.arccos(1 - 2 * (indices + 0.5) / n_frames)
    phi   = (2 * math.pi * indices / golden) % (2 * math.pi)
    cam_x = np.sin(theta) * np.cos(phi)
    cam_y = np.sin(theta) * np.sin(phi)
    cam_z = np.cos(theta)
    cam   = np.stack([cam_x, cam_y, cam_z], axis=1).astype(np.float32)
    return theta, phi, cam


# ──────────────────────────────────────────────
# 2. ROTATION MATRICES R: Z -> cam(i)
# ──────────────────────────────────────────────

def rotmats_from_directions(cam):
    """
    For every camera direction d the rotation that maps Z = (0,0,1) to d
    (Rodrigues formula about the axis Z x d).

    Bit-identical to the fibonacci_sphere logic of an earlier scanner version
    (knot_scanner_v3), so the same trajectory reproduces the same 2D projections.
    """
    n = len(cam)
    rotmats = np.zeros((n, 3, 3), dtype=np.float32)
    Z = np.array([0.0, 0.0, 1.0])
    for i in range(n):
        d = cam[i].astype(np.float64)
        cos_a = float(np.clip(np.dot(Z, d), -1, 1))
        if abs(cos_a - 1.0) < 1e-9:
            rotmats[i] = np.eye(3, dtype=np.float32)
        elif abs(cos_a + 1.0) < 1e-9:
            rotmats[i] = np.diag([1.0, -1.0, -1.0]).astype(np.float32)
        else:
            axis = np.cross(Z, d)
            axis = axis / np.linalg.norm(axis)
            sin_a = math.sqrt(1 - cos_a ** 2)
            K = np.array([
                [    0,      -axis[2],  axis[1]],
                [ axis[2],       0,    -axis[0]],
                [-axis[1],  axis[0],       0   ]
            ])
            R = cos_a * np.eye(3) + sin_a * K + (1 - cos_a) * np.outer(axis, axis)
            rotmats[i] = R.astype(np.float32)
    return rotmats


# ──────────────────────────────────────────────
# 3. ENTRY
# ──────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, default=131072,
                        help='number of camera points (default 2^17 = 131072)')
    parser.add_argument('--out', type=str, default=None,
                        help='stem of the two output files. '
                             'Default: ./camera_<N>')
    args = parser.parse_args()

    N = args.n
    stem = args.out or f'camera_{N}'
    csv_path = stem + '.csv'
    npz_path = stem + '.npz'

    print(f"[CAMERA] Building {N} camera points (Fibonacci spiral pole -> pole)")
    theta, phi, cam = fibonacci_trajectory(N)
    print(f"[CAMERA]   theta range: {theta.min():.6f} .. {theta.max():.6f} rad "
          f"({math.degrees(theta.min()):.2f}° .. {math.degrees(theta.max()):.2f}°)")
    print(f"[CAMERA]   phi range  : {phi.min():.6f} .. {phi.max():.6f} rad")

    print(f"[CAMERA] Computing {N} rotation matrices ...")
    rotmats = rotmats_from_directions(cam)

    print(f"[CAMERA] Writing {csv_path}")
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['t_step', 'theta', 'phi', 'cam_x', 'cam_y', 'cam_z'])
        for i in range(N):
            w.writerow([i,
                        f'{theta[i]:.8f}', f'{phi[i]:.8f}',
                        f'{cam[i, 0]:.8f}', f'{cam[i, 1]:.8f}', f'{cam[i, 2]:.8f}'])

    print(f"[CAMERA] Writing {npz_path}")
    np.savez(
        npz_path,
        t_step=np.arange(N, dtype=np.int32),
        theta=theta.astype(np.float64),
        phi=phi.astype(np.float64),
        cam=cam,
        rotmat=rotmats,
    )

    sz_csv = os.path.getsize(csv_path)
    sz_npz = os.path.getsize(npz_path)
    print(f"[CAMERA] Done.")
    print(f"  {csv_path}  ({sz_csv/1e6:.1f} MB, {N} rows)")
    print(f"  {npz_path}  ({sz_npz/1e6:.1f} MB, t_step+theta+phi+cam[{N},3]+rotmat[{N},3,3])")


if __name__ == '__main__':
    main()
