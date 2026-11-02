#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""build_torus_knot.py — initial conformation of the torus knot T(p,q) on a standard torus, as a VECT polygon.

gamma(t) = ((R + r cos(p t)) cos(q t), (R + r cos(p t)) sin(q t), r sin(p t)),  t in [0, 2 pi):
q turns around the z-axis, p turns around the core circle.  A meridian disc is met q times (spacing 2 pi r / q),
so unit-radius tubes need r >= q/pi; the inner equator carries p passes (spacing 2 pi (R - r) / p).
Radii get a safety factor; ridgerunner -a rescales anyway. Thickness is checked afterwards with ropelength_check.py.

Usage: python3 -I build_torus_knot.py p q n_vertices out.vect [margin]
"""
import sys
import numpy as np

p, q, n, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
margin = float(sys.argv[5]) if len(sys.argv) > 5 else 1.3
r = margin * q / np.pi
R = r + margin * p / np.pi
t = np.linspace(0, 2 * np.pi, n, endpoint=False)
P = np.c_[(R + r * np.cos(p * t)) * np.cos(q * t), (R + r * np.cos(p * t)) * np.sin(q * t), r * np.sin(p * t)]
with open(out, 'w') as f:
    f.write(f'VECT\n1 {n} 0\n-{n}\n0\n')
    for x, y, z in P:
        f.write(f'{x:.12f} {y:.12f} {z:.12f}\n')
print(f'T({p},{q}) n={n} R={R:.4f} r={r:.4f} -> {out}')
