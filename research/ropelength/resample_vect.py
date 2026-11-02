#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""Resample every closed component of a VECT polygon to N equally spaced points (arc length, linear).
Usage: python3 -I resample_vect.py in.vect out.vect N_total [spline]   (default linear; spline = periodic cubic)"""
import sys, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__file__))
from ropelength_check import read_vect
comps = read_vect(sys.argv[1]); N = int(sys.argv[3]); SPLINE = len(sys.argv) > 4 and sys.argv[4] == 'spline'
lens = [np.linalg.norm(np.roll(P, -1, 0) - P, axis=1).sum() for P in comps]
out = []
for P, L in zip(comps, lens):
    n = max(8, round(N * L / sum(lens)))
    Q = np.vstack([P, P[:1]]); s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(Q, axis=0), axis=1))]
    u = np.linspace(0, s[-1], n, endpoint=False)
    if SPLINE:
        from scipy.interpolate import CubicSpline
        cs = CubicSpline(s, Q, bc_type='periodic')       # smooth periodic interpolation, no kinks
        out.append(cs(u))
    else:
        out.append(np.c_[[np.interp(u, s, Q[:, k]) for k in range(3)]].T)
with open(sys.argv[2], 'w') as f:
    f.write(f'VECT\n{len(out)} {sum(map(len, out))} 0\n' + ' '.join(f'-{len(o)}' for o in out) + '\n' + ' '.join('0' for _ in out) + '\n')
    for o in out:
        for x, y, z in o: f.write(f'{x:.12f} {y:.12f} {z:.12f}\n')
