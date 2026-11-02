#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""roundout_check.py — upper bound on the SMOOTH ropelength from a closed polygon (one component), GPU (cupy).

Tao's constant 22a is defined with the smooth (CKS2002) ropelength; ridgerunner / ropelength_check.py give the polygonal
(Rawdon) value.  Procedure of Ashton–Cantarella–Piatek–Rawdon, arXiv:1002.1723 §5 ("roundout"): rounding and exact length as there; the distance bound below is our own
sample-pair version instead of their subsquare elimination:
  1. splice into each corner v_i the circle arc of radius r_i = MinRad(v_i) = min(|e_{i-1}|,|e_i|) / (2 tan(theta_i/2)),
     tangent to both edges (tangent points at distance r_i tan(theta_i/2) <= half the shorter edge from v_i);
     the result V(s) is C^{1,1}, piecewise C^2, minimal radius of curvature r_min = min r_i.
  2. exact length: Len(V(s)) = Len(V) - sum_i (2 r_i tan(theta_i/2) - r_i theta_i).
  3. thickness by Theorem 1.3 of 1002.1723 (Litherland–Simon–Durumeric–Rawdon): Thi = min(r_min, min_dcsd d / 2).
     Every doubly-critical pair is joined by arcs of total curvature >= pi (Fenchel argument, proof of Lemma 2.7 there),
     so with curvature <= 1/r_min its arc separation is >= pi r_min.  Hence
        min_dcsd d >= min { d(s,t) : sample pairs with cyclic separation >= pi r_min - eps } - (1 + K) eps^2,
     sample spacing <= eps along arc length, K = 1/r_min (Proposition 5.1 there; needs d > 1/2).
     Distances are computed in float64 on the GPU (centred coordinates); a margin of 1e-9 covers their rounding.
  4. smooth ropelength upper bound = Len(V(s)) / Thi_lower.
Counter-cases: regular N-gon -> exactly the inscribed circle -> 2 pi (a true lower bound, so the result must be >= 2 pi);
3_1 and 10_124 compared with ACPR's published upper bounds 32.7436 and 71.0739 (both are upper bounds; a lower value
would be a new record to be checked, not a contradiction).
Usage: python -I roundout_check.py file.vect [--eps 0.01]
"""
import argparse, os, sys
import numpy as np
import cupy as cp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ropelength_check import read_vect

ap = argparse.ArgumentParser(); ap.add_argument('files', nargs='+'); ap.add_argument('--eps', type=float, default=0.01)
ap.add_argument('--block', type=int, default=4096)
A = ap.parse_args()


def roundout(P, eps):
    n = len(P)
    e = np.roll(P, -1, 0) - P; le = np.linalg.norm(e, axis=1); u = e / le[:, None]      # edge i: v_i -> v_{i+1}
    uin = np.roll(u, 1, 0); lin = np.roll(le, 1)                                          # edge i-1 into v_i
    th = np.arccos(np.clip((uin * u).sum(1), -1, 1))
    tq = np.tan(th / 2)
    r = np.where(th > 1e-15, np.minimum(lin, le) / (2 * np.where(tq > 0, tq, 1)), np.inf)
    a = np.where(th > 1e-15, np.minimum(lin, le) / 2, 0.0)                                # tangent distance r tan(th/2)
    L_exact = le.sum() - np.sum(np.where(th > 1e-15, 2 * a - r * th, 0.0))
    pts, sarc = [], []; s0 = 0.0
    for i in range(n):
        # arc at v_i
        if th[i] > 1e-15:
            nn = u[i] - (u[i] @ uin[i]) * uin[i]; nn /= np.linalg.norm(nn)
            Apt = P[i] - a[i] * uin[i]; C = Apt + r[i] * nn
            m = max(1, int(np.ceil(r[i] * th[i] / eps)))
            ph = np.arange(m) / m * th[i]
            pts.append(C + r[i] * (-np.cos(ph)[:, None] * nn + np.sin(ph)[:, None] * uin[i]))
            sarc.append(s0 + r[i] * ph); s0 += r[i] * th[i]
        # straight part of edge i
        j = (i + 1) % n
        st = P[i] + a[i] * u[i]; ln = le[i] - a[i] - a[j]
        if ln > 0:
            m = max(1, int(np.ceil(ln / eps)))
            tt = np.arange(m) / m * ln
            pts.append(st + tt[:, None] * u[i]); sarc.append(s0 + tt); s0 += ln
    X = np.vstack(pts); S = np.concatenate(sarc)
    return X, S, s0, L_exact, float(r.min())


def far_min(X, S, Ltot, sepmin, block):
    Xc = cp.asarray(X - X.mean(0), dtype=cp.float64); Sg = cp.asarray(S)
    best = np.inf; m = len(X)
    for i0 in range(0, m, block):
        i1 = min(m, i0 + block)
        d = cp.linalg.norm(Xc[i0:i1, None, :] - Xc[None, :, :], axis=-1)
        sep = cp.abs(Sg[i0:i1, None] - Sg[None, :]); sep = cp.minimum(sep, Ltot - sep)
        d = cp.where(sep >= sepmin, d, cp.inf)
        best = min(best, float(d.min()))
    return best


for f in A.files:
    comps = read_vect(f); assert len(comps) == 1, 'one component only'
    P = comps[0].astype(np.float64)
    Lp = np.linalg.norm(np.roll(P, -1, 0) - P, axis=1).sum()
    X, S, Ls, Lx, rmin = roundout(P, A.eps)
    K = 1.0 / rmin
    eps_eff = float(np.max(np.diff(np.r_[S, Ls])))                       # actual largest sample spacing
    sepmin = np.pi * rmin - eps_eff
    dmin = far_min(X, S, Ls, sepmin, A.block)
    dcsd_lo = dmin - (1 + K) * eps_eff ** 2 - 1e-9
    assert dmin - 2 * eps_eff > 0.5 and dcsd_lo > 0, 'Prop 5.1 needs the critical distances > 1/2 (scale!)'
    thi = min(rmin, dcsd_lo / 2)
    print(f'{f}: poly length {Lp:.6f} | rounded length exact {Lx:.6f} (sampled {Ls:.6f}) r_min {rmin:.6f} '
          f'samples {len(X)} eps {eps_eff:.4g} far-min {dmin:.6f} dcsd>= {dcsd_lo:.6f} thickness>= {thi:.6f} '
          f'SMOOTH ropelength <= {Lx / thi:.6f}', flush=True)
