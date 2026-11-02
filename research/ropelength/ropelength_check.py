#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""ropelength_check.py — independent polygonal ropelength of a closed polygon (VECT file).

Rawdon's polygonal thickness (unit-RADIUS convention, as in ridgerunner/KnotInfo):
    thickness = min( min_v minRad(v) , dcsd / 2 )
    minRad(v) = min(|e_in|, |e_out|) / (2 tan(theta_v / 2)),   theta_v = turning angle at v
    dcsd      = doubly critical self distance: min |x - y| over pairs x != y on the polygon
                with (y - x) in the normal cone at x and (x - y) in the normal cone at y.
    ropelength = length / thickness.
Doubly critical pairs are searched among the closest-point pairs of all non-adjacent edge pairs
(and all vertex/edge combinations those produce); a closest pair is kept only if it is critical
at both ends. Pure numpy, float64 — a calibration tool, not interval arithmetic.

Usage: python3 -I ropelength_check.py file.vect [more.vect ...]
"""
import sys
import numpy as np


def read_vect(path):
    tok = ' '.join(line.split('#')[0] for line in open(path)).split()
    assert tok[0] == 'VECT', path
    ncomp, nv, _nc = int(tok[1]), int(tok[2]), int(tok[3])
    counts = [int(t) for t in tok[4:4 + ncomp]]
    off = 4 + 2 * ncomp
    xyz = np.array([float(t) for t in tok[off:off + 3 * nv]]).reshape(nv, 3)
    comps, i = [], 0
    for c in counts:
        n = abs(c)
        assert c < 0, 'only closed components supported'
        comps.append(xyz[i:i + n]); i += n
    return comps


def minrad(P):
    a = P - np.roll(P, 1, axis=0)          # incoming edge at each vertex
    b = np.roll(P, -1, axis=0) - P         # outgoing edge
    la, lb = np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1)
    cosang = np.clip((a * b).sum(1) / (la * lb), -1.0, 1.0)
    theta = np.arccos(cosang)
    with np.errstate(divide='ignore'):
        r = np.minimum(la, lb) / (2.0 * np.tan(theta / 2.0))
    r[theta == 0] = np.inf
    return r.min()


def seg_closest(p0, p1, q0, q1):
    """Vectorised closest points between segments p0p1 and q0q1; returns s, t in [0,1]."""
    d1, d2, r = p1 - p0, q1 - q0, p0 - q0
    a, e = (d1 * d1).sum(-1), (d2 * d2).sum(-1)
    b, c, f = (d1 * d2).sum(-1), (d1 * r).sum(-1), (d2 * r).sum(-1)
    den = a * e - b * b
    s = np.where(den > 1e-30, np.clip((b * f - c * e) / np.where(den > 1e-30, den, 1), 0, 1), 0.0)
    t = (b * s + f) / e
    # clamp t, recompute s
    t1 = np.clip(t, 0, 1)
    s = np.where(t != t1, np.clip((b * t1 - c) / a, 0, 1), s)
    t = np.clip((b * s + f) / e, 0, 1)
    # (near-)parallel segments: the closest set is an interval; take the midpoint of the overlap
    # of the projections, which is the perpendicular (interior-interior) pair when one exists.
    par = den <= 1e-12 * a * e
    if par.any():
        u = d1[par] / np.sqrt(a[par])[:, None]                 # unit direction of p
        lp = np.sqrt(a[par])
        q0p = ((q0[par] - p0[par]) * u).sum(-1)                # q endpoints projected on p (arc length)
        q1p = ((q1[par] - p0[par]) * u).sum(-1)
        lo = np.maximum(0, np.minimum(q0p, q1p)); hi = np.minimum(lp, np.maximum(q0p, q1p))
        has = lo <= hi
        sm = np.where(has, (lo + hi) / 2 / lp, s[par])
        tm = np.clip((b[par] * sm + f[par]) / e[par], 0, 1)
        s[par] = sm; t[par] = np.where(has, tm, t[par])
    return s, t


def critical_at(P, i, s, w, tol):
    """Is direction w (from the point to the partner) in the normal cone at the point
    P[i] + s (P[i+1]-P[i])?  Interior: w . edge = 0.  At a vertex: w . e_in >= 0 >= w . e_out."""
    n = len(P)
    e = P[(i + 1) % n] - P[i]
    wn = w / np.linalg.norm(w, axis=-1, keepdims=True)
    en = e / np.linalg.norm(e, axis=-1, keepdims=True)
    ok = np.abs((wn * en).sum(-1)) <= tol                     # interior perpendicular
    # vertex at s == 0: vertex i, edges e_in = P[i]-P[i-1], e_out = e
    ein = P[i] - P[(i - 1) % n]
    ein = ein / np.linalg.norm(ein, axis=-1, keepdims=True)
    v0 = ((wn * ein).sum(-1) >= -tol) & ((wn * en).sum(-1) <= tol)
    # vertex at s == 1: vertex i+1, e_in = e, e_out = P[i+2]-P[i+1]
    eout = P[(i + 2) % n] - P[(i + 1) % n]
    eout = eout / np.linalg.norm(eout, axis=-1, keepdims=True)
    v1 = ((wn * en).sum(-1) >= -tol) & ((wn * eout).sum(-1) <= tol)
    return np.where(s <= 1e-12, v0, np.where(s >= 1 - 1e-12, v1, ok))


def dcsd(comps, tol=1e-6):
    best = np.inf
    for ci, P in enumerate(comps):
        for cj, Q in enumerate(comps):
            if cj < ci:
                continue
            n, m = len(P), len(Q)
            I, J = np.meshgrid(np.arange(n), np.arange(m), indexing='ij')
            I, J = I.ravel(), J.ravel()
            if ci == cj:
                keep = (J > I) & ((J - I) % n > 1) & ((I - J) % n > 1)   # non-adjacent, each pair once
                I, J = I[keep], J[keep]
            for k in range(0, len(I), 2_000_000):
                ii, jj = I[k:k + 2_000_000], J[k:k + 2_000_000]
                s, t = seg_closest(P[ii], P[(ii + 1) % n], Q[jj], Q[(jj + 1) % m])
                x = P[ii] + s[:, None] * (P[(ii + 1) % n] - P[ii])
                y = Q[jj] + t[:, None] * (Q[(jj + 1) % m] - Q[jj])
                w = y - x
                d = np.linalg.norm(w, axis=1)
                good = d > 1e-12
                crit = critical_at(P, ii, s, w, tol) & critical_at(Q, jj, t, -w, tol) & good
                if crit.any():
                    best = min(best, d[crit].min())
    return best


def report(path):
    comps = read_vect(path)
    L = sum(np.linalg.norm(np.roll(P, -1, 0) - P, axis=1).sum() for P in comps)
    mr = min(minrad(P) for P in comps)
    dc = dcsd(comps)
    thick = min(mr, dc / 2)
    nv = sum(len(P) for P in comps)
    print(f'{path}: components {len(comps)} vertices {nv} length {L:.6f} minRad {mr:.6f} '
          f'dcsd/2 {dc / 2:.6f} thickness {thick:.6f} ropelength {L / thick:.6f}')
    return L / thick


if __name__ == '__main__':
    for p in sys.argv[1:]:
        report(p)
