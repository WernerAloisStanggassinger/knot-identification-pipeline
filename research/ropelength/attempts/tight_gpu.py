#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""tight_gpu.py — SONO-style GPU tightener for a closed one-component polygon (unit tube RADIUS).

Per iteration:
  1. shrink   x_i += a * (x_{i-1} + x_{i+1} - 2 x_i)                (curve shortening)
  2. K sweeps of constraint projection (Jacobi):
       non-local SEGMENT pairs (arc separation >= ARC_MIN = pi * 1.1) closer than D = 2 (1 + eps)
       are pushed apart along their closest-point direction; the push is split onto the four
       endpoints with the barycentric weights (1-s, s), (1-t, t);
       vertices with Rawdon minRad = h / (2 tan(theta/2)) < 1 + eps are pulled toward the neighbour midpoint.
  3. every R iterations: periodic re-equilateralisation (linear arc-length resampling on the GPU;
     periodic cubic spline on the CPU every R_SPLINE iterations).
Scale is fixed by the constraints (thickness ~ 1), so length ~ ropelength.  The in-run estimate
  est = L / min(minRad, dmin / 2)   uses the minimum over ALL non-local segment pairs (>= dcsd-based thickness
  is NOT assumed; dmin <= dcsd, so est >= the Rawdon value up to the arc cutoff).  The final number is only the
  one printed by ropelength_check.py, never est.  Knot type: polygon_to_pd.py + Alexander, after the run.

Writes on arrival: <out>.jsonl (one line per log step, flushed) and <out> (checkpoint VECT, best est so far).
Usage: python -I tight_gpu.py in.vect out.vect [--nv N] [--minutes M] [--shrink a] [--eps e] [--sweeps K]
"""
import argparse, json, os, sys, time
import numpy as np
import cupy as cp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # ropelength_check.py one folder up
from ropelength_check import read_vect

ap = argparse.ArgumentParser()
ap.add_argument('inp'); ap.add_argument('out')
ap.add_argument('--nv', type=int, default=0, help='resample to N vertices at start (0 = keep)')
ap.add_argument('--minutes', type=float, default=30.0)
ap.add_argument('--shrink', type=float, default=0.05)
ap.add_argument('--eps', type=float, default=0.002)
ap.add_argument('--sweeps', type=int, default=4)
ap.add_argument('--omega', type=float, default=0.8)
ap.add_argument('--resample-every', type=int, default=20)
ap.add_argument('--spline-every', type=int, default=500)
ap.add_argument('--log-every', type=int, default=200)
ap.add_argument('--max-iter', type=int, default=10**9)
ap.add_argument('--block', type=int, default=1024, help='row block for the pair computation')
A = ap.parse_args()

ARC_MIN = np.pi * 1.1
D = 2.0 * (1.0 + A.eps)
RMIN = 1.0 + A.eps


def resample_linear(P, n):
    Q = cp.concatenate([P, P[:1]])
    s = cp.concatenate([cp.zeros(1), cp.cumsum(cp.linalg.norm(cp.diff(Q, axis=0), axis=1))])
    u = cp.linspace(0, float(s[-1]), n, endpoint=False)
    return cp.stack([cp.interp(u, s, Q[:, k]) for k in range(3)], axis=1)


def resample_spline(P, n):
    from scipy.interpolate import CubicSpline
    Pc = cp.asnumpy(P).astype(np.float64)
    Q = np.vstack([Pc, Pc[:1]])
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(Q, axis=0), axis=1))]
    u = np.linspace(0, s[-1], n, endpoint=False)
    return cp.asarray(CubicSpline(s, Q, bc_type='periodic')(u))


def minrad_vec(P):
    a = P - cp.roll(P, 1, axis=0); b = cp.roll(P, -1, axis=0) - P
    la, lb = cp.linalg.norm(a, axis=1), cp.linalg.norm(b, axis=1)
    th = cp.arccos(cp.clip((a * b).sum(1) / (la * lb), -1.0, 1.0))
    r = cp.minimum(la, lb) / (2.0 * cp.tan(cp.maximum(th, 1e-12) / 2.0))
    return r


def seg_pairs(P, i0, i1, kmin):
    """Closest points between segments i in [i0,i1) and all segments j, non-local only."""
    n = len(P)
    I = cp.arange(i0, i1)[:, None]; J = cp.arange(n)[None, :]
    sep = cp.abs(I - J); sep = cp.minimum(sep, n - sep)
    mask = sep >= kmin
    p0 = P[i0:i1][:, None, :]; p1 = P[(cp.arange(i0, i1) + 1) % n][:, None, :]
    q0 = P[None, :, :]; q1 = cp.roll(P, -1, axis=0)[None, :, :]
    d1, d2, r = p1 - p0, q1 - q0, p0 - q0
    a = (d1 * d1).sum(-1); e = (d2 * d2).sum(-1)
    b = (d1 * d2).sum(-1); c = (d1 * r).sum(-1); f = (d2 * r).sum(-1)
    den = a * e - b * b
    s = cp.where(den > 1e-14 * a * e, cp.clip((b * f - c * e) / cp.where(den > 0, den, 1), 0, 1), 0.5)
    t = cp.clip((b * s + f) / e, 0, 1)
    s = cp.clip((b * t - c) / a, 0, 1)
    t = cp.clip((b * s + f) / e, 0, 1)
    x = p0 + s[..., None] * d1; y = q0 + t[..., None] * d2
    w = x - y; d = cp.linalg.norm(w, axis=-1)
    d = cp.where(mask, d, cp.inf)
    return s, t, w, d


def project_pairs(P, kmin):
    n = len(P)
    disp = cp.zeros_like(P); cnt = cp.zeros(n)
    dmin = cp.inf
    for i0 in range(0, n, A.block):
        i1 = min(n, i0 + A.block)
        s, t, w, d = seg_pairs(P, i0, i1, kmin)
        dmin = cp.minimum(dmin, d.min())
        ov = cp.maximum(D - d, 0.0)
        hit = ov > 0
        if not bool(hit.any()):
            continue
        ii, jj = cp.nonzero(hit)
        dd = d[ii, jj]; u = w[ii, jj] / dd[:, None]
        push = 0.5 * ov[ii, jj][:, None] * u          # each side moves half the overlap
        si, tj = s[ii, jj][:, None], t[ii, jj][:, None]
        gi = ii + i0
        for idx, wt, sg in ((gi, 1 - si, 1), ((gi + 1) % n, si, 1), (jj, 1 - tj, -1), ((jj + 1) % n, tj, -1)):
            cp.add.at(disp, (idx[:, None], cp.arange(3)[None, :]), sg * wt * push)
            cp.add.at(cnt, idx, wt[:, 0])
    # each pair is seen twice (i,j) and (j,i) -> halve; normalise crowded vertices
    disp = 0.5 * disp / cp.maximum(1.0, 0.5 * cnt)[:, None]
    return disp, float(dmin)


def project_curv(P):
    r = minrad_vec(P)
    bad = r < RMIN
    if not bool(bad.any()):
        return cp.zeros_like(P), float(r.min())
    mid = 0.5 * (cp.roll(P, 1, 0) + cp.roll(P, -1, 0))
    frac = cp.clip(1.0 - r / RMIN, 0, 1)[:, None] * bad[:, None]
    return frac * (mid - P), float(r.min())


def length(P):
    return float(cp.linalg.norm(cp.roll(P, -1, 0) - P, axis=1).sum())


def write_vect(path, P):
    Pc = cp.asnumpy(P)
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        f.write(f'VECT\n1 {len(Pc)} 0\n-{len(Pc)}\n0\n')
        for x, y, z in Pc:
            f.write(f'{x:.12f} {y:.12f} {z:.12f}\n')
    os.replace(tmp, path)


comps = read_vect(A.inp)
assert len(comps) == 1, 'one component only'
P = cp.asarray(comps[0], dtype=cp.float64)
N = A.nv or len(P)
P = resample_spline(P, N)
# scale so that the in-run thickness estimate is 1 (constraints just satisfied)
h = length(P) / N
kmin = int(np.ceil(ARC_MIN / h))
_, dmin = project_pairs(P, kmin); mr = float(minrad_vec(P).min())
th0 = min(mr, dmin / 2)
P = P / th0
L0 = length(P)
print(f'start {A.inp}: N {N} L/th_est {L0:.4f} (minRad {mr/th0:.4f} dmin/2 {dmin/2/th0:.4f})', flush=True)

logf = open(A.out + '.jsonl', 'a')
best = np.inf; t0 = time.time(); it = 0
while it < A.max_iter and time.time() - t0 < 60 * A.minutes:
    it += 1
    P = P + A.shrink * (cp.roll(P, 1, 0) + cp.roll(P, -1, 0) - 2 * P)
    h = length(P) / N
    kmin = int(np.ceil(ARC_MIN / h))
    for _ in range(A.sweeps):
        dp, dmin = project_pairs(P, kmin)
        P = P + A.omega * dp
        dc, mr = project_curv(P)
        P = P + A.omega * dc
    if it % A.resample_every == 0:
        P = resample_spline(P, N) if it % A.spline_every == 0 else resample_linear(P, N)
    if it % A.log_every == 0:
        _, dmin = project_pairs(P, kmin); mr = float(minrad_vec(P).min())
        L = length(P); th = min(mr, dmin / 2); est = L / th
        rec = dict(it=it, sec=round(time.time() - t0, 1), L=L, minRad=mr, dmin_half=dmin / 2, est=est, h=L / N)
        logf.write(json.dumps(rec) + '\n'); logf.flush()
        print(rec, flush=True)
        if est < best:
            best = est; write_vect(A.out, P)
        if not np.isfinite(est) or th < 0.5:
            print('ABORT: thickness collapsed (possible strand passage)', flush=True); break
print(f'done it {it} best est {best:.4f} -> {A.out} (certify with ropelength_check.py)', flush=True)
