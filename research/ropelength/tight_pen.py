#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""tight_pen.py — penalty / continuation tightener on the GPU (cupy), successor of tight_gpu.py.

tight_gpu.py (SONO-style shrink + Jacobi projection) stopped at 38.06 on a loose trefoil (ridgerunner reaches
32.80 from its own start at 165 vertices).  Here the length is minimised directly:

    E(X) = L(X) + mu * [ sum_{non-local segment pairs} relu(2 - d_ij)^2  +  sum_v relu(c_max - cos theta_v)^2 ]

  d_ij      closest distance of segments i, j with arc separation >= pi + 2h (chord of a unit circle arc of length pi is 2);
            gradient by the envelope rule: d d / d p0 = (1-s) u, d d / d p1 = s u, (and -(1-t) u, -t u on the partner)
  theta_v   turning angle; Rawdon minRad >= 1 with edge h  <=>  theta <= 2 atan(h/2)  =>  c_max = cos(2 atan(h/2))
  mu        raised geometrically per stage (continuation). Steps: --opt gd = backtracking gradient descent (E must drop; per-vertex
            step cap lr*h; tangential gradient part removed), --opt lbfgs = SciPy L-BFGS-B on the exact gradient; periodic cubic-spline
            re-equilateralisation in both. (An Adam variant produced kinks and was dropped.)
Every gradient is analytic; --gradcheck compares it with central finite differences (counter-case before use).
In-run est = L / min(minRad, dmin/2) with dmin over ALL non-local pairs (conservative); the number that counts is only
ropelength_check.py.  Knot type after the run: polygon_to_pd.py + Alexander.
Writes on arrival: <out>.jsonl (flushed per log line) and <out> (best-est checkpoint).

Usage: python -I tight_pen.py in.vect out.vect [--nv N] [--minutes M] [--stages S] [--steps K] [--mu0 m] [--mu-fac f]
       python -I tight_pen.py in.vect - --gradcheck
"""
import argparse, json, os, sys, time
import numpy as np
import cupy as cp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ropelength_check import read_vect

ap = argparse.ArgumentParser()
ap.add_argument('inp'); ap.add_argument('out')
ap.add_argument('--nv', type=int, default=0)
ap.add_argument('--minutes', type=float, default=30.0)
ap.add_argument('--stages', type=int, default=14)
ap.add_argument('--steps', type=int, default=3000)
ap.add_argument('--mu0', type=float, default=10.0)
ap.add_argument('--mu-fac', type=float, default=2.0)
ap.add_argument('--lr0', type=float, default=0.05, help='--opt gd: per-vertex step cap in units of the mean edge length')
ap.add_argument('--lr-fac', type=float, default=0.85)
ap.add_argument('--resample-every', type=int, default=100)
ap.add_argument('--log-every', type=int, default=100)
ap.add_argument('--block', type=int, default=1024)
ap.add_argument('--gradcheck', action='store_true')
ap.add_argument('--opt', choices=['gd', 'lbfgs'], default='gd')
A = ap.parse_args()


def length(P):
    return cp.linalg.norm(cp.roll(P, -1, 0) - P, axis=1).sum()


def grad_length(P):
    e = cp.roll(P, -1, 0) - P
    u = e / cp.linalg.norm(e, axis=1, keepdims=True)
    return cp.roll(u, 1, 0) - u                     # dL/dx_i = u_{i-1} - u_i


def pair_terms(P, kmin, want_grad=True):
    """sum relu(2-d)^2 over ordered non-local segment pairs (each unordered pair twice), its gradient, and dmin."""
    n = len(P); G = cp.zeros_like(P); S = 0.0; dmin = cp.inf
    Q1 = cp.roll(P, -1, 0)
    for i0 in range(0, n, A.block):
        i1 = min(n, i0 + A.block)
        I = cp.arange(i0, i1)[:, None]; J = cp.arange(n)[None, :]
        sep = cp.abs(I - J); sep = cp.minimum(sep, n - sep)
        p0 = P[i0:i1][:, None, :]; p1 = Q1[i0:i1][:, None, :]
        q0 = P[None]; q1 = Q1[None]
        d1, d2, r = p1 - p0, q1 - q0, p0 - q0
        a = (d1 * d1).sum(-1); e = (d2 * d2).sum(-1)
        b = (d1 * d2).sum(-1); c = (d1 * r).sum(-1); f = (d2 * r).sum(-1)
        den = a * e - b * b
        s = cp.where(den > 1e-14 * a * e, cp.clip((b * f - c * e) / cp.where(den > 0, den, 1), 0, 1), 0.5)
        t = cp.clip((b * s + f) / e, 0, 1)
        s = cp.clip((b * t - c) / a, 0, 1)
        t = cp.clip((b * s + f) / e, 0, 1)
        w = (p0 + s[..., None] * d1) - (q0 + t[..., None] * d2)
        d = cp.linalg.norm(w, axis=-1)
        d = cp.where(sep >= kmin, d, cp.inf)
        dmin = cp.minimum(dmin, d.min())
        ov = cp.maximum(2.0 - d, 0.0)
        S = S + (ov * ov).sum()
        if not want_grad:
            continue
        hit = ov > 0
        if not bool(hit.any()):
            continue
        ii, jj = cp.nonzero(hit)
        g = (-2.0 * ov[ii, jj] / d[ii, jj])[:, None] * w[ii, jj]      # dE/dd * u * d ... = -2 ov * w/d
        si, tj = s[ii, jj][:, None], t[ii, jj][:, None]
        gi = ii + i0
        for idx, wt in ((gi, 1 - si), ((gi + 1) % n, si), (jj, -(1 - tj)), ((jj + 1) % n, -tj)):
            cp.add.at(G, (idx[:, None], cp.arange(3)[None, :]), wt * g)
    return S, G, float(dmin)


def curv_terms(P, cmax, want_grad=True):
    a = P - cp.roll(P, 1, 0); b = cp.roll(P, -1, 0) - P
    la = cp.linalg.norm(a, axis=1); lb = cp.linalg.norm(b, axis=1)
    g = (a * b).sum(1) / (la * lb)
    ov = cp.maximum(cmax - g, 0.0)
    S = (ov * ov).sum()
    if not want_grad:
        return S, None
    dEg = (-2.0 * ov)[:, None]
    dga = b / (la * lb)[:, None] - g[:, None] * a / (la * la)[:, None]
    dgb = a / (la * lb)[:, None] - g[:, None] * b / (lb * lb)[:, None]
    # a_i = x_i - x_{i-1}, b_i = x_{i+1} - x_i
    G = dEg * (dga - dgb) + cp.roll(dEg * dgb, 1, 0) - cp.roll(dEg * dga, -1, 0)
    return S, G


def minrad(P):
    a = P - cp.roll(P, 1, 0); b = cp.roll(P, -1, 0) - P
    la, lb = cp.linalg.norm(a, axis=1), cp.linalg.norm(b, axis=1)
    th = cp.arccos(cp.clip((a * b).sum(1) / (la * lb), -1, 1))
    return float((cp.minimum(la, lb) / (2 * cp.tan(cp.maximum(th, 1e-12) / 2))).min())


def params(P):
    n = len(P); h = float(length(P)) / n
    kmin = int(np.ceil((np.pi + 2 * h) / h))
    cmax = float(np.cos(2 * np.arctan(h / 2)))
    return h, kmin, cmax


def energy_grad(P, mu, kmin, cmax, want_grad=True):
    Sp, Gp, dmin = pair_terms(P, kmin, want_grad)
    Sc, Gc = curv_terms(P, cmax, want_grad)
    E = float(length(P)) + mu * float(Sp + Sc)
    G = grad_length(P) + mu * (Gp + Gc) if want_grad else None
    return E, G, dmin


def resample_spline(P, n):
    from scipy.interpolate import CubicSpline
    Pc = cp.asnumpy(P); Q = np.vstack([Pc, Pc[:1]])
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(Q, axis=0), axis=1))]
    return cp.asarray(CubicSpline(s, Q, bc_type='periodic')(np.linspace(0, s[-1], n, endpoint=False)))


def write_vect(path, P):
    Pc = cp.asnumpy(P); tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        f.write(f'VECT\n1 {len(Pc)} 0\n-{len(Pc)}\n0\n')
        for x, y, z in Pc:
            f.write(f'{x:.12f} {y:.12f} {z:.12f}\n')
    os.replace(tmp, path)


comps = read_vect(A.inp); assert len(comps) == 1
P = cp.asarray(comps[0], dtype=cp.float64)
N = A.nv or len(P)
P = resample_spline(P, N)
h, kmin, cmax = params(P)
_, _, dmin = pair_terms(P, kmin, False)
th0 = min(minrad(P), dmin / 2); P = P / th0          # start with est-thickness 1

if A.gradcheck:
    # counter-case for the analytic gradient: overlap forced by shrinking 8 %, so every term is active
    X = P * 0.92; h, kmin, cmax = params(X); mu = 7.0
    cmax = 0.9999                 # forces the curvature term active at most vertices (gradient test only)
    _, G, _ = energy_grad(X, mu, kmin, cmax)
    rng = np.random.default_rng(3); worst = 0.0
    for _ in range(40):
        i, k = int(rng.integers(N)), int(rng.integers(3)); eps = 1e-6
        Xp = X.copy(); Xp[i, k] += eps; Xm = X.copy(); Xm[i, k] -= eps
        fd = (energy_grad(Xp, mu, kmin, cmax, False)[0] - energy_grad(Xm, mu, kmin, cmax, False)[0]) / (2 * eps)
        an = float(G[i, k]); rel = abs(fd - an) / max(1e-8, abs(fd) + abs(an))
        worst = max(worst, rel)
    pair_active = float(pair_terms(X, kmin, False)[0]) > 0; curv_active = float(curv_terms(X, cmax, False)[0]) > 0
    print(f'gradcheck: 40 coords, worst rel err {worst:.2e}, pair term active {pair_active}, curv term active {curv_active}',
          'PASS' if worst < 1e-4 else 'FAIL')
    sys.exit(0)

logf = open(A.out + '.jsonl', 'a')
print(f'start {A.inp}: N {N} est {float(length(P)):.4f}', flush=True)
best = np.inf; t0 = time.time(); it = 0; stop = False


def log_point(P, it, stage, mu, E):
    global best
    h, kmin, cmax = params(P)
    _, _, dmin = pair_terms(P, kmin, False); mr = minrad(P); L = float(length(P))
    est = L / min(mr, dmin / 2)
    rec = dict(it=it, stage=stage, mu=mu, sec=round(time.time() - t0, 1), E=E, L=L, minRad=mr, dmin_half=dmin / 2, est=est)
    logf.write(json.dumps(rec) + '\n'); logf.flush(); print(rec, flush=True)
    if est < best:
        best = est; write_vect(A.out, P)
    return np.isfinite(est) and min(mr, dmin / 2) >= 0.3


if A.opt == 'lbfgs':
    from scipy.optimize import minimize
    for stage in range(A.stages):
        mu = A.mu0 * A.mu_fac ** stage
        for chunk in range(max(1, A.steps // A.resample_every)):
            if time.time() - t0 > 60 * A.minutes:
                stop = True; break
            h, kmin, cmax = params(P)
            def fg(x):
                X = cp.asarray(x.reshape(-1, 3)); E, G, _ = energy_grad(X, mu, kmin, cmax)
                return E, cp.asnumpy(G).ravel()
            r = minimize(fg, cp.asnumpy(P).ravel(), jac=True, method='L-BFGS-B',
                         options=dict(maxiter=A.resample_every, maxcor=30, ftol=1e-15, gtol=1e-12))
            it += int(r.nit)
            P = resample_spline(cp.asarray(r.x.reshape(-1, 3)), N)
            if not log_point(P, it, stage, mu, float(r.fun)):
                print('ABORT: thickness collapsed (possible strand passage)', flush=True); stop = True; break
        if stop:
            break
    print(f'done it {it} best est {best:.4f} -> {A.out} (certify with ropelength_check.py)', flush=True)
    sys.exit(0)
for stage in range(A.stages):
    mu = A.mu0 * A.mu_fac ** stage; lr_rel = A.lr0 * A.lr_fac ** stage
    eta = None
    for step in range(A.steps):
        if time.time() - t0 > 60 * A.minutes:
            stop = True; break
        it += 1
        h, kmin, cmax = params(P)
        E, G, dmin = energy_grad(P, mu, kmin, cmax)
        T = cp.roll(P, -1, 0) - cp.roll(P, 1, 0); T = T / cp.linalg.norm(T, axis=1, keepdims=True)
        G = G - (G * T).sum(1, keepdims=True) * T          # normal part only (still a descent direction)
        gmax = float(cp.linalg.norm(G, axis=1).max())
        cap = lr_rel * h                                    # no vertex moves more than lr_rel * h per step
        if eta is None: eta = cap / max(gmax, 1e-30)
        for _bt in range(30):                               # backtracking: accept only if E decreases
            e_ = min(eta, cap / max(gmax, 1e-30))
            Pn = P - e_ * G
            En = energy_grad(Pn, mu, kmin, cmax, False)[0]
            if En < E:
                P = Pn; eta = e_ * 1.3; break
            eta = e_ * 0.5
        if it % A.resample_every == 0:
            P = resample_spline(P, N)
        if it % A.log_every == 0:
            h, kmin, cmax = params(P)
            _, _, dmin = pair_terms(P, kmin, False); mr = minrad(P); L = float(length(P))
            est = L / min(mr, dmin / 2)
            rec = dict(it=it, stage=stage, mu=mu, sec=round(time.time() - t0, 1), E=E, L=L, minRad=mr, dmin_half=dmin / 2, est=est)
            logf.write(json.dumps(rec) + '\n'); logf.flush(); print(rec, flush=True)
            if est < best:
                best = est; write_vect(A.out, P)
            if not np.isfinite(est) or min(mr, dmin / 2) < 0.3:
                print('ABORT: thickness collapsed (possible strand passage)', flush=True); stop = True; break
    if stop:
        break
print(f'done it {it} best est {best:.4f} -> {A.out} (certify with ropelength_check.py)', flush=True)
