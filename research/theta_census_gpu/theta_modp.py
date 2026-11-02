# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""theta_modp.py — Bar-Natan–van der Veen θ at a point (t1, t2), exactly modulo a prime p.

Same formula as Bar-Natan and van der Veen's Sage implementation of θ (arXiv:2509.18456, Implementation section; our
numerical version Theta_numeric.sage, Implementation.tex line 78), re-implemented here over GF(p):
three matrices A(t) of size 2n+1, G = A^-1, θ_raw = Σ F1 + ΣΣ F2 + Σ Γ1, Δ_k = t_k^e · det A(t_k).
Instead of over ℚ, computed over GF(p). Reduction mod p is a ring homomorphism as long as no denominator is divisible
by p ⟹ GF(p) value = (ℚ value mod p). Unequal mod p ⟹ the ℚ values are provably unequal.

Two routes, same output:
  ref(Cs, phi, t1, t2, p)          — pure Python, Python int, one knot (reference)
  gpu_batch(liste, t1, t2, p)      — cupy, all knots of equal crossing number n at once
Output per knot: dict(det1, det2, det3, theta_roh, e2) with e2 = 2e = -Σphi - Σs (an integer); theta_roh = θ_raw (key name kept);
  None fields if A(t_k) mod p is singular (then take another p — never silently continue computing).
Env: gnottheorie_env (CuPy). Memory latch: cupy pool limit THETA_GPU_MEM_GB (default 3.5).
"""
import os

P_DEFAULT = 2147483629          # prime < 2^31, products < 2^62 fit into int64


def _inv(a, p):
    a %= p
    if a == 0:
        raise ZeroDivisionError('not invertible mod p')
    return pow(a, p - 2, p)


def frac_mod(num, den, p):
    return num % p * _inv(den, p) % p


# ---------------------------------------------------------------- reference (Python int)
def _A(Cs, n, t, p):
    N = 2 * n + 1
    A = [[1 if r == c else 0 for c in range(N)] for r in range(N)]
    for s, i, j in Cs:
        ts = t if s == 1 else _inv(t, p)
        A[i][i + 1] = (-ts) % p
        A[i][j + 1] = (ts - 1) % p
        A[j][j + 1] = p - 1
    return A


def _det_inv(A, p):
    N = len(A)
    M = [row[:] + [1 if r == c else 0 for c in range(N)] for r, row in enumerate(A)]
    det = 1
    for k in range(N):
        piv = next((r for r in range(k, N) if M[r][k] % p), None)
        if piv is None:
            return 0, None
        if piv != k:
            M[k], M[piv] = M[piv], M[k]
            det = -det
        a = M[k][k] % p
        det = det * a % p
        ai = _inv(a, p)
        M[k] = [x * ai % p for x in M[k]]
        for r in range(N):
            if r != k and M[r][k] % p:
                f = M[r][k]
                M[r] = [(x - f * y) % p for x, y in zip(M[r], M[k])]
    return det % p, [row[N:] for row in M]


def ref(Cs, phi, t1, t2, p=P_DEFAULT):
    n = len(Cs)
    t3 = t1 * t2 % p
    d1, G1 = _det_inv(_A(Cs, n, t1, p), p)
    d2, G2 = _det_inv(_A(Cs, n, t2, p), p)
    d3, G3 = _det_inv(_A(Cs, n, t3, p), p)
    e2 = -sum(phi) - sum(c[0] for c in Cs)
    if G1 is None or G2 is None or G3 is None:
        return dict(det1=d1, det2=d2, det3=d3, theta_roh=None, e2=e2)
    inv2 = _inv(2, p)

    def pw(t, s):
        return t if s == 1 else _inv(t, p)

    th = 0
    for s, i, j in Cs:
        t1s, t2s = pw(t1, s), pw(t2, s)
        t12s = t1s * t2s % p
        t2_2s = t2s * t2s % p
        dn = _inv(t2s - 1, p)
        v = (inv2
             + t2s * G1[i][i] * G2[j][i]
             + (t1s - 1) * t2_2s * G1[j][i] % p * G2[j][i] * dn
             - G1[i][i] * G2[j][j]
             - (t1s - 1) * t2s * G1[j][i] % p * G2[j][j] * dn
             - G3[i][i]
             - (t2s - 1) * G2[j][i] * G3[i][i]
             + 2 * G2[j][j] * G3[i][i]
             + (t12s - 1) * G3[j][i] * dn
             - t2s * (t12s - 1) % p * G1[i][i] * G3[j][i] % p * dn
             - (t1s - 1) * (t2s + 1) % p * (t12s - 1) % p * G1[j][i] % p * G3[j][i] % p * dn
             + (t12s - 1) * G2[i][j] % p * G3[j][i] * dn
             + (t12s - 1) * G2[j][i] % p * G3[j][i]
             + (t2s - 2) * (t12s - 1) % p * G2[j][j] % p * G3[j][i] * dn
             + G1[i][i] * G3[j][j]
             + (t1s - 1) * t2s * G1[j][i] % p * G3[j][j] * dn
             - G2[i][i] * G3[j][j]
             - t2s * G2[j][i] * G3[j][j])
        th = (th + s * v) % p
    for s0, i0, j0 in Cs:
        t1s0, t2s0 = pw(t1, s0), pw(t2, s0)
        for s1, i1, j1 in Cs:
            t1s1, t2s1 = pw(t1, s1), pw(t2, s1)
            t12s1 = t1s1 * t2s1 % p
            dn = _inv(t2s1 - 1, p)
            a = (t1s0 - 1) * (t12s1 - 1) % p
            v = (a * t2s0 % p * G1[j1][i0] % p * G2[i1][i0] % p * G3[j0][i1]
                 - a * G1[j1][i0] % p * G2[i1][j0] % p * G3[j0][i1]
                 - a * t2s0 % p * G1[j1][i0] % p * G2[j1][i0] % p * G3[j0][i1]
                 + a * G1[j1][i0] % p * G2[j1][j0] % p * G3[j0][i1]) % p
            th = (th + s1 * v * dn) % p
    for k, ph in enumerate(phi):
        th = (th + ph * G3[k][k] - ph * inv2) % p
    return dict(det1=d1, det2=d2, det3=d3, theta_roh=th % p, e2=e2)


# ---------------------------------------------------------------- GPU (cupy), bundled per n
def gpu_batch(items, t1, t2, p=P_DEFAULT):
    """items: list of (Cs, phi) with equal n. Returns: list of dicts like ref()."""
    import cupy as cp
    cp.get_default_memory_pool().set_limit(size=int(float(os.environ.get('THETA_GPU_MEM_GB', '3.5')) * 2**30))
    B = len(items)
    n = len(items[0][0])
    assert all(len(c) == n for c, _ in items)
    N = 2 * n + 1
    P = cp.int64(p)
    S = cp.asarray([[c[0] for c in Cs] for Cs, _ in items], dtype=cp.int64)      # (B,n)
    I = cp.asarray([[c[1] for c in Cs] for Cs, _ in items], dtype=cp.int64)
    J = cp.asarray([[c[2] for c in Cs] for Cs, _ in items], dtype=cp.int64)
    PHI = cp.asarray([phi for _, phi in items], dtype=cp.int64)                  # (B,N)
    bidx = cp.arange(B)[:, None]

    def tpow(t, Sx):   # t^s for s=±1, element-wise
        return cp.where(Sx == 1, cp.int64(t % p), cp.int64(_inv(t, p)))

    def build(t):
        A = cp.zeros((B, N, N), dtype=cp.int64)
        A[:, cp.arange(N), cp.arange(N)] = 1
        ts = tpow(t, S)
        # Order as in the reference: crossing by crossing, later ones overwrite earlier ones
        for k in range(n):
            i, j, tk = I[:, k], J[:, k], ts[:, k]
            A[cp.arange(B), i, i + 1] = (-tk) % P
            A[cp.arange(B), i, j + 1] = (tk - 1) % P
            A[cp.arange(B), j, j + 1] = P - 1
        return A

    def modpow(a, e):
        r = cp.ones_like(a)
        b = a % P
        while e:
            if e & 1:
                r = r * b % P
            b = b * b % P
            e >>= 1
        return r

    def det_inv(A):
        M = cp.concatenate([A, cp.broadcast_to(cp.eye(N, dtype=cp.int64), (B, N, N))], axis=2).copy()
        det = cp.ones(B, dtype=cp.int64)
        ok = cp.ones(B, dtype=cp.bool_)
        ar = cp.arange(B)
        for k in range(N):
            col = M[:, k:, k] % P                                    # (B, N-k)
            nz = col != 0
            has = nz.any(axis=1)
            ok &= has
            piv = k + cp.argmax(nz, axis=1)                          # first non-zero row
            swap = piv != k
            rk = M[ar, k, :].copy()
            rp = M[ar, piv, :].copy()
            M[ar, k, :] = rp
            M[ar, piv, :] = rk
            det = cp.where(swap, (P - det) % P, det)
            a = M[ar, k, k] % P
            a = cp.where(has, a, cp.int64(1))
            det = det * a % P
            ai = modpow(a, p - 2)
            M[:, k, :] = M[:, k, :] * ai[:, None] % P
            f = M[:, :, k].copy()                                    # (B,N)
            f[:, k] = 0
            M = (M - (f[:, :, None] * M[:, k, :][:, None, :]) % P) % P
        det = cp.where(ok, det, cp.int64(0))
        return det, M[:, :, N:], ok

    t3 = t1 * t2 % p
    d1, G1, o1 = det_inv(build(t1))
    d2, G2, o2 = det_inv(build(t2))
    d3, G3, o3 = det_inv(build(t3))
    ok = o1 & o2 & o3
    inv2 = cp.int64(_inv(2, p))

    def g(G, r, c):    # G[b, r[b,k], c[b,k]]  -> (B,n) or broadcast
        return G[bidx if r.ndim == 2 else bidx[:, :, None], r, c]

    t1s, t2s = tpow(t1, S), tpow(t2, S)
    t12s = t1s * t2s % P
    t2_2s = t2s * t2s % P
    dn = cp.where(S == 1, cp.int64(_inv(t2 - 1, p)), cp.int64(_inv(_inv(t2, p) - 1, p)))
    G1ii, G1ji = g(G1, I, I), g(G1, J, I)
    G2ji, G2jj, G2ii, G2ij = g(G2, J, I), g(G2, J, J), g(G2, I, I), g(G2, I, J)
    G3ii, G3ji, G3jj = g(G3, I, I), g(G3, J, I), g(G3, J, J)
    m = lambda *xs: _mulmod(xs, P)
    v = (inv2
         + m(t2s, G1ii, G2ji)
         + m(t1s - 1 + P, t2_2s, G1ji, G2ji, dn)
         - m(G1ii, G2jj)
         - m(t1s - 1 + P, t2s, G1ji, G2jj, dn)
         - G3ii
         - m(t2s - 1 + P, G2ji, G3ii)
         + m(2, G2jj, G3ii)
         + m(t12s - 1 + P, G3ji, dn)
         - m(t2s, t12s - 1 + P, G1ii, G3ji, dn)
         - m(t1s - 1 + P, t2s + 1, t12s - 1 + P, G1ji, G3ji, dn)
         + m(t12s - 1 + P, G2ij, G3ji, dn)
         + m(t12s - 1 + P, G2ji, G3ji)
         + m(t2s - 2 + P, t12s - 1 + P, G2jj, G3ji, dn)
         + m(G1ii, G3jj)
         + m(t1s - 1 + P, t2s, G1ji, G3jj, dn)
         - m(G2ii, G3jj)
         - m(t2s, G2ji, G3jj)) % P
    th = (m(S % P, v) % P).sum(axis=1) % P                          # Σ F1

    # Σ F2 over pairs (k0 = index of c0, k1 = index of c1)
    S0, I0, J0 = S[:, :, None], I[:, :, None], J[:, :, None]
    S1, I1, J1 = S[:, None, :], I[:, None, :], J[:, None, :]
    b3 = bidx[:, :, None]
    t1s0, t2s0 = tpow(t1, S0), tpow(t2, S0)
    t12s1 = tpow(t1, S1) * tpow(t2, S1) % P
    dn1 = cp.where(S1 == 1, cp.int64(_inv(t2 - 1, p)), cp.int64(_inv(_inv(t2, p) - 1, p)))
    a = m(t1s0 - 1 + P, t12s1 - 1 + P)
    G1j1i0 = G1[b3, J1, I0]
    G3j0i1 = G3[b3, J0, I1]
    inner = (m(t2s0, G2[b3, I1, I0]) - G2[b3, I1, J0] - m(t2s0, G2[b3, J1, I0]) + G2[b3, J1, J0]) % P
    v2 = m(a, G1j1i0, G3j0i1, inner, dn1, S1 % P)
    th = (th + v2.reshape(B, -1).sum(axis=1)) % P
    # Σ Γ1
    G3kk = G3[bidx, cp.arange(N)[None, :], cp.arange(N)[None, :]]
    th = (th + (m(PHI % P, G3kk) - m(PHI % P, inv2)).sum(axis=1)) % P

    e2 = (-PHI.sum(axis=1) - S.sum(axis=1)).get()
    d1, d2, d3, th, ok = d1.get(), d2.get(), d3.get(), th.get(), ok.get()
    return [dict(det1=int(d1[b]), det2=int(d2[b]), det3=int(d3[b]),
                 theta_roh=int(th[b]) if ok[b] else None, e2=int(e2[b])) for b in range(B)]


def _mulmod(xs, P):
    import cupy as cp
    r = None
    for x in xs:
        x = x % P if not isinstance(x, int) else cp.int64(x % int(P))
        r = x if r is None else r * x % P
    return r
