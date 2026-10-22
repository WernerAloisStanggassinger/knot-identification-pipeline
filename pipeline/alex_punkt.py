# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""Alexander polynomial from the PD — computed AT A POINT, not symbolically.

WHY (26.07.2026, after "sage is dead"):
  Sage computes the Alexander determinant symbolically over Z[t,t^-1]. On the smallest
  of the 102 open candidates (c=27, 36 crossings) `alexander_polynomial()` ran
  for over 11 minutes without a result, RSS 739 MB and rising; the set goes up to c=135.
  It was not the time limit that was wrong, it is the method.
  House doctrine (Footgun 4 · Disc #248 · project-theta-throughput-by-interpolation ·
  feedback-point-eval-beats-symbolic-form): mass sweeps in integers AT A POINT.

METHOD
  1. Build the Alexander matrix A(t) from the PD (one row per crossing, one column
     per arc; one row and one column deleted).
  2. Evaluate det A(t) at enough sample points t — in F_p, i.e. pure machine
     arithmetic (Gauss mod p, int64, p < 2^31 so that a*b fits in int64).
  3. Lagrange interpolation mod p gives the coefficients mod p.
  4. Several primes + CRT + symmetric lift -> the integer coefficients.
  No polynomial ring, no term explosion. Cost O(n^3) per (prime, sample point).

⚠️ THE CONVENTION IS NOT GUESSED BUT CALIBRATED.
  For the sign/role assignment at a crossing there are several common
  conventions. Instead of asserting one from memory, ALL variants are
  tried and calibrated on known cases (trefoil t^2-t+1, figure-eight knot
  -t^2+3t-1, plus the 62 candidates on which regina AND Sage agree).
  Only when a variant passes the calibration is it let loose on the 96 open ones.
  Comparison always up to a unit +-t^k (Footgun 12) — never by absolute equality.
"""
from __future__ import annotations

import itertools
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PRIMES = [2147483647, 2147483629, 2147483587, 2147483579]      # < 2^31, a*b fits in int64


# ---------------------------------------------------------------- read PD
def parse_pd(s: str):
    """'X[1,5,2,4] X[3,1,4,6] …' -> [[1,5,2,4], …]"""
    return [[int(x) for x in m.split(',')]
            for m in re.findall(r'X\[([-\d,\s]+)\]', s.replace(' ', ''))]


# ------------------------------------------------------- determinant mod p
def det_mod_p(M: np.ndarray, p: int) -> int:
    """Gaussian elimination in F_p. M is copied. Returns det mod p."""
    A = (M % p).astype(np.int64)
    n = A.shape[0]
    det = 1
    for i in range(n):
        piv = -1
        for r in range(i, n):
            if A[r, i] % p:
                piv = r
                break
        if piv < 0:
            return 0
        if piv != i:
            A[[i, piv]] = A[[piv, i]]
            det = (-det) % p
        a = int(A[i, i]) % p
        det = det * a % p
        inv = pow(a, p - 2, p)
        A[i, i:] = A[i, i:] * inv % p
        col = A[i + 1:, i].copy()
        nz = np.nonzero(col)[0]
        if nz.size:
            A[i + 1 + nz, i:] = (A[i + 1 + nz, i:] - np.outer(col[nz], A[i, i:])) % p
    return det % p


def lagrange_mod_p(xs, ys, p: int):
    """Coefficients of the interpolating polynomial mod p (ascending)."""
    n = len(xs)
    coeff = [0] * n
    for i in range(n):
        # basis polynomial prod_{j!=i} (x - xj) / (xi - xj)
        num = [1] + [0] * (n - 1)
        deg = 0
        den = 1
        for j in range(n):
            if j == i:
                continue
            for k in range(deg + 1, 0, -1):
                num[k] = (num[k] + num[k - 1]) % p
                num[k - 1] = num[k - 1] * (-xs[j]) % p
            # loop above: num *= (x - xj)
            deg += 1
            den = den * (xs[i] - xs[j]) % p
        f = ys[i] * pow(den % p, p - 2, p) % p
        for k in range(n):
            coeff[k] = (coeff[k] + f * num[k]) % p
    return coeff


def crt_lift(residues, primes):
    """CRT + symmetric lift -> integer coefficient."""
    M = 1
    x = 0
    for r, p in zip(residues, primes):
        # x = x mod M, r mod p  ->  combined
        g = pow(M % p, p - 2, p)
        t = (r - x) % p * g % p
        x += M * t
        M *= p
    return x - M if x > M // 2 else x


# --------------------------------------------------- Alexander matrix from PD
def boegen(pd):
    """Merge PD EDGES into ARCS (Wirtinger generators).

    🔑 The mistake of the first draft (26.07., caught by the calibration): taking the 2n
    PD edges as columns. A diagram with n crossings has 2n edges,
    but only n ARCS — and only the arcs are the generators of the Wirtinger presentation.
    At X[a,b,c,d] the strand b->d passes OVER: b and d belong to the SAME
    arc and are therefore merged. The strand a->c passes UNDER: there one
    arc ends and a new one begins, a and c stay separate.
    Returns: dict edge -> arc index, and the number of arcs.
    """
    par = {}

    def find(x):
        par.setdefault(x, x)
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            par[rx] = ry

    for a, b, c, d in pd:
        find(a); find(b); find(c); find(d)
        union(b, d)                      # over-strand: one edge, one arc
    wurzeln = sorted({find(e) for e in par})
    idx = {w: i for i, w in enumerate(wurzeln)}
    return {e: idx[find(e)] for e in par}, len(wurzeln)


def vorzeichen(pd):
    """Sign per crossing from the PD.

    In X[a,b,c,d] the under-strand runs a->c, the over-strand connects b and d.
    The edges are numbered consecutively along the knot (1..2n, cyclic), so
    the direction of the over-strand gives the sign: if it runs b->d (d = b+1), it is
    one sign, if it runs d->b (b = d+1), the other.
    Special case wrap-around: 2n and 1 are adjacent.
    """
    N = 2 * len(pd)
    sg = []
    for a, b, c, d in pd:
        if (d - b) % N == 1:
            sg.append(+1)
        elif (b - d) % N == 1:
            sg.append(-1)
        else:
            sg.append(+1)          # degenerate (loop) — the calibration must catch it
    return sg


def alexander_matrix(pd, t: int, variante: int, p: int):
    """One row per crossing, one column per ARC; last row + column deleted.

    `variante` selects the role assignment at the crossing. Which one is right is decided
    by the calibration — NOTHING is asserted here.
    """
    n = len(pd)
    idx, m = boegen(pd)
    sg = vorzeichen(pd)
    A = np.zeros((n, m), dtype=np.int64)
    for row, (a, b, c, d) in enumerate(pd):
        over = idx[b]                       # == idx[d], merged
        u_in, u_out = idx[a], idx[c]
        # 🔑 Second mistake caught by the calibration (26.07.): uniform weights.
        # On the trefoil this goes unnoticed — it is alternating, all crossings
        # have the same sign. The figure-eight knot has 2 positive and 2 negative ones
        # and thus came out as t^2-t+1 (det 3) instead of t^2-3t+1 (det 5).
        # At a NEGATIVE crossing the roles of incoming and outgoing swap.
        if sg[row] < 0 and (variante & 2):
            u_in, u_out = u_out, u_in
        if variante & 4:
            wo, wi, wu = (1 - t) % p, t % p, (-1) % p
        else:
            wo, wi, wu = (t - 1) % p, t % p, (-1) % p
        if variante & 1:
            wi, wu = wu, wi
        A[row, over] = (A[row, over] + wo) % p
        A[row, u_in] = (A[row, u_in] + wi) % p
        A[row, u_out] = (A[row, u_out] + wu) % p
    return A[:-1, :-1]


def alexander_at(pd, t: int, variante: int, p: int) -> int:
    return det_mod_p(alexander_matrix(pd, t, variante, p), p)


def alexander_poly(pd, variante: int):
    """Integer coefficients (ascending) of det A(t), up to a unit."""
    n = len(pd)
    deg = n + 2                      # safe bound for the degree
    xs = list(range(1, deg + 2))     # sample points
    per_prime = []
    for p in PRIMES:
        ys = [alexander_at(pd, x, variante, p) for x in xs]
        per_prime.append(lagrange_mod_p([x % p for x in xs], ys, p))
    L = len(per_prime[0])
    return [crt_lift([per_prime[i][k] for i in range(len(PRIMES))], PRIMES) for k in range(L)]


# ------------------------------------------------------------ comparison
def normiere(co):
    """Make comparable up to a unit +-t^k (Footgun 12: NEVER compare absolutely)."""
    c = list(co)
    while c and c[0] == 0:
        c.pop(0)
    while c and c[-1] == 0:
        c.pop()
    if not c:
        return ()
    if c[0] < 0 or (c[0] == 0 and sum(c) < 0):
        c = [-x for x in c]
    return tuple(c) if c[0] > 0 else tuple(-x for x in c)


def parse_regina(s: str):
    """regina writes 'x^18 - 5 x^17 + … + 1' — convert to coefficients (ascending)."""
    s = s.replace('-', '+-').replace(' ', '')
    terms = [t for t in s.split('+') if t]
    d = {}
    for t in terms:
        neg = t.startswith('-')
        t = t.lstrip('-')
        if 'x' in t:
            co, _, ex = t.partition('x')
            co = int(co) if co else 1
            e = int(ex[1:]) if ex.startswith('^') else 1
        else:
            co, e = int(t), 0
        d[e] = d.get(e, 0) + (-co if neg else co)
    if not d:
        return ()
    return tuple(d.get(i, 0) for i in range(max(d) + 1))


# ------------------------------------------------------------------ calibration
TREFOIL = parse_pd('X[1,4,2,5] X[3,6,4,1] X[5,2,6,3]')
FIG8 = parse_pd('X[4,2,5,1] X[8,6,1,5] X[6,3,7,2] X[2,7,3,8]')
# ⚠️ The FIG8 PD is a candidate, not a dogma — the calibration checks it too.

if __name__ == '__main__':
    print('=== CALIBRATION: which crossing convention reproduces known results? ===', flush=True)
    soll = {'Kleeblatt': (TREFOIL, normiere(parse_regina('x^2 - x + 1'))),
            'Achterknoten': (FIG8, normiere(parse_regina('x^2 - 3 x + 1')))}
    gut = []
    for v in range(8):
        zeile, ok = [], True
        for name, (pd, s_) in soll.items():
            try:
                got = normiere(alexander_poly(pd, v))
            except Exception as e:
                got = f'Fehler {type(e).__name__}'
            zeile.append(f'{name} {got}')
            ok &= (got == s_)
        print(f'  variant {v}: ' + ' · '.join(zeile) + ('   ✓✓' if ok else ''), flush=True)
        if ok: gut.append(v)
    print(f'\npassed: {gut}', flush=True)
    print(f'arcs in the trefoil: {boegen(TREFOIL)[1]} (expected 3) · '
          f'in the figure-eight knot: {boegen(FIG8)[1]} (expected 4)', flush=True)
