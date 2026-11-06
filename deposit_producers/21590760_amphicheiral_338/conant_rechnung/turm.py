#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Quotient tower after Conant §5 — fetch nu_k from the curve.

WHAT THE TEXT SAYS (2607.18655, l. 363-371):
  (S^3, K_a, g_a) = (S^3, K, g);  (S^3, K_{k-1}, g_{k-1}) = quotient of (S^3, K_k, g_k)
  to r_k := g_k^{2^{k-1}}  — that is the pi ROTATION about the axis C_2, an involution.
  nu_k := Nabla_{K_k}(w),  nu_a = Nabla_K,  nu_1 = phi^2 (Lemma 5.3, Hartley-Kawauchi).

IN COORDINATES: g_a acts on the parameter as t -> t + 2pi/2^a, hence
  r_a = g_a^{2^{a-1}}:  t -> t + pi   and in space  theta -> theta + pi.
Quotient then: fundamental domain t in [0,pi), quotient map (r,theta,z) -> (r, 2theta, z).
One step halves the tower height. a=2: one step. a=3: two.

THE CHECKS that make this worth something in the first place:
  (P1) nu_1 must be a SQUARE in Z[w]            (Lemma 5.3 — can fail)
  (P2) nu_k divides nu_{k+1} EXACTLY in Z[w], R_k(0)=1   (Thm 7.1, INTEGRAL — can fail)
  (P3) Nabla_K = phi^2 * prod R_k                    (Thm 1.3)
"""
import site_config
import sys, os, json, importlib.util
import numpy as np

HIER = os.path.dirname(os.path.abspath(__file__))
SIEB = (site_config.ROOT + '/konversationen/kawauchi_20260725/sieb')
sys.path.insert(0, SIEB)
sp_ = importlib.util.spec_from_file_location('a', os.path.join(SIEB, 'stufe_a_pd.py'))
A = importlib.util.module_from_spec(sp_); sp_.loader.exec_module(A)
sp2 = importlib.util.spec_from_file_location('f', os.path.join(SIEB, 'forge_c.py'))
F = importlib.util.module_from_spec(sp2); sp2.loader.exec_module(F)


def quotient(P):
    """One tower step: take half the period, double the angle."""
    n = len(P) // 2
    Q = P[:n]
    r = np.hypot(Q[:, 0], Q[:, 1])
    th = np.arctan2(Q[:, 1], Q[:, 0])
    return np.column_stack([r * np.cos(2 * th), r * np.sin(2 * th), Q[:, 2]])


def alexander_von(P, NS=1440):
    import regina
    idx = np.linspace(0, len(P), NS, endpoint=False).astype(int)
    sg = A.signed_gauss(P[idx])
    if sg is None:
        return None, 'tangenz'
    L = regina.Link.fromSignedGauss(sg)
    L.simplify()
    if L.size() == 0:
        return '1', 'unknoten'
    return str(L.alexander()), f'c={L.size()}'


if __name__ == '__main__':
    import time
    from mod4_test import to_w
    Q = (site_config.ROOT + '/PRIORITAET_DEPOSIT_20260729/counterexamples_data.json')
    K = json.load(open(Q))['knoten']
    wahl = sys.argv[1] if len(sys.argv) > 1 else None
    proben = [k for k in K if k['id'] == wahl] if wahl else \
             [next(k for k in K if k['a_tower'] == 2), next(k for k in K if k['a_tower'] == 3)]
    for k in proben:
        g = k['erzeugung']; a = k['a_tower']
        t0 = time.time()
        P = F.forge_c(g['k'], g['seed'], g['n_xy'], g['n_z'], g['decay'], c=g['c_rot'])
        print(f"\n=== {k['id']}  a={a}  l={k['achsen_verschlingung_ell']}  c={g['c_rot']} ===")
        print(f"  nu_{a} = Nabla_K (from the deposit): {k['conway_nabla'][:60]}")
        cur = P
        for schritt in range(a - 1):
            cur = quotient(cur)
            alex, info = alexander_von(cur)
            print(f"  stage {a-1-schritt}: {info:14s} alex = {str(alex)[:60]}")
        print(f"  total time: {time.time()-t0:.1f} s")


# ---------------------------------------------------------------------------
# Delta(t) -> Nabla(w):  z = t^{1/2} - t^{-1/2},  w = z^2 = t - 2 + 1/t.
# Expand palindromic Delta in u = t + 1/t (t^j + t^-j = Chebyshev in u),
# then w = u - 2.  Normalization via Nabla(0) = 1 (sign/unit drops out).
# ---------------------------------------------------------------------------
def regina_poly(s):
    """regina writes '17 x^6 - 147 x^5 + ... + 17' — WITHOUT multiplication signs, with spaces.
    Disc #568 already failed on exactly this once (prose parser read t notation).
    Hence a dedicated parser here instead of sympify, with self-test below."""
    import re
    s = s.replace('-', '+-').replace('e+-', 'e-')
    ko = {}
    for term in [x.replace(' ', '') for x in s.split('+') if x.strip()]:
        m = re.fullmatch(r'(-?\d*)x\^(\d+)', term)
        if m:
            c = m.group(1); ko[int(m.group(2))] = int(c if c not in ('', '-') else c + '1')
            continue
        m = re.fullmatch(r'(-?\d*)x', term)
        if m:
            c = m.group(1); ko[1] = int(c if c not in ('', '-') else c + '1')
            continue
        m = re.fullmatch(r'(-?\d+)', term)
        if m:
            ko[0] = int(m.group(1)); continue
        return None
    if not ko:
        return None
    d = max(ko)
    return [ko.get(i, 0) for i in range(d + 1)]


def delta_zu_nabla(alex_str):
    """Delta(t) -> Nabla(w). z = t^(1/2)-t^(-1/2), w = z^2 = t-2+1/t.
    Expand palindromic Delta in u = t+1/t, then w = u-2. Normalized via Nabla(0)=1."""
    import sympy as sp
    c = regina_poly(alex_str)
    if c is None:
        return None
    d = len(c) - 1
    if d % 2:
        return None
    n = d // 2
    for i in range(n + 1):
        if c[n + i] != c[n - i]:
            return None
    u, wv = sp.symbols('u w')
    S = [sp.Integer(2), u]
    for j in range(1, n + 1):
        S.append(sp.expand(u * S[j] - S[j - 1]))
    ausdruck = sp.Integer(c[n])
    for j in range(1, n + 1):
        ausdruck += sp.Integer(c[n + j]) * S[j]
    q = sp.Poly(sp.expand(ausdruck.subs(u, wv + 2)), wv)
    ko = [int(v) for v in q.all_coeffs()[::-1]]
    if not ko or ko[0] == 0:
        return None
    if ko[0] < 0:
        ko = [-v for v in ko]
    return ko if ko[0] == 1 else None


def teilt_exakt(a, b):
    """b / a in Z[w] exactly? returns quotient or None (a,b ascending)."""
    a = a[:]; b = b[:]
    while len(a) > 1 and a[-1] == 0: a.pop()
    while len(b) > 1 and b[-1] == 0: b.pop()
    if len(a) > len(b): return None
    q = [0] * (len(b) - len(a) + 1)
    r = b[:]
    for i in range(len(r) - 1, len(a) - 2, -1):
        k = i - (len(a) - 1)
        if r[i] % a[-1] != 0:
            return None
        f = r[i] // a[-1]
        q[k] = f
        for j in range(len(a)):
            r[k + j] -= f * a[j]
    if any(r):
        return None
    return q
