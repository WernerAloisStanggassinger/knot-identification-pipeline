#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Conant Thm 1.5 / Conjecture 1.1 on our 250 knots — and Thm 7.2 in the form
that does without the intermediate stages nu_k.

WHY THIS WORKS although phi is not in the deposit:
  Modulo 4 a square depends only on the base MOD 2:
      (g + 2g1)^2 = g^2 + 4 g g1 + 4 g1^2  ==  g^2   (mod 4).
  Thus the candidate is UNIQUE and there is nothing to search.

TEST A (Conjecture 1.1 = his main theorem):  Nabla_K = f(z) f(-z)  (mod 4)?
  f = g(w) + z h(w), w = z^2  =>  f(z) f(-z) = g(w)^2 - w h(w)^2.
  Mod 2 the decomposition is forced: even w places -> g, odd -> h.
  So ONE candidate, exactly checkable. The test can fail (see counter-cases).

TEST B (Thm 7.2, last stage): Nabla_K = eps_l [l] * p^2 (mod 4) for some p in Z[w]?
  From nu_{k+1} = eps_l [l] nu_k^2 (mod 4) at k = a-1, with nu_a = Nabla_K.
  p is again uniquely determined mod 2, provided eps_l [l] is invertible mod 2.
  [1]=1, [3]=w+3, [n+2] = (w+2)[n] - [n-2];  [-l] = -[l];  eps_l = +1 if |l|=1 mod 4 else -1.
"""
import site_config
import json, sys
import sympy as sp

w, z = sp.symbols('w z')


def bracket(l):
    """balanced quantum integer [|l|] as a polynomial in w, l odd.

    CONVENTION, calibrated on a single case (29.07.): the absolute value. I had first set [-l] = -[l]
    set AND eps_l taken from |l| — so the sign is in twice, and the
    test broke on exactly the 36 knots with l < 0. Manual check on SK2_1_29_43_07
    (Nabla = 1 - 21 z^2, l = -3): only eps_l[l] = -(w+3) hits Nabla mod 4 = [1,3], so
    eps from |l| and [] from |l|. Footgun 15: the red finding was the checker.
    """
    s = 1
    n = abs(l)
    a, b = sp.Integer(1), sp.Poly(w + 3, w).as_expr()   # [1], [3]
    if n == 1:
        return s * a
    if n == 3:
        return s * b
    cur, prev = b, a
    k = 3
    while k < n:
        cur, prev = sp.expand((w + 2) * cur - prev), cur
        k += 2
    return s * cur


def eps(l):
    return 1 if abs(l) % 4 == 1 else -1


def to_w(nabla_z):
    """Nabla(z) -> coefficient list in w = z^2 (Nabla is even)."""
    p = sp.Poly(sp.expand(sp.sympify(nabla_z)), z)
    c = p.all_coeffs()[::-1]          # ascending by z degree
    for i in range(1, len(c), 2):
        if c[i] != 0:
            raise ValueError('odd z power in Nabla')
    return [int(x) for x in c[0::2]]  # coefficients in w


def mul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    return r


def mod(a, m):
    return [x % m for x in a]


def trim(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def test_A(cw):
    """Nabla = g^2 - w h^2 (mod 4)?  Candidate unique from Nabla mod 2."""
    g = [cw[i] % 2 for i in range(0, len(cw), 2)]     # w^{2i} -> g_i
    h = [cw[i] % 2 for i in range(1, len(cw), 2)]     # w^{2j+1} -> h_j
    lhs = mod(mul(g, g), 4)
    if h:
        wh2 = [0] + mul(h, h)
        n = max(len(lhs), len(wh2))
        lhs = [ (lhs[i] if i < len(lhs) else 0) - (wh2[i] if i < len(wh2) else 0) for i in range(n) ]
        lhs = mod(lhs, 4)
    ziel = mod(cw, 4)
    n = max(len(lhs), len(ziel))
    lhs += [0] * (n - len(lhs)); ziel += [0] * (n - len(ziel))
    return lhs == ziel


def test_B(cw, l):
    """Nabla = eps_l [l] p^2 (mod 4)?  p unique mod 2, provided [l] mod 2 is invertible."""
    bl = sp.Poly(sp.expand(eps(l) * bracket(l)), w).all_coeffs()[::-1]
    bl = [int(x) for x in bl]
    b2 = trim(mod(bl, 2))
    n2 = trim(mod(cw, 2))
    if b2 == [0]:
        return None
    # p^2 mod 2 = p(w^2) (Frobenius) -> solve  b2 * p(w^2) = n2  (mod 2) by polynomial division
    q, r = polydiv2(n2, b2)
    if r != [0]:
        return False
    # q must have the form p(w^2): even exponents only
    if any(q[i] for i in range(1, len(q), 2)):
        return False
    p = q[0::2]
    kand = mod(mul(bl, mul(p, p)), 4)
    ziel = mod(cw, 4)
    n = max(len(kand), len(ziel))
    kand += [0] * (n - len(kand)); ziel += [0] * (n - len(ziel))
    return kand == ziel


def polydiv2(a, b):
    """Division in F_2[w]; a,b ascending, b != 0."""
    a = a[:]; q = [0] * max(1, len(a) - len(b) + 1)
    db = len(b) - 1
    for i in range(len(a) - 1, db - 1, -1):
        if a[i]:
            q[i - db] = 1
            for j in range(len(b)):
                a[i - db + j] ^= b[j]
    return trim(q), trim(a)


if __name__ == '__main__':
    quelle = sys.argv[1] if len(sys.argv) > 1 else \
        (site_config.ROOT + '/PRIORITAET_DEPOSIT_20260729/counterexamples_data.json')
    K = json.load(open(quelle))['knoten']
    stat = {'A_ok': 0, 'A_bad': 0, 'B_ok': 0, 'B_bad': 0, 'B_none': 0}
    bad = []
    per_a = {}
    for k in K:
        cw = to_w(k['conway_nabla'])
        l = k['achsen_verschlingung_ell']; a = k['a_tower']
        A = test_A(cw); B = test_B(cw, l)
        per_a.setdefault(a, [0, 0, 0, 0])
        stat['A_ok' if A else 'A_bad'] += 1
        per_a[a][0 if A else 1] += 1
        if B is None: stat['B_none'] += 1
        else:
            stat['B_ok' if B else 'B_bad'] += 1
            per_a[a][2 if B else 3] += 1
        if not A or B is False:
            bad.append((k['id'], a, l, A, B, k['conway_nabla'][:50]))
    print(f"knots: {len(K)}")
    print(f"TEST A  Nabla = f(z)f(-z) (mod 4)        : {stat['A_ok']}/{len(K)}   breaks {stat['A_bad']}")
    print(f"TEST B  Nabla = eps_l [l] p^2 (mod 4)    : {stat['B_ok']}/{stat['B_ok']+stat['B_bad']}   breaks {stat['B_bad']}  (undecidable {stat['B_none']})")
    for a in sorted(per_a):
        v = per_a[a]
        print(f"   a={a}: A {v[0]}/{v[0]+v[1]}   B {v[2]}/{v[2]+v[3]}")
    for b in bad[:10]:
        print("   BREAK", b)
