#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""M and E per Conant, Corollary 7.3 (case a=2, rotary reflection of order 4).

Statement:  eps_l * [l] * nabla_K = M^2 - 4 E^2   with M,E in Z[w], w = z^2.
We need NO quotient chain: this is a factorization
    (M - 2E)(M + 2E) = eps_l [l] nabla,
so we enumerate divisor pairs A*B and set M=(A+B)/2, E=(B-A)/4 (integrality required).

[l] = balanced quantum number (t^{l/2}-t^{-l/2})/(t^{1/2}-t^{-1/2}), as a polynomial in w.
CALIBRATION against the paper: Conant names [3] == 1+w (mod 2) and [9] == (1+w)(w^3+w+1) (mod 2).
eps_l = +1 if l = 1 (mod 4), else -1  (Lemma 6.13: R_k(0)=eps_l and R_k(0) = l mod 4).

FULL CALIBRATION: the EHW knot. Conant prints M and E explicitly for it (section 8.2).
Our solver must reproduce them - otherwise nothing is computed.
"""
import json, sys
from sympy import Poly, symbols, expand, factor_list, sympify
from itertools import product
w = symbols('w')

def bracket(l):
    """[l] in Z[w]: 1 + sum_{i=1..(l-1)/2} (t^i + t^-i), expressed in w = t+1/t-2."""
    l = abs(l)
    assert l % 2 == 1, 'l muss ungerade sein'  # message: 'l must be odd'
    m = (l-1)//2
    v = Poly(w+2, w)
    p = [Poly(2, w), v]
    for i in range(2, m+1): p.append(Poly(expand(v.as_expr()*p[-1].as_expr()-p[-2].as_expr()), w))
    out = Poly(1, w)
    for i in range(1, m+1): out = Poly(expand(out.as_expr()+p[i].as_expr()), w)
    return out

def eps(l):
    return 1 if abs(l) % 4 == 1 else -1

def loese_ME(P):
    """all (M,E) with M^2-4E^2 = P, from divisor pairs. Return: list (M,E) up to sign."""
    cont, facs = factor_list(P.as_expr(), w)
    facs = [(Poly(p, w), e) for p, e in facs]
    # sign/content as a 'factor' of its own
    loesungen = set()
    ranges = [range(e+1) for _, e in facs]
    for mult in product(*ranges):
        A = Poly(1, w)
        for (p, _), mm in zip(facs, mult): A *= p**mm
        for cA in divisors_int(int(cont)):
            for sA in (1, -1):
                Aa = Poly(expand(sA*cA*A.as_expr()), w)
                if Aa.is_zero: continue
                q, r = divmod_poly(P, Aa)
                if q is None: continue
                B = q
                S = Poly(expand(Aa.as_expr()+B.as_expr()), w)
                D = Poly(expand(B.as_expr()-Aa.as_expr()), w)
                if any(int(c) % 2 for c in S.all_coeffs()): continue
                if any(int(c) % 4 for c in D.all_coeffs()): continue
                M = Poly([int(c)//2 for c in S.all_coeffs()], w)
                E = Poly([int(c)//4 for c in D.all_coeffs()], w)
                if Poly(expand(M.as_expr()**2-4*E.as_expr()**2), w) == P:
                    key = (tuple(M.all_coeffs()), tuple(abs(int(c)) for c in E.all_coeffs()))
                    loesungen.add((tuple(int(c) for c in M.all_coeffs()), tuple(int(c) for c in E.all_coeffs())))
    return loesungen

def divisors_int(n):
    n = abs(n) or 1
    return [d for d in range(1, n+1) if n % d == 0]

def divmod_poly(P, A):
    from sympy import div
    q, r = div(P.as_expr(), A.as_expr(), w)
    qp = Poly(q, w)
    if not Poly(r, w).is_zero: return None, None
    if any(c.q != 1 for c in qp.all_coeffs() if hasattr(c, 'q')): return None, None
    try: return Poly([int(c) for c in qp.all_coeffs()], w), None
    except Exception: return None, None

def eichung():
    print('=== CALIBRATION 1: [l] against the values printed in the paper (mod 2) ===')
    b3, b9 = bracket(3), bracket(9)
    print(f'  [3] = {b3.as_expr()}   mod 2 -> {Poly([int(c)%2 for c in b3.all_coeffs()], w).as_expr()}   (paper: 1+w)')
    print(f'  [9] = {b9.as_expr()}')
    z9 = Poly(expand((1+w)*(w**3+w+1)), w)
    ok9 = [int(c) % 2 for c in b9.all_coeffs()] == [int(c) % 2 for c in z9.all_coeffs()]
    print(f'  [9] mod 2 == (1+w)(w^3+w+1)?  {ok9}   (paper: yes)')
    ok3 = [int(c) % 2 for c in b3.all_coeffs()] == [1, 1]
    print('\n=== CALIBRATION 2: EHW knot, M and E are in the paper (section 8.2) ===')
    nab = Poly(expand((1-w)*(2*w**2-1)**2*(4*w**4+16*w**3+12*w**2-16*w+1)), w)
    M_paper = Poly(expand((2*w**2-1)*(2*w**5+14*w**4+30*w**3+10*w**2-23*w+1)), w)
    E_paper = Poly(expand((2*w**2-1)*(w**5+7*w**4+15*w**3+5*w**2-12*w+1)), w)
    lhs = Poly(expand(eps(3)*bracket(3).as_expr()*nab.as_expr()), w)
    rhs = Poly(expand(M_paper.as_expr()**2-4*E_paper.as_expr()**2), w)
    print(f'  eps_3 = {eps(3)} (paper writes -[3] nabla)   identity -[3]nabla == M^2-4E^2 : {lhs == rhs}')
    sols = loese_ME(lhs)
    getroffen = any(tuple(int(c) for c in M_paper.all_coeffs()) == M or
                    tuple(-int(c) for c in M_paper.all_coeffs()) == M for M, E in sols)
    print(f'  our solver finds {len(sols)} solution(s); among them the paper\'s: {getroffen}')
    for M, E in sorted(sols)[:6]:
        Mp, Ep = Poly(list(M), w), Poly(list(E), w)
        print(f'     M(0)={int(Mp.TC()):5d}  E(0)={int(Ep.TC()):5d}  degM={Mp.degree()} degE={Ep.degree()}')
    print(f'  Conant\'s bound E(0) = (l^2-1)/8 mod 2 for l=3: {(9-1)//8 % 2} (odd expected)')
    return ok3 and ok9 and lhs == rhs and getroffen

if __name__ == '__main__':
    if not eichung():
        print('\n🔴 CALIBRATION FAILED — no computation on our knots.'); sys.exit(1)
    print('\n✓ Calibration passed.')

def nabla_w(nab_z_str):
    from sympy import symbols as S
    z = S('z')
    p = Poly(sympify(nab_z_str), z)
    co = [int(p.coeff_monomial(z**(2*i))) for i in range(p.degree()//2+1)]
    return Poly(list(reversed(co)), w)

def analyse(nab, ell, maxsol=400):
    """Corollary 7.3 on one knot: which sign convention carries a solution?"""
    res = {}
    for lsig in {ell, abs(ell)}:
        for e in (1, -1):
            P = Poly(expand(e*bracket(lsig).as_expr()*nab.as_expr()), w)
            try: sols = loese_ME(P)
            except Exception as ex: res[(lsig, e)] = ('FEHLER', str(ex)[:40]); continue
            if sols:
                e0 = {abs(int(Poly(list(E), w).TC())) % 2 for _, E in sols}
                degs = {Poly(list(E), w).degree() for _, E in sols}
                res[(lsig, e)] = (len(sols), sorted(e0), sorted(degs))
            else:
                res[(lsig, e)] = (0, None, None)
    return res
