#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""fox_milnor.py — the FULL Fox-Milnor criterion, not only the determinant.

Until now (slice_verdict.py, Disc #555) we only check: is det = |Delta(-1)| a square number?
That is a POINT evaluation of the Fox-Milnor theorem. The theorem itself says more:

    K slice  =>  Delta_K(t) = f(t) * f(1/t)   (up to units +- t^k)

As a factorisation condition over Q[t]: write Delta = prod q_i^{e_i} with irreducible q_i,
and let q* be the reciprocal of q (i.e. t^deg(q) * q(1/t), normalised). Then the condition is exactly:
    - q self-reciprocal   ->  e(q) MUST BE EVEN
    - otherwise           ->  e(q) == e(q*)
Because f = prod q^{a} implies f(1/t) = prod (q*)^{a}, hence e(q) = a(q) + a(q*).

This is genuinely stronger than the determinant test: the trefoil has det 3 (not a square,
fails already before), but the figure-eight would have det 5 -- and knots like 5_2 with det 7 too.
What matters are the cases WITH a square determinant: there only the factorisation separates.
Examples: 6_1 and 9_46 (both slice) have Delta = 2t^2-5t+2 = (2t-1)(t-2), and the two
factors are reciprocal to each other -> pass. knots-of-love has (t^2-t+1)^6: self-reciprocal
with even multiplicity -> passes (and is indeed slice, Disc #557).
"""
import re
import sympy as sp

t = sp.Symbol('t')


def parse_alexander(s):
    """Record string (variable x or t) -> sympy polynomial in t, normalised to a constant term != 0."""
    if not s or not str(s).strip():
        return None
    s = str(s).replace('^', '**')
    s = re.sub(r'(?<=[\d\)])\s+(?=[xt])', '*', s)
    s = re.sub(r'(?<=[xt])\s+(?=[xt])', '*', s)
    s = re.sub(r'(\*\*-?\d+)\s+(?=[xt])', r'\1*', s)
    try:
        e = sp.sympify(s, locals={'x': t, 't': t})
    except Exception:
        return None
    e = sp.together(sp.expand(e))
    num, den = sp.fraction(sp.cancel(e))
    try:
        P = sp.Poly(sp.expand(num), t)
    except Exception:
        return None
    if P.total_degree() == 0 and P.LC() == 0:
        return None
    # pull out t^k so that the constant term is != 0 (Delta is only determined up to +- t^k)
    c = P.all_coeffs()
    while len(c) > 1 and c[-1] == 0:
        c = c[:-1]
    return sp.Poly(c, t)


def _reziprok(q):
    """q*(t) = t^deg(q) * q(1/t), primitive and normalised to a positive leading coefficient."""
    c = q.all_coeffs()
    r = sp.Poly(list(reversed(c)), t)
    return _norm(r)


def _norm(q):
    q = q.primitive()[1]
    if q.LC() < 0:
        q = -q
    # pull out t^k
    c = q.all_coeffs()
    while len(c) > 1 and c[-1] == 0:
        c = c[:-1]
    return sp.Poly(c, t)


def fox_milnor(alex_str):
    """-> (besteht, grund, details) = (passes, reason, details). besteht=False means PROVEN not slice.
    besteht=None means: not decidable (no/unreadable Delta).
    Note: the reason strings and the 'faktoren' key are returned data and stay German in this release."""
    P = parse_alexander(alex_str)
    if P is None:
        return None, 'kein lesbares Alexander-Polynom', {}
    if P.total_degree() == 0:
        return True, 'Delta ist konstant (= 1 bis auf Einheiten): Fox-Milnor trivial erfuellt', {}
    try:
        _, faktoren = sp.factor_list(P.as_expr(), t)
    except Exception as e:
        return None, 'Faktorisierung fehlgeschlagen: %s' % str(e)[:80], {}
    mult = {}
    for f, e in faktoren:
        try:
            q = _norm(sp.Poly(f, t))
        except Exception:
            continue
        if q.total_degree() == 0:
            continue
        k = q.as_expr()
        mult[k] = mult.get(k, 0) + e
    details = {'faktoren': {str(k): v for k, v in mult.items()}}
    for k, e in mult.items():
        q = sp.Poly(k, t)
        qs = _reziprok(q).as_expr()
        if sp.simplify(qs - k) == 0:
            if e % 2 == 1:
                return False, ('selbst-reziproker Faktor (%s) mit UNGERADER Vielfachheit %d'
                               % (sp.sstr(k), e)), details
        else:
            if mult.get(qs, 0) != e:
                return False, ('Faktor (%s) hat Vielfachheit %d, sein Reziprokes (%s) aber %d'
                               % (sp.sstr(k), e, sp.sstr(qs), mult.get(qs, 0))), details
    return True, 'Fox-Milnor erfuellt: Delta = f(t)*f(1/t) moeglich', details
