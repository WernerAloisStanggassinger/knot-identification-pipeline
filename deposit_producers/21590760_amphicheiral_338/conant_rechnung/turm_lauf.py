#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Conant Thm 7.1 INTEGRAL on all 250 — not the mod-4 shadow, the exact statement.

  nu_a = Nabla_K,  nu_{k+1} = nu_k * R_k in Z[w] with R_k(0)=1,  nu_1 = phi^2.

What is checked is what CAN fail:
  P1  nu_1 is a square in Z[w]                    (Lemma 5.3, Hartley-Kawauchi)
  P2  nu_k divides nu_{k+1} exactly, quotient(0) = 1      (Thm 7.1, integral)
  P3  product formula hits Nabla_K on the coefficients
"""
import site_config
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from turm import quotient, alexander_von, delta_zu_nabla, teilt_exakt
from mod4_test import to_w
import importlib.util
SIEB = (site_config.ROOT + '/konversationen/kawauchi_20260725/sieb')
sp2 = importlib.util.spec_from_file_location('f', os.path.join(SIEB, 'forge_c.py'))
F = importlib.util.module_from_spec(sp2); sp2.loader.exec_module(F)


def istquadrat(c):
    """c in Z[w] a square? (c[0]=1 assumed) — root by comparison of coefficients."""
    n = len(c) - 1
    if n % 2:
        return None
    m = n // 2
    r = [1] + [0] * m
    for i in range(1, m + 1):
        s = sum(r[j] * r[i - j] for j in range(1, i))
        if (c[i] - s) % 2:
            return None
        r[i] = (c[i] - s) // 2
    prod = [0] * (2 * m + 1)
    for i, x in enumerate(r):
        for j, y in enumerate(r):
            prod[i + j] += x * y
    return r if prod == c else None


def mul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return r


if __name__ == '__main__':
    K = json.load(open((site_config.ROOT + '/PRIORITAET_DEPOSIT_20260729/counterexamples_data.json')))['knoten']
    t0 = time.time()
    st = dict(n=0, P1=0, P1_bad=0, P2=0, P2_bad=0, P3=0, P3_bad=0,
              nu1_unknoten=0, kein_diagramm=0)
    details = []
    for k in K:
        g = k['erzeugung']; a = k['a_tower']
        P = F.forge_c(g['k'], g['seed'], g['n_xy'], g['n_z'], g['decay'], c=g['c_rot'])
        nus = {a: to_w(k['conway_nabla'])}
        cur, fehl = P, None
        # NYQUIST: the angle doubling compresses the curve, 1440 points are not enough.
        # Measured on SK3_1_471_42_055: at 1440 came c=12 and a wrong nu_2, from 2880
        # stable c~26 and the correct one. Hence TWO samplings that must agree.
        for lvl in range(a - 1, 0, -1):
            cur = quotient(cur)
            nu = None
            a1, i1 = alexander_von(cur, 2880)
            a2, i2 = alexander_von(cur, 4320)
            if a1 is None or a2 is None:
                fehl = f'stufe{lvl}:{i1}/{i2}'; break
            n1, n2 = delta_zu_nabla(a1), delta_zu_nabla(a2)
            if n1 is None or n2 is None:
                fehl = f'stufe{lvl}:konversion'; break
            if n1 != n2:
                fehl = f'stufe{lvl}:uneinig'; break
            nus[lvl] = n1
        st['n'] += 1
        if fehl:
            st['kein_diagramm'] += 1
            details.append(dict(id=k['id'], a=a, fehler=fehl)); continue
        if nus[1] == [1]:
            st['nu1_unknoten'] += 1
        phi = istquadrat(nus[1])
        st['P1' if phi is not None else 'P1_bad'] += 1
        Rs, p2ok = [], True
        for lvl in range(1, a):
            q = teilt_exakt(nus[lvl], nus[lvl + 1])
            if q is None or q[0] != 1:
                p2ok = False; break
            Rs.append(q)
        st['P2' if p2ok else 'P2_bad'] += 1
        p3ok = False
        if phi is not None and p2ok:
            prod = mul(phi, phi)
            for R in Rs:
                prod = mul(prod, R)
            p3ok = prod == nus[a]
        st['P3' if p3ok else 'P3_bad'] += 1
        details.append(dict(id=k['id'], a=a, ell=k['achsen_verschlingung_ell'],
                            nu1=nus[1], phi=phi, R=Rs, P1=phi is not None, P2=p2ok, P3=p3ok))
    print(f"knots           : {st['n']}   ({time.time()-t0:.0f} s)")
    print(f"no diagram       : {st['kein_diagramm']}")
    n_ok = st['n'] - st['kein_diagramm']
    print(f"P1 nu_1 square  : {st['P1']}/{n_ok}   breaks {st['P1_bad']}   (of these nu_1=1: {st['nu1_unknoten']})")
    print(f"P2 nu_k | nu_k+1: {st['P2']}/{n_ok}   breaks {st['P2_bad']}")
    print(f"P3 product form.: {st['P3']}/{n_ok}   breaks {st['P3_bad']}")
    json.dump(dict(stat=st, details=details), open('ergebnis_turm_20260729.json', 'w'), indent=1)
    for d in details:
        if d.get('fehler') or d.get('P2') is False or d.get('P1') is False:
            print("   ", d.get('id'), d.get('a'), d.get('fehler') or f"P1={d.get('P1')} P2={d.get('P2')}")
