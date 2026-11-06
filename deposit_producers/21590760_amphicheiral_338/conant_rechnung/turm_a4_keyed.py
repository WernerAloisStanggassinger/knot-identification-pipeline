#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Conant Thm 7.1 INTEGRAL on the order-16 knots (a=4) — KEYED.

WHY THIS FILE EXISTS (lesson 2026-07-29):
`ergebnis_turm_a4_20260729.json` carries 48 entries with ONLY `seed` as identifier. The key
of a candidate is, however, `(k, seed, n_xy, n_z, decay, c_rot)` — 48 candidates are distributed over
27 seeds. Thus that result cannot be assigned to the DISTINCT knots, and the number
"28/28" belongs to a different inventory set than the balance number "7".
A letter that puts both side by side mixes two sets.

This file computes the tower on exactly the distinct a=4 knots from `BILANZ_DATEN.json`,
with full key, and writes it out again with key.

What is checked is what CAN fail (identical to turm_lauf.py):
  P1  nu_1 is a square in Z[w]                    (Lemma 5.3)
  P2  nu_k divides nu_{k+1} exactly, quotient(0) = 1      (Thm 7.1, integral)
  P3  product formula hits Nabla_K on the coefficients
For a=4 the tower has THREE stages (two for a=3).

Werner Stanggassinger + Claude Code (Anthropic), 2026-07-29.
"""
import site_config
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from turm import quotient, alexander_von, delta_zu_nabla, teilt_exakt
from mod4_test import to_w
from turm_lauf import istquadrat, mul
import importlib.util

SIEB = (site_config.ROOT + '/konversationen/kawauchi_20260725/sieb')
sp2 = importlib.util.spec_from_file_location('f', os.path.join(SIEB, 'forge_c.py'))
F = importlib.util.module_from_spec(sp2); sp2.loader.exec_module(F)

BILANZ = (site_config.ROOT + '/konversationen/kawauchi_20260725/BILANZ_DATEN.json')
A_MIN = int(os.environ.get('A_MIN', '4'))


def schluessel(r):
    return (r['k'], r['seed'], r['nxy'], r['nz'], r['dc'], r['crot'])


if __name__ == '__main__':
    zeilen = [r for r in json.load(open(BILANZ)) if (r.get('a_tower') or 0) >= A_MIN]
    print(f"distinct knots with a>={A_MIN}: {len(zeilen)}")
    t0 = time.time()
    st = dict(n=0, P1=0, P1_bad=0, P2=0, P2_bad=0, P3=0, P3_bad=0, kein_diagramm=0)
    details = []

    for r in zeilen:
        a = int(r['a_tower'])
        k, seed, nxy, nz, dc, crot = schluessel(r)
        P = F.forge_c(k, seed, nxy, nz, dc, c=crot)
        nus = {a: to_w(r['nabla'])}
        cur, fehl = P, None
        # TWO samplings must agree — the angle doubling halves the
        # effective resolution, and three agreeing sampling rates can see the same false
        # shadow (finding 29.07., d_min_rel below the margin).
        for lvl in range(a - 1, 0, -1):
            cur = quotient(cur)
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
        eintrag = dict(k=k, seed=seed, nxy=nxy, nz=nz, dc=dc, crot=crot,
                       a=a, ell=r.get('ell'), deg=r.get('deg'), lc=r.get('lc'),
                       vol=r.get('vol'))
        if fehl:
            st['kein_diagramm'] += 1
            eintrag['fehler'] = fehl
            details.append(eintrag)
            print(f"  k={k} seed={seed:<4} crot={crot} 🔴 {fehl}")
            continue
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
        eintrag.update(nu1=nus[1], phi=phi, R=Rs, stufen=len(Rs),
                       P1=phi is not None, P2=p2ok, P3=p3ok)
        details.append(eintrag)
        print(f"  k={k} seed={seed:<4} crot={crot} nxy={nxy} nz={nz} dc={dc} "
              f"stages={len(Rs)} P1={'✓' if phi is not None else '🔴'} "
              f"P2={'✓' if p2ok else '🔴'} P3={'✓' if p3ok else '🔴'}")

    n_ok = st['n'] - st['kein_diagramm']
    print(f"\nknots           : {st['n']}   ({time.time()-t0:.0f} s)")
    print(f"no diagram       : {st['kein_diagramm']}")
    print(f"P1 nu_1 square  : {st['P1']}/{n_ok}   breaks {st['P1_bad']}")
    print(f"P2 nu_k | nu_k+1: {st['P2']}/{n_ok}   breaks {st['P2_bad']}")
    print(f"P3 product form.: {st['P3']}/{n_ok}   breaks {st['P3_bad']}")
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       f'ergebnis_turm_a{A_MIN}_keyed_20260729.json')
    json.dump(dict(stat=st, details=details), open(out, 'w'), indent=1)
    print(f"-> {out}")
