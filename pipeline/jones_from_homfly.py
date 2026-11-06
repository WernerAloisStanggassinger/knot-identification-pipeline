#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Jones from HOMFLY — by MEASURING the substitution, not by setting it.

WHY THIS EXISTS
    On 24.08.2026 there were 759 records with HOMFLY, but only 553 with Jones.
    Jones is a specialisation of HOMFLY: the 206 missing polynomials
    were already on disk, just not computed out. No knot has to be touched
    anew for this.

THE SUBSTITUTION IS MEASURED, NOT SET
    752 records carry BOTH. On them the candidate space
    x -> s1*X^a,  y -> s2*(X^b - X^-b)     (X = t^(1/2))
    was tried through completely for a in {+-1,+-2}, b in {+-1}, s in {+-1}.
    Exactly one class hits, and bit-exactly:
        (x,y) convention (regina):  x -> X^-2,  y -> (X - X^-1)
        749 of 752 records identical term for term.
    The signs s1,s2 are FREE, because in this corpus all exponents
    are even — that is measured (all eight sign variants give
    the same hit count), not assumed.

    The three deviants were NOT an error: they are in the
    Lickorish-Millett convention (L,M). There the same substitution hits
    with the additional factor (-1)^((i+j)/2) — bit-exactly on all three.
    Related: feedback-calibrate-each-convention-alone-not-the-product.

⚠️ COMPUTE INTEGER-ONLY
    During calibration my own checker swallowed a bit: `s1**i` with
    negative i gives 1.0 in Python — a float — and with that the whole
    substitution computed in floating point. With coefficients beyond 2^53 that loses
    exact integers and reported 747 instead of 749 hits. The same error class
    as Disc #688 on the same day, an hour later, in the tool next door.
    Hence: NOTHING runs through float here. Test 5 in the selftest pins this down.

Author: Werner Stanggassinger + Claude Code, 2026-08-24.
"""
import site_config
import argparse
import glob
import json
import os
import re
import sys
from math import comb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from homfly_conventions import parse, variablen          # noqa: E402

# Measured on 24.08.2026 (see header). Numbers, not belief.
SUBST = {('x', 'y'): dict(a=-2, b=1, vz_parität=0),
         ('L', 'M'): dict(a=-2, b=1, vz_parität=1)}
EICHSTAND = {'paare': 752, 'bitgenau_xy': 749, 'bitgenau_LM': 3}


def parse_jones(s, u='x'):
    """'x^10 - 6 x^8 + 11 - 48 x^-2' -> {exponent: coefficient}, integer."""
    t = s.strip().replace('**', '^').replace('*', ' ')
    # Split only at binary signs — '-' after '^' belongs to the
    # exponent (footgun 15, happened on 23.08. with two lj_ records).
    t = re.sub(r'(?<!\^)([+-])', r'|\1', t)
    out = {}
    for st in t.split('|'):
        st = st.strip()
        if not st:
            continue
        vz = -1 if st.startswith('-') else 1
        st = st.lstrip('+-').strip()
        koef, e, gefunden = 1, 0, False
        for stueck in st.split():
            m = re.fullmatch(rf'{u}(?:\^(-?\d+))?', stueck)
            if m:
                e = int(m.group(1)) if m.group(1) else 1
                gefunden = True
                continue
            m = re.fullmatch(r'(\d+)', stueck)
            if m:
                koef *= int(m.group(1))
                gefunden = True
                continue
            return None
        if not gefunden:
            return None
        out[e] = out.get(e, 0) + vz * koef
    return {k: v for k, v in out.items() if v} or None


def jones_aus(H, konvention=('x', 'y')):
    """HOMFLY dict {(i,j): c} -> Jones dict {exponent of X: c}. INTEGER."""
    s = SUBST.get(tuple(konvention))
    if s is None:
        return None
    a, b, vzp = s['a'], s['b'], s['vz_parität']
    erg = {}
    for (i, j), c in H.items():
        for k in range(j + 1):
            ex = a * i + b * (j - 2 * k)
            f = c * comb(j, k) * ((-1) ** k)
            if vzp:
                f *= (-1) ** (((i + j) // 2) * vzp)
            erg[ex] = erg.get(ex, 0) + f
    return {k: v for k, v in erg.items() if v}


def als_string(P):
    """Back into the notation found in the corpus: 'x^10 - 6 x^8 + 11'."""
    if not P:
        return '0'
    teile = []
    for e in sorted(P, reverse=True):
        c = P[e]
        vz = '-' if c < 0 else '+'
        b = abs(c)
        if e == 0:
            t = f'{b}'
        else:
            t = ('' if b == 1 else f'{b} ') + ('x' if e == 1 else f'x^{e}')
        teile.append((vz, t))
    s = ('-' if teile[0][0] == '-' else '') + teile[0][1]
    for vz, t in teile[1:]:
        s += f' {vz} {t}'
    return s


def auswerten(P, t):
    """Evaluate P at an integer X — exactly, via fractions."""
    from fractions import Fraction
    return sum(Fraction(c) * Fraction(t) ** e for e, c in P.items())


def bei_t_minus_1(P):
    """V(t=-1) — and t=-1 means X=i, NOT X=-1.

    The polynomial is in X = t^(1/2). Whoever substitutes X=-1 here evaluates at
    t=+1 and always gets 1: a witness that never contradicts
    looks like a witness that always agrees. (Built exactly like that on 24.08.
    and refuted at once by the trefoil.)
    For even exponents i^e = (-1)^(e/2) — integer, no complex numbers.
    ⚠️ Do NOT write it as (-1)**(e//2): for negative e//2 Python gives 1.0,
    i.e. float, and the whole witness computes in floating point. Happened
    exactly here on 24.08. — the third manifestation of the same error class on one
    day (Disc #688, then in the calibration checker, then here). Hence parity, not
    power.
    """
    if any(e % 2 for e in P):
        return None
    return sum(c if (e // 2) % 2 == 0 else -c for e, c in P.items())


def zeugen(P, det=None):
    """What must hold for every KNOT? Three cheap, hard checks."""
    gruende = []
    if any(e % 2 for e in P):
        ung = [e for e in P if e % 2][:3]
        gruende.append(f'odd exponents in X=t^(1/2): {ung} — '
                       'a knot has only integer t-powers '
                       '(odd ⟹ link, not knot)')
    v1 = auswerten(P, 1)
    if v1 != 1:
        gruende.append(f'V(1) = {v1}, must be 1')
    if det is not None:
        vm1 = bei_t_minus_1(P)
        if vm1 is None:
            pass                      # already reported above as odd
        elif abs(vm1) != det:
            gruende.append(f'|V(-1)| = {abs(vm1)} != determinant {det}')
    return gruende


# ── collect records ────────────────────────────────────────────────────
def records(wurzeln=None):
    wurzeln = wurzeln or [(site_config.GENERATED_KNOTS_DIR),
                          (site_config.ROOT + '/own_knots')]
    gesehen = set()
    for r in wurzeln:
        for f in glob.glob(r + '/**/*_complete.json', recursive=True):
            if 'pipeline_failures' in f:
                continue
            rp = os.path.realpath(f)
            if rp in gesehen:
                continue
            gesehen.add(rp)
            yield f


def det_aus(A):
    for k in ('determinant', 'det'):
        if isinstance(A.get(k), int):
            return A[k]
    return None


def bearbeite(f):
    rec = json.load(open(f))
    A = rec['pflichtliste']['ALGEBRA']
    h, j = A.get('homfly_polynomial'), A.get('jones_polynomial')
    if not h:
        return None
    vp = variablen(h)
    if not vp:
        return {'datei': f, 'ok': False, 'hat_jones': bool(j),
                'gruende': ['variable pair in the HOMFLY string not recognised']}
    H = parse(h, *vp)
    if not H:
        return {'datei': f, 'ok': False, 'hat_jones': bool(j),
                'gruende': [f'HOMFLY in ({vp[0]},{vp[1]}) not parsable']}
    P = jones_aus(H, vp)
    if P is None:
        return {'datei': f, 'ok': False, 'hat_jones': bool(j),
                'gruende': [f'no substitution calibrated for convention {vp} '
                            '— measure first, then enter']}
    g = zeugen(P, det_aus(A))
    aus = {'datei': f, 'konvention': vp, 'hat_jones': bool(j),
           'jones': als_string(P), 'terme': len(P), 'gruende': g,
           'ok': not g}
    if j:
        gespeichert = parse_jones(j)
        aus['reproduziert'] = (gespeichert == P) if gespeichert else None
    return aus


def schreibe(r, datum):
    f = r['datei']
    rec = json.load(open(f))
    rec['pflichtliste']['ALGEBRA']['jones_polynomial'] = r['jones']
    rec['jones_herkunft'] = {
        'datum': datum,
        'weg': ('specialisation of the stored HOMFLY polynomial: '
                f"convention {r['konvention'][0]},{r['konvention'][1]} → "
                'x → X^-2, y → (X − X^-1) with X = t^(1/2), integer '
                '(tools/jones_aus_homfly.py)'),
        'eichung': (f"substitution measured on {EICHSTAND['paare']} records "
                    f"carrying BOTH polynomials; {EICHSTAND['bitgenau_xy']} "
                    'of them reproduced identically term for term'),
        'zeugen': ['only even exponents in X (knot, not link)',
                   'V(1) = 1',
                   '|V(-1)| = determinant in the record'],
        'kein_neuer_knoten_gerechnet': True}
    json.dump(rec, open(f, 'w'), indent=1, ensure_ascii=False)


# ── Selftest: six cases, four of them counter-cases ──────────────────────
def selbsttest():
    ok = True

    # (1) Trefoil: HOMFLY in (x,y) as regina writes it
    H3 = parse('x^-2 y^2 + 2 x^-2 - x^-4', 'x', 'y')
    P3 = jones_aus(H3)
    g = zeugen(P3, 3)
    gut = not g
    ok &= gut
    print(f'  3_1  from HOMFLY : {als_string(P3)}')
    print(f'       V(1)={auswerten(P3,1)} · |V(-1)|={abs(bei_t_minus_1(P3))} '
          f'(det 3) {"✅" if gut else "🔴 " + "; ".join(g)}')

    # (2) Figure-eight knot
    H4 = parse('x^2 - y^2 - 1 + x^-2', 'x', 'y')
    P4 = jones_aus(H4)
    g = zeugen(P4, 5)
    gut = not g
    ok &= gut
    print(f'  4_1  from HOMFLY : {als_string(P4)}')
    print(f'       V(1)={auswerten(P4,1)} · |V(-1)|={abs(bei_t_minus_1(P4))} '
          f'(det 5) {"✅" if gut else "🔴 " + "; ".join(g)}')

    # (3) COUNTER-CASE: one coefficient bent — the witnesses must bite
    kaputt = dict(H4)
    schluessel = sorted(kaputt)[0]
    kaputt[schluessel] += 1
    g = zeugen(jones_aus(kaputt), 5)
    gut = bool(g)
    ok &= gut
    print(f'  COUNTER-CASE coefficient +1 : '
          f'{"caught ✅ (" + g[0][:46] + "…)" if gut else "🔴 slips through"}')

    # (4) COUNTER-CASE: wrong determinant slipped in
    g = zeugen(P4, 7)
    gut = bool(g)
    ok &= gut
    print(f'  COUNTER-CASE det 5 as 7    : '
          f'{"caught ✅" if gut else "🔴 slips through"}')

    # (5) COUNTER-CASE floating point: the result must be an INTEGER.
    #     Exactly here the checker swallowed a bit itself on 24.08.
    gross = {(-2, 0): 9007199254740993}          # 2^53 + 1
    P = jones_aus(gross)
    gut = all(isinstance(v, int) for v in P.values()) \
        and 9007199254740993 in P.values()
    # and the WITNESS itself — it has to process NEGATIVE exponents,
    # and exactly there (-1)**(e//2) gave a float. Without this case the
    # counter-case is blind at the place where the error sat.
    w = bei_t_minus_1({-4: 9007199254740993, 0: 1})
    gut &= isinstance(w, int)
    ok &= gut
    print(f'  COUNTER-CASE 2^53+1 exact  : '
          f'{"substitution AND witness stay integer ✅" if gut else "🔴 float swallows a bit"}')

    # (6) THE LOAD-BEARING CASE: against the records that carry BOTH.
    #     Does the substitution reproduce what is stored independently?
    treffer = fehl = 0
    kaputte = []
    for f in records():
        try:
            A = json.load(open(f))['pflichtliste']['ALGEBRA']
        except Exception:
            continue
        if not (A.get('homfly_polynomial') and A.get('jones_polynomial')):
            continue
        r = bearbeite(f)
        if not r or r.get('reproduziert') is None:
            continue
        if r['reproduziert']:
            treffer += 1
        else:
            fehl += 1
            kaputte.append(os.path.basename(f)[:48])
    gut = (fehl == 0)
    ok &= gut
    print(f'  against {treffer + fehl} records with BOTH polynomials: '
          f'{treffer} equal term for term, {fehl} off '
          f'{"✅" if gut else "🔴"}')
    for k in kaputte[:5]:
        print(f'      🔴 {k}')

    # (7) COUNTER-CASE ON THE WITNESS ITSELF: at X=-1 every knot would give 1.
    #     The witness must DISTINGUISH between 3_1 (det 3) and 4_1 (det 5).
    blind = (abs(auswerten(P3, -1)) == abs(auswerten(P4, -1)))
    scharf = (abs(bei_t_minus_1(P3)) != abs(bei_t_minus_1(P4)))
    gut = blind and scharf
    ok &= gut
    print(f'  COUNTER-CASE witness at X=-1 : '
          f'{"would be blind (both 1), sharp at X=i (3 vs 5) ✅" if gut else "🔴"}')

    print('  ⇒', '✅ Jones from HOMFLY bites' if ok else '🔴 selftest failed')
    return 0 if ok else 1


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--selbsttest', action='store_true')
    ap.add_argument('--schreiben', action='store_true',
                    help='without this flag: dry run, nothing is touched')
    ap.add_argument('--datum', default='2026-08-24')
    ap.add_argument('--limit', type=int)
    ap.add_argument('--out')
    a = ap.parse_args()

    if a.selbsttest:
        sys.exit(selbsttest())

    ziele = []
    for f in records():
        try:
            A = json.load(open(f))['pflichtliste']['ALGEBRA']
        except Exception:
            continue
        if A.get('homfly_polynomial') and not A.get('jones_polynomial'):
            ziele.append(f)
    if a.limit:
        ziele = ziele[:a.limit]
    print(f'{len(ziele)} records with HOMFLY, without Jones', flush=True)

    erg, gut, schlecht = [], 0, 0
    for i, f in enumerate(ziele, 1):
        r = bearbeite(f)
        if r is None:
            continue
        erg.append(r)
        nm = os.path.basename(f).replace('_pflicht_invarianten_complete.json', '')
        if r['ok']:
            gut += 1
            if a.schreiben:
                schreibe(r, a.datum)
            print(f'[{i}/{len(ziele)}] ✅ {nm[:52]:<52} {r["terme"]:>3} terms',
                  flush=True)
        else:
            schlecht += 1
            print(f'[{i}/{len(ziele)}] 🔴 {nm[:52]:<52} '
                  f'{"; ".join(r["gruende"])[:70]}', flush=True)
    if a.out:
        json.dump(erg, open(a.out, 'w'), indent=1, ensure_ascii=False)
    print(f'\n{gut} through all three witnesses, {schlecht} rejected'
          f'{" — WRITTEN" if a.schreiben else " (dry run, nothing touched)"}')
    sys.exit(0 if schlecht == 0 else 1)
