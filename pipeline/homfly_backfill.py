#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""homfly_backfill.py — add missing HOMFLY polynomials, only with agreeing witnesses.

OCCASION (2026-08-23): 229 records carry no HOMFLY. The reason was NOT that
regina failed, but that the run ran in point_eval mode and did not even send
the exponential lines along
(`pflicht_invarianten.py:116`, `_exp_lines = [] if point_eval else [...]`).
A skipped computation is not a finding.

Sibling tool of det_backfill.py, same doctrine:
TWO AGREEING WITNESSES, otherwise nothing is written.
  Witness 1  det from the NEW HOMFLY: |P(a=1, z=2i)|   (HOMFLY -> Alexander -> t=-1)
  Witness 2  the determinant LYING in the record (from regina-Alexander / Seifert),
             i.e. an entirely different computation.
  Witness 3  |Alexander(-1)| from the stored Alexander string, where present.

CALIBRATION is mandatory and built in (--eichung): on records that ALREADY
carry a HOMFLY, the same route (signed_gauss -> regina.Link -> simplify -> homfly)
must reproduce the same value. Without a passed calibration nothing is written —
green over an empty base set is not green.
"""
import sys, os, json, argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from homfly_calibration import parse_homfly, eval_homfly
from homfly_crosscheck import check_single
from homfly_conventions import variablen, gleich_ueber_konventionen


def poly_gleich(s1, s2, tol=1e-9):
    """Two HOMFLY strings for equality — as POLYNOMIALS, not as text, and
    across CONVENTION BOUNDARIES.

    The corpus carries HOMFLY in three notations: regina's AZ (x,y), LM (L,M)
    and SymPy's vz (v,z, partly with fractions). A text comparison then reports
    deviations that are none — happened exactly so on 23.08. on 4 of 89 calibration records.
    The substitutions are MEASURED, not set:
        L,M  →  x = i*u,   y = i*w      (calibrated on regina's own AZ/LM pair)
        v,z  →  x = 1/u,   y = w
    (homfly_conventions.py, selftest with two counter-cases)."""
    vp = variablen(s1)
    if vp and vp != ('x', 'y'):
        g, _ = gleich_ueber_konventionen(s1, s2)
        return g
    p1, p2 = parse_homfly(s1), parse_homfly(s2)
    if p1 is None or p2 is None:
        return None
    for k in set(p1) | set(p2):
        if abs(p1.get(k, 0) - p2.get(k, 0)) > tol:
            return False
    return True


def eichung(wanddatei, laut=True):
    """Records that already carry a HOMFLY: does the route reproduce it?"""
    erg = json.load(open(wanddatei))
    ok = fehl = 0
    for e in erg:
        if e.get('status') != 'OK' or not e.get('homfly'):
            continue
        rec = json.load(open(e['record']))
        alt = rec.get('pflichtliste', {}).get('ALGEBRA', {}).get('homfly_polynomial')
        if not alt:
            continue
        g = poly_gleich(alt, e['homfly'])
        if g:
            ok += 1
        else:
            fehl += 1
            if laut:
                print(f"  🔴 CALIBRATION OFF {e['name'][:38]}")
                print(f"     Record: {str(alt)[:110]}")
                print(f"     new   : {str(e['homfly'])[:110]}")
    return ok, fehl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--wand', required=True,
                    help='result of homfly_wand.py (records WITHOUT HOMFLY)')
    ap.add_argument('--eich-wand', required=True,
                    help='result of homfly_wand.py on records WITH HOMFLY')
    ap.add_argument('--datum', required=True)
    ap.add_argument('--min-eich', type=int, default=20,
                    help='at least this many known values must be hit')
    ap.add_argument('--schreiben', action='store_true')
    a = ap.parse_args()

    print('── CALIBRATION on records that already carry a HOMFLY ──')
    e_ok, e_fehl = eichung(a.eich_wand)
    print(f'  {e_ok} reproduced, {e_fehl} off')
    if e_fehl:
        print('\n🔴 Calibration off ⟹ NOTHING will be written.')
        return 1
    if e_ok < a.min_eich:
        print(f'\n🔴 only {e_ok} known values calibrated (< {a.min_eich}) ⟹ '
              'NOTHING will be written.')
        return 2

    print('\n── WITNESSES per candidate ──')
    erg = json.load(open(a.wand))
    kand, ab = [], []
    for e in erg:
        if e.get('status') != 'OK' or not e.get('homfly'):
            continue
        rec = json.load(open(e['record']))
        if rec.get('pflichtliste', {}).get('ALGEBRA', {}).get('homfly_polynomial'):
            continue                       # already carries one, do not touch
        r = check_single(e['record'], e['homfly'])
        if not r['ok']:
            ab.append((e['name'], r)); continue
        kand.append((e, r))
    print(f'  {len(kand)} with agreeing witnesses, {len(ab)} rejected')
    for nm, r in ab[:10]:
        print(f'   🔴 {nm[:38]:<38} {r.get("grund","det deviation")} '
              f'rec={r.get("det_record")} homfly={r.get("det_homfly")}')

    drei = sum(1 for _, r in kand if r.get('zeugen') == 3)
    print(f'  of which with THREE witnesses: {drei}')

    for e, r in kand:
        print(f"  {'+' if a.schreiben else '·'} {e['name'][:38]:<38} "
              f"n={e.get('n'):>3} c={e.get('c'):>4} det={r['det_record']}")
        if not a.schreiben:
            continue
        f = e['record']
        rec = json.load(open(f))
        rec['pflichtliste']['ALGEBRA']['homfly_polynomial'] = e['homfly']
        # Origin comes from the entry, where it brings one along. The default text
        # below describes the regina route; a polynomial that came via another
        # route (e.g. schmaler.py+mirror) carries its own — otherwise
        # the record contains a computation that never took place.
        rec['homfly_herkunft'] = {
            'datum': a.datum,
            'weg': e.get('weg') or
                   'signed_gauss -> regina.Link.fromSignedGauss -> simplify() -> homfly()',
            'konvention': 'regina AZ, Variablen gedruckt als x=a, y=z',
            'warum_fehlte_es': e.get('warum_fehlte_es') or
                               ('Lauf im point_eval-Modus: die exponentiellen '
                                'regina-Zeilen (homfly/jones/bracket) wurden gar '
                                'not sent along in the first place — skipped, not '
                                'gescheitert.'),
            'zeuge_1': f"det aus HOMFLY |P(1,2i)| = {r['det_homfly']:.0f}",
            'zeuge_2': f"Determinante im Record = {r['det_record']}",
            'zeuge_3': (f"|Alexander(-1)| = {r['det_alexander']}"
                        if r.get('det_alexander') is not None else None),
            'zeugen_einig': True,
            'eichung': (e.get('eichung_text') or
                        f'{e_ok} records with existing HOMFLY on the same route '
                        f'bitgleich reproduziert, 0 daneben'),
            'wandzeit_s': e.get('t_homfly'),
            'simplified_size': e.get('size_simplified'),
        }
        tmp = f + '.tmp'
        json.dump(rec, open(tmp, 'w'), indent=1)
        os.replace(tmp, f)                  # atomic, never a half file
    print(f"\n{'WRITTEN' if a.schreiben else 'DRY RUN'}: {len(kand)} records")
    return 0


if __name__ == '__main__':
    sys.exit(main())
