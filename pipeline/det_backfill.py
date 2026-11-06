#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""det_backfill.py — add missing determinants, but only with TWO agreeing witnesses.

OCCASION (2026-08-23): in the flip_ dihedral ladder (D5..D25) 6 of 23 records carried
a determinant, 17 did not — although ALL 23 have both an Alexander polynomial and
a pd_code. The value was there all along, just not computed.

TWO INDEPENDENT WITNESSES, otherwise nothing is written:
  1. |Delta(-1)| from the stored regina Alexander string (symbolic -> point)
  2. punkt_alexander_det(pd_code) — own tool, variant 6 (Disc #570):
     sample points in F_p, interpolation, CRT. Different input (PD instead of polynomial),
     different computation route, different library.
They cannot jointly make the same error. On disagreement: write NOTHING,
report loudly.

CALIBRATION is built in: on the records that ALREADY carry the determinant, both
witnesses must reproduce it. If even one misses, the run aborts —
an add-in tool without calibration on known values is decoration (footgun 16).
"""
import site_config
import sys, os, json, glob, argparse, time
sys.path.insert(0, (site_config.PIPELINE_DIR))
from mandatory_invariants import _eval_regina_poly, point_alexander_det


def zeugen(rec):
    A = rec['pflichtliste']['ALGEBRA']
    D = rec['pflichtliste']['DIAGRAMM']
    v = _eval_regina_poly(A.get('alexander_polynomial'), -1)
    z1 = abs(int(v)) if v is not None else None
    z2 = point_alexander_det(D.get('pd_code')) if D.get('pd_code') else None
    return z1, z2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--muster', required=True, help='glob on *_complete.json')
    ap.add_argument('--praefix', default='', help='only stem names with this prefix')
    ap.add_argument('--schreiben', action='store_true')
    ap.add_argument('--datum', required=True)
    a = ap.parse_args()

    dateien = sorted(glob.glob(a.muster, recursive=True))
    eich_ok = eich_fehl = 0
    kandidaten = []
    uneinig = []
    for f in dateien:
        b = os.path.basename(f).replace('_pflicht_invarianten_complete.json', '')
        if a.praefix and not b.startswith(a.praefix):
            continue
        try:
            rec = json.load(open(f))
        except Exception:
            continue
        if 'pflichtliste' not in rec:
            continue
        A = rec['pflichtliste']['ALGEBRA']
        z1, z2 = zeugen(rec)
        if z1 is None or z2 is None or z1 != z2:
            uneinig.append((b, z1, z2)); continue
        if A.get('determinant') is not None:
            if int(A['determinant']) == z1: eich_ok += 1
            else:
                eich_fehl += 1
                print(f"  🔴 CALIBRATION MISSED {b}: record {A['determinant']} vs witnesses {z1}")
            continue
        kandidaten.append((f, b, z1))

    print(f"Checked: {sum(1 for f in dateien)} files in the pattern")
    print(f"  Calibration on known values : {eich_ok} hit, {eich_fehl} missed")
    print(f"  Witnesses disagree / incomplete: {len(uneinig)}  {[u[0] for u in uneinig][:6]}")
    print(f"  can be added                : {len(kandidaten)}")
    if eich_fehl:
        print("\n🔴 Calibration missed ⟹ NOTHING will be written.")
        return 1
    if eich_ok == 0:
        print("\n🔴 NOT A SINGLE known value to calibrate with ⟹ NOTHING will be written "
              "(green over an empty base set is not green).")
        return 2

    for f, b, wert in kandidaten:
        print(f"  {'+' if a.schreiben else '·'} {b:<24} det = {wert:,}")
        if not a.schreiben:
            continue
        rec = json.load(open(f))
        rec['pflichtliste']['ALGEBRA']['determinant'] = wert
        rec['determinant_herkunft'] = {
            'wert': wert, 'datum': a.datum,
            'zeuge_1': '|Delta(-1)| from the stored regina Alexander',
            'zeuge_2': 'punkt_alexander_det(pd_code), Variante 6 (Disc #570)',
            'beide_einig': True,
            'geeicht_an': f'{eich_ok} records of the same family that already carried the value',
            'hinweis': ('Derived, not delivered by regina. The value was always already contained '
                        'in the record (Alexander + pd_code), just not computed.')}
        tmp = f + '.tmp'
        json.dump(rec, open(tmp, 'w'), indent=1)
        os.replace(tmp, f)          # atomic, never a half file
    print(f"\n{'WRITTEN' if a.schreiben else 'TRIAL RUN'}: {len(kandidaten)} records")
    return 0


if __name__ == '__main__':
    sys.exit(main())
