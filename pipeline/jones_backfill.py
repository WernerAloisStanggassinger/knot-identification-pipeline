#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""jones_backfill.py — Jones from the skein run (#722) into the records.

── THE CALIBRATION IS THE PURPOSE, NOT THE EXTRA ────────────────────────────
So far the run was calibrated against `regina.bracket()` on sample knots
and against two witnesses (V(1)=1, |V(-1)|=det). What was missing: the comparison with
what is in the CORPUS. `--eiche N` recomputes N records with an EXISTING
`jones_polynomial` and compares **the finished regina string**, character
by character. Only afterwards is anything written.

── THE FORMAT IS NOT REBUILT ────────────────────────────────────────────────
The record holds Jones as a regina string in x with t = x^2. Our result
is in A with t = A^-4, i.e. x^k <-> A^-2k. The string is printed by
`regina.Laurent` itself — a hand-rebuilt format would be
exactly the kind of silent deviation that nobody finds later.

── PROVENANCE ───────────────────────────────────────────────────────────────
Written with source `skein_durchlauf` and date. Disc #719: the
provenance belongs to the date, not in the writer.
"""
import site_config
import sys, os, json, time, argparse, glob, shutil
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, (site_config.FREEZE_DIR))
from pipeline_seal import seal_all
seal_all(verbose=False)
import regina
import skein_pass as SD

ROOT = site_config.ROOT


def jones_A_zu_laurent(jones_A):
    """A exponent e  ->  x exponent -e/2. Only even e are permitted."""
    p = regina.Laurent()
    for e, k in jones_A.items():
        if e % 2:
            raise ValueError(f'odd A exponent {e} — not a knot Jones')
        p.set(-e // 2, k)
    return p


def record_datei(folder):
    for muster in ('*_pflicht_invarianten_complete.json', '*_pflicht_complete.json'):
        for wurzel in ('erzeugte_knoten', 'own_knots'):
            t = glob.glob(f'{ROOT}/{wurzel}/**/{folder}/{muster}', recursive=True)
            if t:
                return t[0]
    return None


def rechne_einen(arg):
    folder, nur_pruefen = arg
    f = record_datei(folder)
    if not f:
        return {'folder': folder, 'status': 'KEINE DATEI'}
    try:
        d = json.load(open(f))
        alg = d.get('pflichtliste', {}).get('ALGEBRA', {})
        sg = d.get('signed_gauss') or (d.get('raw_subprocess_outputs', {})
                                       .get('regina', {}).get('signed_gauss_reduced'))
        if not sg:
            return {'folder': folder, 'status': 'KEIN SIGNED_GAUSS'}
        L = regina.Link.fromSignedGauss(sg)
        L.simplify()
        t0 = time.time()
        r = SD.rechne(L, SD.WAHL, neustarts=20)
        jA = SD.jones_aus_klammer(r['bracket'], int(L.writhe()))
        ok, wie = SD.zeugen_pruefe(jA, alg.get('determinant'))
        neu = str(jones_A_zu_laurent(jA))
        alt = alg.get('jones_polynomial')
        return {'folder': folder, 'status': 'OK', 'datei': f,
                'c': r['c'], 'breite': r['schnittbreite'],
                'zustaende': r['max_zustaende'], 'sek': round(time.time() - t0, 2),
                'neu': neu, 'alt': alt, 'gleich': (alt == neu) if alt else None,
                'zeugen_ok': ok, 'zeugen': wie}
    except Exception as ex:
        return {'folder': folder, 'status': f'{type(ex).__name__}: {ex}'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--eiche', type=int, default=0,
                    help='recompute and compare N records WITH jones')
    ap.add_argument('--fehlende', action='store_true',
                    help='compute records WITHOUT jones')
    ap.add_argument('--schreiben', action='store_true')
    ap.add_argument('--kerne', type=int, default=10)
    ap.add_argument('--out')
    a = ap.parse_args()

    db = json.load(open(site_config.PIPELINE_DIR + '/corpus_db.json'))
    mit = [r['folder'] for r in db['knots'] if r.get('jones') and r.get('folder')]
    ohne = [r['folder'] for r in db['knots'] if not r.get('jones') and r.get('folder')]
    print(f'  Corpus: {len(mit)} with Jones, {len(ohne)} without')

    if a.eiche:
        schritt = max(1, len(mit) // a.eiche)
        ziele = mit[::schritt][:a.eiche]          # spread over the corpus
        print(f'  CALIBRATION on {len(ziele)} records with existing Jones\n')
    elif a.fehlende:
        ziele = ohne
        print(f'  {len(ziele)} records WITHOUT Jones\n')
    else:
        sys.exit('🔴 --eiche N or --fehlende')

    erg = []
    with ProcessPoolExecutor(max_workers=a.kerne) as ex:
        futs = {ex.submit(rechne_einen, (z, bool(a.eiche))): z for z in ziele}
        for fu in as_completed(futs):
            r = fu.result()
            erg.append(r)
            if r['status'] != 'OK':
                print(f"  {r['folder'][:30]:32} 🔴 {r['status']}")
            elif a.eiche:
                m = '✅' if r['gleich'] else '🔴 DEVIATION'
                print(f"  {r['folder'][:30]:32} c={r['c']:<4} width {r['breite']:<3} "
                      f"{r['sek']:>6}s  {m}")
            else:
                print(f"  {r['folder'][:30]:32} c={r['c']:<4} width {r['breite']:<3} "
                      f"{r['sek']:>6}s  {'✅' if r['zeugen_ok'] else '🔴'} {r['zeugen']}")

    gut = [r for r in erg if r['status'] == 'OK']
    if a.eiche:
        gl = sum(1 for r in gut if r['gleich'])
        print(f"\n  ⇒ CALIBRATION: {gl}/{len(gut)} character-identical, "
              f"{len(erg)-len(gut)} not computable")
        abweich = [r for r in gut if not r['gleich']]
        for r in abweich[:3]:
            print(f"     🔴 {r['folder']}\n        old: {str(r['alt'])[:100]}"
                  f"\n        new: {str(r['neu'])[:100]}")
    else:
        z = sum(1 for r in gut if r['zeugen_ok'])
        print(f"\n  ⇒ {len(gut)} computed, {z} with both witnesses")
        if a.schreiben:
            if z != len(gut):
                sys.exit('🔴 not all witnesses green — nothing is written')
            sich = f'{ROOT}/claude_setup_v2/tmp_durable/backup_jones_skein_20260826'
            os.makedirs(sich, exist_ok=True)
            heute = time.strftime('%Y-%m-%d')
            for r in gut:
                shutil.copy2(r['datei'], f"{sich}/{os.path.basename(r['datei'])}")
                d = json.load(open(r['datei']))
                d['pflichtliste']['ALGEBRA']['jones_polynomial'] = r['neu']
                hk = d.setdefault('herkunft', {})
                hk[f'jones_polynomial@{heute}'] = 'skein_durchlauf (Disc #722)'
                json.dump(d, open(r['datei'], 'w'), indent=1)
                print(f"  written: {r['folder']}")
            print(f'  Backups: {sich}')
    if a.out:
        json.dump(erg, open(a.out, 'w'), indent=1)
        print(f'  written: {a.out}')


if __name__ == '__main__':
    main()
