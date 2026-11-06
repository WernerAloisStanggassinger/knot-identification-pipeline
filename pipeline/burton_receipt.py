#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""burton_receipt.py — fetch the receipt before 59 records are rewritten.

THE OCCASION (2026-09-04): The Burton run of 10.07.2026 FOUND 17 of our targets in the
census and booked them as `matched_c18_downgrade`. **14 of these still stand until
today as `novel_hyperbolic` in the corpus_db** — among them `sundown`(=trainstation,
one of the five sculptures built in steel), `spiral` and the 25-member cluster
S_18_d243_D3_chiral. The result lay unbooked on the disk for seven weeks.

WHY NOT SIMPLY ADOPT IT: `stage_b_result.json` stores per hit only
`{id, det, crossing}` — not WHICH census knot was hit. So the
finding is not re-checkable, and the Maguire experience says that errors sit
exactly here: there 2 of 42 retractions were exposed as (det,vol) near-collisions
(delta ~1e-5) that failed `is_isometric_to`. Retracting a novelty
is as much a statement of truth as asserting one.

WHAT THIS TOOL DOES: for each affected target, search for the census candidate with
matching determinant, compute its volume, and on volume proximity run
`is_isometric_to` — and then write NAME, DT code and both volumes into a
receipt file that a human can re-check.

--eichen runs three counter-cases, among them the decisive one: a target that according to
10.07. is NOT in the census must not get a hit here either. A tool
that only confirms cannot be distinguished from one that confirms everything.
"""
import site_config
import json, os, sys, time, argparse
import multiprocessing as mp

WORK = (site_config.MACHINE_YARD + '/corpus_vs_burton/work')
sys.path.insert(0, (site_config.MACHINE_YARD + '/corpus_vs_maguire/work'))
DB = (site_config.PIPELINE_DIR + '/corpus_db.json')
GEOM = 'all tetrahedra positively oriented'
PREFILTER = 1e-3


def kandidaten_streamen(ziel_dets):
    """Filter 745 MB of candidates line by line instead of loading them at once (kf320 is running alongside)."""
    puffer, tief, n, raus = [], 0, 0, []
    with open(WORK + '/candidates_stageA.json') as f:
        for zeile in f:
            s = zeile.strip()
            if s.startswith('{'):
                tief, puffer = 1, [s]
                continue
            if tief:
                puffer.append(s)
                if s.startswith('}'):
                    tief = 0
                    n += 1
                    txt = ''.join(puffer).rstrip(',')
                    try:
                        o = json.loads(txt)
                    except Exception:
                        continue
                    if o.get('det') in ziel_dets:
                        raus.append(o)
    return raus, n


def _vol(cand):
    import snappy, knotlib as K
    try:
        M = snappy.Manifold('DT:' + repr([tuple(K.decode_dt(cand['dt']))]))
        if str(M.solution_type()) != GEOM:
            return (cand, None)
        return (cand, float(M.volume()))
    except Exception:
        return (cand, None)


def _isometrisch(unser_isosig, burton_dt):
    import snappy, knotlib as K
    try:
        a = snappy.Manifold(unser_isosig)
        b = snappy.Manifold('DT:' + repr([tuple(K.decode_dt(burton_dt))]))
        if a.is_isometric_to(b):
            return True
        for _ in range(4):
            b.randomize()
            if a.is_isometric_to(b):
                return True
        return False
    except Exception:
        return None


def ziele_laden(nur_ids=None):
    """our targets with isosig + vol, from the corpus_db (not from the old target file)."""
    db = json.load(open(DB))
    ks = db['knots'] if 'knots' in db else db['records']
    ks = list(ks.values()) if isinstance(ks, dict) else ks
    ziele = {}
    for r in ks:
        if not r.get('is_distinct_representative'):
            continue
        i = r.get('distinct_knot_id')
        if nur_ids and i not in nur_ids:
            continue
        if not r.get('isosig') or r.get('vol') in (None, ''):
            continue
        ziele[i] = {'id': i, 'det': int(r['det']) if r.get('det') is not None else None,
                    'vol': float(r['vol']), 'isosig': r['isosig'],
                    'c_ub': r.get('c'), 'klasse': r.get('realizes_class'),
                    'n_real': r.get('n_realisierungen'), 'stem': r.get('stem')}
    return ziele


def lauf(ids, ausgabe, gegenprobe_ids=()):
    alle = list(ids) + list(gegenprobe_ids)
    ziele = ziele_laden(set(alle))
    fehlend = [i for i in alle if i not in ziele]
    if fehlend:
        print(f"⚠️ without isosig/vol in the corpus_db, not checkable: {fehlend}", flush=True)
    ziel_dets = {z['det'] for z in ziele.values() if z['det'] is not None}
    print(f"Targets: {len(ziele)}  (of which counter-check: {len(gegenprobe_ids)})  Determinants: {sorted(ziel_dets)}", flush=True)

    import snappy
    cal = float(snappy.Manifold('4_1').volume())
    assert abs(cal - 2.029883) < 1e-3, f"CALIBRATION FAILED 4_1 vol={cal}"
    print(f"CALIBRATION OK 4_1 vol={cal:.6f}", flush=True)

    t0 = time.time()
    kand, gesamt = kandidaten_streamen(ziel_dets)
    print(f"Candidates with matching det: {len(kand)} of {gesamt} ({time.time()-t0:.0f}s)", flush=True)

    treffer = {}
    with mp.Pool(int(os.environ.get('NPROC', os.cpu_count() or 24))) as pool:
        n = 0
        for cand, v in pool.imap_unordered(_vol, kand, chunksize=50):
            n += 1
            if v is None:
                continue
            for z in ziele.values():
                if z['det'] != cand['det']:
                    continue
                dv = abs(v - z['vol'])
                if dv < PREFILTER and z['id'] not in treffer:
                    iso = _isometrisch(z['isosig'], cand['dt'])
                    if iso is True:
                        treffer[z['id']] = {
                            'unser_id': z['id'], 'stem': z['stem'], 'det': z['det'],
                            'unser_vol': z['vol'], 'unser_isosig': z['isosig'],
                            'c_ub_bisher': z['c_ub'], 'klasse_bisher': z['klasse'],
                            'n_realisierungen': z['n_real'],
                            'burton_name': cand['name'], 'burton_tabelle': cand['crossing'],
                            'burton_dt': cand['dt'], 'burton_sig': cand.get('sig'),
                            'burton_vol': v, 'vol_differenz': dv,
                            'bestaetigt_durch': 'is_isometric_to'}
                        print(f"  ✅ {z['id']:<32} == {cand['name']} ({cand['crossing']})  dvol={dv:.2e}", flush=True)
            if n % 20000 == 0:
                print(f"    {n}/{len(kand)} candidates, {len(treffer)} confirmed, {time.time()-t0:.0f}s", flush=True)

    ohne = [i for i in ids if i not in treffer]
    gegen_getroffen = [i for i in gegenprobe_ids if i in treffer]
    erg = {'gefahren_am': time.strftime('%Y-%m-%d %H:%M'),
           'quelle': 'candidates_stageA.json (Stage A 10.07.2026, 17a/17n/18a/18n)',
           'bestaetigt': list(treffer.values()),
           'nicht_bestaetigt': ohne,
           'gegenprobe_ids': list(gegenprobe_ids),
           'gegenprobe_faelschlich_getroffen': gegen_getroffen}
    json.dump(erg, open(ausgabe, 'w'), indent=1)
    print(f"\n=== RESULT ===")
    print(f"  confirmed (census hit, isometry-receipted): {len(treffer)}/{len(ids)}")
    print(f"  NOT confirmed: {len(ohne)}  {ohne}")
    print(f"  Counter-check ({len(gegenprobe_ids)} targets that must NOT be in the census): "
          f"{'🔴 ' + str(gegen_getroffen) if gegen_getroffen else '✅ none hit'}")
    print(f"  written: {ausgabe}")
    return erg


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--ids', nargs='+', required=True)
    ap.add_argument('--gegenprobe', nargs='*', default=[])
    ap.add_argument('--ausgabe', default=WORK + '/burton_beleg_20260904.json')
    a = ap.parse_args()
    lauf(a.ids, a.ausgabe, a.gegenprobe)
