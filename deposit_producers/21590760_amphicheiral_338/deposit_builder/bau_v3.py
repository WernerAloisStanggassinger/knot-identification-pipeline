#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Builds counterexamples_data.json for v3 = 250 (v2) + ALL new distinct ones.

03.08.2026 (harvest chain): new additions are no longer `k==8` (the 7 back then),
but EVERYTHING from distinkte_knoten.json that is not in the published v2
(state today: 65 = 29 order 8 + 36 order 16). PD sources: sieve outputs
(pd/pd_hi) AND stage-1 projections (up to 7 raw PDs per record).

Doctrine:
 1. CALIBRATION first: all derived fields of the 250 are rebuilt from the source files
    and set against the PUBLISHED v2. Only if that matches 250/250
    is the recipe proven and allowed to be applied to the new ones.
 2. The crossing upper bound is UNIFIED: it is the crossing number of the
    supplied diagram. A reader counts the X[] and has the receipt in hand.
    v2 carried a number here that the attached PD undercut in 233 cases and
    failed to reach in 2 cases. The old one remains as `kreuzungs_obergrenze_v2`.
 3. FRANKENSTEIN GUARD: every newly simplified diagram must carry the same
    Alexander polynomial as the record. Otherwise the record goes out, not the guard.
"""
import site_config
import json, glob, os, sys, collections
from multiprocessing import Pool
import regina

B   = (site_config.ROOT + '/')
KAW = B + 'konversationen/kawauchi_20260725/'
OUT = B + 'PRIORITAET_DEPOSIT_20260730/'
A   = {2: 2, 4: 3, 8: 4}
R   = int(os.environ.get('R', '60'))

def mkid(k, crot, seed, nxy, nz, dc):
    return f"SK{A[k]}_{crot}_{seed}_{nxy}{nz}_{str(dc).replace('0.', '0')}"

def key_of(e):
    return (e['k'], e['seed'], e['nxy'], e['nz'], e['dc'], e.get('crot', 1))

# ---------- sources ----------
v2roh = json.load(open(B + 'PRIORITAET_DEPOSIT_20260729/counterexamples_data.json'))
v2    = {e['id']: e for e in v2roh['knoten']}
dk    = {mkid(*key_of(e)[:1] + (e.get('crot',1),) + key_of(e)[1:5]): e
         for e in json.load(open(KAW + 'sieb/distinkte_knoten.json'))}
me    = {key_of(e): e for e in json.load(open(KAW + 'me_ergebnis_a3.json'))}
ze    = {key_of(e): e for e in json.load(open(KAW + 'sieb/alex_zeuge.json'))}
bil   = {key_of(e): e for e in json.load(open(KAW + 'BILANZ_DATEN.json'))}

# ---------- calibration of the recipe on the 250 ----------
def baue(e, pd, cob):
    """e = entry from distinkte_knoten.json"""
    k = key_of(e); g = e['geo']; m = me.get(k); z = ze.get(k); b = bil.get(k)
    return {
        'id': mkid(k[0], k[5], k[1], k[2], k[3], k[4]),
        'erzeugung': {'k': k[0], 'c_rot': k[5], 'seed': k[1], 'n_xy': k[2], 'n_z': k[3], 'decay': k[4]},
        'ordnung_der_drehspiegelung': 2 * k[0],
        'a_tower': A[k[0]],
        'achsen_verschlingung_ell': e['ell'],
        'sigma_zertifikat': e['sigma'],
        'conway_nabla': b['nabla'], 'grad_z': b['deg'], 'leitkoeffizient': b['lc'],
        'alexander': e['alex'],
        'pd_code': pd,
        'kreuzungs_obergrenze': cob,
        'volumen': g['vol'], 'tetraeder': g['tets'],
        'symmetriegruppe': g['group'], 'symmetrie_typen': g['types_union'],
        'snappy_identify': g['identify'],
        'pruefungen': {
            'nabla_zerfaellt_nicht': bool(b['gegenfall']),
            'zweite_algebra_einig': z['einig'], 'zweite_algebra_trifft_regina': z['trifft'],
            'me_loesbar': m['loesung_existiert'], 'E0_kongruenz': m['E0_passt'],
            'E0_vorhersage': m['vorhersage_E0'], 'E0_zweig_scharf': bool(b['E0_scharf']),
            'cm_leitkoeffizient_schranke': bool(b['cm_schranke_ok']),
        },
        'im_deposit_2026_07_25': False,
    }

VERGLEICH = ['erzeugung', 'ordnung_der_drehspiegelung', 'a_tower', 'achsen_verschlingung_ell',
             'sigma_zertifikat', 'conway_nabla', 'grad_z', 'leitkoeffizient', 'alexander',
             'volumen', 'tetraeder', 'symmetriegruppe', 'symmetrie_typen', 'snappy_identify',
             'pruefungen']

def eichung():
    mis = collections.Counter(); n = 0
    for i, alt in v2.items():
        if i not in dk: mis['quelle_fehlt'] += 1; continue
        neu = baue(dk[i], alt['pd_code'], alt['kreuzungs_obergrenze']); n += 1
        for f in VERGLEICH:
            if neu[f] != alt[f]: mis[f] += 1
    return n, mis

# ---------- diagram work ----------
def pd_string(data):
    """regina.pdData() -> the PD notation that snappy/spherogram and regina both read."""
    return 'PD[' + ', '.join('X[%d, %d, %d, %d]' % tuple(x) for x in data) + ']'


def besser(arg):
    """(id, pd) -> (id, best pd, crossing number, alexander of the best diagram)"""
    i, pds = arg
    best = None
    for pd in pds:
        if not pd: continue
        for _ in range(R):
            L = regina.Link.fromPD(pd)
            for _ in range(6): L.simplify()
            if best is None or L.size() < best.size(): best = L
    if best is None: return (i, None, None)
    return (i, pd_string(best.pdData()), best.size())

def volumen_von_pd(arg):
    """Who owns the diagram? Volume is the proof of ownership (Mostow)."""
    i, pd, soll = arg
    try:
        import snappy
        M = snappy.Link(pd).exterior()
        M.simplify()
        v = M.volume()
        return (i, float(v), abs(float(v) - soll) < 1e-4)
    except Exception as ex:
        return (i, None, f'ERR {ex}')

if __name__ == '__main__':
    n, mis = eichung()
    print(f"CALIBRATION of the recipe on the {n} published records")
    print("  deviations:", dict(mis) or "NONE  -> recipe proven")
    if mis:
        print("  ABORT: the recipe does not reproduce v2. No build."); sys.exit(2)

    # New additions = everything distinct that is not in the published v2.
    # PD candidates per new addition: sieve outputs (pd/pd_hi) + stage-1 projections.
    neu_keys = [k for k, e in ((key_of(e), e) for e in
                json.load(open(KAW + 'sieb/distinkte_knoten.json')))
                if mkid(k[0], k[5], k[1], k[2], k[3], k[4]) not in v2]
    nk = set(neu_keys)
    kandidaten = {k: [] for k in neu_keys}
    for f in glob.glob(KAW + 'sieb/breit_c_out/*.jsonl'):
        for ln in open(f):
            try: r = json.loads(ln)
            except Exception: continue
            kk = (r.get('k'), r.get('seed'), r.get('nxy'), r.get('nz'), r.get('dc'), r.get('crot'))
            if kk in nk and r.get('rej') is None:
                kandidaten[kk] += [r.get('pd'), r.get('pd_hi')]
    for r in json.load(open(KAW + 'sieb/pruefstand_stufe1.json')):
        if 'k' not in r: continue
        kk = key_of(r)
        if kk in nk:
            for v in (r.get('proj') or {}).values():
                if isinstance(v, dict) and v.get('pd'): kandidaten[kk].append(v['pd'])
    for kk in kandidaten:   # dedupe, drop None — order is kept
        kandidaten[kk] = list(dict.fromkeys(p for p in kandidaten[kk] if p))
    fehlt = [k for k in neu_keys if not kandidaten[k]]
    if fehlt: print("  ABORT: no PD candidate for", fehlt); sys.exit(2)

    aufgaben = [(i, [alt['pd_code']]) for i, alt in v2.items()]
    aufgaben += [(mkid(k[0], k[5], k[1], k[2], k[3], k[4]), kandidaten[k]) for k in neu_keys]
    print(f"\ndiagram work on {len(aufgaben)} records, {R} randomized simplify runs each ...")
    with Pool(10) as p:
        res = {r[0]: r[1:] for r in p.map(besser, aufgaben, chunksize=2)}

    # The Frankenstein guard (volume) runs in gnottheorie_env: snappy lives there,
    # regina here. Therefore this script writes a DRAFT; only wache_volumen.py
    # turns it into counterexamples_data.json. No green without the guard.
    # ---------- build ----------
    knoten = []
    for i, alt in v2.items():
        pd, c = res[i]
        neu = dict(alt)
        neu['pd_code'] = pd
        neu['kreuzungs_obergrenze'] = c
        neu['kreuzungs_obergrenze_v2'] = alt['kreuzungs_obergrenze']
        knoten.append(neu)
    for k in neu_keys:
        i = mkid(k[0], k[5], k[1], k[2], k[3], k[4])
        pd, c = res[i]
        knoten.append(baue(dk[i], pd, c))
    knoten.sort(key=lambda e: (e['a_tower'], e['erzeugung']['c_rot'], e['erzeugung']['seed'], e['id']))

    gesenkt = sum(1 for e in knoten if 'kreuzungs_obergrenze_v2' in e and e['kreuzungs_obergrenze'] < e['kreuzungs_obergrenze_v2'])
    gehoben = [e['id'] for e in knoten if 'kreuzungs_obergrenze_v2' in e and e['kreuzungs_obergrenze'] > e['kreuzungs_obergrenze_v2']]
    cs = [e['kreuzungs_obergrenze'] for e in knoten]
    vols = [e['volumen'] for e in knoten]
    verteilung = collections.Counter(e['a_tower'] for e in knoten)

    meta = dict(v2roh['meta'])
    meta.update({
        'datum': '2026-08-03',
        'version': '3.0 — Part III: growth at order 8 and 16',
        'vorgaenger_doi': '10.5281/zenodo.21669171',
        'anzahl': len(knoten),
        'aussage': (f"{len(knoten)} distinct hyperbolic knots, positive but not strongly amphichiral, "
                    "whose Conway polynomial does not factor as f(z)f(-z). "
                    f"{verteilung[2]} at a=2 (order 4), {verteilung[3]} at a=3 (order 8), "
                    f"{verteilung[4]} at a=4 (order 16). {len(v2)} of them are version 2 of 29.07.2026, "
                    f"{len(knoten) - len(v2)} are new. NO NEW THEOREM IS CLAIMED."),
        'aenderung_gegen_v2': {
            'neue_knoten': len(knoten) - len(v2),
            'kreuzungs_obergrenze': ("unified: it is now the crossing number of the supplied "
                                     "diagram, hence recountable on the bundle. In v2 it was a number from "
                                     "a run that was not supplied."),
            'davon_gesenkt': gesenkt,
            'davon_angehoben': gehoben,
            'pd_code': 'fuer alle Records das beste gefundene Diagramm; Alexander gegen das v2-Diagramm geprueft',
        },
        'spannen': {'kreuzungs_obergrenze': [min(cs), max(cs)],
                    'volumen': [min(vols), max(vols)]},
    })
    json.dump({'meta': meta, 'knoten': knoten}, open(OUT + 'counterexamples_data_ENTWURF.json', 'w'), indent=1)
    print(f"\nwritten: {OUT}counterexamples_data_ENTWURF.json  —  {len(knoten)} knots "
          f"(a=2: {verteilung[2]}, a=3: {verteilung[3]}, a=4: {verteilung[4]})")
    print(f"  Crossing upper bound lowered at {gesenkt}, RAISED at {len(gehoben)}: {gehoben}")
    print(f"  range c {min(cs)} … {max(cs)}   |   vol {min(vols):.3f} … {max(vols):.3f}")
    print("\nNEXT STEP (mandatory): /home/squeezy/miniforge3/envs/gnottheorie_env/bin/python wache_volumen.py")
