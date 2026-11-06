#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""burton_receipt_c19.py — fetch the receipt for the 7 stage-B hits of 06.09.2026,
BEFORE burton_retractions.json is extended.

OCCASION: the c19 run (work_c19_20260904, stage B v4, EXIT=0) found 7 targets in the
Burton census (1x 18n + 6x 19a/19n). `stage_b_result.json` stores per
hit only {id, det, crossing} — NOT which census knot was hit.
A retraction is a statement of truth; it needs census name + DT + both
volumes + is_isometric_to, re-checkable by a human (pattern: burton_receipt.py
of 04.09., there for the July run in the old JSON format; here: JSONL stream).

WHAT IT DOES: stream candidates_stageA.jsonl once, collect candidates of the target dets,
compute volumes on all cores, volume proximity (1e-3) -> MANDATORY
is_isometric_to against the target isosig. Output: one receipt file per target with
census name, DT, sig, both volumes and the difference.

--eichen runs the counter-cases (footgun 16: a tool that only confirms
is worthless):
  G1: counter-check T_11_3_2_10 (span>=21 PROVEN, det=19) — must give 0 hits,
      although the det-19 base set in the census is large.
  G2: pz_3_5_11 with deliberately falsified volume (+0.5) — must give 0 hits.
  G3: an empty base set is reported as an ERROR (exit 3), never as "no hit".

Call (gnottheorie_env directly, NEVER conda run):
  python -u burton_receipt_c19.py [--eichen]      (an env with snappy + regina, see pipeline/site_config.py)
"""
import site_config
import json, os, sys, time
import multiprocessing as mp

WORK = (site_config.MACHINE_YARD + '/corpus_vs_burton/work_c19_20260904')
OUT = os.path.join(WORK, 'beleg_tighten3_20260906.json' if '--tighten' in sys.argv
                   else 'beleg_7_treffer_20260906.json')
PREFILTER = 1e-3
GEOM = 'all tetrahedra positively oriented'
sys.path.insert(0, (site_config.MACHINE_YARD + '/corpus_vs_maguire/work'))
import knotlib as K  # noqa: E402


def lade_ziele(eichen):
    res = json.load(open(os.path.join(WORK, 'stage_b_result.json')))
    tgt = {t['id']: t for t in json.load(open(os.path.join(WORK, 'beyond17_targets.json')))}
    ziele = []
    if '--tighten' in sys.argv:
        # receipt for the downgrades of the tighten run (06.09.) instead of the stage-B hits
        tr = json.load(open(os.path.join(WORK, 'tighten_65_result_20260906.json')))
        for e in tr['detail']['now_downgrade']:
            t = tgt[e['id']]
            ziele.append({'id': t['id'], 'det': t['det'], 'vol': t['vol'],
                          'isosig': t['isosig'], 'erwartet_crossing': e['crossing']})
        return ziele
    for bucket in ('matched_c18_downgrade', 'matched_c19_downgrade'):
        for e in res['detail'][bucket]:
            t = tgt[e['id']]
            ziele.append({'id': t['id'], 'det': t['det'], 'vol': t['vol'],
                          'isosig': t['isosig'], 'erwartet_crossing': e['crossing']})
    if eichen:
        d2t = json.load(open(os.path.join(WORK, 'det_to_targets.json')))
        g = next(t for t in d2t['19'] if t['id'] == 'T_11_3_2_10')
        ziele = [
            dict(id='EICH_G1_T_11_3_2_10', det=19, vol=g['vol'], isosig=g['isosig'],
                 erwartet_crossing=None),
            dict(id='EICH_G2_pz_3_5_11_vol_verfaelscht', det=103,
                 vol=tgt['pz_3_5_11']['vol'] + 0.5, isosig=tgt['pz_3_5_11']['isosig'],
                 erwartet_crossing=None),
        ]
    return ziele


def _vol(cand):
    import snappy
    try:
        M = snappy.Manifold('DT:' + repr([tuple(K.decode_dt(cand['dt']))]))
        if str(M.solution_type()) == GEOM:
            return (cand, float(M.volume()))
        # as in tighten_65: rescue a non-geometric solution by retriangulation,
        # otherwise the receipt misses exactly the candidates that tighten found
        for _ in range(8):
            M.randomize()
            if str(M.solution_type()) == GEOM:
                return (cand, float(M.volume()))
        return (cand, None)
    except Exception:
        return (cand, None)


def _isometrisch(ziel_isosig, cand_dt):
    """MANDATORY confirmation; retriangulation up to 6 attempts, None = undecidable."""
    import snappy
    try:
        A = snappy.Manifold(ziel_isosig)
        B = snappy.Manifold('DT:' + repr([tuple(K.decode_dt(cand_dt))]))
        for _ in range(6):
            try:
                return A.is_isometric_to(B)
            except RuntimeError:
                A.randomize(); B.randomize()
    except Exception:
        pass
    return None


def main():
    eichen = '--eichen' in sys.argv
    ziele = lade_ziele(eichen)
    dets = {z['det'] for z in ziele}
    print(f"{len(ziele)} targets, dets={sorted(dets)}", flush=True)

    t0 = time.time()
    kands = {d: [] for d in dets}
    with open(os.path.join(WORK, 'candidates_stageA.jsonl')) as f:
        for line in f:
            o = json.loads(line)
            if o['det'] in dets:
                kands[o['det']].append(o)
    for d in sorted(dets):
        print(f"  det={d}: {len(kands[d])} census candidates", flush=True)
        if not kands[d]:
            print(f"🔴 BASE SET EMPTY for det={d} — this is an error, not a finding (G3).",
                  flush=True)
            sys.exit(3)

    alle = [c for d in dets for c in kands[d]]
    vols = {}
    n = 0
    with mp.Pool(24) as pool:
        for cand, v in pool.imap_unordered(_vol, alle, chunksize=200):
            n += 1
            if v is not None:
                vols[(cand['det'], cand['dt'])] = (cand, v)
            if n % 200000 == 0:
                print(f"  vol {n}/{len(alle)} {time.time()-t0:.0f}s", flush=True)
    print(f"Volumes: {len(vols)}/{len(alle)} geometric, {time.time()-t0:.0f}s", flush=True)

    beleg, problems = {}, []
    for z in ziele:
        treffer = []
        for (det, dt), (cand, v) in vols.items():
            if det != z['det'] or abs(v - z['vol']) >= PREFILTER:
                continue
            iso = _isometrisch(z['isosig'], dt)
            if iso is True:
                treffer.append({'burton_name': cand['name'], 'burton_dt': dt,
                                'burton_sig': cand['sig'], 'det': det,
                                'ziel_vol': z['vol'], 'burton_vol': v,
                                'vol_diff': abs(v - z['vol']),
                                'crossing': cand['crossing'],
                                'confirmed_by': 'is_isometric_to'})
            elif iso is None:
                problems.append({'ziel': z['id'], 'burton_name': cand['name'],
                                 'note': 'vol-close, isometry undecidable'})
        beleg[z['id']] = treffer
        print(f"  {z['id']}: {len(treffer)} confirmed census hits "
              f"{[t['burton_name'] for t in treffer]}", flush=True)

    if eichen:
        ok = all(len(v) == 0 for v in beleg.values())
        print(("✅ CALIBRATION PASSED: both counter-cases 0 hits (G3 implicitly green)"
               if ok else "🔴 CALIBRATION FAILED: counter-case has hits!"), flush=True)
        sys.exit(0 if ok else 1)

    fehl = [z['id'] for z in ziele if not beleg[z['id']]]
    json.dump({'_meta': {'date': '2026-09-06', 'lauf': 'work_c19_20260904 stage B v4',
                         'quelle_kandidaten': 'candidates_stageA.jsonl (Stage A 04.-06.09., '
                         '350.449.955 knots, all 6 categories exact, detfails=0)',
                         'methode': 'det filter -> volume 1e-3 -> MANDATORY is_isometric_to'},
               'beleg': beleg, 'unentscheidbar': problems},
              open(OUT, 'w'), indent=1)
    print(f"-> {OUT}", flush=True)
    if fehl:
        print(f"🔴 REMAINED WITHOUT RECEIPT: {fehl} — do NOT adopt into retractions!", flush=True)
        sys.exit(2)
    print("✅ all targets receipted", flush=True)


if __name__ == '__main__':
    main()
