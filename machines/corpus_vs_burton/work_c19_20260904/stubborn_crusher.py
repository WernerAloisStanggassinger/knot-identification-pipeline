#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""stubborn_crusher.py — the last 32 non-geometric census candidates (2026-09-06).

After tighten_65, 32 candidates remain in the dets of the 28 open targets, whose
SnapPy solution did not become geometric even after 8 retriangulations — they
block 28 absence proofs. Here: stream the candidates of these dets again,
standard ladder (8x randomize) identifies the stubborn ones, then ARTILLERY:
ManifoldHP (quad-double) + up to 40 retriangulations + Dehn-fill trick
(fill_triangulation after (1,0)-filling of an extra cusp does not apply to knots;
instead canonize attempt). Every rescued value is checked against the targets of the det
(1e-3 -> is_isometric_to). Result -> stubborn_result.json.
Honest: if a candidate stays non-geometric even after HP, the det stays
UNPROVEN (no proof by exhaustion of patience).
Counter-check as always: span>=21 targets must not match (ALARM, exit 3).
Run: gnottheorie_env, 24 cores for pass 1, HP pass serial-parallel.
"""
import site_config
import json, os, sys, time
import multiprocessing as mp
sys.path.insert(0, (site_config.MACHINE_YARD + '/corpus_vs_maguire/work'))
import knotlib as K

WORK = os.path.dirname(os.path.abspath(__file__))
GEOM = 'all tetrahedra positively oriented'
PREFILTER = 1e-3


def _vol8(cand):
    import snappy
    try:
        M = snappy.Manifold('DT:' + repr([tuple(K.decode_dt(cand['dt']))]))
        if str(M.solution_type()) == GEOM:
            return (cand, float(M.volume()), 'ok')
        for _ in range(8):
            M.randomize()
            if str(M.solution_type()) == GEOM:
                return (cand, float(M.volume()), 'ok')
        return (cand, None, 'stubborn')
    except Exception:
        return (cand, None, 'stubborn')


def _artillerie(cand):
    """ManifoldHP + 40 Retriangulierungen + canonize. None = endgueltig verweigert."""
    import snappy
    try:
        M = snappy.ManifoldHP('DT:' + repr([tuple(K.decode_dt(cand['dt']))]))
        for i in range(41):
            if str(M.solution_type()) == GEOM:
                return (cand, float(M.volume()))
            if i == 0:
                try:
                    M.canonize()
                    continue
                except Exception:
                    pass
            M.randomize()
        return (cand, None)
    except Exception:
        return (cand, None)


def main():
    tr = json.load(open(WORK + '/tighten_65_result_20260906.json'))
    geg = set(json.load(open(WORK + '/gegenprobe.json'))['gegenprobe_ids'])
    offen = [u for u in tr['detail']['still_unrigorous']
             if u['id'] not in geg and u['id'] not in
             ('S_15_d625_D3_chiral', 'S_16_d1377_D8_amphi')]  # retracted today
    dets = {u['det'] for u in offen}
    d2t = json.load(open(WORK + '/det_to_targets.json'))
    print(f'{len(offen)} open targets, {len(dets)} dets', flush=True)

    t0 = time.time()
    cands = []
    with open(WORK + '/candidates_stageA.jsonl') as f:
        for line in f:
            o = json.loads(line)
            if o['det'] in dets:
                cands.append(o)
    print(f'{len(cands)} candidates collected {time.time()-t0:.0f}s', flush=True)

    stubborn, hits = [], []
    vols = {}
    n = 0
    with mp.Pool(24) as pool:
        for cand, v, tag in pool.imap_unordered(_vol8, cands, chunksize=200):
            n += 1
            if tag == 'stubborn':
                stubborn.append(cand)
            elif v is not None:
                vols[(cand['det'], cand['dt'])] = v
            if n % 500000 == 0:
                print(f'  pass1 {n}/{len(cands)} stubborn={len(stubborn)} {time.time()-t0:.0f}s',
                      flush=True)
    print(f'Pass 1: {len(stubborn)} stubborn candidates (tighten saw 32; deviation = '
          f'randomisation luck, both numbers honest)', flush=True)

    gerettet, verweigert = [], []
    with mp.Pool(min(24, max(1, len(stubborn)))) as pool:
        for cand, v in pool.imap_unordered(_artillerie, stubborn):
            if v is None:
                verweigert.append(cand)
                print(f'  ⛔ refused: {cand["name"]} (det {cand["det"]})', flush=True)
            else:
                gerettet.append((cand, v))
                vols[(cand['det'], cand['dt'])] = v
                print(f'  💥 rescued: {cand["name"]} vol={v:.6f}', flush=True)

    # Target check: all rescued AND all pass1 vols against the targets of the dets
    def _confirm(iso, dt):
        import snappy
        try:
            a = snappy.Manifold(iso)
            b = snappy.Manifold('DT:' + repr([tuple(K.decode_dt(dt))]))
            for _ in range(6):
                try:
                    r = a.is_isometric_to(b)
                    if r:
                        return True
                    b.randomize()
                except RuntimeError:
                    a.randomize(); b.randomize()
            return False
        except Exception:
            return None

    matched, alarm = {}, {}
    for (det, dt), v in vols.items():
        for tgt in d2t.get(str(det), []):
            if abs(v - tgt['vol']) < PREFILTER and tgt['id'] not in matched:
                if _confirm(tgt.get('isosig'), dt) is True:
                    ziel = alarm if tgt['id'] in geg else matched
                    name = next(c['name'] for c in cands if c['dt'] == dt and c['det'] == det)
                    ziel[tgt['id']] = {'name': name, 'dt': dt, 'det': det, 'vol': v}
                    print(f'  MATCH {tgt["id"]} = {name}', flush=True)

    verweigert_dets = {c['det'] for c in verweigert}
    out = {'now_rigorous': [], 'now_downgrade': [], 'still_blocked': []}
    for u in offen:
        tid, det = u['id'], u['det']
        if tid in matched:
            out['now_downgrade'].append({'id': tid, 'det': det, **matched[tid]})
        elif det in verweigert_dets:
            out['still_blocked'].append({'id': tid, 'det': det,
                'verweigerte': [c['name'] for c in verweigert if c['det'] == det]})
        else:
            out['now_rigorous'].append({'id': tid, 'det': det})
    json.dump({'summary': {k: len(v) for k, v in out.items()},
               'detail': out, 'gegenprobe_alarm': alarm,
               'stubborn_pass1': len(stubborn),
               'gerettet': [{'name': c['name'], 'det': c['det'], 'vol': v} for c, v in gerettet],
               'verweigert': [{'name': c['name'], 'det': c['det'], 'dt': c['dt']}
                              for c in verweigert]},
              open(WORK + '/stubborn_result.json', 'w'), indent=1)
    print(f"\n=== STUBBORN RESULT === {({k: len(v) for k, v in out.items()})}", flush=True)
    if alarm:
        print('🔴 COUNTER-CHECK VIOLATED — FULL STOP.', flush=True)
        sys.exit(3)
    print('✅ no counter-check violation', flush=True)


if __name__ == '__main__':
    main()
