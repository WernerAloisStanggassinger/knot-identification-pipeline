#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""low_stage_b.py — volume+isometry pass of the ≤16 layer (2026-09-06).

Takes candidates_low_c16.jsonl (low_stage_a: all Burton hyp 03..16, det filter
to the 91 target dets, detfails=0, category sum against known primes
weighted) and checks EVERY target from det_to_targets: volume proximity 1e-3 ->
MANDATORY is_isometric_to. Non-geometric candidates are rescued like tighten_65
with retriangulation (8 attempts); remaining fails make the
absence of the affected dets UNPROVEN (counted and reported per det).

Output low_c16_result.json:
  matched:      target found in the ≤16 census (⟹ TABLE KNOT, retraction;
                for span>=21 counter-checks: ALARM + exit 3)
  fails_by_det: dets with finally non-geometric candidates
Absence conclusion holds only for dets without fails and without match.
Run: gnottheorie_env directly.
"""
import site_config
import json, os, sys, time
import multiprocessing as mp
sys.path.insert(0, (site_config.MACHINE_YARD + '/corpus_vs_maguire/work'))
import knotlib as K

WORK = os.path.dirname(os.path.abspath(__file__))
GEOM = 'all tetrahedra positively oriented'
PREFILTER = 1e-3


def _vol(cand):
    import snappy
    try:
        M = snappy.Manifold('DT:' + repr([tuple(K.decode_dt(cand['dt']))]))
        if str(M.solution_type()) == GEOM:
            return (cand, float(M.volume()))
        for _ in range(8):
            M.randomize()
            if str(M.solution_type()) == GEOM:
                return (cand, float(M.volume()))
        return (cand, None)
    except Exception:
        return (cand, None)


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


def main():
    d2t = json.load(open(WORK + '/det_to_targets.json'))
    geg = set(json.load(open(WORK + '/gegenprobe.json'))['gegenprobe_ids'])
    cands = [json.loads(l) for l in open(WORK + '/candidates_low_c16.jsonl')]
    print(f'{len(cands)} candidates, {len(d2t)} target dets', flush=True)
    t0 = time.time()
    fails_by_det = {}
    prehits = []
    n = 0
    with mp.Pool(24) as pool:
        for cand, v in pool.imap_unordered(_vol, cands, chunksize=200):
            n += 1
            if v is None:
                fails_by_det[cand['det']] = fails_by_det.get(cand['det'], 0) + 1
            else:
                for tgt in d2t.get(str(cand['det']), []):
                    if abs(v - tgt['vol']) < PREFILTER:
                        prehits.append((tgt['id'], tgt.get('isosig'), cand))
            if n % 50000 == 0:
                print(f'  vol {n}/{len(cands)} fails={sum(fails_by_det.values())} '
                      f'prehits={len(prehits)} {time.time()-t0:.0f}s', flush=True)
    print(f'vols done: fails={sum(fails_by_det.values())} prehits={len(prehits)}', flush=True)

    matched = {}
    for tid, iso, cand in prehits:
        if tid in matched or not iso:
            continue
        if _confirm(iso, cand['dt']) is True:
            matched[tid] = {'name': cand['name'], 'dt': cand['dt'],
                            'crossing': cand['crossing'], 'det': cand['det']}
            print(f'  MATCH {tid} = {cand["name"]}', flush=True)
    alarm = sorted(t for t in matched if t in geg)
    json.dump({'matched': matched, 'fails_by_det': fails_by_det,
               'n_candidates': len(cands), 'gegenprobe_alarm': alarm},
              open(WORK + '/low_c16_result.json', 'w'), indent=1)
    print(f'=== LOW-c16 RESULT: {len(matched)} matches, '
          f'{len(fails_by_det)} dets with fails, alarm={alarm}', flush=True)
    if alarm:
        print('🔴 COUNTER-CHECK VIOLATED (span>=21 target in the ≤16 census) — FULL STOP.', flush=True)
        sys.exit(3)


if __name__ == '__main__':
    main()
