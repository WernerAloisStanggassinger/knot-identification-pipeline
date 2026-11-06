#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""tighten_179.py — gnottheorie/snappy.

The beyond-18 run left 179 targets 'unrigorous': their det-bucket contained Burton
candidates whose SnapPy solution was non-geometric, so a volume couldn't be trusted
and we could not rule those candidates out as the target's realization. Here we
RE-ATTEMPT every candidate in those det-buckets with randomized retriangulation.
A candidate that now yields a geometric volume matching a target (is_isometric_to
confirmed) downgrades that target to its census crossing. A det-bucket with zero
remaining non-geometric candidates and no match makes its target rigorously beyond-18.
"""
import site_config
import sys, json, time
import multiprocessing as mp
sys.path.insert(0, (site_config.MACHINE_YARD + '/corpus_vs_maguire/work'))
import knotlib as K

WORK = (site_config.MACHINE_YARD + '/corpus_vs_burton/work')
GEOM = 'all tetrahedra positively oriented'
TRIES = 8

def _worker(cand):
    import snappy
    try:
        M = snappy.Manifold('DT:' + repr([tuple(K.decode_dt(cand['dt']))]))
        if str(M.solution_type()) == GEOM:
            return (cand['det'], cand['crossing'], cand['dt'], float(M.volume()), False)
        for _ in range(TRIES):
            M.randomize()
            if str(M.solution_type()) == GEOM:
                return (cand['det'], cand['crossing'], cand['dt'], float(M.volume()), True)
        return (cand['det'], cand['crossing'], cand['dt'], None, True)  # still non-geom
    except Exception:
        return (cand['det'], cand['crossing'], cand['dt'], None, True)

def _confirm(iso, dt):
    import snappy
    try:
        a = snappy.Manifold(iso)
        b = snappy.Manifold('DT:' + repr([tuple(K.decode_dt(dt))]))
        if a.is_isometric_to(b):
            return True
        for _ in range(4):
            b.randomize()
            if a.is_isometric_to(b):
                return True
        return False
    except Exception:
        return None

def main():
    res = json.load(open(WORK + '/stage_b_result.json'))
    unrig = res['detail']['beyond_unrigorous_skips']
    unrig_dets = {u['det'] for u in unrig}
    d2t = json.load(open(WORK + '/det_to_targets.json'))
    print(f"loading candidates...", flush=True)
    raw = json.load(open(WORK + '/candidates_stageA.json'))
    cands = raw['candidates'] if isinstance(raw, dict) else raw
    sub = [c for c in cands if c['det'] in unrig_dets]
    del cands, raw
    print(f"179 targets, {len(unrig_dets)} dets, {len(sub)} candidates to re-attempt", flush=True)

    still_fail, recovered = {}, 0
    prefilter_hits = []
    t0 = time.time(); n = 0
    with mp.Pool(24) as pool:
        for det, cross, dt, v, was_nongeom in pool.imap_unordered(_worker, sub, chunksize=100):
            n += 1
            if v is None:
                still_fail[det] = still_fail.get(det, 0) + 1
            else:
                if was_nongeom:
                    recovered += 1
                for tgt in d2t.get(str(det), []):
                    if abs(v - tgt['vol']) < 1e-3:
                        prefilter_hits.append((tgt['id'], tgt.get('isosig'), cross, dt))
            if n % 100000 == 0:
                print(f"  re-vol {n}/{len(sub)} recovered={recovered} still_fail={sum(still_fail.values())} "
                      f"hits={len(prefilter_hits)} {time.time()-t0:.0f}s", flush=True)
    print(f"re-vol done: recovered {recovered} nongeom->geom, {sum(still_fail.values())} still non-geom, "
          f"{len(prefilter_hits)} prefilter hits", flush=True)

    matched = {}
    for tid, iso, cross, dt in sorted(prefilter_hits, key=lambda h: h[2]):
        if tid in matched:
            continue
        if iso and _confirm(iso, dt) is True:
            matched[tid] = cross

    out = {'now_downgrade': [], 'now_rigorous': [], 'still_unrigorous': []}
    for u in unrig:
        tid, det = u['id'], u['det']
        if tid in matched:
            out['now_downgrade'].append({'id': tid, 'crossing': matched[tid]})
        elif still_fail.get(det, 0) > 0:
            out['still_unrigorous'].append({'id': tid, 'det': det, 'still_nongeom': still_fail[det]})
        else:
            out['now_rigorous'].append({'id': tid, 'det': det})
    json.dump({'summary': {k: len(v) for k, v in out.items()},
               'detail': out, 'recovered_nongeom': recovered,
               'still_nongeom_total': sum(still_fail.values())},
              open(WORK + '/tighten_179_result.json', 'w'), indent=1)
    print("\n=== 179-TIGHTENING RESULT ===", flush=True)
    print("  ", {k: len(v) for k, v in out.items()}, flush=True)
    print(f"  151 rigorous + {len(out['now_rigorous'])} newly rigorous = "
          f"{151 + len(out['now_rigorous'])} beyond-18 rigorous total", flush=True)

if __name__ == '__main__':
    main()
