#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""tighten_65.py — gnottheorie/snappy. (2026-09-06, pattern: tighten_179.py of 10.07.)

The c19 run left 65 UNRESOLVED targets 'unrigorous': their det buckets contain
Census candidates with non-geometric SnapPy solution (35,392 fail volumes)
— there the absence is NOT proven. Here EVERY candidate of these 54 dets
is touched again, fails with randomized retriangulation (8 attempts).
  - new geometric vol + target match (is_isometric_to) => retraction.
  - det bucket without remaining fails and without match => target RIGOROUSLY beyond-19 (c>=20).
Deviations from the July pattern:
  - candidates_stageA.jsonl is STREAMED (2.7 GB; footgun 17 — no full load).
  - BASE-SET LATCH: counted number of candidates must match the independent count
    of 14:0x (14,281,808) — otherwise the run checks something other than intended.
  - COUNTER-CHECK ALARM: if a span>=21 target matches (shares a det), that is a
    proof contradiction => loud block in the result, exit 3 (as stage_b_v4).
Call: gnottheorie_env python directly, NEVER conda run.
"""
import site_config
import sys, json, time
import multiprocessing as mp
sys.path.insert(0, (site_config.MACHINE_YARD + '/corpus_vs_maguire/work'))
import knotlib as K

WORK = (site_config.MACHINE_YARD + '/corpus_vs_burton/work_c19_20260904')
GEOM = 'all tetrahedra positively oriented'
TRIES = 8
ERWARTETE_GRUNDMENGE = 14281808  # counted independently 06.09. 14:0x (jsonl stream, 54 dets)


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
        return (cand['det'], cand['crossing'], cand['dt'], None, True)
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


def cands_stream(dets):
    with open(WORK + '/candidates_stageA.jsonl') as f:
        for line in f:
            o = json.loads(line)
            if o['det'] in dets:
                yield o


def main():
    res = json.load(open(WORK + '/stage_b_result.json'))
    geg = set(json.load(open(WORK + '/gegenprobe.json'))['gegenprobe_ids'])
    unrig = [u for u in res['detail']['beyond_unrigorous_skips'] if u['id'] not in geg]
    unrig_dets = {u['det'] for u in unrig}
    d2t = json.load(open(WORK + '/det_to_targets.json'))
    print(f"{len(unrig)} unresolved targets, {len(unrig_dets)} dets, "
          f"expected base set {ERWARTETE_GRUNDMENGE}", flush=True)

    still_fail, recovered = {}, 0
    prefilter_hits = []
    t0 = time.time(); n = 0
    with mp.Pool(24) as pool:
        for det, cross, dt, v, was_nongeom in pool.imap_unordered(
                _worker, cands_stream(unrig_dets), chunksize=200):
            n += 1
            if v is None:
                still_fail[det] = still_fail.get(det, 0) + 1
            else:
                if was_nongeom:
                    recovered += 1
                for tgt in d2t.get(str(det), []):
                    if abs(v - tgt['vol']) < 1e-3:
                        prefilter_hits.append((tgt['id'], tgt.get('isosig'), cross, dt))
            if n % 200000 == 0:
                print(f"  re-vol {n} recovered={recovered} "
                      f"still_fail={sum(still_fail.values())} hits={len(prefilter_hits)} "
                      f"{time.time()-t0:.0f}s", flush=True)
    assert n == ERWARTETE_GRUNDMENGE, (
        f"GRUNDMENGEN-RIEGEL: gestreamt {n} != unabhaengig gezaehlt {ERWARTETE_GRUNDMENGE}")
    print(f"re-vol done: {n} candidates, recovered {recovered} nongeom->geom, "
          f"{sum(still_fail.values())} still non-geom, {len(prefilter_hits)} prefilter hits",
          flush=True)

    matched, geg_verletzt = {}, {}
    for tid, iso, cross, dt in sorted(prefilter_hits, key=lambda h: h[2]):
        ziel_map = geg_verletzt if tid in geg else matched
        if tid in ziel_map:
            continue
        if iso and _confirm(iso, dt) is True:
            ziel_map[tid] = {'crossing': cross, 'dt': dt}

    out = {'now_downgrade': [], 'now_rigorous': [], 'still_unrigorous': []}
    for u in unrig:
        tid, det = u['id'], u['det']
        if tid in matched:
            out['now_downgrade'].append({'id': tid, 'det': det, **matched[tid]})
        elif still_fail.get(det, 0) > 0:
            out['still_unrigorous'].append({'id': tid, 'det': det,
                                            'still_nongeom': still_fail[det]})
        else:
            out['now_rigorous'].append({'id': tid, 'det': det})
    json.dump({'summary': {k: len(v) for k, v in out.items()},
               'detail': out,
               'gegenprobe_verletzt': geg_verletzt,
               'recovered_nongeom': recovered,
               'still_nongeom_total': sum(still_fail.values()),
               'still_fail_by_det': still_fail,
               'grundmenge': n},
              open(WORK + '/tighten_65_result_20260906.json', 'w'), indent=1)
    print("\n=== 65-TIGHTENING RESULT ===", flush=True)
    print("  ", {k: len(v) for k, v in out.items()}, flush=True)
    if geg_verletzt:
        print(f"🔴 COUNTER-CHECK VIOLATED — span>=21 target(s) confirmed in the census: "
              f"{sorted(geg_verletzt)} — this refutes a PROOF. FULL STOP.", flush=True)
        sys.exit(3)
    print("✅ no counter-check violation", flush=True)


if __name__ == '__main__':
    main()
