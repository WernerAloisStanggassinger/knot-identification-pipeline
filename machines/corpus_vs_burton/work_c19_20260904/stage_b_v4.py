#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Stage B (v4, 2026-09-05): like v3 (identity-safe, Disc #516: vol prefilter ->
is_isometric_to confirmation), adapted to the JSONL stage A of 04.09.:

  1. STREAMS candidates_stageA.jsonl (no 10-GB list load; footgun 17).
  2. COMPLETENESS GATE from the stage-A log (CVB_ALOG): all six
     "DONE <cat>: N knots (expected E)" lines must carry N==E and
     "STAGE A DONE" must be present — otherwise abort BEFORE the first volume
     (feedback: a result file is not a corpus; v3 checked this via a
     provenance manifest that does not exist in the JSONL format).
  3. COUNTER-CHECK LATCH: the span>=21 targets from gegenprobe.json MUST NOT
     match. A match there refutes not the knot, but one of our
     proofs -> loud 🔴 + own result field + exit 3.

Call (gnottheorie_env, directly -- NEVER conda run, error log 04.09.):
  CVB_WORK=<this folder> CVB_ALOG=<stage_a-log> \
    .../envs/gnottheorie_env/bin/python -u stage_b_v4.py
"""
import site_config
import sys, os, json, time, re
import multiprocessing as mp
sys.path.insert(0, os.environ.get('CVB_KNOTLIB', (site_config.MACHINE_YARD + '/corpus_vs_maguire/work')))
import knotlib as K

WORK = os.environ.get('CVB_WORK', (site_config.MACHINE_YARD + '/corpus_vs_burton/work_c19_20260904'))
ALOG = os.environ.get('CVB_ALOG', (site_config.ROOT + '/claude_setup_v2/jetzt/stage_a_c19_20260904.log'))
PREFILTER = 1e-3
HP_TOL = 1e-9
GEOM = 'all tetrahedra positively oriented'
CROSS_ORDER = {'17a': 17.0, '17n': 17.1, '18a': 18.0, '18n': 18.1, '19a': 19.0, '19n': 19.1}
CENSUS_SIZE = {'17a': 1769978, '17n': 6283385, '18a': 8400285, '18n': 39866095,
               '19a': 40619384, '19n': 253510828}


def _vol_worker(cand):
    import snappy
    try:
        M = snappy.Manifold('DT:' + repr([tuple(K.decode_dt(cand['dt']))]))
        if str(M.solution_type()) != GEOM:      # Finding 5: no degenerate volumes
            return (cand['det'], cand['crossing'], cand['dt'], None)
        return (cand['det'], cand['crossing'], cand['dt'], float(M.volume()))
    except Exception:
        return (cand['det'], cand['crossing'], cand['dt'], None)


def _confirm_isometric(target_isosig, burton_dt):
    import snappy
    try:
        ours = snappy.Manifold(target_isosig)
        burton = snappy.Manifold('DT:' + repr([tuple(K.decode_dt(burton_dt))]))
        if ours.is_isometric_to(burton):
            return True
        for _ in range(4):
            burton.randomize()
            if ours.is_isometric_to(burton):
                return True
        return False
    except Exception:
        return None


def stage_a_vollstaendig():
    """Gate: stage-A log must carry all six categories exactly + completion."""
    txt = open(ALOG, errors='replace').read()
    done = dict(re.findall(r"^DONE (\d\d[an]): (\d+) knots", txt, re.M))
    fehler = []
    for cat, soll in CENSUS_SIZE.items():
        ist = int(done.get(cat, -1))
        if ist != soll:
            fehler.append(f"{cat}: {ist} != expected {soll}")
    if 'STAGE A DONE' not in txt:
        fehler.append("Schlusszeile 'STAGE A DONE' fehlt")
    if 'COMPLETENESS FAIL' in txt:
        fehler.append("Log traegt COMPLETENESS FAIL")
    m = re.search(r"STAGE A DONE: \d+ knots in \d+s, \d+ det-matches, detfails=(\d+)", txt)
    if m and int(m.group(1)) != 0:
        fehler.append(f"detfails={m.group(1)} != 0")
    return fehler


def cands_stream(pfad):
    with open(pfad) as f:
        for line in f:
            line = line.strip()
            if line:
                yield json.loads(line)


def main():
    fehler = stage_a_vollstaendig()
    if fehler:
        raise SystemExit("🔴 Stage A NOT complete — Stage B refuses:\n  " + "\n  ".join(fehler))
    print(f"Stage A completeness OK: all {len(CENSUS_SIZE)} categories exact, detfails=0", flush=True)

    cand_path = WORK + '/candidates_stageA.jsonl'
    n_cands = sum(1 for _ in open(cand_path))
    det_to_targets = json.load(open(WORK + '/det_to_targets.json'))
    geg = json.load(open(WORK + '/gegenprobe.json'))
    geg_ids = set(geg.get('gegenprobe_ids', []))
    print(f"candidates={n_cands}  target-dets={len(det_to_targets)}  counter_checks={len(geg_ids)}", flush=True)

    import snappy
    cal = float(snappy.Manifold('4_1').volume())
    assert abs(cal - 2.029883) < 1e-3, f"CALIB FAIL 4_1 vol={cal}"
    print(f"CALIB OK 4_1 vol={cal:.5f}", flush=True)

    t0 = time.time(); n = 0
    fails_by_det, succ_by_det = {}, {}
    prefilter_hits = []
    with mp.Pool(int(os.environ.get('CVB_NPROC', os.cpu_count() or 24))) as pool:
        for det, crossing, dt, v in pool.imap_unordered(_vol_worker, cands_stream(cand_path), chunksize=200):
            n += 1
            if v is None:
                fails_by_det[det] = fails_by_det.get(det, 0) + 1
            else:
                succ_by_det[det] = succ_by_det.get(det, 0) + 1
                for tgt in det_to_targets.get(str(det), []):
                    dv = abs(v - tgt['vol'])
                    if dv < PREFILTER:
                        prefilter_hits.append((tgt['id'], tgt.get('isosig'), tgt['vol'],
                                               crossing, dt, v, dv))
            if n % 200000 == 0:
                print(f"  vol {n}/{n_cands} prefilter_hits={len(prefilter_hits)} "
                      f"fails={sum(fails_by_det.values())} {time.time()-t0:.0f}s", flush=True)
    assert n == n_cands, f"gestreamt {n} != gezaehlt {n_cands}"
    print(f"PASS 1 done: {n} vols, {len(prefilter_hits)} prefilter hits, "
          f"{sum(fails_by_det.values())} non-geometric/failed", flush=True)

    confirmed, vol_only, ambiguous = {}, {}, []
    for tid, isosig, tvol, crossing, bdt, bvol, dv in sorted(prefilter_hits,
                                                             key=lambda h: CROSS_ORDER[h[3]]):
        if tid in confirmed:
            continue
        if isosig:
            r = _confirm_isometric(isosig, bdt)
            if r is True:
                confirmed[tid] = crossing
            else:
                ambiguous.append({'id': tid, 'crossing': crossing, 'dvol': dv,
                                  'note': 'vol-near but NOT isometric' if r is False else 'unverifiable'})
        elif dv < HP_TOL and tid not in vol_only:
            vol_only[tid] = crossing

    beyond17 = json.load(open(WORK + '/beyond17_targets.json'))
    out = {'matched_c17_ALARM': [], 'matched_c18_downgrade': [], 'matched_c19_downgrade': [],
           'beyond_rigorous': [], 'beyond_unrigorous_skips': [], 'needs_isosig_review': []}
    for b in beyond17:
        tid, det = b['id'], b['det']
        if tid in confirmed:
            c = confirmed[tid]
            bucket = ('matched_c17_ALARM' if c.startswith('17')
                      else 'matched_c18_downgrade' if c.startswith('18')
                      else 'matched_c19_downgrade')
            out[bucket].append({'id': tid, 'det': det, 'crossing': c})
        elif tid in vol_only:
            out['needs_isosig_review'].append({'id': tid, 'det': det, 'crossing': vol_only[tid]})
        elif fails_by_det.get(det, 0) > 0:
            out['beyond_unrigorous_skips'].append({'id': tid, 'det': det,
                                                    'vol_fails_in_det': fails_by_det[det]})
        else:
            out['beyond_rigorous'].append({'id': tid, 'det': det})

    # COUNTER-CHECK LATCH (Werner 04./05.09.: the span>=21 targets MUST NOT show up)
    verletzt = sorted(g for g in geg_ids
                      if g in confirmed or g in vol_only)
    json.dump({'summary': {k: len(v) for k, v in out.items()},
               'detail': out, 'ambiguous_prefilter_hits': ambiguous,
               'gegenprobe_verletzt': verletzt,
               'gegenprobe_n': len(geg_ids),
               'total_candidates': n_cands,
               'total_nongeom_or_fail': sum(fails_by_det.values()),
               'fails_by_det': fails_by_det},
              open(WORK + '/stage_b_result.json', 'w'), indent=1)
    print("\n=== STAGE B RESULT (identity-safe, v4) ===", flush=True)
    for k, v in out.items():
        print(f"  {k}: {len(v)}", flush=True)
    print(f"  ambiguous (vol-near, not isometric/unverifiable): {len(ambiguous)}", flush=True)
    print(f"  non-geometric/failed vols: {sum(fails_by_det.values())} / {n_cands}", flush=True)
    if out['matched_c17_ALARM']:
        print("  ⚠️ c17 MATCHES (Maguire-hole/bug!):", [m['id'] for m in out['matched_c17_ALARM']], flush=True)
    if verletzt:
        print(f"🔴 COUNTER-CHECK VIOLATED — {len(verletzt)} span>=21 targets match in the c<=19 census: "
              f"{verletzt}\n   This refutes a PROOF (span or genus route), not the assignment. STOP.", flush=True)
        sys.exit(3)
    print(f"✅ counter-check holds: 0 of {len(geg_ids)} span>=21 targets found in the census.", flush=True)


if __name__ == '__main__':
    main()
