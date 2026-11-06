#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Stage B (v3, identity-safe + audit-hardened): volume-confirm the det-matching
Burton candidates from Stage A against the 347 beyond-17 targets, then CONFIRM
every volume near-match with `is_isometric_to` — NOT vol±tolerance (Disc #516).
Runs in gnottheorie_env (snappy).

Hardening:
  - Disc #516: loose vol prefilter (1e-3) -> is_isometric_to confirmation.
  - Audit Finding 5: reject non-geometric SnapPy solutions (solution_type) so a
    degenerate volume can't spuriously match — routes to failure/rigor accounting.
  - Accepts both candidate formats (bare list = old Stage A; {'candidates',
    'provenance'} = hardened Stage A) and asserts provenance completeness if present.

Classification per target: matched_c17_ALARM, matched_c18_downgrade,
matched_c19_downgrade, beyond_rigorous, beyond_unrigorous_skips, needs_isosig_review.
"""
import site_config
import sys, os, json, time
import multiprocessing as mp
sys.path.insert(0, os.environ.get('CVB_KNOTLIB', (site_config.MACHINE_YARD + '/corpus_vs_maguire/work')))
import knotlib as K

WORK = os.environ.get('CVB_WORK', (site_config.MACHINE_YARD + '/corpus_vs_burton/work'))
PREFILTER = 1e-3
HP_TOL = 1e-9
GEOM = 'all tetrahedra positively oriented'
CROSS_ORDER = {'17a': 17.0, '17n': 17.1, '18a': 18.0, '18n': 18.1, '19a': 19.0, '19n': 19.1}

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

def main():
    raw = json.load(open(WORK + '/candidates_stageA.json'))
    if isinstance(raw, dict):
        cands = raw['candidates']
        prov = raw.get('provenance', {})
        bad = [c for c, p in prov.items() if not p.get('complete', True)]
        if bad:
            raise SystemExit(f"🔴 Stage A completeness FAIL for {bad} — refusing to run Stage B")
        print(f"Stage A provenance OK: { {c: p['lines'] for c, p in prov.items()} }", flush=True)
    else:
        cands = raw
        print("⚠️ candidates in legacy bare-list format (no provenance manifest)", flush=True)
    det_to_targets = json.load(open(WORK + '/det_to_targets.json'))
    print(f"candidates={len(cands)}  target-dets={len(det_to_targets)}", flush=True)

    import snappy
    cal = float(snappy.Manifold('4_1').volume())
    assert abs(cal - 2.029883) < 1e-3, f"CALIB FAIL 4_1 vol={cal}"
    print(f"CALIB OK 4_1 vol={cal:.5f}", flush=True)

    t0 = time.time(); n = 0
    fails_by_det, succ_by_det = {}, {}
    prefilter_hits = []
    with mp.Pool(int(os.environ.get('CVB_NPROC', os.cpu_count() or 24))) as pool:
        for det, crossing, dt, v in pool.imap_unordered(_vol_worker, cands, chunksize=200):
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
                print(f"  vol {n}/{len(cands)} prefilter_hits={len(prefilter_hits)} "
                      f"fails={sum(fails_by_det.values())} {time.time()-t0:.0f}s", flush=True)
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

    json.dump({'summary': {k: len(v) for k, v in out.items()},
               'detail': out, 'ambiguous_prefilter_hits': ambiguous,
               'total_candidates': len(cands),
               'total_nongeom_or_fail': sum(fails_by_det.values()),
               'fails_by_det': fails_by_det},
              open(WORK + '/stage_b_result.json', 'w'), indent=1)
    print("\n=== STAGE B RESULT (identity-safe) ===", flush=True)
    for k, v in out.items():
        print(f"  {k}: {len(v)}", flush=True)
    print(f"  ambiguous (vol-near, not isometric/unverifiable): {len(ambiguous)}", flush=True)
    print(f"  non-geometric/failed vols: {sum(fails_by_det.values())} / {len(cands)}", flush=True)
    if out['matched_c17_ALARM']:
        print("  ⚠️ c17 MATCHES (Maguire-hole/bug!):", [m['id'] for m in out['matched_c17_ALARM']], flush=True)

if __name__ == '__main__':
    main()
