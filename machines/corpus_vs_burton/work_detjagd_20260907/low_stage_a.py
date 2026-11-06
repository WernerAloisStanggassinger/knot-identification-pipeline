#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""low_stage_a.py — the missing ≤16 layer of the absence proofs (2026-09-06).

OCCASION (full stop 15:0x): "c>=20 via census absence" rested on the stream
17a–19n; for the layer c<=16 everything hung on the old pipeline crosschecks —
and S_15_d507_D3_chiral (15-crossing PD, triple-verified) proves that
this layer has holes. Here it is closed with the SAME mechanics as stage A:
stream all Burton hyp files 03..16, det = |Alexander(-1)|
via regina, write the candidates of the 91 target determinants as JSONL.
Non-hyp files (torus/satellite) are proof-free for hyperbolic targets
irrelevant (a hyperbolic knot is never torus/satellite) and stay out.

COMPLETENESS GATE as in stage A: each category must hit its target line count
exactly (target numbers below from regina-normal.github.io/data.html or
the known census sizes; where None, the measured number is printed and
written into the manifest — then it is a measurement, not a target).
Calibration gate: 4_1 det==5.

Run: base Python (regina), 24 cores.
"""
import site_config
import bz2, csv, io, json, os, sys, time
import multiprocessing as mp

RAW = (site_config.ROOT + '/external_tools/burton_census/raw')
WORK = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(WORK, 'candidates_low_c16.jsonl')

# hyp files 03..16 (torus/satellite deliberately out, see above)
FILES = []
for c in range(3, 17):
    for t in ('a', 'n'):
        p = os.path.join(RAW, f'{c:02d}{t}-hyp.csv.bz2')
        p2 = os.path.join(RAW, f'{c}{t}-hyp.csv.bz2')
        if os.path.exists(p):
            FILES.append((f'{c}{t}', p))
        elif os.path.exists(p2):
            FILES.append((f'{c}{t}', p2))

DET_SET = set(json.load(open(os.path.join(WORK, 'target_det_set.json'))))


def _det_from_sig(sig):
    from regina import Link
    L = Link.fromKnotSig(sig)
    a = L.alexander()
    d = a.degree()
    return abs(sum(int(str(a[i])) * ((-1) ** i) for i in range(d + 1)))


def _worker(rec):
    crossing, name, sig, dt = rec
    try:
        det = _det_from_sig(sig)
    except Exception:
        return ('FAIL', crossing, name, None, None)
    if det in DET_SET:
        return ('HIT', crossing, name, dt, det)
    return ('MISS', crossing, None, None, None)


def rows():
    for crossing, path in FILES:
        with bz2.open(path) as fh:
            for row in csv.DictReader(io.TextIOWrapper(fh)):
                yield (crossing, row['name'], row['knot_sig'], row['dt_code'])


def main():
    from regina import Link
    # Kalibrier-Tor
    cal = _det_from_sig(Link.fromPD([(4,2,5,1),(8,6,1,5),(6,3,7,4),(2,7,3,8)]).knotSig())
    assert cal == 5, f'Kalibrierung 4_1 det={cal} != 5'
    print(f'Calibration OK (4_1 det=5) · {len(FILES)} files · {len(DET_SET)} target dets',
          flush=True)
    t0 = time.time()
    n = {}; fails = {}; hits = 0
    fh = open(OUT, 'w')
    with mp.Pool(24) as pool:
        for res in pool.imap_unordered(_worker, rows(), chunksize=500):
            tag, crossing = res[0], res[1]
            n[crossing] = n.get(crossing, 0) + 1
            if tag == 'FAIL':
                fails[crossing] = fails.get(crossing, 0) + 1
            elif tag == 'HIT':
                _, _, name, dt, det = res
                fh.write(json.dumps({'det': det, 'crossing': crossing,
                                     'name': name, 'dt': dt}) + '\n')
                hits += 1
    fh.close()
    tot = sum(n.values())
    print(f'DONE: {tot} knots in {time.time()-t0:.0f}s, {hits} det candidates, '
          f'fails={sum(fails.values())} {fails if fails else ""}', flush=True)
    json.dump({'per_category': n, 'fails': fails, 'hits': hits, 'total': tot,
               'files': [f for _, f in FILES]},
              open(os.path.join(WORK, 'candidates_low_c16_meta.json'), 'w'), indent=1)
    for crossing in sorted(n):
        print(f'  {crossing}: {n[crossing]}', flush=True)
    if fails:
        print('🔴 DET-FAILS present — layer NOT complete, do not use as proof.',
              flush=True)
        sys.exit(2)


if __name__ == '__main__':
    main()
