#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""scan_triage.py — splits failed Stage-2 runs into TWO kinds,
instead of calling them all "tracker failed".

The tracker aborts when no cell in the (ti,tj) grid reaches the threshold
--min-cell-count. That has two completely different causes:

  UNKNOTEN   (unknot) There are projection directions WITHOUT a crossing. A closed
             curve with a crossing-free projection is embedded in the plane
             ⟹ unknot. The tracker rightly refuses, there is nothing to
             identify. NOT a parameter problem.

  VERSCHMIERT (smeared) Every projection shows crossings (min>0), but none sits
             still in arc-parameter space. The tracker's model assumption
             ("real crossing = compact cluster") is violated. THAT is the
             real gap.

⚠️ THE TEST IS ONE-SIDED — and that is intentional.
  "There is a crossing-free projection  ⟹  unknot"   is a theorem.
  The CONVERSE does NOT hold: a curve can show crossings from EVERY direction
  and still be the unknot (a kink that no perspective
  turns away). VERSCHMIERT therefore does NOT mean "is a knot", but
  "the cheap unknot proof did not succeed".
  Shown on 23.08.: of 18 curves classified as VERSCHMIERT,
  36__tg_s1 (min/frame=1) and 46__la_crisscross_n5 (min/frame=3) resolved as K0_1
  in the full pipeline.
  For the purpose, exactly this direction is the right one: as a gate before
  lowering the tracker threshold, the test may at most be TOO CAUTIOUS.
  A false "unknot" would throw away a real knot; a false
  "knot" costs only one pipeline run, at the end of which the truth stands.

Read-only. Touches neither the pipeline nor any record.
"""
import sys, glob, os, json, argparse
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor

def eine(pq):
    try:
        df = pd.read_parquet(pq, columns=['t_step'])
    except Exception as e:
        return dict(pq=pq, fehler=str(e)[:80])
    if len(df) == 0:
        return dict(pq=pq, n=0, frames=0, minc=0, medc=0, nullanteil=1.0, sorte='UNKNOTEN')
    ts = df['t_step'].to_numpy()
    n_frames = int(ts.max()) + 1
    voll = np.bincount(ts, minlength=n_frames)
    null = float((voll == 0).mean())
    minc = int(voll.min())
    sorte = 'UNKNOTEN' if minc == 0 else 'VERSCHMIERT'
    return dict(pq=pq, n=int(len(df)), frames=n_frames, minc=minc,
                medc=int(np.median(voll)), nullanteil=null, sorte=sorte)

def selbsttest():
    """Counter-case built in: the test must catch both kinds AND a mix-up."""
    import tempfile
    ok = 0; total = 0
    def pruef(name, ist, soll):
        nonlocal ok, total
        total += 1
        gut = ist == soll
        ok += gut
        print(f"  {'✅' if gut else '🔴'} {name}: {ist!r} (expected {soll!r})")
    with tempfile.TemporaryDirectory() as d:
        # (1) curve with a crossing-free projection ⟹ UNKNOTEN
        p1 = os.path.join(d, 'a.parquet')
        pd.DataFrame({'t_step': [0,0,0, 2,2, 3]}).to_parquet(p1)   # frame 1 empty
        pruef('leerer Frame ⟹ UNKNOTEN', eine(p1)['sorte'], 'UNKNOTEN')
        # (2) every projection carries crossings ⟹ VERSCHMIERT
        p2 = os.path.join(d, 'b.parquet')
        pd.DataFrame({'t_step': [0,0, 1,1, 2,2, 3,3]}).to_parquet(p2)
        pruef('kein leerer Frame ⟹ VERSCHMIERT', eine(p2)['sorte'], 'VERSCHMIERT')
        # (3) COUNTER-CASE: a single gap must be enough to flip (2) —
        #     otherwise the test measures the amount instead of the minimum.
        p3 = os.path.join(d, 'c.parquet')
        pd.DataFrame({'t_step': [0]*500 + [1]*500 + [3]*500}).to_parquet(p3)
        pruef('EINE Luecke unter 1500 Punkten kippt es', eine(p3)['sorte'], 'UNKNOTEN')
        # (4) empty file ⟹ UNKNOTEN, no crash
        p4 = os.path.join(d, 'e.parquet')
        pd.DataFrame({'t_step': pd.Series([], dtype='int64')}).to_parquet(p4)
        pruef('leere Datei ⟹ UNKNOTEN', eine(p4)['sorte'], 'UNKNOTEN')
        # (5) canary: a result that must NOT be true
        pruef('Kanarie (muss rot sein, wenn der Test blind ist)',
              eine(p2)['sorte'] == 'UNKNOTEN', False)
    print(f"\n  {ok}/{total} probes")
    return 0 if ok == total else 1

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('wurzel', nargs='?')
    ap.add_argument('--selbsttest', action='store_true')
    ap.add_argument('--gate', metavar='CROSSINGS.PARQUET',
                    help='gate mode for the pipeline: exit 0 = UNKNOTEN (unknot) '
                         '(there is a crossing-free projection), '
                         'exit 1 = knot, exit 2 = unreadable/empty.')
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--json', help='write the result here')
    a = ap.parse_args()
    if a.selbsttest:
        sys.exit(selbsttest())
    if a.gate:
        r = eine(a.gate)
        if 'fehler' in r:
            print(f"[GATE] unreadable: {r['fehler']}", file=sys.stderr)
            sys.exit(2)
        if r['sorte'] == 'UNKNOTEN':
            print(f"[GATE] UNKNOTEN — {r['nullanteil']*100:.1f}% of the {r['frames']:,} "
                  f"projections are crossing-free (min/frame = {r['minc']}). "
                  f"A closed curve with a crossing-free projection is embedded in the "
                  f"plane.")
            sys.exit(0)
        print(f"[GATE] NO unknot proof — each of the {r['frames']:,} projections "
              f"carries crossings (min/frame = {r['minc']}, median {r['medc']}). "
              f"This is NO proof that it is a knot (the test is one-sided); "
              f"it only means: lowering the threshold is allowed, the pipeline decides.")
        sys.exit(1)
    if not a.wurzel:
        ap.error('wurzel, --gate or --selbsttest is required')

    pqs = []
    for d in sorted(glob.glob(os.path.join(a.wurzel, '*/'))):
        if glob.glob(os.path.join(d, '*_pflicht_invarianten_complete.json')):
            continue
        g = glob.glob(os.path.join(d, '*_crossings.parquet'))
        if g:
            pqs.append(g[0])
    print(f"{len(pqs)} candidates without a record, with scanner output")
    with ProcessPoolExecutor(a.workers) as ex:
        res = list(ex.map(eine, pqs))
    res = [r for r in res if 'fehler' not in r]
    for r in res:
        r['stamm'] = os.path.basename(os.path.dirname(r['pq']))
    for sorte in ('UNKNOTEN', 'VERSCHMIERT'):
        teil = sorted([r for r in res if r['sorte'] == sorte],
                      key=lambda r: -r['n'])
        print(f"\n=== {sorte}: {len(teil)} ===")
        print(f"  {'stem':<44}{'crossings':>11}{'min/frame':>10}{'med':>6}{'zero frames':>13}")
        for r in teil:
            print(f"  {r['stamm']:<44}{r['n']:>11,}{r['minc']:>10}{r['medc']:>6}"
                  f"{r['nullanteil']*100:>12.1f}%")
    if a.json:
        json.dump(res, open(a.json, 'w'), indent=1)
        print(f"\n→ {a.json}")
