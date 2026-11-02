#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""jones_census.py — Jones polynomial at a point and at its inverse, exactly mod p, for every knot of Burton census files (GPU).

CPU pool: DT code -> spherogram PD -> jones_gpu.prepare_fast (one crossing order and one exponent form, evaluated at x and 1/x).
GPU: knots bucketed by (crossings, frontier width F); a bucket is contracted when it holds MAXCELLS / 2^F knots (memory latch).
Output per knot (written on arrival, bz2): name, crossings, F, V(t), V(1/t) with t = x^-8 mod p.
Reading: given the calibrated kernel and a PD of the census knot, V(t) != V(1/t) mod p ==> V_K != V_mirror(K) ==> K is chiral
         (the mirror's Jones polynomial is V(1/t)); the stats key 'chiral_proved' means exactly this.
         V(t) == V(1/t) mod p  ==>  only "Jones symmetric at the point" (flagged).
Counter-cases first: 3_1 (chiral, must differ), 4_1 (amphichiral, must agree); GPU (fast path, as used here) == jones_tn.jones_at
is test_gpu_vs_ref.py.
Usage: python -I jones_census.py <workers> <out_prefix> <census.csv.bz2> ...   (env JONES_MAXCELLS, default 2^23)
"""
import sys, os, bz2, csv, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multiprocessing import Pool
from jones_tn import P, inv, jones_at

X0 = 987654321 % P
MAXCELLS = int(os.environ.get('JONES_MAXCELLS', str(1 << 23)))
def dt(s): return [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]


def init():
    global spherogram, prepare_fast
    import spherogram
    from jones_gpu import prepare_fast


def conv(row):
    name, code = row
    try:
        pd = [tuple(X) for X in spherogram.Link('DT:' + str([tuple(dt(code))])).PD_code()]
        steps, F, w, ((Ta, Aa), (Tb, Ab)) = prepare_fast(pd, (X0, inv(X0)))
        return name, len(pd), (Ta, steps, F, Aa, w), (Tb, steps, F, Ab, w)
    except Exception as e:
        return name, -1, 'ERR:' + repr(e)[:100], None


def rows(files):
    for f in files:
        with bz2.open(f, 'rt') as fh:
            r = csv.reader(fh); next(r)
            for x in r: yield (x[0], x[2])


if __name__ == '__main__':
    import spherogram as sg
    from jones_gpu import prepare as prep0, contract_batch, jones_from_bracket
    workers, pre, files = int(sys.argv[1]), sys.argv[2], sys.argv[3:]
    print(json.dumps({'workers': workers, 'files': [os.path.basename(f) for f in files], 'maxcells': MAXCELLS}), flush=True)
    # counter-cases
    for nm, must_differ in (('3_1', True), ('4_1', False)):
        pd = [tuple(X) for X in sg.Link(nm).PD_code()]
        va = jones_at(pd, X0)[0]; vb = jones_at(pd, inv(X0))[0]
        assert (va != vb) == must_differ, (nm, va, vb)
    print('counter-cases OK (3_1 differs, 4_1 agrees)', flush=True)
    out = bz2.open(pre + '_jones.tsv.bz2', 'wt'); err = open(pre + '_errors.tsv', 'w')
    stats = dict(knots=0, symmetric=0, chiral_proved=0, errors=0); t0 = time.time()
    buckets = {}; checked = 0

    def flush(key):
        global checked
        items = buckets.pop(key)
        F = key[1]
        ba = contract_batch([(it[2][0], it[2][1]) for it in items], F)
        bb = contract_batch([(it[3][0], it[3][1]) for it in items], F)
        for it, x, y in zip(items, ba, bb):
            name, c, a, b = it
            va = jones_from_bracket(x, a[3], a[4]); vb = jones_from_bracket(y, b[3], b[4])
            stats['knots'] += 1
            if va == vb: stats['symmetric'] += 1
            else: stats['chiral_proved'] += 1
            out.write(f'{name}\t{c}\t{F}\t{va}\t{vb}\n')
    with Pool(workers, initializer=init) as pool:
        for name, c, a, b in pool.imap_unordered(conv, rows(files), chunksize=200):
            if c < 0:
                stats['errors'] += 1; err.write(f'{name}\t{a}\n'); continue
            F = max(a[2], b[2]); key = (c, F)
            buckets.setdefault(key, []).append((name, c, a, b))
            if len(buckets[key]) >= max(1, MAXCELLS >> F):
                flush(key); out.flush()
                print(json.dumps(dict(stats, sec=round(time.time() - t0))), flush=True)
        for key in list(buckets): flush(key)
    out.close(); err.close(); stats['sec'] = round(time.time() - t0); stats['p'] = P; stats['x'] = X0
    json.dump(stats, open(pre + '_stats.json', 'w')); print('DONE', json.dumps(stats), flush=True)
