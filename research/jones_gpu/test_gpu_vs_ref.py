# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Counter-case for jones_gpu.py: GPU batch result == jones_tn.jones_at (calibrated against Sage) knot by knot, at x and at 1/x,
# for both the plain path (prepare) and the fast path (prepare_fast, used by jones_census.py). Usage: python -I test_gpu_vs_ref.py census.csv.bz2 N
import sys, bz2, csv, itertools, random, time
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import spherogram
from jones_tn import jones_at, P, inv
from jones_gpu import prepare, prepare_fast, contract_batch, jones_from_bracket
dt = lambda s: [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]
f, n = sys.argv[1], int(sys.argv[2])
random.seed(3)
rows = list(itertools.islice(csv.DictReader(bz2.open(f, 'rt')), 20000))
pds = [[tuple(X) for X in spherogram.Link('DT:' + str([tuple(dt(r['dt_code']))])).PD_code()] for r in random.sample(rows, n)]
x = 123456789 % P
bad = 0; byc = {}
for pd in pds: byc.setdefault(len(pd), []).append(pd)
t0 = time.time(); nk = 0; maxF = 0
for c, group in byc.items():
    fast = [prepare_fast(pd, (x, inv(x))) for pd in group]           # the path jones_census.py uses
    for xi, xx in enumerate((x, inv(x))):
        prep = [prepare(pd, xx) for pd in group]                      # the plain path
        for path, items in (('plain', [(pp[0], pp[1], pp[2], pp[3], pp[4]) for pp in prep]),
                            ('fast', [(fp[3][xi][0], fp[0], fp[1], fp[3][xi][1], fp[2]) for fp in fast])):
            F = max(it[2] for it in items); maxF = max(maxF, F)
            br = contract_batch([(it[0], it[1]) for it in items], F)
            for pd, it, b in zip(group, items, br):
                vg = jones_from_bracket(b, it[3], it[4]); vr, _ = jones_at(pd, xx)
                nk += 1; bad += (vg != vr)
print(f'knot x point x path checks {nk}, mismatches {bad}, crossings {sorted(byc)}, max F {maxF}, {time.time()-t0:.1f} s')
assert bad == 0
print('GPU == reference')
