#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""calibrate_jones_tn.py — fix the two conventions of jones_tn.py (which smoothing is A, sign of the writhe) and the orientation of
t against Sage's Jones polynomial (spherogram Link.jones_polynomial(new_convention=False), a Laurent polynomial in q).
Every combination is tested on every knot at 3 random points; exactly the mirror-dual pair of combinations must agree on ALL knots
(the chiral knots separate q = t from q = 1/t).
Usage (Python with SageMath + spherogram): python -I calibrate_jones_tn.py [n_census_knots]
"""
import sys, random, itertools, bz2, csv, json
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import spherogram
from jones_tn import jones_at, P, inv

random.seed(7)
def sage_eval(J, q):
    tot = 0
    for e, c in J.dict().items():
        e = int(e); c = int(c)
        tot = (tot + c * (pow(q, e, P) if e >= 0 else pow(inv(q), -e, P))) % P
    return tot

names = ['3_1', '4_1', '5_1', '5_2', '6_1', '6_2', '6_3', '7_4', '8_19', '8_20', '9_42', '10_124', 'K11n34', 'K11n42']
links = [(n, spherogram.Link(n)) for n in names]
ncen = int(sys.argv[1]) if len(sys.argv) > 1 else 40
dt = lambda s: [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]
RAW = __import__('os').environ.get('BURTON_CENSUS_DIR', 'burton_census/raw')   # folder with Burton's 14n-hyp.csv.bz2
rows = list(itertools.islice(csv.DictReader(bz2.open(RAW + '/14n-hyp.csv.bz2', 'rt')), 5000))
for r in random.sample(rows, ncen):
    links.append((r['name'], spherogram.Link('DT:' + str([tuple(dt(r['dt_code']))]))))
combos = list(itertools.product([True, False], [1, -1], ['t', '1/t']))
ok = {c: 0 for c in combos}; total = 0
for name, L in links:
    pd = [tuple(X) for X in L.PD_code()]
    J = L.jones_polynomial(new_convention=False)
    total += 1
    for x in [random.randrange(2, P - 1) for _ in range(3)]:
        res = {}
        for A02, ws in itertools.product([True, False], [1, -1]):
            res[(A02, ws)] = jones_at(pd, x, A02, ws)
        for (A02, ws, tq) in combos:
            v, t = res[(A02, ws)]
            q = t if tq == 't' else inv(t)
            if v != sage_eval(J, q): ok[(A02, ws, tq)] = ok.get((A02, ws, tq), 0) - 10**6
    for c in combos: ok[c] += 1
good = [c for c in combos if ok[c] == total]
print(json.dumps({'knots': total, 'per_combo_agreeing_knots': {str(c): max(v, -1) for c, v in ok.items()}, 'consistent': [str(c) for c in good]}))
# Mirroring swaps the two smoothings, flips the writhe sign and sends t -> 1/t, so the conventions come in dual pairs:
# exactly the pair {(True, 1, 't'), (False, -1, '1/t')} must agree on every knot (both describe the same evaluation).
assert set(good) == {(True, 1, 't'), (False, -1, '1/t')}, good
print('CALIBRATION OK: A-smoothing = slots (0,1)(2,3), writhe sign +1, jones_at returns V(t) at t = x^-8 (Sage old convention q = t)')
