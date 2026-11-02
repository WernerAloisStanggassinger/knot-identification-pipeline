# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""Calibration: (1) Python reference mod p == reduction of the Q value from Sage (400 knots, 3 to 11 crossings); (2) GPU == reference.
Counter-case: at a wrong point (T1 + 1) no knot with theta_raw != 0 may agree."""
import json, sys, time
from fractions import Fraction
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from theta_modp import ref, gpu_batch, frac_mod, P_DEFAULT as p
R = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ref_sage_22_7_21_13.json')))
t1, t2 = frac_mod(22, 7, p), frac_mod(21, 13, p)
def red(s):
    f = Fraction(s); return frac_mod(f.numerator, f.denominator, p)
bad_ref = 0
for r in R:
    o = ref(r['Cs'], r['phi'], t1, t2, p)
    expected = dict(det1=red(r['d1']), det2=red(r['d2']), det3=red(r['d3']), theta_roh=red(r['theta']))
    got = {k: o[k] for k in expected}
    if got != expected:
        bad_ref += 1
        if bad_ref <= 3: print('REFERENCE DEVIATES', r['name'], got, expected)
print(f'(1) reference vs Sage-Q mod p: {len(R)-bad_ref}/{len(R)} equal')
# Counter-case: another point must (almost always) give another theta
g = sum(ref(r['Cs'], r['phi'], t1 + 1, t2, p)['theta_roh'] == red(r['theta']) for r in R)
g_nz = sum(ref(r['Cs'], r['phi'], t1 + 1, t2, p)['theta_roh'] == red(r['theta']) for r in R if red(r['theta']) != 0)
n_nz = sum(red(r['theta']) != 0 for r in R)
print(f'    counter-case t1+1: {g}/{len(R)} equal; among the {n_nz} with theta != 0 at the point: {g_nz} equal (expected 0)')
# (2) GPU vs reference, bundled per n
from collections import defaultdict
grp = defaultdict(list)
for r in R: grp[len(r['Cs'])].append(r)
bad_gpu = 0; t0 = time.time()
for n, rs in sorted(grp.items()):
    out = gpu_batch([(r['Cs'], r['phi']) for r in rs], t1, t2, p)
    for r, o in zip(rs, out):
        expected = dict(det1=red(r['d1']), det2=red(r['d2']), det3=red(r['d3']), theta_roh=red(r['theta']))
        if {k: o[k] for k in expected} != expected:
            bad_gpu += 1
            if bad_gpu <= 3: print('GPU DEVIATES', r['name'], {k: o[k] for k in expected}, expected)
print(f'(2) GPU vs Sage-Q mod p: {len(R)-bad_gpu}/{len(R)} equal  ({time.time()-t0:.1f} s incl. warm-up)')
