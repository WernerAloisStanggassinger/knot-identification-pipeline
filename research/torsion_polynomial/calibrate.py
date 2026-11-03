#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Calibration of the DFJ comparator inside the range DFJ checked (<= 15 crossings), and counter-cases that must not 'agree'.
# Usage: python -I calibrate.py out.jsonl [N_RANDOM=400]   (SageMath's Python; Burton's census files in BURTON_CENSUS_DIR)
# g, fibered from SnapPy knot_floer_homology (Szabo's program) on the same diagram -- independent of T_K.
import sys, json, bz2, random, time, csv
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import dfj_core as D
D.set_memlimit(6)
import snappy

RAW = __import__('os').environ.get('BURTON_CENSUS_DIR', 'burton_census/raw').rstrip('/') + '/'   # Burton's <c><a|n>-<hyp|torus|satellite>.csv.bz2
out = open(sys.argv[1], 'w')
n_rand = int(sys.argv[2]) if len(sys.argv) > 2 else 400


def M_from_dt(dt):
    c = chr(96 + len(dt))
    return snappy.Manifold('DT:' + c + 'a' + c + dt)


def emit(r):
    out.write(json.dumps(r) + '\n'); out.flush()
    print(r.get('name'), r.get('status'), r.get('deg'), r.get('four_g_minus_2'), r.get('monic'), r.get('fibered'),
          '%.4g%+.4gi' % (r.get('lc_re', 0), r.get('lc_im', 0)), r.get('notes'), flush=True)


def run(name, M, kind, g_override=None, fib_override=None):
    try:
        h = M.link().knot_floer_homology()
        g, fib = h['seifert_genus'], h['fibered']
        if g_override is not None: g = g_override(g)
        if fib_override is not None: fib = fib_override(fib)
        r = D.check(M, g, fib)
        r.update(name=name, kind=kind)
    except Exception as e:  # noqa: BLE001
        r = dict(name=name, kind=kind, status='ERROR', err=repr(e)[:300])
    emit(r)


# 1. named knots in DFJ range
for nm in ['4_1', '5_2', '6_1', '6_2', '6_3', '7_4', '8_20', '8_21', '9_42', '10_132', 'K11n42', 'K11n34', 'K12n242',
           '10_153', '10_157', 'K14a506', '8_17', 'K12a1202', 'K12n706']:
    try:
        M = snappy.Manifold(nm)
    except Exception as e:  # noqa: BLE001
        emit(dict(name=nm, kind='named', status='ERROR', err=repr(e)[:200])); continue
    run(nm, M, 'named')

# 2. counter-cases that MUST NOT come out as 'agree'
run('4_1 g+1', snappy.Manifold('4_1'), 'counter_wrong_genus', g_override=lambda g: g + 1)
run('5_2 fib flipped', snappy.Manifold('5_2'), 'counter_wrong_fibered', fib_override=lambda f: not f)
run('8_20 fib flipped', snappy.Manifold('8_20'), 'counter_wrong_fibered', fib_override=lambda f: not f)
run('K11n42 g-1', snappy.Manifold('K11n42'), 'counter_wrong_genus', g_override=lambda g: g - 1)
for nm in ['5_1', '7_1', 'K10n21']:  # torus knots: T(2,5), T(2,7), T(3,5)
    run(nm, snappy.Manifold(nm), 'counter_nonhyperbolic')
for f in ['16n-torus.csv.bz2', '16n-satellite.csv.bz2']:
    with bz2.open(RAW + f, 'rt') as fh:
        rd = csv.reader(fh); next(rd)
        for row in rd:
            name, sig, dt = row[:3]
            run(name, M_from_dt(dt), 'counter_nonhyperbolic')

# 3. random sample of Burton census <= 15 crossings (DFJ's range): 15n, 15a, 14n, 13n, 12n
random.seed(20261016)
for f, k in [('15n-hyp.csv.bz2', n_rand), ('15a-hyp.csv.bz2', n_rand // 4), ('14n-hyp.csv.bz2', n_rand // 4),
             ('13n-hyp.csv.bz2', n_rand // 8), ('12n-hyp.csv.bz2', n_rand // 8)]:
    with bz2.open(RAW + f, 'rt') as fh:
        next(fh)
        rows = [l.strip().split(',') for l in fh]
    for row in random.sample(rows, min(k, len(rows))):
        name, sig, dt = row[:3]
        run(name, M_from_dt(dt), 'sample_le15')
