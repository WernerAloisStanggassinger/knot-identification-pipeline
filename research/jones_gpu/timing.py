#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""timing.py — per-knot CPU times on the first N rows of a census file, one core: spherogram DT->PD, prepare_fast (2 points),
and, for comparison, Regina Link.fromDT + Link.jones() (the exact Jones polynomial). Usage: python -I timing.py census.csv.bz2 N"""
import sys, os, bz2, csv, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spherogram, regina
from jones_gpu import prepare_fast
from jones_tn import P, inv
dt = lambda s: [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]
f, n = sys.argv[1], int(sys.argv[2])
rows = [r['dt_code'] for r in itertools.islice(csv.DictReader(bz2.open(f, 'rt')), n)]
t = time.time(); pds = [[tuple(X) for X in spherogram.Link('DT:' + str([tuple(dt(c))])).PD_code()] for c in rows]; a = (time.time() - t) / n
x = 987654321 % P
t = time.time(); [prepare_fast(pd, (x, inv(x))) for pd in pds]; b = (time.time() - t) / n
t = time.time(); Ls = [regina.Link.fromDT(' '.join(map(str, dt(c)))) for c in rows]; c1 = (time.time() - t) / n
t = time.time(); [L.jones() for L in Ls]; c2 = (time.time() - t) / n
print(f'{os.path.basename(f)} first {n} rows, ms/knot, one core: spherogram DT->PD {1000*a:.2f} | prepare_fast (2 points) {1000*b:.2f} | '
      f'Regina fromDT {1000*c1:.3f} + jones {1000*c2:.3f}')
