# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# PD codes (SnapPy) of K and of its mirror for khovanov_witness.py, which runs in an environment with khoca but without SnapPy.
# From the authors' make_pd_kh.py (same construction). Usage: python -I make_kh_pd.py names.txt|counter out.json
import sys, json, os, snappy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from census_codes import codes, DT
out = []
if sys.argv[1] == 'counter':
    for nm in ['3_1', '4_1', '9_42', '10_125']: L = snappy.Link(nm); out.append([nm, 'counter', L.PD_code(), L.mirror().PD_code()])
else:
    names = [l.split()[0] for l in open(sys.argv[1]) if l.strip()]; cd = codes(names)
    for nm in names: L = snappy.Link('DT:' + str([tuple(DT(cd[nm]))])); out.append([nm, 'target', L.PD_code(), L.mirror().PD_code()])
json.dump(out, open(sys.argv[2], 'w')); print(len(out), 'knots')
