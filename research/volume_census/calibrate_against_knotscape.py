# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""Calibration of volume_census.py: the multiset of volumes of Burton's 15n-hyp (output of volume_census.py) is compared with the
multiset of volumes of KnotScape's 15n hyperbolic knots (a table name<TAB>volume<TAB>… with a header line, e.g. volumes computed from another knot table, 12 digits),
up to rounding. Names differ between the tables, so only the multisets are compared."""
import sys, collections
ours = [float(l.split('\t')[1]) for l in open(sys.argv[1]) if not l.startswith('#') and l.split('\t')[1] not in ('None', '')]
flags = collections.Counter(l.split('\t')[2] for l in open(sys.argv[1]) if not l.startswith('#'))
ks = []
for i, l in enumerate(open(sys.argv[2])):
    if i == 0: continue
    v = l.split('\t')[1]
    try:
        x = float(v)
    except ValueError:
        continue
    if x > 0.5: ks.append(x)
print('ours', len(ours), 'flags', dict(flags), '| knotscape hyperbolic', len(ks))
for d in (6, 8, 9):
    a = collections.Counter(round(x, d) for x in ours); b = collections.Counter(round(x, d) for x in ks)
    print(f'digits {d}: multisets equal = {a == b}; only ours {sum((a - b).values())}, only knotscape {sum((b - a).values())}')
