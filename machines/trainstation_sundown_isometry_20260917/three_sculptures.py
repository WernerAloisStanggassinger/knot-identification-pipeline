# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Double-3 / Always-a-Pleasure / Proxima B: same knot 12a_1019? Isometry from the own .fseries (17.09.2026, for Ode speaking text)
import sys, json
sys.argv=['x']
exec(open('isometry.py').read().split('T = aus_kurve("trainstation")')[0])
R = {}
X = {s: next(x for x in aus_kurve(s) if 'E' in x) for s in ('double-3','always-a-pleasure','proxima-b')}
K = snappy.Manifold('K12a1019')
for s, x in X.items(): print(s, {k: v for k, v in x.items() if k not in ('E','L')})
R['double-3~always-a-pleasure'] = iso(X['double-3']['E'], X['always-a-pleasure']['E'])
R['double-3~proxima-b'] = iso(X['double-3']['E'], X['proxima-b']['E'])
for s in X: R[s+'~K12a1019'] = iso(X[s]['E'], K)
R['GEGENFALL double-3~K12a1018'] = iso(X['double-3']['E'], snappy.Manifold('K12a1018'))
R['GEGENFALL double-3~trainstation'] = iso(X['double-3']['E'], next(x for x in aus_kurve('trainstation', seeds=(0,)) if 'E' in x)['E'])
print('ERGEBNIS', json.dumps(R))
