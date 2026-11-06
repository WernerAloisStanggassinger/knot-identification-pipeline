# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
import random, signal, snappy
random.seed(1)
H = snappy.HTLinkExteriors
M = H['K13n4587']; print('K13n4587 vol', M.volume(), 'solution', M.solution_type())
ok = miss = nonhyp = 0; misses = []
knoten = [n for n in (H[i].name() for i in random.sample(range(len(H)), 3000)) if n.startswith('K')][:300]
for n in knoten:
    signal.alarm(60)
    M = H[n]
    if M.solution_type() != 'all tetrahedra positively oriented' or M.volume() < 1: nonhyp += 1; continue
    for _ in range(3): M.randomize()
    r = H.identify(M)
    if r and any(x.name().startswith(n) for x in [r] if hasattr(x, 'name')) or (r and n in str(r)): ok += 1
    else: miss += 1; misses.append((n, str(r), M.solution_type()))
print('HT knots checked', ok + miss, 'found', ok, 'missed', miss, '(non-hyp/degenerate skipped', nonhyp, ')')
print(misses[:10])
