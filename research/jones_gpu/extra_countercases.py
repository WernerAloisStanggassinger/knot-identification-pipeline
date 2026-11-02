#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""extra_countercases.py — V(t) = V(1/t) at the census point vs Sage's exact test J(q) = J(1/q), on knots whose answer is known:
chiral but Jones-symmetric (9_42, 10_48, 10_71), amphichiral (8_17, 4_1, 6_3), chiral (3_1).
Usage (Python with SageMath + spherogram): python -I extra_countercases.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spherogram
from jones_tn import jones_at, inv, P
x = 987654321 % P                     # the point of jones_census.py
for nm in ['9_42', '10_48', '10_71', '3_1', '8_17', '4_1', '6_3']:
    L = spherogram.Link(nm); pd = [tuple(X) for X in L.PD_code()]
    a = jones_at(pd, x)[0]; b = jones_at(pd, inv(x))[0]
    J = L.jones_polynomial(new_convention=False); q = J.parent().gen()
    print(nm, 'V(t)==V(1/t) at the point:', a == b, '| Sage J(q)==J(1/q):', J == J.subs({q: 1 / q}))
