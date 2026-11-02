# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Usage (sage_stable env): python -I check_torus_alexander.py p q polygon_to_pd_output.log
# Compare a printed Alexander polynomial with the torus-knot formula (t^(pq)-1)(t-1)/((t^p-1)(t^q-1)).
import sys
from sage.all import PolynomialRing, ZZ
R = PolynomialRing(ZZ, 't'); t = R.gen()
p, q, log = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
s = [l for l in open(log) if l.startswith('alexander ')][-1].split(' ', 1)[1].strip()
A = R(s)
F = R((t**(p*q) - 1) * (t - 1) // ((t**p - 1) * (t**q - 1)))
print('degree', A.degree(), 'formula degree', F.degree(), 'equal:', A == F or A == -F)
