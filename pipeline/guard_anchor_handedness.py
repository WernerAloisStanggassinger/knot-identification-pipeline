#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""guard_anchor_handedness.py — check anchor 1 CHIRALITY-SENSITIVELY.

Why this guard exists: det = 3^6 and Alexander (t^2-t+1)^6 — the two bit-exact
receipts of the anchor — are CHIRALITY-BLIND. Exactly for that reason knots-of-love stood for two months
wrongly as 3_1#3_1#3_1#3_1#3_1#3_1 in CLAUDE.md, in the record AND publicly on the site,
while the existing Alexander guard glowed green the whole time.

The check is bit-exact: Jones is multiplicative under connected sum and distinguishes
mirror images. V(knots-of-love) must be V(3_1)^3 * V(m3_1)^3.

COUNTER-CASE BUILT IN (footgun 16): the guard ALSO computes V(3_1)^6 and must
REJECT this variant. If it does not, its green is worthless and it reports red.
Everything fresh from the PD of the record with regina — no stored field, no table.
"""
import site_config
import glob, re, sys
from regina import Link
import sympy as sp

x = sp.symbols('x')


def jones_expr(PD):
    s = Link.fromPD(PD).jones().str().replace('^', '**')
    s = re.sub(r'(\d)\s+x', r'\1*x', s)
    s = re.sub(r'\bx\*\*-(\d+)', r'x**(-\1)', s)
    return sp.expand(sp.sympify(s, locals={'x': x}))


cands = glob.glob((site_config.ROOT + '/own_knots/knots-of-love/*pflicht*complete.json'))
if not cands:
    print("knots-of-love pflicht missing"); sys.exit(4)
m = re.search(r'"pd_list_str"\s*:\s*"((?:[^"\\]|\\.)*)"', open(cands[0]).read())
if not m:
    print("no pd_list_str -> BLIND"); sys.exit(5)
PD = [tuple(int(q) for q in g) for g in
      re.findall(r'[\(\[]\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*[\)\]]', m.group(1))]

V = jones_expr(PD)
VR = jones_expr([(1, 4, 2, 5), (3, 6, 4, 1), (5, 2, 6, 3)])   # trefoil, fresh
VL = sp.expand(VR.subs(x, 1 / x))                              # mirror image
drei_und_drei = sp.expand(VR**3 * VL**3)
sechs_gleich = sp.expand(VR**6)

ok = sp.simplify(sp.expand(V - drei_und_drei)) == 0
gegenfall_abgelehnt = sp.simplify(sp.expand(V - sechs_gleich)) != 0

if not gegenfall_abgelehnt:
    print("🔴 COUNTER-CASE NOT REJECTED — the guard cannot separate 3+3 from 6+0, green would be worthless.")
    sys.exit(8)
if not ok:
    print(f"🔴 DRIFT: Jones != V(3_1)^3 * V(m3_1)^3  (got {str(V)[:70]})")
    sys.exit(7)
print("MATCH 3_1#3_1#3_1#m3_1#m3_1#m3_1 (Jones bit-exact) · counter-case 6x same-handed rejected")
sys.exit(0)
