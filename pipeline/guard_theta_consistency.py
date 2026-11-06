#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""guard_theta_consistency.py — every stored theta_full against the numerical route.

Why: on 26.07. three records carried a WRONG theta_full (among other reasons because it had been computed from a different
diagram). No guard saw that — the blind-spot register named exactly
this gap ("no auto-re-derive on PD change"). This guard closes it.

Relation: theta_full == ±theta_bare * D(t1)*D(t2)*D(t1*t2), where theta_at_22_7__21_13 is the
bare theta and D the balanced Alexander polynomial. theta_full comes from Theta.sage,
theta_bare from Theta_numeric — two computation routes, one target value.

COUNTER-CASE BUILT IN (footgun 16): a deliberately falsified theta (x2) MUST be noticed.
If it is not noticed, the guard reports red instead of green — a test that rejects nothing checks nothing.
"""
import site_config
import json, os, re, sys
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from large_expr_eval import evaluiere            # noqa: E402

# 🔴 CATCH-UP 2026-08-14 — the coverage had collapsed (20 "not checkable"), and the
# cause was NOT the one assumed here:
#
#  (a) `sys.setrecursionlimit(50000)` never helped. The big theta expressions die in the
#      CPython *COMPILER* ("maximum recursion depth exceeded during compilation"), and that
#      depends on the C stack, not on the interpreter limit. `a + b + c + …` with 20,000 summands
#      becomes a left-deep tree. The fix of 26.07. only moved the limit —
#      and because the theta night run keeps delivering LARGER knots, the number of
#      silently unchecked ones crept up again (17 cases, all T_11_9_* / T_13_*).
#      Cure: `gross_eval.evaluiere` — tokenizer + shunting-yard, iterative, exact in Fraction,
#      no compile(). Calibrated: 755 records that `eval` manages come out BIT-IDENTICAL,
#      0 deviations, plus the 17 that `eval` cannot manage. 36 s for all 772.
#
#  (b) The three Delta=1 knots (`conway_kt`, `conway_cw`, `W_trefoil_untw`) failed on
#      `max()` of an empty list: their Alexander polynomial is "1", so there is no `x^n`.
#      The guard was thus blind precisely on the mathematically most interesting cases.
#      Cure: no x^ term ⟹ g = 0 (constant polynomial, no balancing necessary).
#      Counter-check in the corpus: the highest exponent is even in ALL 785 remaining cases,
#      so `//2` is exact and not rounded.
ROOT = site_config.ROOT
P1, P2 = F(22, 7), F(21, 13)


def alex_bal(a, t):
    s = re.sub(r'(\d)\s+([a-z])', r'\1*\2', a)
    exp = [int(m) for m in re.findall(r'x\^(\d+)', a)]
    g = (max(exp) // 2) if exp else 0        # Delta = 1 has no x^ term
    return evaluiere(s, {"x": t}) / (t ** g)


def am_punkt(s):
    return evaluiere(str(s), {"t1": P1, "t2": P2})


def urteil(rec, grund=None):
    """True = consistent, False = broken, None = not checkable (reason is appended)."""
    try:
        if rec.get('theta_at_22_7__21_13') is None or not rec.get('alex'):
            if grund is not None:
                grund.append((rec.get('stem'), "Field missing (theta_at_22_7__21_13 or alex)"))
            return None
        D = alex_bal(rec['alex'], P1) * alex_bal(rec['alex'], P2) * alex_bal(rec['alex'], P1 * P2)
        soll = F(str(rec['theta_at_22_7__21_13'])) * D
        return abs(am_punkt(rec['theta_full'])) == abs(soll)
    except Exception as e:
        if grund is not None:
            grund.append((rec.get('stem'), "%s: %s" % (type(e).__name__, str(e)[:80])))
        return None


db = json.load(open(site_config.PIPELINE_DIR + "/corpus_db.json"))
recs = [r for r in db["knots"] if r.get('theta_full') is not None]

kaputt, ok, blind, gruende = [], 0, 0, []
for r in recs:
    u = urteil(r, gruende)
    if u is True: ok += 1
    elif u is False: kaputt.append(r.get('stem'))
    else: blind += 1

# --- COUNTER-CASE: a falsified theta must be noticed ---
probe = next((r for r in recs if urteil(r) is True and am_punkt(r['theta_full']) != 0), None)
if probe is None:
    print("🔴 no probe record with theta != 0 — counter-case cannot be run, green would be worthless.")
    sys.exit(8)
gefaelscht = dict(probe); gefaelscht['theta_full'] = "2*(" + str(probe['theta_full']) + ")"
if urteil(gefaelscht) is not False:
    print("🔴 COUNTER-CASE SLIPPED THROUGH — falsified theta was reported as consistent.")
    sys.exit(9)

print(f"{ok} consistent · {len(kaputt)} broken · {blind} not checkable · counter-case rejected ✓")
# Threshold since 14.08. at ZERO (before 10): with the iterative evaluator and the Delta=1 fix
# every record of the corpus is checkable. Every non-check is from now on a finding, not
# noise — and it says right away WHY, instead of staying silent.
if blind > 0:
    print(f"🔴 COVERAGE: {blind} records not checkable (expected 0) — "
          f"silent non-checking is not green. Reasons:")
    for stem, g in gruende[:12]:
        print(f"     {stem}: {g}")
    sys.exit(4)
if kaputt:
    print("BROKEN:", kaputt[:10])
    sys.exit(3)
sys.exit(0)
