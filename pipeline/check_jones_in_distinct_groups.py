#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Sweep: do all members of a distinct_knot_id carry the same Jones polynomial?

OCCASION (19.09.2026, sector Z of the attack swarm): su_4_4_tm1 and su_8_2_tm1 carry
the same distinct_knot_id, but have Jones SPAN 32 against 40 — two different knots.
Both non-hyperbolic, therefore the isometry cure of 18.09. (Disc #798, hyperbolic
groups only) could not see them.

WHY THE SPAN AND NOT THE POLYNOMIAL:
There are two string conventions in the corpus ('x^k' and 'v**k*z**m'), and one of them also mirrors. A comparison of the polynomial STRINGS therefore reports false red. The SPAN
(max exponent - min exponent) is immune:
  * independent of the variable naming,
  * independent of mirroring (V(Kbar)(x) = V(K)(1/x) merely reverses the exponents),
  * independent of a global shift.
Different span ==> different knots. That is a proof, not a heuristic.
Equal span proves NOTHING — therefore the polynomial is additionally compared afterwards,
but a difference THERE counts only as suspicion (convention!), not as verdict.
"""
import site_config
import json, re, sys
from collections import defaultdict

DB = (site_config.PIPELINE_DIR + '/corpus_db.json')

def exponenten(s):
    """All exponents of the first variable from a polynomial string.

    🔴 TWO SPELLINGS IN THE CORPUS, and the second reverses the sign:
          'x^-12'   negative exponent
          '1/x^12'  fraction  ->  likewise exponent -12
    Whoever searches only for 'x^(-?\\d+)' reads +12 instead of -12 in the fraction case and measures a
    span that is too small. That is exactly where the first run of this checker failed on 19.09.2026
    (two identical lj records reported as "different").
    Therefore it is read TERM BY TERM and every denominator part counted negative.
    """
    if not isinstance(s, str) or not s.strip():
        return None
    t = s.replace("**", "^").replace(" ", "")
    m = re.search(r"([A-Za-z])", t)
    if m is None:
        return [0]                      # constant polynomial
    var = m.group(1)
    ex = []
    # split into summands: at + / -, except in the exponent (x^-12)
    for stueck in re.split(r"(?<![\^\*(])[+-]", t):
        if not stueck:
            continue
        zaehler, _, nenner = stueck.partition("/")
        e = 0
        for teil, vz in ((zaehler, +1), (nenner, -1)):
            if not teil:
                continue
            tr = re.search(re.escape(var) + r"\^(-?\d+)", teil)
            if tr:
                e += vz * int(tr.group(1))
            elif re.search(re.escape(var) + r"(?![\^\w])", teil):
                e += vz * 1
        ex.append(e)
    return ex or None

def spanne(s):
    e = exponenten(s)
    return (max(e) - min(e)) if e else None

def selbsttest():
    print("=" * 78); print("SELF-TEST — counter-cases before the sweep"); print("=" * 78)
    faelle = [
        ("x^16 - 2 x^14 + 39 - 2 x^-14 + x^-16", 32, "symmetric, as in the corpus"),
        ("x^20 + 1 - x^-20",                     40, "different span"),
        # mirroring must NOT change the span:
        ("x^4 - x^2 + 1",                         4, "chiral, original"),
        ("x^-4 - x^-2 + 1",                       4, "same, mirrored -> same span"),
        # other variable, other spelling -> same span
        ("v**4-v**2+1",                           4, "v instead of x, ** instead of ^"),
        ("1",                                     0, "constant polynomial"),
    ]
    ok = True
    for s, soll, was in faelle:
        ist = spanne(s)
        gut = ist == soll
        ok &= gut
        print(f"[{'OK ' if gut else 'ROT'}] {was:<38} span {ist} (expected {soll})")
    # falsified probe: two polynomials of equal span, different coefficients
    gleich = spanne("x^4-x^2+1") == spanne("x^4-7x^2+1")
    print(f"[{'OK ' if gleich else 'ROT'}] equal span with different coefficients -> "
          f"span alone does NOT decide (hence second stage)")
    ok &= gleich
    print(f"Self-test: {'ALL GREEN' if ok else '🔴 RED'}\n")
    return ok

def sweep():
    db = json.load(open(DB))
    grp = defaultdict(list)
    for k in db["knots"]:
        d = k.get("distinct_knot_id")
        if d: grp[d].append(k)
    mehr = {d: v for d, v in grp.items() if len(v) > 1}
    print("=" * 78)
    print(f"SWEEP over {len(grp)} distinct groups, of which {len(mehr)} with more than one record")
    print("=" * 78)
    harte, weiche, ohne = [], [], 0
    for d, v in sorted(mehr.items()):
        sp = {}
        for k in v:
            s = spanne(k.get("jones"))
            sp.setdefault(s, []).append(k["stem"])
        if None in sp:
            ohne += 1
        echte = {s: n for s, n in sp.items() if s is not None}
        if len(echte) > 1:
            harte.append((d, echte, v))
        else:
            pol = {}
            for k in v:
                pol.setdefault(str(k.get("jones")), []).append(k["stem"])
            if len(pol) > 1:
                weiche.append((d, pol, v))
    print(f"\n🔴 HARD (different Jones SPAN in a group = provably different knots): {len(harte)}")
    for d, e, v in harte:
        hyp = {k["stem"]: k.get("is_hyperbolic") for k in v}
        print(f"  {d}")
        for s, namen in sorted(e.items()):
            for n in namen:
                k = [x for x in v if x["stem"] == n][0]
                print(f"     span {s:>4}  c={k.get('c')}  alt={k.get('alternating')} "
                      f"hyp={k.get('is_hyperbolic')}  det={k.get('det')} sig={k.get('sig')}  {n}")
    print(f"\n🟡 SOFT (same span, but polynomial string different — may be convention): {len(weiche)}")
    for d, p, v in weiche[:15]:
        print(f"  {d}: {[n for ns in p.values() for n in ns]}")
    print(f"\nGroups with at least one record without Jones: {ohne}")
    return harte, weiche

if __name__ == "__main__":
    if not selbsttest(): sys.exit(1)
    sweep()
