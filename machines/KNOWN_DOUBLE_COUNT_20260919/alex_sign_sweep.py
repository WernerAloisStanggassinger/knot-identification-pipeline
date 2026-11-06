#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""SWEEP: does the OVERALL SIGN of the Alexander polynomial separate distinct knots? (19.09.2026)

OCCASION: poly_s08 and kf332_polys08_sh+0.5_tw-1 stand as TWO distinct novel-hyperbolic
knots (S_24_d25721_D8_amphi_sig0 and ...#2), although they agree in EVERYTHING:
det 25721, sig 0, Jones bit-identical, HOMFLY bit-identical (only another term order), genus 9,
fibered, theta = 0, amphi, vol to 14 digits, complements ISOMETRIC (is_isometric_to).
The only difference: `alex` of A is exactly the NEGATIVE of B.

The Alexander polynomial is only defined up to ±t^k. -Δ and +Δ are the SAME polynomial.
But ident() in label_distinct_knots.py:149 compares the STRING.

COUNTER-CASE (Hard-Rule 3a):
  (a) POSITIVE: poly_s08 / kf332_polys08 MUST be found.
  (b) NEGATIVE: the sweep must not find everything — knots with different det/sig/vol
      must never come together, and a knot never with itself.
  (c) A found pair must additionally agree in Jones OR HOMFLY, otherwise
      the negated alex is a REAL difference and not a convention artifact.
"""
import site_config
import json
import re
import sys
from collections import defaultdict

DB = (site_config.PIPELINE_DIR + '/corpus_db.json')


def polyterme(s):
    """'-x^18 + 13 x^17 - ... - 1' -> {monomial: coefficient}. Format without '*'."""
    if not s or not isinstance(s, str):
        return None
    out = {}
    for t in re.findall(r'[+-]?\s*[^+-]+', s):
        t = t.strip()
        if not t:
            continue
        vz = -1 if t.startswith('-') else 1
        t = t.lstrip('+- ').strip()
        m = re.match(r'^(\d+)?\s*(.*)$', t)
        if not m:
            return None
        k = int(m.group(1)) if m.group(1) else 1
        out[m.group(2).strip()] = out.get(m.group(2).strip(), 0) + vz * k
    return {k: v for k, v in out.items() if v != 0}


def ist_negiert(a, b):
    ta, tb = polyterme(a), polyterme(b)
    if not ta or not tb:
        return False
    return ta == {k: -v for k, v in tb.items()}


def gleich_als_menge(a, b):
    """Same polynomial, only another term order?"""
    ta, tb = polyterme(a), polyterme(b)
    if ta is None or tb is None:
        return a == b
    return ta == tb


def main():
    db = json.load(open(DB))
    reps = [m for m in db['knots'] if m.get('is_distinct_representative')]
    print(f"distinct representatives: {len(reps)}")

    # bucket = ident() key WITHOUT alex
    buckets = defaultdict(list)
    for m in reps:
        v = m.get('vol')
        volk = round(v, 6) if isinstance(v, (int, float)) and v > 0.5 else None
        buckets[(str(m.get('det')), str(m.get('sig')),
                 str(m.get('theta_at_22_7__21_13')), volk)].append(m)

    mehr = {k: v for k, v in buckets.items() if len(v) > 1}
    print(f"Buckets (det, sig, theta, vol@6) with >1 distinct knot: {len(mehr)}")

    treffer = []
    for k, ms in sorted(mehr.items()):
        for i, a in enumerate(ms):
            for b in ms[i + 1:]:
                if ist_negiert(a.get('alex'), b.get('alex')):
                    treffer.append((a, b, 'alex NEGIERT'))
                elif gleich_als_menge(a.get('alex'), b.get('alex')) \
                        and str(a.get('alex')) != str(b.get('alex')):
                    treffer.append((a, b, 'alex nur umsortiert'))

    print("\n" + "=" * 78)
    print(f"HITS: {len(treffer)} pairs of distinct knots that differ ONLY in the "
          f"overall Alexander sign (or in the spelling)")
    print(f"  (since the fix of 19.09.2026 0 is the RIGHT value — regression guard)")
    print("=" * 78)
    for a, b, art in treffer:
        j = gleich_als_menge(a.get('jones'), b.get('jones'))
        h = gleich_als_menge(a.get('homfly'), b.get('homfly'))
        print(f"\n[{art}]  class: {a.get('realizes_class')} / {b.get('realizes_class')}")
        print(f"  A {a.get('distinct_knot_id'):40s} stem={a.get('stem')}")
        print(f"  B {b.get('distinct_knot_id'):40s} stem={b.get('stem')}")
        print(f"  det={a.get('det')} sig={a.get('sig')} vol={a.get('vol')} "
              f"theta={str(a.get('theta_at_22_7__21_13'))[:20]}")
        print(f"  Jones equal: {j}   HOMFLY equal: {h}   "
              f"{'✅ convention artifact' if (j and h) else '⚠️ NOT only convention — real difference'}")
        print(f"  n_realisierungen: {a.get('n_realisierungen')} / {b.get('n_realisierungen')}")

    # ---- COUNTER-CASES
    print("\n" + "=" * 78)
    print("COUNTER-CASES")
    print("=" * 78)
    # ⚠️ SINCE THE REPAIR (19.09.2026, Disc #809) THIS SCRIPT IS A REGRESSION
    # GUARD, NO LONGER A SEARCH RUN. Before the fix it had to FIND poly_s08 / kf332_polys08
    # (positive probe). Since ident() compares the alex normal form, the two are ONE
    # knot — the right finding is now ZERO hits. Whoever finds something here again
    # has either lost the normal form or has a new case of the same class.
    # The positive probe therefore runs against the NORMAL FORM instead of the sweep.
    import sys as _s
    _s.path.insert(0, (site_config.TOOLS_DIR))
    from alexander_normal_form import canonical as _kan
    _A = [m for m in reps if m.get('stem') == 'poly_s08']
    _B = [m for m in reps if m.get('stem') == 'kf332_polys08_sh+0.5_tw-1']
    _alle = json.load(open(DB))['knots']
    _A = _A or [m for m in _alle if m.get('stem') == 'poly_s08']
    _B = _B or [m for m in _alle if m.get('stem') == 'kf332_polys08_sh+0.5_tw-1']
    a_ok = bool(_A and _B) and _kan(_A[0].get('alex')) == _kan(_B[0].get('alex')) \
           and _A[0].get('distinct_knot_id') == _B[0].get('distinct_knot_id')
    print(f"  (a) POSITIVE poly_s08 / kf332_polys08 carry ONE distinct_knot_id "
          f"and the same alex normal form: {a_ok} "
          f"{'✅ (the fix holds)' if a_ok else '🔴 REGRESSION — the fix is gone'}")
    selbst = any(x is y for x, y, _ in treffer)
    print(f"  (b1) NEGATIVE no knot with itself: {not selbst} {'✅' if not selbst else '🔴'}")
    versch = [(x, y) for x, y, _ in treffer
              if (x.get('det'), x.get('sig')) != (y.get('det'), y.get('sig'))]
    print(f"  (b2) NEGATIVE no pair with different det/sig: {not versch} "
          f"{'✅' if not versch else '🔴'}")
    print(f"  (b3) NEGATIVE sweep does not find everything: {len(treffer)} of "
          f"{len(reps)*(len(reps)-1)//2} possible pairs "
          f"{'✅' if len(treffer) < 20 else '🔴 implausibly many'}")
    nicht_konv = [(x, y) for x, y, _ in treffer
                  if not (gleich_als_menge(x.get('jones'), y.get('jones'))
                          and gleich_als_menge(x.get('homfly'), y.get('homfly')))]
    print(f"  (c) every pair also equal in Jones AND HOMFLY: {not nicht_konv} "
          f"{'✅' if not nicht_konv else '⚠️ '+str(len(nicht_konv))+' pair(s) NOT'}")

    print("\n" + "=" * 78)
    from collections import Counter
    kl = Counter(a.get('realizes_class') for a, b, _ in treffer)
    print(f"AFFECTED CLASSES: {dict(kl)}")
    print("=" * 78)
    return treffer


if __name__ == '__main__':
    main()
