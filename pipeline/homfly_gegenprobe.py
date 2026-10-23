#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""homfly_gegenprobe.py — checks a freshly computed HOMFLY against the record
BEFORE it is written into it.

Two witnesses, otherwise nothing is written (the same doctrine as
det_nachtragen.py on 23.08.). The second witness is the **determinant** already
stored in the record — it comes from regina-Alexander or Sage, i.e.
from a different computation than the HOMFLY.

Specialisation (regina prints the AZ convention, the variables are called x=a, y=z):
    Delta(t) = P(a=1, z = t^(1/2) - t^(-1/2))
    t = -1  =>  z = i - (-i) = 2i  =>  det = |P(1, 2i)|
Additionally, where the record holds an Alexander polynomial: |Alexander(-1)| must give
the same value — then there are three witnesses.
"""
import sys, os, json, argparse, cmath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from homfly_eichung import parse_homfly, eval_homfly


def det_from_homfly(homfly_str):
    """det = |P(a=1, z=2i)| — EXACT in integers, never in floating point.

    Floating point is not a convenience here but an error: at y-degree 20
    the terms run up to 2^20 ~ 1e6 and cancel down to a determinant of
    15 — pure cancellation. On 23.08. the floating-point version reported
    8 records as "det mismatch", among them flip_n15_p3_s3 with 265365 instead of
    265335 (off by 0.011 %). Not a single one really deviated.
    Since x=1 and y=2i, P(1,2i) = Σ c_ij · 2^j · i^j, and i^j cycles in
    four steps — that is integral and exact."""
    p = parse_homfly(homfly_str)
    if p is None:
        return None
    re_, im_ = 0, 0
    for (i, j), c in p.items():
        if j < 0:                      # negative z power: not integral
            return None
        w = c * (2 ** j)               # x=1 contributes nothing
        r = j % 4
        if r == 0:   re_ += w
        elif r == 1: im_ += w
        elif r == 2: re_ -= w
        else:        im_ -= w
    from math import isqrt
    q = re_ * re_ + im_ * im_
    w = isqrt(q)
    return float(w) if w * w == q else (q ** 0.5)


def _alex_terms(s):
    """regina writes 'x', Sage 't' — both as {exponent: coefficient}.
    Footgun 15: whoever reads only one notation measures himself.
    In regina, coefficient and variable are separated by a SPACE
    ('3 x^20') — whoever tokenises at spaces tears the term apart
    and again measures himself (happened exactly so in the self-test on 23.08.)."""
    import re
    s = s.strip().replace('*', ' ')
    if not s:
        return None
    var = 'x' if 'x' in s else ('t' if 't' in s else None)
    # split terms ONLY at binary signs, not at an exponent minus
    s = re.sub(r'(?<!\^)([+-])', r'|\1', s)   # '^-1' belongs to the exponent
    out = {}
    for st in s.split('|'):
        st = re.sub(r'\s+', '', st)
        if not st:
            continue
        vz = -1 if st.startswith('-') else 1
        st = st.lstrip('+-')
        if not st:
            return None
        if var and var in st:
            m = re.match(rf'^(\d+)?{var}(?:\^(-?\d+))?$', st)
            if not m:
                return None
            coef_k = int(m.group(1)) if m.group(1) else 1
            e = int(m.group(2)) if m.group(2) else 1
        else:
            if not st.isdigit():
                return None
            coef_k, e = int(st), 0
        out[e] = out.get(e, 0) + vz * coef_k
    return out


def det_from_alexander(alex_str):
    t = _alex_terms(alex_str)
    if not t:
        return None
    return abs(sum(c * ((-1) ** e) for e, c in t.items()))


def record_det(d):
    """The determinant is NOT always in `ALGEBRA.determinant` — there many
    records hold None, while `determinant_alexander` or the
    consistency quartet carries the value (seen on 46 records on 23.08.: the too
    narrow searcher reported 'no det', it was not the mathematics that was red, Footgun 15).
    Order = reliability, first hit wins."""
    pf = d.get('pflichtliste', {})
    raw = d.get('raw_subprocess_outputs', {})
    cand_k = [
        pf.get('ALGEBRA', {}).get('determinant'),
        pf.get('KLASSISCH', {}).get('determinant'),
        pf.get('ALGEBRA', {}).get('determinant_alexander'),
        pf.get('ALGEBRA', {}).get('determinant_sage'),
        (d.get('konsistenz_quartett') or {}).get('det_value'),
        raw.get('sage', {}).get('determinant'),
        raw.get('seifert_zeuge', {}).get('determinant_seifert'),
    ]
    for v in cand_k:
        if v is not None:
            return int(v)
    return None


def record_alexander(d):
    pf = d.get('pflichtliste', {})
    v = pf.get('ALGEBRA', {}).get('alexander_polynomial')
    if v:
        return v
    return d.get('raw_subprocess_outputs', {}).get('regina', {}).get('alexander')


def check_single(rec_file, homfly_str, tolerance=1e-6):
    d = json.load(open(rec_file))
    det_rec = record_det(d)
    dh = det_from_homfly(homfly_str)
    alex = record_alexander(d)
    da = det_from_alexander(alex) if alex else None
    if dh is None:
        return {'ok': False, 'grund': 'homfly_unparsbar'}
    if det_rec is None:
        return {'ok': False, 'grund': 'kein_det_im_record', 'det_homfly': dh}
    ab = abs(dh - det_rec) / max(det_rec, 1)
    ok = ab < tolerance
    r = {'ok': ok, 'det_record': det_rec, 'det_homfly': dh, 'rel': ab,
         'zeugen': 2}
    if da is not None:
        r['det_alexander'] = da
        r['zeugen'] = 3 if da == det_rec else 2
        if da != det_rec:
            r['warnung'] = 'Alexander witness differs from the record det'
    return r


def selbsttest():
    """Counter-case duty (Footgun 16): a checker that can only confirm
    is decoration. Four probes, two of them MUST be rejected."""
    ok = True
    # 1) trefoil: P = 2a^-2 - a^-4 + a^-2 z^2, det must be 3
    d1 = det_from_homfly('x^-2 y^2 + 2 x^-2 - x^-4')
    print(f'  probe 1 trefoil     det={d1:.6f}  expected 3   ',
          '✅' if abs(d1 - 3) < 1e-9 else '🔴')
    ok &= abs(d1 - 3) < 1e-9
    # 2) figure-eight 4_1: P = -a^-2 + 3 - a^2 - z^2  (det 5)
    d2 = det_from_homfly('- x^-2 + 3 - x^2 - y^2')
    print(f'  probe 2 figure-8    det={d2:.6f}  expected 5   ',
          '✅' if abs(d2 - 5) < 1e-9 else '🔴')
    ok &= abs(d2 - 5) < 1e-9
    # 3) COUNTER-CASE: a falsified trefoil must NOT give 3
    d3 = det_from_homfly('x^-2 y^2 + 3 x^-2 - x^-4')
    print(f'  probe 3 COUNTER-CASE det={d3:.6f}  must NOT be 3 ',
          '✅' if abs(d3 - 3) > 1e-9 else '🔴')
    ok &= abs(d3 - 3) > 1e-9
    # 4) Alexander witness, both notations, plus counter-case
    a1 = det_from_alexander('x^2 - x + 1')          # trefoil -> 3
    a2 = det_from_alexander('t^2 - 3 t + 1')        # 4_1       -> 5
    a3 = det_from_alexander('x^2 - x + 2')          # counter-case -> 4, not 3
    print(f'  probe 4 Alexander   x:{a1} t:{a2} counter-case:{a3} ',
          '✅' if (a1 == 3 and a2 == 5 and a3 == 4) else '🔴')
    ok &= (a1 == 3 and a2 == 5 and a3 == 4)
    print('  ⇒', '✅ the cross-check bites' if ok else '🔴 self-test failed')
    return 0 if ok else 1


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--wand', help='result JSON of homfly_wand.py')
    ap.add_argument('--out')
    ap.add_argument('--selbsttest', action='store_true')
    a = ap.parse_args()
    if a.selbsttest:
        sys.exit(selbsttest())
    erg = json.load(open(a.wand))
    good = bad = 0
    out_list = []
    for e in erg:
        if e.get('status') != 'OK' or not e.get('homfly'):
            continue
        r = check_single(e['record'], e['homfly'])
        r.update({'name': e['name'], 'n': e.get('n'), 'c': e.get('c'),
                  'record': e['record'], 'homfly': e['homfly']})
        out_list.append(r)
        if r['ok']:
            good += 1
        else:
            bad += 1
            print(f"  🔴 {e['name'][:40]:<40} {r.get('grund','det-mismatch')} "
                  f"rec={r.get('det_record')} homfly={r.get('det_homfly')}")
    print(f"\n{good} passed, {bad} failed "
          f"(of {len(out_list)} checked).")
    three = sum(1 for r in out_list if r.get('zeugen') == 3)
    print(f"of these with THREE witnesses (det + Alexander + HOMFLY): {three}")
    if a.out:
        json.dump(out_list, open(a.out, 'w'), indent=1)
    sys.exit(0 if bad == 0 else 1)
