#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Knot-Identification Pipeline — Werner Alois Stanggassinger & Claude Code (Anthropic).
"""span_lower_bound.py — c >= span V(K), without any table.

PURPOSE (2026-09-04): the prime knots are enumerated and classified up to AND INCLUDING 20 crossings
(Thistlethwaite, AGT 25 (2025) 329-344: 1,847,319,428 items; the tabulation
was made in 2018, independently twice by Thistlethwaite and Burton). The
DATA is public only up to c=19 (regina-normal.github.io/data.html, 352,152,252
knots). For c=20 the table exists, we do not have it — so we cannot
compare against it.

THE WAY OUT that needs no table: Kauffman-Murasugi-Thistlethwaite gives
    span V(K) <= c(K)
for EVERY knot, with equality exactly for reduced alternating diagrams.
So span V is a LOWER bound for the crossing number.

🔴 CORRECTION 2026-09-20 (Disc #814): "span V >= 21 ==> beyond EVERY enumeration" was
TOO STRONG. There are TWO limits, not one:
   all prime knots          complete up to c = 20  (Thistlethwaite, AGT 25 (2025))
   ALTERNATING prime knots  complete up to c = 23 published (Rankin/Flint/Schermann,
                            JKTR 13 (2004) 57-150), computed up to c = 24 (Fontaine/
                            Rankin/Flint 2007, OEIS A002864 a(24) = 128,564,665,125)
Verbatim from Thistlethwaite himself, AGT 25 (2025) p. 2: "This has indeed been
accomplished quite dramatically for n <= 23 [Rankin et al. 2004]."
A knot is therefore beyond if and only if
   c >= 25   OR   c >= 21 AND K proven NON-alternating.
Receipt: konversationen/alternierende_enumeration_20260920/

⚠️ The bound is one-sided: span V < 21 proves NOTHING. It then only says that
this route does not suffice for this knot, not that the knot is small.

CALIBRATION (--calibrate-census): on records with a PROVEN crossing number (reduced alternating,
c_exact) span == c must come out. If the convention does not fit, the number is misinterpreted
and everything built on it is worthless (Footgun 15/12).
"""
import json, re, sys, argparse
import site_config
from collections import Counter

DB = '@PIPELINE_DIR@/corpus_db.json'.replace('@PIPELINE_DIR@', site_config.PIPELINE_DIR)


def exponents(poly):
    """All x exponents of a Regina polynomial string. 'x^10 - 6 x^8 + ... + 11 - x^-22'."""
    if not isinstance(poly, str) or not poly.strip():
        return None
    s = poly.replace(' ', '')
    ex = []
    # terms with x^k (also negative), terms with bare x, and constant terms
    for m in re.finditer(r'x\^(-?\d+)', s):
        ex.append(int(m.group(1)))
    for m in re.finditer(r'(?<!\^)(?<!\d)x(?!\^)', s):
        ex.append(1)
    # a constant term (no x) counts as exponent 0
    rest = re.sub(r'[+-]?\d*x(\^-?\d+)?', '', s)
    if re.search(r'\d', rest):
        ex.append(0)
    return ex or None


def span_in_t(poly):
    """span V(K) in the variable t. Regina writes Jones in x with x^2 = t."""
    ex = exponents(poly)
    if not ex:
        return None
    span_x = max(ex) - min(ex)
    if span_x % 2:
        return None          # odd exponents => different convention, do not guess
    return span_x // 2


# The two enumeration limits (Disc #814). Numbers, not feelings.
ENUM_ALL = 20            # Thistlethwaite, AGT 25 (2025) 329-344
ENUM_ALTERNATING = 24    # Rankin/Flint/Schermann 2004 (c<=23) + Fontaine 2007 (c=24)


def alex_coefficients(poly):
    """Coefficients of a Regina Alexander polynomial as {exponent: coefficient}.

    Format: '4 x^14 - 27 x^13 + 49 x^12 + 158 x^11 - ... + 4'. Returns None
    if the string cannot be parsed completely — better to say nothing than to guess.
    """
    if not isinstance(poly, str) or not poly.strip():
        return None
    t = poly.replace(' ', '')
    if not re.fullmatch(r'[-+0-9x^]*', t):
        return None
    coeff_ = {}
    for m in re.finditer(r'([+-]?)(\d*)(?:x(?:\^(-?\d+))?)?', t):
        vz, number, exp = m.group(1), m.group(2), m.group(3)
        if not (vz or number or m.group(0)):
            continue
        if m.group(0) == '':
            continue
        hat_x = 'x' in m.group(0)
        if not hat_x and not number:
            continue
        c = int(number) if number else 1
        if vz == '-':
            c = -c
        e = int(exp) if exp is not None else (1 if hat_x else 0)
        coeff_[e] = coeff_.get(e, 0) + c
    return coeff_ or None


def alexander_is_alternating(poly):
    """Crowell-Murasugi test: does Delta satisfy the condition (-1)^(i+j) c_i c_j >= 0?

    Verbatim (Boden-Karimi, arXiv:2204.09767 p.1, here on disk):
      "In [Cro59a, Mur58], Crowell and Murasugi independently proved that the Alexander
       polynomial of an alternating link is alternating. Here, a Laurent polynomial
       Delta_L(t) = sum c_i t^i is said to be alternating if its coefficients satisfy
       (-1)^(i+j) c_i c_j >= 0."
    The condition allows ZERO coefficients (">= 0"); there are alternating knots with
    gaps. It is violated only when two coefficients carry the wrong relative sign.
    True/False/None (None = not parsable).
    """
    k = alex_coefficients(poly)
    if not k:
        return None
    nz = {e: c for e, c in k.items() if c}
    if not nz:
        return None
    e0 = min(nz)
    eps = 1 if (nz[e0] > 0) == (e0 % 2 == 0) else -1
    for e, c in nz.items():
        expected = eps if e % 2 == 0 else -eps
        if (c > 0) != (expected > 0):
            return False
    return True


def jones_is_alternating(poly):
    """Thistlethwaite test: does V(K) satisfy the sign condition (-1)^(i+j) c_i c_j >= 0?

    Verbatim (Boden-Karimi, arXiv:2204.09767 p.11, here on disk):
      "In [Thi87], Thistlethwaite proved that the Jones polynomial V_L(t) of any non-split,
       alternating classical link L is alternating."
    Regina writes V in x with t = x^2; half-integer t powers ==> None instead of guessing.
    Pure EXCLUSION: an alternating V proves nothing in the opposite direction.
    """
    k = alex_coefficients(poly)
    if not k:
        return None
    nz = {}
    for e, c in k.items():
        if not c:
            continue
        if e % 2:
            return None
        nz[e // 2] = c
    if not nz:
        return None
    e0 = min(nz)
    eps = 1 if (nz[e0] > 0) == (e0 % 2 == 0) else -1
    for e, c in nz.items():
        expected = eps if e % 2 == 0 else -eps
        if (c > 0) != (expected > 0):
            return False
    return True


def _proven_alternating():
    """Knots for which an alternating diagram has been CONSTRUCTED and the identity is receipted via
    the isometry of the complements (Gordon-Luecke). File instead of hand fix,
    so that the provenance stays readable — like proven_distinctions.json."""
    import json as _j, os as _os
    p = ('@PIPELINE_DIR@/'
         'proven_alternating.json'.replace('@PIPELINE_DIR@', site_config.PIPELINE_DIR))
    if not _os.path.exists(p):
        return {}
    try:
        return _j.load(open(p))
    except Exception:
        return {}


def alternating_status(r, c_lower=None):
    """Is the knot alternating?  True / False / None (= undecided).

    Both inferences from Kauffman-Murasugi-Thistlethwaite, valid for PRIME knots
    (secured in the corpus via hyperbolicity: hyperbolic ==> prime):

      (1) A record carries `alternating=True`  ==>  an alternating diagram
          EXISTS  ==>  K is alternating. Holds at any diagram size.
      (2) span V(K) == c_ub (size of a REALISED diagram)  ==>  c = c_ub exactly
          and K alternating, because a prime NON-alternating diagram would have
          span <= c-1.
      (3) span V(K) < c_unten (proven lower bound from another route)  ==>
          K is NOT alternating, because otherwise span V = c >= c_unten.

    ⚠️ `alternating=False` alone proves NOTHING AT ALL: Regina measures the flag after
    L.simplify() on the diagram then at hand, while `crossing_number_min` can stem from an
    earlier, better run through the latch "never raise c" (Disc #562).
    Flag and bound then belong to DIFFERENT diagrams. This is exactly what fourteen cases
    slipped past on 19.09.2026 (Disc #814).
    """
    sv = span_in_t(r.get('jones')) if r.get('jones') else None
    c_ub = r.get('c')
    _proven = _proven_alternating().get(str(r.get('stem') or ''))
    if _proven:
        return bool(_proven.get('alternierend', True))
    if r.get('alternating') is True:
        # Latch: alternating AND span < lower bound would be a contradiction in the record.
        if sv is not None and isinstance(c_lower, int) and sv < c_lower:
            raise AssertionError(
                f"{r.get('stem')}: alternating=True, but span V={sv} < proven bound "
                f"c>={c_lower} — for alternating knots span V = c holds. Check the record.")
        return True
    if sv is not None and isinstance(c_ub, int) and sv == c_ub:
        return True
    if sv is not None and isinstance(c_lower, int) and sv < c_lower:
        return False
    # (3b) Thistlethwaite: the JONES polynomial of an alternating link is alternating.
    if jones_is_alternating(r.get('jones')) is False:
        return False
    # (4) Crowell-Murasugi: Delta of an alternating knot is alternating.
    #     If it violates the sign condition, K is not alternating. Pure
    #     EXCLUSION — an alternating Delta proves nothing in the opposite direction.
    if alexander_is_alternating(r.get('alex')) is False:
        return False
    # (5) Bankwitz 1930: det(L) >= c(L) for every non-split alternating link
    #     (likewise verbatim in Boden-Karimi p.1). det < lower c bound ==> not alternating.
    try:
        _det = abs(int(r.get('det')))
    except (TypeError, ValueError):
        _det = None
    if _det is not None and isinstance(c_lower, int) and _det < c_lower:
        return False
    return None


def beyond_every_enumeration(r, c_lower):
    """Does K lie provably beyond BOTH enumeration limits?  True/False.

    c >= 25                                   ==> yes
    c >= 21 AND proven non-alternating        ==> yes
    otherwise                                 ==> no (enumerated range or open)
    """
    if not isinstance(c_lower, int) or c_lower <= ENUM_ALL:
        return False
    if c_lower > ENUM_ALTERNATING:
        return True
    return alternating_status(r, c_lower) is False


def laden():
    db = json.load(open(DB))
    ks = db['knots'] if 'knots' in db else db['records']
    return list(ks.values()) if isinstance(ks, dict) else ks


def calibrate_census():
    """span == c must hold for a proven crossing number (reduced alternating)."""
    ks = laden()
    print('CALIBRATION span_untere_schranke — the convention on proven cases\n')
    errors = 0

    print('1) REDUCED ALTERNATING: span V must be EXACTLY c (Kauffman-Murasugi-Thistlethwaite)')
    hits = misses = 0
    examples = []
    for r in ks:
        if not r.get('alternating') or not isinstance(r.get('c'), int) or not r.get('jones'):
            continue
        s = span_in_t(r['jones'])
        if s is None:
            continue
        if s == r['c']:
            hits += 1
        else:
            misses += 1
            if len(examples) < 6:
                examples.append((r.get('stem'), r['c'], s))
    ok = misses == 0 and hits > 20
    print(f'   {"OK " if ok else "ERROR"} {hits} records hit span == c, {misses} deviate')
    for st, c, s in examples:
        print(f'        {st[:44]:<46} c={c:<5} span={s}')
    errors += (not ok)

    print('\n2) THE BOUND MUST HOLD — span <= c for ALL records with c')
    violated = []
    for r in ks:
        if not isinstance(r.get('c'), int) or not r.get('jones'):
            continue
        s = span_in_t(r['jones'])
        if s is not None and s > r['c']:
            violated.append((r.get('stem'), r['c'], s))
    print(f'   {"OK " if not violated else "ERROR"} {len(violated)} violations of span <= c')
    for st, c, s in violated[:8]:
        print(f'        {st[:44]:<46} c={c:<5} span={s}   ← span GREATER than c')
    errors += bool(violated)

    print('\n3) COUNTER-CASE — the parser must not silently deliver a number on nonsense')
    for junk in ('', 'no formula', None, 'x^1 + x^2'):
        s = span_in_t(junk)
        mark = 'None' if s is None else str(s)
        expected = s is None
        print(f'   {"OK " if expected else "ERROR"} {str(junk)[:20]:<22} -> {mark}')
        errors += (not expected)

    print('\n' + ('ALL COUNTER-CASES PASSED' if errors == 0 else f'{errors} ERRORS'))
    return errors


def balance():
    ks = laden()
    nov = [r for r in ks if r.get('is_distinct_representative')
           and r.get('realizes_class') == 'novel_hyperbolic']
    print(f'distinct novel-hyperbolic knots: {len(nov)}\n')
    beyond20, beyond19, open_, without = [], [], [], []
    for r in nov:
        s = span_in_t(r.get('jones')) if r.get('jones') else None
        if s is None:
            without.append(r); continue
        r['_span'] = s
        if s >= 21:
            beyond20.append(r)
        elif s >= 20:
            beyond19.append(r)
        else:
            open_.append(r)
    # Disc #814: span >= 21 means c >= 21 — that is NOT the same as "beyond every
    # enumeration". For alternating knots Rankin/Flint/Schermann reaches up to c = 24.
    alt_inside = [r for r in beyond20
                if alternating_status(r, r.get('c_lower_bound_proven')) is not False
                and (r.get('c_lower_bound_proven') or 0) <= ENUM_ALTERNATING]
    print(f'  span >= 21  ⟹ c >= 21 PROVEN                                : {len(beyond20)}')
    print(f'     of these alternating/open at c <= 24 (Rankin/Flint 2004)  : {len(alt_inside)}'
          f'   ← NOT beyond')
    print(f'     ⟹ of these really beyond BOTH limits                     : {len(beyond20) - len(alt_inside)}')
    print(f'  span == 20  ⟹ c >= 20, only the c=20 table could have it    : {len(beyond19)}')
    print(f'  span <= 19  ⟹ this route does not suffice (says NOTHING on c): {len(open_)}')
    print(f'  without Jones                                                : {len(without)}')
    print(f'\n  ⟹ with certificate beyond the enumeration: {len(beyond20)}/{len(nov)}')
    print(f'  ⟹ reliant on the census comparison        : {len(beyond19)+len(open_)+len(without)}')
    if open_:
        print('\n  the smallest span values (there the comparison is most urgent):')
        for r in sorted(open_, key=lambda r: r['_span'])[:15]:
            print(f'      span={r["_span"]:<4} c_ub={str(r.get("c")):<5} det={str(r.get("det"))[:12]:<13} {r.get("distinct_knot_id")}')
    return beyond20, beyond19, open_, without


def calibrate_alternating():
    """Counter-cases for alternierend_status()/jenseits_jeder_enumeration() (Disc #814).

    The test must BITE: every case has a known answer, and three of them would have been
    answered wrongly under the old logic.
    """
    ks = laden()
    by = {r.get('stem'): r for r in ks if r.get('is_distinct_representative')}
    errors = 0

    def probe(name, actual, expected_, why):
        nonlocal errors
        ok = actual == expected_
        errors += (not ok)
        print(f"   {'OK ' if ok else 'ERROR'} {name:<34} -> {str(actual):<5} (expected {expected_})  {why}")

    # (1) Flag True on the c diagram
    r = by.get('pz_3_5_13')
    probe('pz_3_5_13 (Flag True, span==c)', alternating_status(r, 21), True,
          'Pretzel P(3,5,13), c=21=span')
    # (2) Flag False on a LARGER diagram, span == c_ub -> alternating anyway
    r = by.get('fleet_n5_3_3_5_5_5')
    probe('fleet_n5_3_3_5_5_5 (Flag False!)', alternating_status(r, 21), True,
          'span 21 == c_ub 21, flag belongs to the 25-crossing diagram')
    # (3) span < lower bound -> proven NOT alternating
    r = by.get('S_24_d2197_Z3_chiral__kf65_helix_twist_t1.25')
    probe('kf65_helix (span 16 < L 24)', alternating_status(r, 24), False,
          'were it alternating, span = c >= 24 would hold')
    # (4) span == L < c_ub -> undecidable
    r = by.get('tw_k4')
    probe('tw_k4 (span 22 == L, c_ub 26)', alternating_status(r, 22), None,
          'no contradiction in either direction')
    # (5) Contradiction latch must strike
    synthetic = {'stem': 'kunstfall', 'alternating': True,
             'jones': (by.get('pz_3_5_13') or {}).get('jones'), 'c': 21}
    try:
        alternating_status(synthetic, 99)
        probe('contradiction alternating+span<L', 'no abort', 'AssertionError', 'latch')
    except AssertionError:
        probe('contradiction alternating+span<L', 'AssertionError', 'AssertionError',
              'alternating ==> span = c >= L')
    # (6) the enumeration limits themselves
    probe('beyond: L=25, alternating', beyond_every_enumeration(by.get('pz_3_5_13'), 25), True,
          'c>=25 lies above BOTH limits')
    probe('beyond: L=21, alternating', beyond_every_enumeration(by.get('pz_3_5_13'), 21), False,
          'alternating and c<=24 ==> Rankin/Flint has it')
    probe('beyond: L=24, non-alt.',
          beyond_every_enumeration(by.get('S_24_d2197_Z3_chiral__kf65_helix_twist_t1.25'), 24), True,
          'non-alternating, c>=21 ==> only the c<=20 limit counts')
    probe('beyond: L=20', beyond_every_enumeration(by.get('pz_3_5_13'), 20), False,
          'in the completely enumerated range')
    # (7) Crowell-Murasugi + Bankwitz as EXCLUSION routes, calibrated on the whole corpus
    alt_recs = [r for r in ks if r.get('is_distinct_representative') and r.get('alternating')]
    cm_bad = [r for r in alt_recs if alexander_is_alternating(r.get('alex')) is False]
    cm_none = [r for r in alt_recs if alexander_is_alternating(r.get('alex')) is None]
    bw_bad = [r for r in alt_recs
              if r.get('det') and isinstance(r.get('c'), int) and abs(int(r['det'])) < r['c']]
    probe('Crowell-Murasugi, all alternating', (len(cm_bad), len(cm_none)), (0, 0),
          f'{len(alt_recs)} records: Delta MUST be alternating')
    probe('Bankwitz det>=c, the same', len(bw_bad), 0, 'det(L) >= c(L) for alternating')
    j_bad = [r for r in alt_recs if jones_is_alternating(r.get('jones')) is False]
    j_none = [r for r in alt_recs if jones_is_alternating(r.get('jones')) is None]
    probe('Thistlethwaite, all alternating', (len(j_bad), len(j_none)), (0, 0),
          'Jones MUST be alternating')
    # Counter-case to the Jones route: it must bite on a DEMONSTRABLY non-alternating one
    _kf = next((r for r in ks if r.get('stem') == 'kf354_spiral_amp0.4_freq4'), None)
    if _kf:
        probe('Jones route bites (kf354)', jones_is_alternating(_kf.get('jones')), False,
              'HFK rank 69147 != det 11163 confirms it independently')
    # The proven file must arrive, otherwise it is a silent no-op
    _bw = _proven_alternating()
    probe('proven_alternating read', sorted(k for k in _bw if not k.startswith('_')),
          ['fleet_n5_2_3_5_7_7', 'st118_k2'], 'file reaches the status function')
    _st = next((r for r in ks if r.get('stem') == 'st118_k2'), None)
    if _st:
        probe('proven overrules heuristic', alternating_status(_st, 22), True,
              'isometric to the closure of an alternating 3-braid')
    for txt, expected_ in [('x^4 - x^3 + x^2 - x + 1', True), ('x^4 + x^3', False),
                      ('x^4 + x^2 + 1', True), ('no polynomial', None)]:
        probe(f'CM-Parser {txt[:20]!r}', alexander_is_alternating(txt), expected_, 'hand case')

    print('\n' + ('ALL COUNTER-CASES PASSED' if errors == 0 else f'{errors} ERRORS'))
    return errors


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--calibrate-census', action='store_true')
    ap.add_argument('--balance', action='store_true')
    ap.add_argument('--calibrate-alternating', action='store_true')
    a = ap.parse_args()
    if a.calibrate_census:
        sys.exit(1 if calibrate_census() else 0)
    elif a.calibrate_alternating:
        sys.exit(1 if calibrate_alternating() else 0)
    elif a.balance:
        balance()
    else:
        ap.print_help()
