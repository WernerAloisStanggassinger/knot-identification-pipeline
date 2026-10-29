#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""yamada_c_bound.py — c >= 2g + n - 1, the sharp version of the genus route.

WHAT FOR (2026-09-09): label_distinct_knots.c_untere_schranke_mit_via uses
    c >= 2g + s - 1   with s >= 2
and plugs in s = 2. But the comment there itself already cites Yamada 1987:
the MINIMAL number of Seifert circles over all diagrams IS the braid index n.
So s >= n holds, not just s >= 2 — and n is typically 4..9 for us, not 2.

THE INEQUALITY, step by step (none of it guessed):
  Let D be a diagram with minimal crossing number, c(D) = c(K), and s = s(D) the
  number of its Seifert circles. The Seifert algorithm yields a surface F from it
  with  2 g(F) = c - s + 1.
  The Seifert genus is the minimum over all surfaces:  g(K) <= g(F).
  Yamada 1987 (Invent. Math. 89, 347-356):  min_D s(D) = n(K), the braid index,
  hence  s(D) >= n(K)  for EVERY diagram, in particular for this one.
  ==>  2 g(K) <= c - s + 1 <= c - n + 1
  ==>  c(K) >= 2 g(K) + n(K) - 1.
  With n >= n_MFW (Morton-Franks-Williams lower bound from HOMFLY) the
  inequality stays valid: c >= 2g + n_MFW - 1.

TWO HARD PREREQUISITES, both checked instead of believed:
  * g must be the TRUE Seifert genus, not a bound-as-value.
    Only what knot_floer_homology.seifert_genus certifies is taken over
    (Ozsvath-Szabo; the same test as label_distinct_knots.hfk_exaktes_genus).
  * n_MFW must be a LOWER bound. braid_index_lower_bound is exactly that
    (enrich_braid_index_mfw.py); braid_index/upper_bound is NOT and is used
    nowhere here.

TRANSFER WITHIN A CLUSTER: members with the same distinct_knot_id are
the same knot, hence the same braid index and the same genus. Ingredients may
come from a sibling — but ONLY if the siblings do not contradict each other.
Contradiction => the cluster is skipped and reported loudly, never smoothed over.

--calibrate-census runs four counter-cases. Without a green calibration every number here is worthless.
"""
import json, os, re, sys, argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from span_lower_bound import span_in_t, laden, exponents

ENUMERATION_LIMIT = 20          # up to c=20 enumerated completely


def genus_lower_from_alexander(a):
    """g >= span Δ(t) / 2 (Seifert: deg Δ <= 2g). LOWER bound, no HFK needed.

    Why that suffices: c >= 2g + n - 1 holds with the TRUE genus; whoever estimates g from below
    obtains a valid, only weaker bound. The Alexander polynomial is present for
    almost every record — so even the knots for which HFK dies at the girth limit get a bound
    (09.09.2026: 3 more EXACT crossing numbers up to c=146,
    without a single new computation).
    CALIBRATED 09.09.2026 (probe 5 in --calibrate-census, runs along on every call): on the
    646 records that carry genus AND Alexander, span/2 <= g holds, 0 overshoots.
    The 775 named here earlier was a different denominator (without the genus condition) —
    corrected, because a number in a docstring is read exactly like a measured one.
    """
    if not isinstance(a, str) or not a.strip():
        return None
    raw = a.replace('**', '^').replace(' ', '')
    if raw in ('1', '-1', '1.0'):
        return None                      # trivial Alexander => no statement
    # WHICH variable is there? Do not guess and do not query "empty":
    # exponenten() returns the list [0] for EVERY string with a digit and without x
    # (constant branch). That looks like a result and is not one — the same silent
    # trap as Disc #766 (MFW blind to v,z). So first establish the convention.
    hat_x = re.search(r'(?<![A-Za-z])x(?![A-Za-z])', raw) is not None
    hat_t = re.search(r'(?<![A-Za-z])t(?![A-Za-z])', raw) is not None
    if hat_x and hat_t:
        return None                      # two variables => not responsible
    if hat_x:
        ex = exponents(a)               # Regina convention (829/829 in the corpus)
    elif hat_t:
        ex = [int(m.group(1)) for m in re.finditer(r'(?<![A-Za-z0-9])t\^(-?\d+)', raw)]
        if re.search(r'(?<![A-Za-z0-9])t(?!\^)(?![A-Za-z0-9])', raw):
            ex.append(1)
        if re.search(r'(?:^|[+-])\d+(?:$|[+-])', raw):
            ex.append(0)
    else:
        return None                      # no variable => no polynomial
    if not ex:
        return None
    sp = max(ex) - min(ex)
    if sp <= 0:
        return None
    return sp // 2


def mfw_from_homfly(h):
    """MFW lower bound b(K) >= (v-span)/2 + 1 from HOMFLY P(v,z). None in case of doubt.

    Knows all three conventions occurring in the corpus: Regina 'x,y', Sage 'L,M'
    and 'v,z'. The pipeline parser (pflicht_invarianten_extra.mfw_braid_index_bound)
    does NOT know the third and silently returns lower_bound=1 there — a number made of
    nonsense is exactly the trap against which None stands here (finding 09.09.2026).
    """
    if not isinstance(h, str) or not h.strip():
        return None
    s = h.strip()
    if s in ('1', '1.0'):        # unknot: HOMFLY == 1, no braid argument
        return None
    for var in ('L', 'v', 'x'):
        if re.search(rf'(?<![A-Za-z0-9]){var}(?![A-Za-z0-9])', s):
            break
    else:
        return None
    t = s.replace('**', '^').replace('*', ' ')
    t = re.sub(r'\s+-\s+', ' + -', t)
    if not t.startswith(('+', '-')):
        t = '+ ' + t
    pat = re.compile(rf'(?<![A-Za-z0-9]){var}(?:\s*\^\s*(-?\d+))?')
    ex = []
    for part in [p.strip() for p in re.split(r'\s+\+\s+', t) if p.strip()]:
        m = pat.search(part)
        ex.append(int(m.group(1)) if (m and m.group(1) is not None) else (1 if m else 0))
    if not ex:
        return None
    span_ = max(ex) - min(ex)
    if span_ % 2:               # odd span => convention not understood
        return None
    return span_ // 2 + 1


def ingredients_of(r):
    """MFW lower bound and HFK-exact genus from the mandatory record of this member."""
    out = {'stem': r.get('stem'), 'did': r.get('distinct_knot_id'),
           'c': r.get('c'), 'genus': r.get('genus'), 'alt': r.get('alternating'),
           'rc': r.get('realizes_class'), 'rep': bool(r.get('is_distinct_representative')),
           'clb': r.get('c_lower_bound_proven'), 'jones': r.get('jones'),
           'bi_l': None, 'hfk_g': None, 'g_alex': None, 'bi_konflikt': None, 'g_konflikt': None}
    # MFW directly from the corpus HOMFLY — covers the records whose mandatory record
    # never received the field. Own parser: the pipeline parser is blind to the
    # v,z convention and silently returns 1 there instead of None (finding 09.09.).
    m = mfw_from_homfly(r.get('homfly'))
    if m is not None:
        out['bi_l'] = m
    out['g_alex'] = genus_lower_from_alexander(r.get('alex'))
    pf = r.get('pflicht_path')
    if pf and os.path.exists(pf):
        try:
            d = json.load(open(pf))
            diag_block = d.get('pflichtliste', {}).get('DIAGRAMM', {})
            bl = diag_block.get('braid_index_lower_bound')
            # A field counts only if its OWN note covers it. su_4_4_tm1 carried
            # field=11 with note "b(K) >= 9" — upward, i.e. overclaim direction.
            note = diag_block.get('braid_index_mfw_note')
            mn = re.search(r'b\(K\)\s*>=\s*(\d+)', str(note)) if note else None
            if isinstance(bl, int) and mn and int(mn.group(1)) != bl:
                out['bi_konflikt'] = (bl, int(mn.group(1)))
                bl = None
            if isinstance(bl, int):
                out['bi_l'] = bl if out['bi_l'] is None else max(out['bi_l'], bl)
            sg = (d.get('pflichtliste', {}).get('ALGEBRA', {})
                  .get('knot_floer_homology') or {}).get('seifert_genus')
            # HFK IS the Seifert genus (Ozsvath-Szabo). The corpus field is a
            # cross-check, not a replacement: if it CONTRADICTS, it is discarded and reported
            # (then one of the two was a bound-as-value); if it is MISSING, there is
            # nothing to contradict and HFK counts. Measured 09.09.: 591 records
            # agree, 109 have no corpus genus, 0 contradict.
            if isinstance(sg, int):
                g_corp = r.get('genus')
                if isinstance(g_corp, int) and g_corp != sg:
                    out['g_konflikt'] = (g_corp, sg)
                else:
                    out['hfk_g'] = sg
        except Exception:
            pass
    return out


def collect(ks, workers=12):
    with ProcessPoolExecutor(max_workers=workers) as ex:
        return list(ex.map(ingredients_of, ks, chunksize=20))


def je_cluster(rows):
    """Bundle ingredients per distinct_knot_id. Returns (table, contradictions)."""
    tab, conflicts = {}, []
    for r in rows:
        d = r['did']
        if not d:
            continue
        t = tab.setdefault(d, {'bi_l': None, 'hfk_g': None, 'g_alex': None, 'members': 0})
        t['members'] += 1
        if isinstance(r['bi_l'], int):
            # MFW lower bounds may scatter (different diagrams) -> the
            # largest is the sharpest and remains a valid lower bound.
            t['bi_l'] = r['bi_l'] if t['bi_l'] is None else max(t['bi_l'], r['bi_l'])
        if isinstance(r.get('g_alex'), int):
            t['g_alex'] = r['g_alex'] if t.get('g_alex') is None else max(t['g_alex'], r['g_alex'])
        if isinstance(r['hfk_g'], int):
            if t['hfk_g'] is not None and t['hfk_g'] != r['hfk_g']:
                conflicts.append((d, t['hfk_g'], r['hfk_g'], r['stem']))
            t['hfk_g'] = r['hfk_g']
    return tab, conflicts


def bound(g, n):
    """c >= 2g + n - 1, or None if an ingredient is missing."""
    if not isinstance(g, int) or not isinstance(n, int) or g < 0 or n < 2:
        return None
    return 2 * g + n - 1


# ------------------------------------------------------------------ Calibration
def calibrate_census():
    print('CALIBRATION yamada_c_schranke — four counter-cases\n')
    errors = 0
    ks = laden()
    rows = collect(ks)
    tab, conflicts = je_cluster(rows)

    print('1) THE BOUND MUST HOLD ON PROVEN CROSSING NUMBERS')
    print('   (reduced alternating: span V == c, so c is known EXACTLY)')
    checked = violated = 0
    samples = []
    for r in rows:
        if not r['alt'] or not isinstance(r['c'], int) or not r['jones']:
            continue
        if span_in_t(r['jones']) != r['c']:
            continue                      # c not proven exact here
        t = tab.get(r['did']) or {}
        s = bound(t.get('hfk_g'), t.get('bi_l'))
        if s is None:
            continue
        checked += 1
        if s > r['c']:
            violated += 1
            if len(samples) < 6:
                samples.append((r['stem'], r['c'], t['hfk_g'], t['bi_l'], s))
    ok = violated == 0 and checked > 10
    print(f'   {"OK " if ok else "ERROR"} {checked} cases checked, {violated} violations of 2g+n-1 <= c')
    for st, c, g, n, s in samples:
        print(f'        {str(st)[:40]:<42} c={c} g={g} n={n} -> bound {s}  ← LARGER than c')
    if checked <= 10:
        print(f'        ⚠ only {checked} calibration cases — too thin to establish the convention')
    errors += (not ok)

    print('\n2) CANARY trainstation — c=18 exact (Tait), g=6, n=7: the bound MUST hit 18')
    tr = [r for r in rows if r['stem'] and 'trainstation' in str(r['stem']).lower()]
    hit = False
    for r in tr:
        t = tab.get(r['did']) or {}
        s = bound(t.get('hfk_g'), t.get('bi_l'))
        print(f'        {str(r["stem"])[:40]:<42} c={r["c"]} g={t.get("hfk_g")} n={t.get("bi_l")} -> {s}')
        if s == 18 and r['c'] == 18:
            hit = True
    print(f'   {"OK " if hit else "ERROR"} Canary {"hits c=18 exactly" if hit else "does NOT hit — ingredients or formula wrong"}')
    errors += (not hit)

    print('\n3) COUNTER-CASE — a missing or nonsensical ingredient must deliver NO number')
    for g, n, label in ((None, 5, 'g missing'), (10, None, 'n missing'), (10, 1, 'n=1 (no knot braid)'),
                        ('10', 5, 'g as string'), (-1, 5, 'g negative')):
        s = bound(g, n)
        gut = s is None
        print(f'   {"OK " if gut else "ERROR"} {label:<24} -> {s}')
        errors += (not gut)

    print('\n4) CLUSTER CONTRADICTIONS — same knot, different genus?')
    print(f'   {"OK " if not conflicts else "ERROR"} {len(conflicts)} contradictions')
    for d, a, b, st in conflicts[:6]:
        print(f'        {str(d)[:34]:<36} g={a} vs g={b}  ({st})')
    errors += bool(conflicts)

    print('\n5) ALEXANDER ROUTE — span(Delta)/2 must never overtake the TRUE genus')
    ks = laden()
    checked = over = 0
    worst = []
    for r in ks:
        # The yardstick is the corpus genus: it agrees with HFK on 591
        # records and contradicts on none (measurement 09.09.) — exactly right as an UPPER
        # comparison quantity for a lower bound.
        g_true = r.get('genus')
        g_alex = genus_lower_from_alexander(r.get('alex'))
        if not isinstance(g_true, int) or g_alex is None:
            continue
        checked += 1
        if g_alex > g_true:
            over += 1
            if len(worst) < 6:
                worst.append((r.get('stem'), g_alex, g_true))
    ok = over == 0 and checked > 100
    print(f'   {"OK " if ok else "ERROR"} {checked} records with both, {over} overshoots')
    for st, a, e in worst:
        print(f'        {str(st)[:44]:<46} alex={a} > g={e}')
    errors += (not ok)

    print('\n6) COUNTER-CASE ALEXANDER — wrong variable, no polynomial, two variables')
    for a, expected_, label in (('', None, 'empty'), ('1', None, 'trivial'),
                           ('nonsense 5', None, 'no variable'),
                           ('x + t', None, 'two variables'),
                           ('2 x^8 - 4 x^6 + 9 - 4 x^-6', 7, 'Regina with coefficient'),
                           ('t^2 - t + 1', 1, 't convention')):
        got = genus_lower_from_alexander(a)
        gut = got == expected_
        print(f'   {"OK " if gut else "ERROR"} {label:<26} -> {got} (expected {expected_})')
        errors += (not gut)

    print('\n' + ('ALL COUNTER-CASES PASSED' if errors == 0 else f'{errors} ERRORS — numbers unusable'))
    return errors


# ------------------------------------------------------------------ Balance
def balance(only_novel=True):
    ks = laden()
    rows = collect(ks)
    tab, conflicts = je_cluster(rows)
    if conflicts:
        print(f'⚠ {len(conflicts)} cluster contradictions — these are SKIPPED:')
        for d, a, b, st in conflicts[:5]:
            print(f'    {d}: g={a} vs {b} ({st})')
        print()
    bad = {d for d, *_ in conflicts}

    reps = [r for r in rows if r['rep'] and (r['rc'] == 'novel_hyperbolic' or not only_novel)]
    better, exact_val, without = [], [], 0
    for r in reps:
        if r['did'] in bad:
            continue
        t = tab.get(r['did']) or {}
        # HFK is sharp; where it is missing, the Alexander lower bound carries (weaker,
        # but valid — and it needs no computation).
        g = t.get('hfk_g')
        source = 'HFK'
        if g is None:
            g, source = t.get('g_alex'), 'Alex'
        s = bound(g, t.get('bi_l'))
        if s is None:
            without += 1
            continue
        r['_y'] = s
        r['_g'], r['_n'], r['_q'] = g, t['bi_l'], source
        alt = r['clb'] if isinstance(r['clb'], int) else 0
        if s > alt:
            better.append(r)
        if isinstance(r['c'], int) and s == r['c']:
            exact_val.append(r)
        if isinstance(r['c'], int) and s > r['c']:
            print(f'🔴 CONTRADICTION {r["stem"]}: bound {s} > c upper bound {r["c"]}')

    print(f'Representatives considered: {len(reps)}   (without both ingredients: {without})')
    print(f'\n  Yamada bound SHARPER than the one stored so far       : {len(better)}')
    print(f'  Bound == c upper bound    ⟹ CROSSING NUMBER EXACT     : {len(exact_val)}')
    if exact_val:
        print('\n  EXACT crossing numbers (c = 2g + n - 1 = c_ub):')
        for r in sorted(exact_val, key=lambda x: -x['c']):
            print(f'      c={r["c"]:<4} g={r["_g"]:<3} n={r["_n"]:<3} [{r["_q"]:<4}] {str(r["did"])[:34]:<36} {str(r["stem"])[:26]}')
    if better:
        print('\n  Improvements (old -> new):')
        for r in sorted(better, key=lambda x: -(x['_y'] - (x['clb'] or 0)))[:25]:
            print(f'      {str(r["clb"]):>4} -> {r["_y"]:<4} (g={r["_g"]}, n={r["_n"]}, c_ub={r["c"]})  {str(r["did"])[:34]}')
    return better, exact_val


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--calibrate-census', action='store_true')
    ap.add_argument('--balance', action='store_true')
    ap.add_argument('--all-classes', action='store_true', help='not only novel_hyperbolic')
    a = ap.parse_args()
    if a.calibrate_census:
        sys.exit(1 if calibrate_census() else 0)
    elif a.balance:
        balance(only_novel=not a.all_classes)
    else:
        ap.print_help()
