#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Knot-Identification Pipeline — Werner Alois Stanggassinger & Claude Code (Anthropic).
"""turaev_bound.py — crossing-number lower bound c(K) >= span V(K) + g_T(K), with g_T(K) >= w_HF(K) - 1. No table.

THEOREMS (receipts, versions read 2026-10-28; PDFs + text in quellen_scans/):
  [CK]  Champanerkar–Kofman, "A survey on the Turaev genus of knots", arXiv:1406.1945v1, eq. (1):
        span V_L(t) <= c(L) - g_T(L)        (Turaev; L non-split, here a knot)
        property (c): L is alternating iff g_T(L) = 0.
  [Lo]  Lowrance, "On knot Floer width and Turaev genus", arXiv:0709.0720v1 = AGT 8 (2008) 1141-1162,
        Def. 2.3: w_HF(K) = Delta_K - delta_K + 1, Delta/delta = max/min of A - M over generators of HFK-hat;
        Thm 1.1:  w_HF(K) <= g_T(K) + 1.
  ==> c(K) >= span V(K) + max(0, w_HF(K) - 1).   w_HF >= 2 also PROVES K non-alternating ([Lo] + [CK](c)).
  Beyond every enumeration iff c >= 25, or c >= 21 and K non-alternating (CLAUDE.md, Disc #814).

FRESH MEASUREMENT per knot: Jones via regina (Link.fromPD(...).jones(), Laurent in sqrt(t) -> span/2), HFK via
knot_floer_homology.pd_to_hfk (Szabo's program). The stored record Jones is compared as a second witness.

Usage (sage_stable env):
  turaev_bound.py --calibrate-census N [N ...]      calibration on Burton census knots of crossing number N (bound <= N, must never exceed)
  turaev_bound.py --records [--limit K]   all distinct representatives of class unresolved_c_le_20
Werner Alois Stanggassinger + Claude Code, 2026-10-28.
"""
import argparse, ast, bz2, csv, glob, json, re, signal, sys
import site_config
import regina, knot_floer_homology as kfh

DB = '@PIPELINE_DIR@/corpus_db.json'.replace('@PIPELINE_DIR@', site_config.PIPELINE_DIR)
CENSUS = '@ROOT@/external_tools/burton_census/raw'.replace('@ROOT@', site_config.ROOT)


class Timeout(Exception):
    pass


def _alarm(*_):
    raise Timeout()


def span_jones(pd):
    """span V in powers of t (regina returns V in sqrt(t))."""
    j = regina.Link.fromPD(pd).jones()
    return (j.maxExp() - j.minExp()) / 2, str(j)


def hfk_width(pd):
    r = kfh.pd_to_hfk(pd)
    d = [a - m for (a, m), rk in r['ranks'].items() if rk]
    return max(d) - min(d) + 1, r['total_rank']


def bound(pd, secs):
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(secs)
    try:
        s, jstr = span_jones(pd)
        w, tot = hfk_width(pd)
    finally:
        signal.alarm(0)
    return {'span_V': s, 'w_HF': w, 'hfk_total_rank': tot, 'nonalt_proven': w >= 2,
            'c_lower': int(s + max(0, w - 1)), 'jones_regina': jstr}


def record_span(jones_str):
    """span of a stored record Jones string, written in x = sqrt(t) -> span in t."""
    # two spellings occur in the corpus: 'x^-2' and '197/x^2', 'x^(-14)', '2*x^22' (lj_4_7_11 misread before 28.10. 16:3x)
    s = jones_str.replace('*', '')
    s = re.sub(r'x\^\((-?\d+)\)', r'x^\1', s)
    s = re.sub(r'/x\^(\d+)', r'x^-\1', s)
    s = re.sub(r'/x(?![\^\d])', 'x^-1', s)
    toks = re.split(r'\s(?=[+-]\s)', s.strip())
    ex = []
    for t in toks:
        m = re.search(r'x\^(-?\d+)', t)
        ex.append(int(m.group(1)) if m else (1 if 'x' in t else 0))
    return (max(ex) - min(ex)) / 2


def dt_to_pd(dt_letters):
    import spherogram as S
    code = [(-1 if ch.isupper() else 1) * 2 * (ord(ch.lower()) - 96) for ch in dt_letters]
    # spherogram labels edges 0..2n-1, regina.fromPD wants 1..2n
    return [tuple(e + 1 for e in c) for c in S.Link('DT: [(%s)]' % ','.join(map(str, code))).PD_code()]


def calibrate_census(ns, secs):
    bad = 0
    for n in ns:
        files = sorted(glob.glob(f'{CENSUS}/{n}a-*.csv.bz2') + glob.glob(f'{CENSUS}/{n}n-*.csv.bz2'))
        cnt = {'knots': 0, 'nonalt_census': 0, 'nonalt_proven': 0, 'tight': 0, 'violations': 0, 'alt_but_thick': 0}
        for f in files:
            alt_file = f'/{n}a-' in f
            for row in csv.DictReader(bz2.open(f, 'rt')):
                b = bound(dt_to_pd(row['dt_code']), secs)
                cnt['knots'] += 1
                cnt['nonalt_census'] += not alt_file
                cnt['nonalt_proven'] += b['nonalt_proven']
                cnt['tight'] += b['c_lower'] == n
                if b['c_lower'] > n:
                    cnt['violations'] += 1; bad += 1
                    print('VIOLATION', row['name'], row['dt_code'], b, flush=True)
                if alt_file and b['nonalt_proven']:
                    cnt['alt_but_thick'] += 1; bad += 1
                    print('ALTERNATING BUT HFK-THICK', row['name'], b, flush=True)
        print(f'n={n}', json.dumps(cnt), flush=True)
    print('CALIBRATION', 'RED' if bad else 'GREEN', f'({bad} violations)')
    return 1 if bad else 0


def records(limit, secs):
    db = json.load(open(DB))['knots']
    recs = [x for x in (db.values() if isinstance(db, dict) else db)
            if x.get('is_distinct_representative') and x.get('realizes_class') == 'unresolved_c_le_20']
    recs.sort(key=lambda x: x['c'])
    out = []
    for x in recs[:limit]:
        pd = ast.literal_eval(json.load(open(x['pflicht_path']))['pd_list_str'])
        row = {'id': x['distinct_knot_id'], 'c_ub': x['c'], 'c_lb_stored': x['c_lower_bound_proven'], 'n_pd': len(pd)}
        try:
            row.update(bound(pd, secs))
            row.pop('jones_regina')
            rs = record_span(x['jones']) if x.get('jones') else None
            row['span_V_record'] = rs
            row['span_witness_agree'] = rs == row['span_V']
            row['beyond'] = row['nonalt_proven'] and max(row['c_lower'], x['c_lower_bound_proven']) >= 21
        except Timeout:
            row['error'] = f'timeout {secs}s'
        except Exception as e:
            row['error'] = repr(e)[:200]
        print(json.dumps(row), flush=True)
        out.append(row)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--calibrate-census', type=int, nargs='*')
    ap.add_argument('--records', action='store_true')
    ap.add_argument('--limit', type=int, default=10**9)
    ap.add_argument('--secs', type=int, default=1800)
    a = ap.parse_args()
    if a.calibrate_census:
        sys.exit(calibrate_census(a.calibrate_census, a.secs))
    if a.records:
        records(a.limit, a.secs)
