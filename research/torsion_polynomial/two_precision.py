#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""two_precision.py — stability check of check_census.py results: recompute T_K at 1000 bits on a re-randomized triangulation.

Takes N random knots with status 'agree' from a check_census.py output (seed 7), plus the 10 numerically tightest non-monic ones (the
non-fibered knots with the smallest |lc| or |lc − 1|, i.e. closest to a wrong monic or wrong zero reading), and compares degree, monic class, real class and the leading coefficient (double precision, |Δlc| < 1e-12).

  python -I two_precision.py <census.csv.bz2> <check_census_out.jsonl> N out.jsonl      (SageMath's Python)
"""
import bz2, csv, json, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dfj_core as D

DT = lambda s: [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]


def main():
    census, prod_path, n, out_path = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
    D.set_memlimit(int(os.environ.get('MEM_GB', '6')))
    import snappy
    prod = {}
    for ln in open(prod_path):
        try:
            r = json.loads(ln)
        except Exception:  # noqa: BLE001
            continue
        if r.get('status') == 'agree': prod[r['name']] = r
    names = sorted(prod)
    random.Random(7).shuffle(names)
    # non-monic only: for a monic T_K, d_lc1_log2 is -1e9 and would sort first (the authors' script had that slip, so it picked 10 of
    # the clearest cases instead of the tightest)
    tight = sorted((r for r in prod.values() if r['monic'] != 'zero'), key=lambda r: min(r['lc_abs'], 2 ** r['d_lc1_log2']))[:10]
    pick = set(names[:n]) | {r['name'] for r in tight}
    codes = {r['name']: r['dt_code'] for r in csv.DictReader(bz2.open(census, 'rt')) if r['name'] in pick}
    bad = 0
    with open(out_path, 'w') as out:
        for nm in sorted(pick):
            p = prod[nm]
            M = snappy.Link('DT:' + str([tuple(DT(codes[nm]))])).exterior()
            M.randomize()
            T, _, _ = D.torsion_record(M, 1000)
            a = D.analyse(T)
            dlc = abs(complex(a['lc_re'], a['lc_im']) - complex(p['lc_re'], p['lc_im']))
            same = (a['deg'] == p['deg'] and a['monic'] == p['monic'] and a['real'] == p['real'] and dlc < 1e-12)
            bad += (not same)
            out.write(json.dumps(dict(name=nm, deg500=p['deg'], deg1000=a['deg'], lc500=[p['lc_re'], p['lc_im']],
                                      lc1000=[a['lc_re'], a['lc_im']], dlc=dlc, same=same)) + '\n'); out.flush()
    print('checked', len(pick), 'differences', bad, flush=True)


if __name__ == '__main__':
    main()
