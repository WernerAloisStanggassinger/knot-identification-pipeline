#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""wd_homfly.py — HOMFLY polynomial of the untwisted Whitehead doubles of a knot, and a mirror-blind comparison of two knots by it.

The Whitehead double is built by Regina (Link.whiteheadDouble, positive or negative clasp, untwisted), its HOMFLY polynomial by
SageMath (libhomfly). For a pair A, B the doubles of A are compared with the doubles of B AND of mirror(B), clasp by clasp, so a
mirror image does not count as a difference. Different values for a clasp prove that A and B are different knots up to mirror image;
equal values prove nothing.

Names: a knot of Burton's census ('<c><a|n><h|s|t>_<i>', DT code read from BURTON_CENSUS_DIR, default ./burton_census/raw) or any name
SnapPy's link tables accept (e.g. K14n22185, 8_19).

  python -I wd_homfly.py counter                       counter-cases (below); all are computed, exit code 1 if any fails
  python -I wd_homfly.py pairs.tsv out.jsonl WORKERS   pairs.tsv: two names per line (tab or space). One JSON line per
                                                       (knot, mirror, clasp), written the moment it is computed; a restart skips
                                                       what out.jsonl already holds. Then one verdict line per pair on stdout.
Counter-cases (no census files needed):
  K14n22185 / K14n22589 — Stoimenow's 14-crossing pair (14_41721, 14_42125) (= KnotScape 14n 22185 / 22589): equal volume,
                          Alexander and Jones polynomial; his tables list it as separated by the Whitehead-double skein polynomial:
                          the doubles must DIFFER for both clasps;
  15ah_27440 / 15ah_27441 — Conway mutants (explicit mutation; = KnotScape 15a_13571 / 15a_13629, listed by Stoimenow): the doubles
                          must be EQUAL for both clasps — a regression check (equal on 3 November 2026), not a consequence of a
                          cited theorem. Their two census DT codes are written into this file, so no census file is read.
Run with SageMath's Python; needs regina, snappy. Running times: see README.md.
"""
import bz2, csv, json, os, sys, time
from multiprocessing import Pool

DT = lambda s: [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]
KIND = {'h': 'hyp', 's': 'satellite', 't': 'torus'}
CENSUS_PREFIX = ('a', 'n')


def is_census(name):
    pre = name.split('_')[0]
    return '_' in name and len(pre) >= 3 and pre[:-2].isdigit() and pre[-2] in CENSUS_PREFIX and pre[-1] in KIND


def dt_codes(names):
    """name -> DT code (list of even ints) for census names and SnapPy table names."""
    import snappy
    out, want = {}, {}
    for n in names:
        if is_census(n):
            pre = n.split('_')[0]
            f = os.path.join(os.environ.get('BURTON_CENSUS_DIR', 'burton_census/raw'), f'{pre[:-2]}{pre[-2]}-{KIND[pre[-1]]}.csv.bz2')
            want.setdefault(f, set()).add(n)
        else:
            d = snappy.Link(n).DT_code()
            assert len(d) == 1, (n, 'not a knot')
            out[n] = list(d[0])
    for f, ns in want.items():
        got = 0
        for r in csv.DictReader(bz2.open(f, 'rt')):
            if r['name'] in ns:
                out[r['name']] = DT(r['dt_code']); got += 1
                if got == len(ns): break
    missing = [n for n in names if n not in out]
    assert not missing, ('names not found', missing[:5])
    return out


def wd_job(args):
    import regina, snappy
    import sage.all  # noqa: F401  (loads SageMath for sage_link().homfly_polynomial())
    name, code, mirror, positive = args
    t = time.time()
    L = regina.Link.fromDT(code)
    if mirror: L.reflect()
    W = L.whiteheadDouble(positive)
    S = snappy.Link([list(x) for x in W.pdData()])
    assert len(S.link_components) == 1
    h = str(S.sage_link().homfly_polynomial())
    return dict(name=name, mirror=mirror, positive=positive, homfly=h, sec=round(time.time() - t, 1))


def compute(names, out_path, workers, codes=None):
    C = dict(codes or {}); C.update(dt_codes([n for n in names if n not in C]))
    jobs = [(n, C[n], m, p) for n in names for m in (False, True) for p in (True, False)]
    have = {}
    if out_path and os.path.exists(out_path):
        for ln in open(out_path):
            r = json.loads(ln); have[(r['name'], r['mirror'], r['positive'])] = r['homfly']
    todo = [j for j in jobs if (j[0], j[2], j[3]) not in have]
    todo.sort(key=lambda j: len(j[1]))
    print(f'{len(jobs)} doubles, {len(have)} on disk, {len(todo)} to compute', flush=True)
    f = open(out_path, 'a') if out_path else None
    with Pool(workers, maxtasksperchild=1) as pool:
        for r in pool.imap_unordered(wd_job, todo):
            have[(r['name'], r['mirror'], r['positive'])] = r['homfly']
            if f: f.write(json.dumps(r) + '\n'); f.flush()
            print(f'  {r["name"]} mirror={r["mirror"]} positive={r["positive"]} {r["sec"]} s', flush=True)
    if f: f.close()
    return have


def equal(H, a, b, positive):
    """Mirror-blind: WD(A) == WD(B) or WD(A) == WD(mirror B), for one clasp."""
    return H[(a, False, positive)] in (H[(b, False, positive)], H[(b, True, positive)])


def counter(workers=4):
    cases = [('K14n22185', 'K14n22589', False), ('15ah_27440', '15ah_27441', True)]
    fixed = {'15ah_27440': DT('nejgbaclmfdoihk'), '15ah_27441': DT('negfajclmdboihk')}   # Burton's census, 15a-hyp.csv.bz2
    H = compute(sorted({n for a, b, _ in cases for n in (a, b)}), None, workers, fixed)
    ok = True
    for a, b, expect in cases:
        for p in (True, False):
            got = equal(H, a, b, p)
            print(f'counter-case {a} / {b} clasp={"positive" if p else "negative"}: equal={got} (expect {expect})', flush=True)
            ok &= got == expect
    print('COUNTER-CASES OK' if ok else 'COUNTER-CASE FAILED', flush=True)
    return ok


def main():
    if sys.argv[1:2] == ['counter']:
        sys.exit(0 if counter(int(sys.argv[2]) if len(sys.argv) > 2 else 4) else 1)
    pairs_path, out_path, workers = sys.argv[1], sys.argv[2], int(sys.argv[3])
    pairs = [ln.split()[:2] for ln in open(pairs_path) if ln.strip() and not ln.startswith('#')]
    H = compute(sorted({n for p in pairs for n in p}), out_path, workers)
    for a, b in pairs:
        eq = {p: equal(H, a, b, p) for p in (True, False)}
        verdict = 'different knots (up to mirror image)' if not all(eq.values()) else 'not separated'
        print(f'{a}\t{b}\tpositive_equal={eq[True]}\tnegative_equal={eq[False]}\t{verdict}', flush=True)


if __name__ == '__main__':
    main()
