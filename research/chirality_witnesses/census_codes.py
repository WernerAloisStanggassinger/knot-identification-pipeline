#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""census_codes.py — shared driver for the chirality witnesses of this folder.

DT codes are read from Burton's census files in BURTON_CENSUS_DIR (default ./burton_census/raw): a name '<c><a|n><h|s|t>_<i>' is looked
up in '<c><a|n>-<hyp|satellite|torus>.csv.bz2'. Every witness module calls main(test, counter, from_dt):
  python -I <witness>.py counter                 runs the witness's counter-cases (known knots, expected answers) and stops on a failure
  python -I <witness>.py names.txt out.jsonl     one JSON line per name, written on arrival; names already in out.jsonl are skipped
"""
import bz2, csv, json, os, sys, time

DT = lambda s: [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]
KIND = {'h': 'hyp', 's': 'satellite', 't': 'torus'}


def census_file(name):
    pre = name.split('_')[0]                      # e.g. 15nh
    c, alt, kind = pre[:-2], pre[-2], pre[-1]
    return os.path.join(os.environ.get('BURTON_CENSUS_DIR', 'burton_census/raw'), f'{c}{alt}-{KIND[kind]}.csv.bz2')


def codes(names):
    want = {}
    for n in names: want.setdefault(census_file(n), set()).add(n)
    out = {}
    for f, ns in want.items():
        for r in csv.DictReader(bz2.open(f, 'rt')):
            if r['name'] in ns: out[r['name']] = r['dt_code']
            if len(out) == len(names): break
    missing = [n for n in names if n not in out]
    assert not missing, ('names not found in the census files', missing[:5])
    return out


def main(test, counter, from_dt):
    """test(obj) -> dict with 'verdict'; counter = list of (name, make, check): make() -> obj, check(result) -> bool;
    from_dt(code) -> obj."""
    if len(sys.argv) > 1 and sys.argv[1] == 'counter':
        bad = 0
        for nm, make, check in counter:
            t = time.time(); r = test(make()); ok = bool(check(r))
            print(json.dumps(dict(name=nm, counter_case_ok=ok, sec=round(time.time() - t, 2), **{k: v for k, v in r.items() if k != 'name'}), default=str), flush=True)
            bad += not ok
        print('COUNTER-CASES', 'OK' if not bad else f'{bad} FAILED', flush=True)
        sys.exit(1 if bad else 0)
    names = [l.split()[0] for l in open(sys.argv[1]) if l.strip()]
    out = sys.argv[2]
    done = {json.loads(l)['name'] for l in open(out)} if os.path.exists(out) else set()
    cd = codes([n for n in names if n not in done])
    with open(out, 'a') as fo:
        for n in names:
            if n in done: continue
            t = time.time(); r = dict(name=n)
            try: r.update(test(from_dt(cd[n])))
            except Exception as e: r.update(verdict='ERROR', error=type(e).__name__ + ' ' + str(e)[:160])
            r['sec'] = round(time.time() - t, 2); fo.write(json.dumps(r, default=str) + '\n'); fo.flush()
            print(n, r['verdict'], flush=True)


def gauss_sequence(pd):
    """Crossing indices in the order the knot passes them (each twice), from a one-component PD code with consecutive edge labels
    (KnotTheory convention: slot 0 = incoming under strand; the over strand runs j -> l iff l = j + 1 mod 2n)."""
    n2 = 2 * len(pd); head = {}
    for i, X in enumerate(pd):
        head[X[0]] = i
        head[X[1] if (X[3] - X[1]) % n2 == 1 else X[3]] = i
    assert len(head) == n2, 'PD edge labels not consecutive along one component'
    return [head[e] for e in sorted(head)]


def is_reduced(pd):
    """True iff the diagram has no nugatory crossing: a crossing is nugatory iff its chord in the Gauss word is interlaced with
    no other chord (every other crossing occurs 0 or 2 times between its two occurrences)."""
    seq = gauss_sequence(pd); pos = {}
    for k, c in enumerate(seq): pos.setdefault(c, []).append(k)
    for c, (a, b) in pos.items():
        inside = seq[a + 1:b]
        if not any(inside.count(d) == 1 for d in set(inside)): return False
    return True


SQUARE_KNOT_7 = 'O1+U2+O3+U4+O2+U3+O4+U1+O5-U6-O7-U5-O6-U7-'   # 7-crossing alternating, NOT reduced, diagram of 3_1 # m3_1 (amphichiral)
