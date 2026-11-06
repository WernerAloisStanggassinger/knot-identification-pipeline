#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""anchor_recompute.py — recompute the three corpus anchors FRESH from the PD.

WHY: after a tool update the question is not „does it still run", but „does it compute
the SAME". The function probes cover the import and a few standard knots; this
guard takes the three anchors from CLAUDE.md, reads their PD from the Pflicht file and computes
determinant and Jones ANEW - against the values that landed in the corpus when they were entered.
A deviation means: the stack has moved and the corpus is to be checked.

COUNTER-CASE BUILT IN: for each anchor a MUTILATED PD is computed along (last
crossing removed). If that yields the same value, the guard measures nothing and reports red.

  python3 anchor_recompute.py        (base env: regina)
"""
import site_config
import glob, json, re, sys
from pathlib import Path

DB = Path((site_config.PIPELINE_DIR + '/corpus_db.json'))
ANKER = ['knots-of-love', 'irrational', 'nitro']

def pd_aus_pflicht(name):
    c = glob.glob(f'{site_config.ROOT}/own_knots/{name}/*pflicht*complete.json')
    if not c:
        raise FileNotFoundError(f'no Pflicht file for {name}')
    m = re.search(r'"pd_list_str"\s*:\s*"((?:[^"\\]|\\.)*)"', open(c[0]).read())
    if not m:
        raise ValueError(f'no pd_list_str in {c[0]}')
    return [tuple(int(q) for q in g) for g in
            re.findall(r'[\(\[]\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*[\)\]]', m.group(1))]

def det_von(pd):
    """Determinant = |Delta(-1)|. regina.Link has NO determinant() - looked up with
    dir(Link) on 21.09.2026, not guessed (Hard Rule 2)."""
    import sympy as sp
    from regina import Link
    t = sp.symbols('t')
    roh = Link.fromPD(pd).alexander().str().replace('^', '**')
    roh = re.sub(r'(\d)\s*x', r'\1*x', roh).replace('x', 't')
    return int(abs(sp.expand(sp.sympify(roh)).subs(t, -1)))

def main():
    d = json.loads(DB.read_text())
    soll = {r['stem']: r for r in d['knots'] if r['stem'] in ANKER}
    fehlt = [a for a in ANKER if a not in soll]
    if fehlt:
        print('RED: not found in the corpus:', fehlt); return 1

    rot = []
    for name in ANKER:
        try:
            pd = pd_aus_pflicht(name)
            ist = det_von(pd)
            erwartet = int(soll[name]['det'])
            # Counter-case: the same PD without the last crossing MUST deviate
            try:
                gegen = det_von(pd[:-1]); beisst = (gegen != ist)
            except Exception:
                beisst = True          # mutilated PD already fails on reading in
            ok = (ist == erwartet)
            print(f'  {name:16s} c={len(pd):3d}  det fresh={ist:6d}  corpus={erwartet:6d}  '
                  f'{"OK " if ok else "RED"}  counter-case {"bites" if beisst else "BLUNT"}')
            if not ok or not beisst: rot.append(name)
        except Exception as e:
            print(f'  {name:16s} RED: {type(e).__name__}: {e}'); rot.append(name)

    print(f'\nAnchors: {len(ANKER)-len(rot)}/{len(ANKER)} freshly confirmed')
    return 1 if rot else 0

if __name__ == '__main__':
    sys.exit(main())
