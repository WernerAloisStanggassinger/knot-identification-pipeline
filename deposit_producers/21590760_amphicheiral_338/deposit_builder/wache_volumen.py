#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Frankenstein guard for the v3 draft: does the SUPPLIED diagram carry
the volume OF ITS record? (Mostow: the volume belongs to the knot.)

Full guard on all new ones (the assignment is new there) + sample on the carried-over ones.
Only when it holds does the draft become counterexamples_data.json.
Env: gnottheorie_env (spherogram/snappy).

--selbsttest swaps TWO diagrams. It passes only if
  (a) exactly these two break AND (b) all the others remain ok.
A checker failure breaks everything too — that would otherwise pass as „passed".
Exactly that happened on the first attempt (32/32 break, self-test reported green).
"""
import site_config
import json, os, re, sys
from multiprocessing import Pool
import snappy          # MUST be loaded before the pool fork: otherwise spherogram loads SnapPy
import spherogram     # only afterwards in the worker -- and fails there with "SnapPy doesn.t seem to be available"

OUT = (site_config.ROOT + '/PRIORITAET_DEPOSIT_20260730/')
TOL = 1e-4

def pd_liste(s):
    return [[int(x) for x in m] for m in re.findall(r'X\[\s*(-?\d+),\s*(-?\d+),\s*(-?\d+),\s*(-?\d+)\s*\]', s)]

def vol(arg):
    i, pd, soll = arg
    try:
        M = spherogram.Link(pd_liste(pd)).exterior()
        M.simplify()
        v = float(M.volume())
        return (i, v, abs(v - soll) < TOL, '')
    except Exception as ex:
        return (i, None, False, f'{type(ex).__name__}: {str(ex)[:60]}')

if __name__ == '__main__':
    selbsttest = '--selbsttest' in sys.argv
    d = json.load(open(OUT + 'counterexamples_data_ENTWURF.json'))
    ks = d['knoten']; soll = {e['id']: e['volumen'] for e in ks}
    # 03.08.2026: „new" NO LONGER means a=4 — the discriminator is the origin:
    # records carried over from v2 carry `kreuzungs_obergrenze_v2`, new ones do not.
    neue = [e for e in ks if 'kreuzungs_obergrenze_v2' not in e]
    uebern = [e for e in ks if 'kreuzungs_obergrenze_v2' in e]
    stich = neue + uebern[::10]
    ids_neu = {e['id'] for e in neue}
    auf = [(e['id'], e['pd_code'], e['volumen']) for e in stich]
    vertauscht = set()
    if selbsttest:
        a, b = auf[0], auf[1]
        auf[0], auf[1] = (a[0], b[1], a[2]), (b[0], a[1], b[2])
        vertauscht = {a[0], b[0]}
        print(f"SELF-TEST: diagrams of {a[0]} and {b[0]} swapped — exactly these two must break.\n")
    print(f"volume guard: {len(neue)} new + {len(stich)-len(neue)} sample = {len(stich)} records\n", flush=True)
    with Pool(10) as p:
        res = p.map(vol, auf, chunksize=1)
    bruch = [r for r in res if not r[2]]
    for i, v, ok, err in res:
        if not ok or i in ids_neu:
            print(f"  {'ok' if ok else 'BREAK':<6} {i:<22} PD volume {('—' if v is None else round(v,6)):>12}"
                  f"  expected {round(soll[i],6):>12}  {err}", flush=True)
    print(f"\n  checked {len(res)}  ok {len(res)-len(bruch)}  broken {len(bruch)}")
    if selbsttest:
        gebrochen = {r[0] for r in bruch}
        gut = gebrochen == vertauscht
        print(("SELF-TEST PASSED: exactly the two swapped ones break, all others hold."
               if gut else
               f"SELF-TEST FAILED: broken {sorted(gebrochen)}, expected {sorted(vertauscht)} "
               "— on total failure the checker measures itself."))
        sys.exit(0 if gut else 3)
    if bruch:
        print("🔴 ABORT: the draft stays a draft."); sys.exit(2)
    os.replace(OUT + 'counterexamples_data_ENTWURF.json', OUT + 'counterexamples_data.json')
    print("✅ guard holds → counterexamples_data.json written.")
