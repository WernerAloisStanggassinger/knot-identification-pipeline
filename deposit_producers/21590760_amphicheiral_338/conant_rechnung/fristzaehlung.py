#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""COUNT deadline deaths — by a feature of the RESULT, not by log lines (footgun 23).

Recipe as in `BELEG_zahlen_nachfassen_20260729.md` A1, only re-run here because the road
has kept running since and candidates that once died at the deadline may later
have been evaluated after all. That is exactly why this number must not be carried
forward: it can DROP while everything else grows.

  1) read all jsonl from breit_out/ and breit_c_out/
  2) rows with rej == 'frist' -> distinct keys (k, seed, nxy, nz, dc, crot)
  3) subtract those who have a regular row IN THE SAME jsonl (rej == None):
     the same candidate, computed through after all in another run.

🔴 FIRST DRAFT OF THIS FILE WAS WRONG (01.08.2026) and the error is instructive: step 3
subtracted against `pruefstand_stufe1.json` and arrived at 3702 instead of 3153. Stage 1 sees only the
HITS of the sieves, not every evaluated candidate — the subtraction was thus almost empty (62 instead of
611). It was noticed only because the number ROSE against the previous day's value, where it can at most fall:
evaluated candidates do not become unevaluated again. **The contradiction was the finding,
not the own number** (footgun 29). With the correct subtraction the receipt state is reproduced
bit-exactly, including the distribution per c and per k.

Werner Stanggassinger + Claude Code (Anthropic), 2026-08-01.
"""
import site_config
import json, os, glob, collections

SIEB = (site_config.ROOT + '/konversationen/kawauchi_20260725/sieb')
KEY = lambda r: (r['k'], r['seed'], r['nxy'], r['nz'], r['dc'], r.get('crot', 1))

dateien = sorted(glob.glob(os.path.join(SIEB, 'breit_out', '*.jsonl')) +
                 glob.glob(os.path.join(SIEB, 'breit_c_out', '*.jsonl')))
zeilen = frist_zeilen = 0
frist, ausgewertet = set(), set()
for p in dateien:
    for ln in open(p):
        ln = ln.strip()
        if not ln: continue
        zeilen += 1
        try: r = json.loads(ln)
        except Exception: continue
        try: k = KEY(r)
        except KeyError: continue
        rj = r.get('rej')
        if rj == 'frist':
            frist_zeilen += 1; frist.add(k)
        elif rj is None:
            ausgewertet.add(k)      # regular result row: this candidate was computed

nie = frist - ausgewertet

print(f'=== DEADLINE COUNT {len(dateien)} jsonl, {zeilen} rows ===')
print(f"  rej='frist' in rows          : {frist_zeilen}")
print(f'  distinct keys                : {len(frist)}')
print(f'  of these computed elsewhere  : {len(frist & ausgewertet)}')
print(f'  🔴 NEVER evaluated           : {len(nie)}      (letter state 296 named 3153)')
if len(nie) > 3153:
    print('  🔴 ROSE — that cannot be, evaluated ones do not become unevaluated.')
    print('     First suspect the counter, not the data (footgun 15/29).')

je_c = collections.Counter(k[5] for k in nie)
je_k = collections.Counter(k[0] for k in nie)
print(f'  per c: {dict(sorted(je_c.items()))}')
print(f'  per k: {dict(sorted(je_k.items()))}   (k=8 is order 16)')

s1 = json.load(open(os.path.join(SIEB, 'pruefstand_stufe1.json')))
c13_gerechnet = [r for r in s1 if r.get('k') == 8 and r.get('crot') == 13]
c13_entschieden = [r for r in c13_gerechnet if r.get('ok_stufe1') is not None]
print(f'\n  c=13: dead at the deadline {je_c.get(13, 0)} · reached the evaluation '
      f'{len(c13_gerechnet)} · of these decided {len(c13_entschieden)}')
print(f'        (letter state 296: 604 dead · 6 evaluated · 4 decided)')
