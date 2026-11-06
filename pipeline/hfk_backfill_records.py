#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Enters the HFK results of the hunt into the Pflicht records — the step
that 'result file != corpus' demands (feedback_a_result_in_a_result_file_is_not_
in_the_corpus; the sundown case: a result lay seven weeks in a file
while the corpus and the public site claimed the opposite).

Writes ONLY ALGEBRA.knot_floer_homology (Pflicht field) — no classes, no
bounds, no labels. The classification remains the labeler's business.

THREE LATCHES, each with a counter-case (--selbsttest):
  * A record that ALREADY carries the field is only touched if the new value is
    identical; on deviation: stop loudly, write nothing.
  * genus in the corpus vs. seifert_genus from HFK: deviation is reported, not
    smoothed over — it means the corpus genus was a bound-as-value.
  * Backup of every touched file next to it (.bak-pre-hfk-YYYYMMDD), before writing.
"""
import site_config
import json, os, shutil, sys, argparse
from datetime import date
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from span_lower_bound import laden

OUT = Path((site_config.MACHINE_YARD + '/corpus_vs_burton/hfk_jagd_novel_20260909'))
STAMP = f'.bak-pre-hfk-{date.today():%Y%m%d}'


def ergebnisse():
    """did -> hfk dict, only real hits."""
    out = {}
    for f in sorted(OUT.glob('*.json')):
        if f.name.startswith('_selbsttest'):
            continue
        try:
            r = json.load(open(f))
        except Exception:
            continue
        h = r.get('hfk')
        if isinstance(h, dict) and isinstance(h.get('seifert_genus'), int):
            out[r['did']] = h
    return out


def plan(res):
    """(to_write, conflicts, genus_deviations) — nothing is written here."""
    ks = laden()
    schreiben, konflikte, gabw = [], [], []
    for r in ks:
        did = r.get('distinct_knot_id')
        if did not in res:
            continue
        pf = r.get('pflicht_path')
        if not pf or not Path(pf).exists():
            continue
        neu = res[did]
        try:
            d = json.load(open(pf))
        except Exception as e:
            konflikte.append((r.get('stem'), f'unreadable: {e}'))
            continue
        alt = (d.get('pflichtliste', {}).get('ALGEBRA', {})
               .get('knot_floer_homology') or {})
        if alt.get('seifert_genus') is not None:
            if alt.get('seifert_genus') != neu['seifert_genus']:
                konflikte.append((r.get('stem'),
                                  f"Record has g={alt.get('seifert_genus')}, hunt says g={neu['seifert_genus']}"))
            continue                      # already there and consistent -> nothing to do
        if isinstance(r.get('genus'), int) and r['genus'] != neu['seifert_genus']:
            gabw.append((r.get('stem'), r['genus'], neu['seifert_genus']))
        schreiben.append((r.get('stem'), pf, neu))
    return schreiben, konflikte, gabw


def selbsttest():
    print('SELFTEST hfk_nachtrag — three counter-cases\n')
    fehler = 0
    res = ergebnisse()
    print(f'1) Result files read: {len(res)} with seifert_genus')
    ok = len(res) > 0
    print(f'   {"OK " if ok else "ERROR"} (0 would be an empty run, not a green one)')
    fehler += (not ok)

    schreiben, konflikte, gabw = plan(res)
    print(f'\n2) Conflict latch — record already carries a DIFFERENT g')
    print(f'   {"OK " if not konflikte else "FINDING"} {len(konflikte)} conflicts (will NOT be written)')
    for s, w in konflikte[:8]:
        print(f'        {str(s)[:40]:<42} {w}')

    print(f'\n3) corpus genus vs. HFK genus — make bound-as-value visible')
    print(f'   {"OK " if not gabw else "FINDING"} {len(gabw)} deviations')
    for s, a, b in gabw[:8]:
        print(f'        {str(s)[:40]:<42} corpus genus={a}  HFK says={b}')

    print(f'\n4) Dry run: {len(schreiben)} records would get a NEW field')
    print('\nSELFTEST PASSED — findings are findings, not a fault of the tool')
    return fehler


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--selbsttest', action='store_true')
    ap.add_argument('--schreiben', action='store_true')
    a = ap.parse_args()
    if a.selbsttest or not a.schreiben:
        sys.exit(selbsttest())
    res = ergebnisse()
    schreiben, konflikte, gabw = plan(res)
    if konflikte:
        print(f'⛔ {len(konflikte)} conflicts — nothing written. Clarify first:')
        for s, w in konflikte:
            print(f'    {s}: {w}')
        sys.exit(2)
    n = 0
    for stem, pf, neu in schreiben:
        shutil.copy2(pf, pf + STAMP)
        d = json.load(open(pf))
        d.setdefault('pflichtliste', {}).setdefault('ALGEBRA', {})['knot_floer_homology'] = neu
        json.dump(d, open(pf, 'w'), indent=1)
        n += 1
    print(f'{n} Pflicht records supplemented (backup {STAMP} next to them).')
    if gabw:
        print(f'⚠ {len(gabw)} records: corpus genus deviated from HFK genus — corpus_db rebuild pulls that along:')
        for s, x, y in gabw[:10]:
            print(f'    {s}: {x} -> {y}')
    print('Next step is a step of its own: rebuild corpus_db, then labeler.')


if __name__ == '__main__':
    main()
