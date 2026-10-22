#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
enrich_braid_index_mfw.py
=========================
Standalone enrichment: reads an existing pflicht_invarianten_*.json
and adds the Morton-Franks-Williams lower bound for braid_index(K) to the
DIAGRAMM block, computed from the homfly_polynomial already present.

Background (Discovery #99):
  The mandatory-list invariant 'braid_index' used to be reported as n_strands from
  Sage's K.braid() — that is an UPPER bound, NOT
  guaranteed minimal (Markov reduction is NP-hard, Discovery #56).
  MFW gives a genuine LOWER bound: b(K) >= (max v-exp - min v-exp)/2 + 1
  from the HOMFLY polynomial. This yields a clear interval [lower, upper].

Usage:
  python enrich_braid_index_mfw.py <pflicht_*.json> [--in-place]

Without --in-place: writes to <input>.enriched.json
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from pflicht_invarianten_extra import mfw_braid_index_bound


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('json_path', help='pflicht_invarianten_*.json')
    ap.add_argument('--in-place', action='store_true',
                    help='Write the enrichment directly into the input file')
    args = ap.parse_args()

    p = Path(args.json_path)
    with p.open() as f:
        report = json.load(f)

    a = report.get('pflichtliste', {}).get('ALGEBRA', {})
    d = report.get('pflichtliste', {}).get('DIAGRAMM', {})
    homfly = a.get('homfly_polynomial')
    n_strands = d.get('braid_n_strands') or d.get('braid_index')

    if not homfly:
        print('ERROR: homfly_polynomial is missing from the JSON.', file=sys.stderr)
        sys.exit(2)
    # Graceful skip instead of a hard exit: braid_n_strands/braid_index is LEGITIMATELY
    # missing for plat-closure knots (symmetric unions etc.). The MFW LOWER bound
    # needs only HOMFLY and is set anyway; only the upper bound (Sage strand count)
    # stays open. An optional enrichment must NEVER block stage 8 (summary)
    # (finding 2026-06-13: su_9_2_t1 det=361 dragged the whole pipeline down to rc=3).
    if not n_strands:
        print('NOTE: braid_n_strands / braid_index is missing — upper bound open, '
              'MFW lower bound is set nevertheless.', file=sys.stderr)
        n_strands = None

    mfw = mfw_braid_index_bound(homfly)
    if not mfw:
        print('ERROR: MFW parsing failed on homfly_polynomial.', file=sys.stderr)
        sys.exit(3)

    # Before values for the diff
    old_lower = d.get('braid_index_lower_bound')
    old_interval = d.get('braid_index_interval')

    d['braid_index_lower_bound']    = mfw['lower_bound']
    d['braid_index_upper_bound']    = n_strands
    d['braid_index_interval']       = [mfw['lower_bound'], n_strands]
    d['braid_index_mfw_v_spread']   = mfw['spread']
    if n_strands is not None:
        d['braid_index_mfw_note']   = (
            f"MFW-Theorem: b(K) >= {mfw['lower_bound']} aus HOMFLY-v-Spread "
            f"{mfw['spread']} (v-exp ∈ [{mfw['min_v_exp']}, {mfw['max_v_exp']}]); "
            f"obere Grenze = aktuelle Sage-Braid-Straenge = {n_strands} "
            f"(NICHT garantiert minimal — Markov-Reduktion NP-schwer, "
            f"Discovery #56). Wenn lower_bound == upper_bound, ist b(K) "
            f"exakt bestimmt.")
    else:
        d['braid_index_mfw_note']   = (
            f"MFW-Theorem: b(K) >= {mfw['lower_bound']} aus HOMFLY-v-Spread "
            f"{mfw['spread']} (v-exp ∈ [{mfw['min_v_exp']}, {mfw['max_v_exp']}]); "
            f"obere Grenze offen — kein Braid-Wort im JSON (plat-closure-Knoten).")

    out_path = p if args.in_place else p.with_suffix('.enriched.json')
    with out_path.open('w') as f:
        json.dump(report, f, indent=2, default=str)

    print(f'==== braid_index MFW enrichment ====')
    print(f'Input   : {p}')
    print(f'Output  : {out_path}')
    print()
    print(f'HOMFLY variable detected: {mfw["variable_used"]}')
    print(f'v-exp Range: [{mfw["min_v_exp"]}, {mfw["max_v_exp"]}]  '
          f'(Spread {mfw["spread"]})')
    print()
    print(f'Before: braid_index_lower_bound = {old_lower}')
    print(f'        braid_index_interval    = {old_interval}')
    print(f'After:   braid_index_lower_bound = {mfw["lower_bound"]}')
    print(f'         braid_index_upper_bound = {n_strands}')
    print(f'         braid_index_interval    = [{mfw["lower_bound"]}, {n_strands}]')
    if n_strands is not None and mfw['lower_bound'] == n_strands:
        print(f'\nb(K) = {n_strands} EXACT — MFW bound attained sharply.')
    elif n_strands is not None:
        print(f'\nb(K) lies in the interval — exact value undetermined '
              f'(Markov reduction NP-hard, Discovery #56).')
    else:
        print(f'\nOnly the lower bound b(K) >= {mfw["lower_bound"]} set — '
              f'upper open (no braid word, plat-closure knot).')


if __name__ == '__main__':
    main()
