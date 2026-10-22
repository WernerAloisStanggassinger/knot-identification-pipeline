#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
enrich_bridge_number_bounds.py
==============================
Standalone enrichment: reads an existing pflicht_invarianten_*.json
and adds bridge-number bounds to the DIAGRAMM block, plus the trivial
upper bound and the alternating heuristic. Imports bridge_number_bounds() from
pflicht_invarianten_extra (shared logic).

Background:
  Discovery #102: 'bridge_number' from cyclic O-blocks of the signed Gauss code
    was trivially equal to c for alternating knots (Discovery #93), i.e. not
    informative.
  Discovery #109 (CORRECTION): the Schubert bound 'b(K) >= g(K)+1' quoted
    here earlier is NOT a general theorem (counterexample: the
    (5,2) torus knot, g=2, b=2). Lower bound now: b >= 2 for a
    nontrivial knot. Sharpenings via the Wirtinger meridional_rank in
    enrich_optional_invariants.py, wave 5.

Usage:
  python enrich_bridge_number_bounds.py <pflicht_*.json> [--in-place]
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from pflicht_invarianten_extra import bridge_number_bounds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('json_path')
    ap.add_argument('--in-place', action='store_true')
    args = ap.parse_args()

    p = Path(args.json_path)
    with p.open() as f:
        report = json.load(f)

    a = report['pflichtliste']['ALGEBRA']
    d = report['pflichtliste']['DIAGRAMM']
    g = a.get('genus') or a.get('knot_floer_homology', {}).get('seifert_genus')
    c = d.get('crossing_number_min')
    det = a.get('determinant_sage') or a.get('determinant')
    altern = d.get('alternating', False)

    if g is None or c is None:
        print(f'ERROR: genus={g} or crossing_number_min={c} is missing.', file=sys.stderr)
        sys.exit(2)

    bnb = bridge_number_bounds(genus=g, crossing_number=c, det=det, alternating=altern)

    old_lower = d.get('bridge_number_lower_bound')
    old_interval = d.get('bridge_number_interval')

    d['bridge_number_lower_bound']  = bnb['lower_bound']
    d['bridge_number_upper_bound']  = bnb['upper_bound']
    d['bridge_number_interval']     = bnb['interval']
    d['bridge_number_is_2_excluded'] = bnb['is_2_bridge_excluded']
    d['bridge_number_lower_source'] = bnb['lower_bound_source']
    d['bridge_number_upper_source'] = bnb['upper_bound_source']
    if 'alt_upper_bound_heuristic' in bnb:
        d['bridge_number_alt_heuristic']        = bnb['alt_upper_bound_heuristic']
        d['bridge_number_alt_heuristic_source'] = bnb['alt_upper_bound_heuristic_source']

    out_path = p if args.in_place else p.with_suffix('.bridge.json')
    with out_path.open('w') as f:
        json.dump(report, f, indent=2, default=str)

    print(f'==== Bridge-number bounds ====')
    print(f'Input   : {p}')
    print(f'Output  : {out_path}')
    print()
    print(f'Input values: g={g}, c={c}, det={det}, alternating={altern}')
    print()
    print(f'Before: bridge_number_lower_bound = {old_lower}')
    print(f'        bridge_number_interval    = {old_interval}')
    print(f'After:   bridge_number_interval   = {bnb["interval"]}')
    print(f'         is_2_bridge_excluded     = {bnb["is_2_bridge_excluded"]}')
    print(f'         lower_bound_source       = {bnb["lower_bound_source"]}')
    print(f'         upper_bound_source       = {bnb["upper_bound_source"]}')
    if 'alt_upper_bound_heuristic' in bnb:
        print(f'         alt_upper_bound (heur)   = {bnb["alt_upper_bound_heuristic"]}')


if __name__ == '__main__':
    main()
