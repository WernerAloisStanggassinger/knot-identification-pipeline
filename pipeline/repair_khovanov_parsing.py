#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
repair_khovanov_parsing.py
==========================
Standalone repair/enrichment: reads an existing
mandatory_invariants_*.json with an (old, faulty or incomplete)
Khovanov evaluation, applies the corrected parser to the stored
raw_subprocess_outputs.knotjob.raw_stdout and writes structured
fields back.

Fields that are computed/updated (Discovery #97 + #101):
  - polynomial, terms_t_q, delta_diagonals, n_distinct_diagonals, width, thin
  - torsion (list per order), torsion_orders, has_torsion, z2_thin_only
  - reduced_polynomial, reduced_terms_t_q, reduced_delta_diagonals,
    reduced_n_distinct_diagonals, reduced_width, reduced_thin

Usage:
  python repair_khovanov_parsing.py <pflicht_*.json> [--in-place]
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from mandatory_invariants_extra import parse_knotjob_khovanov_stdout


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('json_path', help='pflicht_*.json (with raw_subprocess_outputs.knotjob.raw_stdout)')
    ap.add_argument('--in-place', action='store_true',
                    help='Write the repair directly into the input file')
    args = ap.parse_args()

    p = Path(args.json_path)
    with p.open() as f:
        report = json.load(f)

    knotjob = report.get('raw_subprocess_outputs', {}).get('knotjob', {})
    raw_stdout = knotjob.get('raw_stdout', '')
    if not raw_stdout:
        # Disc #289: an enrichment without input is NOT an error. For monster knots
        # (c>~40) KnotJob (Khovanov) times out -> no raw_stdout. This repair
        # must NOT kill the harvest (formerly sys.exit(2) -> stage 7 fatal ->
        # rc=3 -> a knot with real core invariants was lost). Skip = ok.
        print('SKIP: no raw_subprocess_outputs.knotjob.raw_stdout — KnotJob '
              'timed out/empty (monster?). Nothing to repair, the knot keeps '
              'its core invariants.', file=sys.stderr)
        sys.exit(0)

    # Old values for the diff
    old_kh = report.get('mandatory_list', {}).get('ALGEBRA', {}).get('khovanov_homology', {})
    old_terms = old_kh.get('terms_t_q', {}) or {}
    old_width = old_kh.get('width')
    old_thin  = old_kh.get('thin')
    old_torsion = old_kh.get('torsion')
    old_reduced = old_kh.get('reduced_polynomial')

    # Reparse via helper (same logic as the main() path in mandatory_invariants_extra.py)
    parsed = parse_knotjob_khovanov_stdout(raw_stdout)
    if not parsed.get('polynomial'):
        # Disc #289: raw_stdout present, but no parsable Khovanov block
        # (truncated/partial output for large knots). That too is only a missing
        # enrichment, not a harvest killer -> skip instead of sys.exit(3).
        print('SKIP: no parsable "Integral unreduced Khovanov Homology" block '
              'in raw_stdout. Nothing to repair, the knot keeps its core invariants.',
              file=sys.stderr)
        sys.exit(0)

    # Persist into ALGEBRA (same fields as main())
    new_kh = {
        'polynomial':                   parsed.get('polynomial'),
        'terms_t_q':                    parsed.get('terms_t_q'),
        'delta_diagonals':              parsed.get('delta_diagonals'),
        'n_distinct_diagonals':         parsed.get('n_distinct_diagonals'),
        'width':                        parsed.get('width'),
        'thin':                         parsed.get('thin'),
        'torsion':                      parsed.get('torsion'),
        'torsion_orders':               parsed.get('torsion_orders'),
        'has_torsion':                  parsed.get('has_torsion'),
        'z2_thin_only':                 parsed.get('z2_thin_only'),
        'reduced_polynomial':           parsed.get('reduced_polynomial'),
        'reduced_terms_t_q':            parsed.get('reduced_terms_t_q'),
        'reduced_delta_diagonals':      parsed.get('reduced_delta_diagonals'),
        'reduced_n_distinct_diagonals': parsed.get('reduced_n_distinct_diagonals'),
        'reduced_width':                parsed.get('reduced_width'),
        'reduced_thin':                 parsed.get('reduced_thin'),
    }
    report['mandatory_list']['ALGEBRA']['khovanov_homology'] = new_kh

    # Synchronise raw_subprocess_outputs.knotjob
    knotjob.update({k: v for k, v in parsed.items() if k != 'raw_stdout'})
    report['raw_subprocess_outputs']['knotjob'] = knotjob

    out_path = p if args.in_place else p.with_suffix('.repaired.json')
    with out_path.open('w') as f:
        json.dump(report, f, indent=2, default=str)

    # Diff report
    print(f'==== Khovanov parsing repair ====')
    print(f'Input  : {p}')
    print(f'Output : {out_path}')
    print()
    print(f'Free-part:')
    print(f'  Number of terms: old={len(old_terms)}, new={len(new_kh["terms_t_q"] or {})}')
    print(f'  Width:         old={old_width}, new={new_kh["width"]}')
    print(f'  Thin:          old={old_thin}, new={new_kh["thin"]}')
    print(f'  δ-diagonals:   {new_kh["delta_diagonals"]}')
    print()
    print(f'Torsion:')
    print(f'  before:                {"present" if old_torsion else "MISSING"}')
    print(f'  has_torsion:           {new_kh["has_torsion"]}')
    print(f'  torsion_orders:        {new_kh["torsion_orders"]}')
    if new_kh['torsion']:
        for t in new_kh['torsion']:
            print(f'  order {t["order"]}: {t["total_torsion_count"]} generators, '
                  f'{len(t["terms_t_q"])} bigradings')
    print(f'  z2_thin_only:          {new_kh["z2_thin_only"]}')
    print()
    print(f'Reduced Khovanov:')
    print(f'  before:                {"present" if old_reduced else "MISSING"}')
    if new_kh.get('reduced_polynomial'):
        print(f'  reduced_width:         {new_kh["reduced_width"]}')
        print(f'  reduced_thin:          {new_kh["reduced_thin"]}')
        print(f'  reduced δ-diagonals:   {new_kh["reduced_delta_diagonals"]}')
        print(f'  reduced n_terms:       {len(new_kh["reduced_terms_t_q"])}')


if __name__ == '__main__':
    main()
