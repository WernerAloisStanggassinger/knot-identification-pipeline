#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
enrich_unknotting_signed_gauss.py
=================================
Standalone enrichment: runs the signed-Gauss crossing-flip brute force
(Discovery #100) on an existing pflicht_invarianten_*.json and
writes unknotting_number_{lower,upper}_bound + the witness subset back.

Background:
  Discovery #93 documented u(K) as not determinable for Saturn,
  because brute force on the non-minimal 9-strand × 40-letter Sage braid
  (Discovery #60) gave only a loose bound. With the signed-Gauss crossing flip
  (Discovery #100) the search on the 12-crossing diagram is tractable in
  seconds.

Usage:
  python enrich_unknotting_signed_gauss.py <pflicht_*.json> [--max-subset 4] [--processes 12] [--in-place]
"""
import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from pflicht_invarianten_extra import unknotting_via_signed_gauss_flip


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('json_path', help='pflicht_invarianten_*.json')
    ap.add_argument('--max-subset', type=int, default=4)
    ap.add_argument('--processes',  type=int, default=12)
    ap.add_argument('--in-place',   action='store_true',
                    help='Write the enrichment directly into the input file')
    args = ap.parse_args()

    p = Path(args.json_path)
    with p.open() as f:
        report = json.load(f)

    sg = report.get('signed_gauss')
    if not sg:
        print('ERROR: signed_gauss is missing from the JSON.', file=sys.stderr)
        sys.exit(2)
    a = report['pflichtliste']['ALGEBRA']
    det = a.get('determinant_sage') or a.get('determinant')
    jones = (a.get('jones_polynomial') or '').strip()

    # Lower bound from det
    if det is not None and abs(int(det)) == 1 and jones in ('1', 't^0'):
        print('Pre-check: |det|=1 ∧ Jones=1 ⇒ trivial knot, u = 0')
        a['unknotting_number_upper_bound'] = 0
        a['unknotting_number_lower_bound'] = 0
        a['unknotting_number'] = 0
        a['unknotting_method'] = 'pre_check_det_jones'
    else:
        u_lower = 1 if det is not None else 0
        print(f'Pre-check: det={det}, Jones polynomial non-trivial ⇒ u >= {u_lower}')
        print(f'Brute-Force signed-gauss-flip max_subset={args.max_subset}, '
              f'processes={args.processes} ...')
        t0 = time.time()
        unk = unknotting_via_signed_gauss_flip(sg,
                                                max_subset=args.max_subset,
                                                processes=args.processes)
        dt = time.time() - t0
        print(f'  ... {dt:.1f}s')
        u_upper = unk.get('unknotting_number_upper_bound')
        a['unknotting_number_upper_bound'] = u_upper
        a['unknotting_number_lower_bound'] = u_lower
        a['unknotting_witness_subset']     = unk.get('witness_subset')
        a['unknotting_method']             = unk.get('method')
        a['unknotting_number_searched_to'] = unk.get('searched_up_to_subset_size')
        # Discovery #114c: unknotting_number = upper_bound for a stage-1 witness
        # (kept in sync with pflicht_invarianten_extra.py, Discovery #114c).
        if u_upper is None:
            a['unknotting_number']              = None
            a['unknotting_number_strict_lower'] = u_lower
            a['unknotting_number_strict_upper'] = None
            a['unknotting_number_note']         = (
                f'u > {args.max_subset} on this diagram')
        else:
            progs = unk.get('progress', []) or []
            ks_searched = sorted({p['k'] for p in progs if 'k' in p})
            ks_required = list(range(u_lower, u_upper))
            ks_with_hits = {p['k'] for p in progs
                            if (p.get('unknot_simplify_zero', 0) > 0
                                or p.get('unknot_jones_only', 0) > 0)}
            exhaustive_below_upper = (
                set(ks_required).issubset(set(ks_searched))
                and not (set(ks_required) & ks_with_hits)
            )

            if u_lower == u_upper:
                a['unknotting_number']              = u_upper
                a['unknotting_number_strict_lower'] = u_lower
                a['unknotting_number_strict_upper'] = u_upper
                a['unknotting_number_assignment_basis'] = (
                    f'lower == upper == {u_upper}: eindeutig.')
            elif exhaustive_below_upper and unk.get('witness_subset'):
                a['unknotting_number']              = u_upper
                a['unknotting_number_strict_lower'] = u_lower
                a['unknotting_number_strict_upper'] = u_upper
                a['unknotting_number_indicated']    = True
                a['unknotting_number_assignment_basis'] = (
                    f'Discovery #114c: Stage-1-Witness {unk.get("witness_subset")} '
                    f'(Regina simplify→0) at k={u_upper}; k={ks_required} all '
                    f'exhaustive without hit. diagram_u(K) = {u_upper} provable; '
                    f'global u(K) could be < {u_upper} on another diagram, '
                    f'hence strict_lower = {u_lower} stays formal.')
            else:
                a['unknotting_number']              = None
                a['unknotting_number_strict_lower'] = u_lower
                a['unknotting_number_strict_upper'] = u_upper
                a['unknotting_number_assignment_basis'] = (
                    f'k={ks_required} not completely tested '
                    f'(searched k={ks_searched}, with hits={sorted(ks_with_hits)}); '
                    f'no Discovery #114c auto-assign.')

            a['unknotting_number_note'] = (
                f'{u_lower} <= u(K) <= {u_upper}; Witness-Subset '
                f'{unk.get("witness_subset")} gives the unknot via crossing flip')
        report['raw_subprocess_outputs']['unknotting'] = unk

        # Progress report
        for prog in unk.get('progress', []):
            if 'k' in prog:
                print(f"  k={prog['k']:>2d}: {prog['n_subsets']:>4d} subsets, "
                      f"simplify→0={prog.get('unknot_simplify_zero', 0)}, "
                      f"jones=1={prog.get('unknot_jones_only', 0)}")
        if u_upper is not None:
            print(f'\nu(K) ∈ [{u_lower}, {u_upper}]; '
                  f'Witness {unk.get("witness_subset")}')

    out_path = p if args.in_place else p.with_suffix('.unknotting.json')
    with out_path.open('w') as f:
        json.dump(report, f, indent=2, default=str)
    print(f'\nWritten: {out_path}')


if __name__ == '__main__':
    main()
