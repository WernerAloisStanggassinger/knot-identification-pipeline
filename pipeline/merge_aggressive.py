#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
merge_aggressive.py — Patches kinematic_result.json with the better diagram
from aggressive_simplify.json, if aggressive_simplify found a diagram with a
smaller crossing_number.

Background: Discovery #91 declares aggressive_simplify MANDATORY after
kinematic_resolve, because single-frame stagnation is typical (~3% of the frames
find the Reidemeister-3 sequence that leads to further reduction).
The pipeline stage pflicht_invarianten.py, however, takes kinematic_result.json
as input, NOT aggressive_simplify.json — so this step must be inserted before
stage 5.

Usage:
  python merge_aggressive.py <kinematic_result.json> <aggressive_simplify.json>
    [--in-place | --output <new_kinematic.json>]

If aggressive.best.size_simp_regina < kinematic.regina.simplified_size:
  - regina.signed_gauss      ← aggressive.best.signed_gauss_simp
  - regina.simplified_size   ← aggressive.best.size_simp_regina
  - regina.homfly            ← aggressive.best.homfly
  - regina.jones             ← aggressive.best.jones
  - chosen.signed_gauss      ← aggressive.best.signed_gauss_simp
  - merge_meta               ← {original_size, aggressive_size, source}
Otherwise: kinematic_result unchanged; merge_meta records "no merge needed".
"""
import argparse
import json
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('kinematic_json')
    ap.add_argument('aggressive_json')
    ap.add_argument('--in-place', action='store_true')
    ap.add_argument('--output', default=None)
    args = ap.parse_args()

    with open(args.kinematic_json) as f:
        kin = json.load(f)
    with open(args.aggressive_json) as f:
        agg = json.load(f)

    kin_size = kin.get('regina', {}).get('simplified_size')
    agg_best = agg.get('best', {}) or {}
    agg_size = agg_best.get('size_simp_regina')
    agg_min  = agg.get('min_size_global')

    # ── Disc #562 (2026-07-25): pass on the SMALLEST code that has a diagram receipt ──
    # Before, this was hard-wired to agg_best.signed_gauss_simp = Regina's simplification.
    # Since Disc #562 aggressive_simplify additionally delivers 'handoff': the smallest
    # code for which a diagram exists (spherogram-reduced, converted back via
    # regina.fromPD, det cross-checked). If the key is missing (old JSONs), everything
    # stays as before — backwards compatible.
    handoff  = agg.get('handoff') or {}
    ho_size  = handoff.get('size')
    ho_sg    = handoff.get('signed_gauss')
    ho_src   = handoff.get('source') or 'regina'
    use_ho   = bool(ho_sg) and ho_size is not None and (agg_size is None or ho_size < agg_size)
    if use_ho:
        print(f'[MERGE] Handoff receipt (Disc #562): {ho_size} crossings from "{ho_src}" '
              f'(regina-best would be {agg_size}) — det cross-checked in stage 4.')
        agg_size = ho_size
        agg_sg   = ho_sg
    else:
        agg_sg   = agg_best.get('signed_gauss_simp')
    if agg_min is not None and agg_size is not None and agg_min < agg_size:
        print(f'[MERGE] ⚠️  measured min_size_global={agg_min} < passed-on '
              f'{agg_size} — the diagram receipt for the optimum is missing '
              f'(crossing_number_min stays an upper bound, Disc #562).')

    if kin_size is None:
        print(f'ERROR: kinematic_result.regina.simplified_size is missing', file=sys.stderr)
        sys.exit(2)

    if agg_size is None:
        # ── Root fix 2026-06-17 (Werner's go) ────────────────────────────────
        # aggressive_simplify delivered NO usable best diagram. Actually
        # measured case (novelty_run 06-16): at high c (c≈54+) the Regina sweep
        # runs into global_timeout -> regina_status.partial=True,
        # global_timed_out=True, the 'best' key is missing ENTIRELY, min_size_global=None,
        # no_classical_found set.
        # BEFORE: sys.exit(2) -> stage 4a dead -> whole knot lost, ALTHOUGH
        # kinematic_result.json carries a fully valid diagram (stage 3).
        # NOW: graceful fallback to the kinematic diagram. Stage 5 computes
        # det/Alex/Jones/HOMFLY/vol/sym/theta on it; crossing_number_min is then an
        # upper bound instead of minimal (corpus convention Disc #231). This fix thereby
        # automatically does what the manual crutch SKIP_AGGRESSIVE_SIMPLIFY=1 did before
        # — but targeted at exactly the knots for which stage 4 found nothing better.
        rstat = agg.get('regina_status') or {}
        reason = []
        if rstat.get('global_timed_out'): reason.append('global_timeout')
        if rstat.get('partial'):          reason.append('regina_partial')
        if agg.get('no_classical_found'): reason.append('no_classical_found')
        if agg.get('stage4_partial'):     reason.append('stage4_partial')
        reason = reason or ['kein_brauchbares_best']
        print(f'[MERGE] aggressive_simplify without usable best ({",".join(reason)}) '
              f'⇒ fallback to kinematic diagram (size={kin_size}); stage 5 continues.')
        if rstat:
            print(f'        regina: {rstat.get("n_codes_processed")}/{rstat.get("n_total_codes")} '
                  f'codes in {rstat.get("duration_s")}s, partial={rstat.get("partial")}')
        kin['merge_meta'] = {
            'merged_from':          None,
            'kinematic_size':       kin_size,
            'aggressive_best_size': None,
            'aggressive_partial':   True,
            'fallback_reason':      reason,
            'regina_status':        rstat or None,
            'note':                 'aggressive_simplify lieferte kein brauchbares '
                                    'Diagramm (Hoch-c/Timeout) — kinematic-Diagramm '
                                    'als Fallback. crossing_number_min = Obergrenze '
                                    '(nicht-minimal, Disc #231).',
        }
        out_path = args.output or (args.kinematic_json if args.in_place
                                    else args.kinematic_json.replace('.json', '_merged.json'))
        with open(out_path, 'w') as f:
            json.dump(kin, f, indent=2)
        print(f'[MERGE] Written (fallback): {out_path}')
        return

    if agg_size < kin_size:
        print(f'[MERGE] aggressive_simplify has a better diagram:')
        print(f'        kinematic regina.simplified_size = {kin_size}')
        print(f'        aggressive best.size_simp_regina = {agg_size}')
        print(f'        aggressive min_size_global       = {agg_min}')
        print(f'        ⇒ patching kinematic_result with aggressive-best.')

        kin['regina']['signed_gauss']      = agg_sg
        kin['regina']['simplified_size']   = agg_size
        # With a spherogram handoff, NO polynomials from the Regina diagram travel along
        # (Frankenstein doctrine Disc #256: no foreign algebra on the diagram).
        # Stage 5 computes the algebra fresh on the passed-on code anyway.
        kin['regina']['homfly']            = None if use_ho else agg_best.get('homfly')
        kin['regina']['jones']             = None if use_ho else agg_best.get('jones')
        kin['chosen']['signed_gauss']      = agg_sg
        kin['merge_meta'] = {
            'merged_from':            'aggressive_simplify (Discovery #91)',
            'original_kinematic_size': kin_size,
            'aggressive_best_size':    agg_size,
            'aggressive_min_global':   agg_min,
            'aggressive_t_step':       handoff.get('t_step') if use_ho
                                       else agg_best.get('t_step'),
            'handoff_source':          ho_src,
            'handoff_det_checked':     handoff.get('det') if use_ho else None,
            'regina_best_size':        agg_best.get('size_simp_regina'),
            'note':                    'Pflichtliste-Werte werden auf das '
                                       'reduzierte Diagramm aus '
                                       'aggressive_simplify gerechnet.',
        }
    else:
        print(f'[MERGE] no merge needed:')
        print(f'        kinematic regina.simplified_size = {kin_size}')
        print(f'        aggressive best.size_simp_regina = {agg_size}')
        print(f'        ⇒ kinematic is equal or better, unchanged.')
        kin['merge_meta'] = {
            'merged_from':           None,
            'kinematic_size':        kin_size,
            'aggressive_best_size':  agg_size,
            'note':                  'kinematic_result war bereits minimal '
                                     'oder besser als aggressive_simplify.',
        }

    out_path = args.output or (args.kinematic_json if args.in_place
                                else args.kinematic_json.replace('.json', '_merged.json'))
    with open(out_path, 'w') as f:
        json.dump(kin, f, indent=2)
    print(f'[MERGE] Written: {out_path}')


if __name__ == '__main__':
    main()
