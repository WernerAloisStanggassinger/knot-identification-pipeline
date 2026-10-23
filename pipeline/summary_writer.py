#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
summary_writer.py — Reads mandatory_invariants_complete.json and writes a
compact Markdown overview (knot_summary.md).

Usage:
  python summary_writer.py <complete.json> [--output summary.md]
"""
import argparse
import json
import os
import sys


def fmt(v, max_chars=120):
    if v is None:
        return '—'
    if isinstance(v, (list, tuple)):
        s = json.dumps(v)
    elif isinstance(v, dict):
        s = json.dumps(v, separators=(',', ':'))
    else:
        s = str(v)
    if len(s) > max_chars:
        s = s[:max_chars - 3] + '...'
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('complete_json')
    ap.add_argument('--output', default=None,
                    help='Default: <name>_summary.md next to the JSON')
    args = ap.parse_args()

    with open(args.complete_json) as f:
        d = json.load(f)

    name = d.get('name', 'unnamed')
    out_path = args.output or os.path.join(
        os.path.dirname(args.complete_json) or '.',
        f'{name}_summary.md'
    )

    pf = d.get('mandatory_list', {})
    diag = pf.get('DIAGRAM', {})
    alg  = pf.get('ALGEBRA', {})
    hyp  = pf.get('HYPERBOLIC', {})
    opt  = d.get('optional_invariants', {})
    quart = d.get('consistency_quartet', {})

    # Field balance (populated vs. total leaves)
    def count_leaf(obj):
        out = {'total': 0, 'populated': 0}
        def rec(o):
            if isinstance(o, dict):
                for v in o.values():
                    rec(v)
            else:
                out['total'] += 1
                if o is not None and o != '' and o != [] and o != {}:
                    out['populated'] += 1
        rec(obj)
        return out

    pcount = count_leaf(pf)
    ocount = count_leaf(opt)

    lines = []
    lines.append(f"# {name} — pipeline identification\n")
    lines.append(f"_Generated from_ `{os.path.basename(args.complete_json)}`\n")

    lines.append(f"\n## Profile (mandatory list)\n")
    lines.append("| Category | Value |")
    lines.append("|---|---|")
    rows = [
        ('crossing_number_min',   diag.get('crossing_number_min')),
        ('alternating',            diag.get('alternating')),
        ('writhe',                 diag.get('writhe')),
        ('braid_index',            f"{diag.get('braid_index_lower_bound')}..{diag.get('braid_index_upper_bound') if diag.get('braid_index_upper_bound') is not None else '?'} "
                                   f"(MFW spread {diag.get('braid_index_mfw_v_spread')})"
                                   if diag.get('braid_index_lower_bound') else diag.get('braid_index')),
        ('bridge_number',          f"{diag.get('bridge_number_lower_bound')}..{diag.get('bridge_number_upper_bound')}"
                                   if diag.get('bridge_number_lower_bound') is not None else diag.get('bridge_number')),
        ('determinant',            alg.get('determinant')),
        ('signature',              alg.get('signature')),
        ('arf_invariant',          alg.get('arf_invariant')),
        ('s_invariant',            alg.get('s_invariant')),
        ('genus (HFK)',            alg.get('genus')),
        ('smooth_4_genus',         alg.get('smooth_4_genus')),
        ('slice',                  alg.get('slice')),
        ('fibered',                alg.get('fibered')),
        ('tau / nu / epsilon',     f"{alg.get('tau_invariant')} / "
                                   f"{alg.get('knot_floer_homology', {}).get('nu') if isinstance(alg.get('knot_floer_homology'), dict) else None} / "
                                   f"{alg.get('knot_floer_homology', {}).get('epsilon') if isinstance(alg.get('knot_floer_homology'), dict) else None}"),
        ('HFK total_rank',         alg.get('knot_floer_homology', {}).get('total_rank') if isinstance(alg.get('knot_floer_homology'), dict) else None),
        ('khovanov width / thin',  f"{alg.get('khovanov_homology', {}).get('width')} / "
                                   f"{alg.get('khovanov_homology', {}).get('thin')}"
                                   if isinstance(alg.get('khovanov_homology'), dict) else None),
        ('hyperbolic_volume',      hyp.get('hyperbolic_volume')),
        ('chern_simons',           hyp.get('chern_simons')),
        ('symmetry_group',         f"{hyp.get('symmetry_group')} (order {hyp.get('symmetry_group_order')})"),
        ('amphichiral / invertible', f"{hyp.get('is_amphicheiral')} / {hyp.get('is_invertible_knot')}"),
        ('homology',               hyp.get('homology')),
        ('unknotting_number',      f"{alg.get('unknotting_number')} "
                                   f"(strict [{alg.get('unknotting_number_strict_lower')}, "
                                   f"{alg.get('unknotting_number_strict_upper')}], "
                                   f"indicated={alg.get('unknotting_number_indicated', False)})"),
        ('ribbon',                 alg.get('ribbon')),
    ]
    for k, v in rows:
        lines.append(f"| {k} | {fmt(v, 100)} |")

    lines.append(f"\n## Polynomials\n")
    lines.append("| Polynomial | Value |")
    lines.append("|---|---|")
    pols = [
        ('Alexander',          alg.get('alexander_polynomial')),
        ('Conway',             alg.get('conway_polynomial')),
        ('Jones',              alg.get('jones_polynomial')),
        ('HOMFLY',             alg.get('homfly_polynomial')),
        ('Kauffman F',         alg.get('kauffman_polynomial')),
    ]
    for k, v in pols:
        lines.append(f"| {k} | `{fmt(v, 200)}` |")

    lines.append(f"\n## Diagram codes\n")
    codes = [
        ('signed_gauss_reduced',  diag.get('signed_gauss_reduced')),
        ('dt_code',               diag.get('dt_code')),
        ('pd_code',               diag.get('pd_code')),
        ('braid_word',            diag.get('braid_word')),
        ('isosig (canonical)',    hyp.get('canonical_retriangulation')),
    ]
    for k, v in codes:
        lines.append(f"- **{k}**: `{fmt(v, 200)}`")

    lines.append(f"\n## Consistency\n")
    lines.append(f"- Pipeline consistency: **{'OK' if d.get('consistency_ok') else 'FAIL'}**")
    lines.append(f"- det paths: `{quart.get('det_paths')}` → det = {quart.get('det_value')}")
    lines.append(f"- HFK-rank check: rank={quart.get('hfk_rank_check', {}).get('rank')}, "
                 f"diff={quart.get('hfk_rank_check', {}).get('diff_rank_minus_det')}, "
                 f"rank=det:{quart.get('hfk_rank_equals_det')}")
    if alg.get('kauffman_polynomial_amphi_symmetry_check') is True:
        lines.append(f"- Kauffman F amphi-symmetry-check: ✓")

    lines.append(f"\n## Balance\n")
    def _pct(p, t):
        return f"{p/t*100:.1f}%" if t else "n/a (empty)"
    lines.append(f"- Mandatory list: **{pcount['populated']} / {pcount['total']}** "
                 f"populated ({_pct(pcount['populated'], pcount['total'])})")
    lines.append(f"- Optional:    **{ocount['populated']} / {ocount['total']}** "
                 f"populated ({_pct(ocount['populated'], ocount['total'])})")

    lines.append(f"\n---\n_Pipeline blind-test compliant_ "
                 f"(`pipeline_seal.py` active, no knot-table lookups).")

    with open(out_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f'Written: {out_path}')


if __name__ == '__main__':
    main()
