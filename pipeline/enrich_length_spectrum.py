#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
enrich_length_spectrum.py
=========================
Standalone enrichment: computes length_spectrum_short via
snappy.Manifold.length_spectrum_alt_gen() and writes it into the HYPERBOLISCH
block. Works for high-symmetry manifolds (e.g. Saturn D_6),
on which the Discovery #95 Dirichlet approach fails.

Background (Discovery #103):
  The mandatory-list invariant 'length_spectrum_short' was documented on Saturn as
  Discovery #95 (Dirichlet failure). With length_spectrum_alt
  (a numerical method independent of the Dirichlet polyhedron) it is
  computable after all. First geodesic for Saturn ~70s; further ones come
  out of the generator in <1s.

Usage:
  python enrich_length_spectrum.py <pflicht_*.json> [--n-geos 10] [--timeout 600] [--in-place]
"""
import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from mandatory_invariants_extra import length_spectrum_via_canonize


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('json_path')
    ap.add_argument('--n-geos',  type=int, default=10)
    ap.add_argument('--timeout', type=int, default=600)
    ap.add_argument('--in-place', action='store_true')
    args = ap.parse_args()

    p = Path(args.json_path)
    with p.open() as f:
        report = json.load(f)

    h = report['mandatory_list']['HYPERBOLIC']
    d = report['mandatory_list']['DIAGRAM']
    pd = d.get('pd_code')
    isosig = h.get('canonical_retriangulation') or \
             report.get('raw_subprocess_outputs', {}).get('snappy', {}).get('isosig')

    if not pd and not isosig:
        print('ERROR: neither pd_code nor isosig in the JSON.', file=sys.stderr)
        sys.exit(2)

    print(f'Input path: {"isosig" if isosig else "pd_code"}')
    if isosig:
        print(f'  isosig: {isosig[:60]}...')
    print(f'  n_geos:  {args.n_geos}')
    print(f'  timeout: {args.timeout}s')
    print(f'  Volume   : {h.get("hyperbolic_volume")}')
    print()

    t0 = time.time()
    print(f'Starting length_spectrum_alt_gen ...')
    lsp = length_spectrum_via_canonize(pd, n_geos=args.n_geos,
                                        timeout=args.timeout, isosig=isosig)
    dt = time.time() - t0

    if 'error' in lsp:
        print(f'\nERROR after {dt:.1f}s:')
        print(f'  {lsp.get("error")}')
        if 'stderr' in lsp:
            print(f'  stderr: {lsp["stderr"][-500:]}')
        sys.exit(3)

    short_geos = lsp.get('length_spectrum_short')
    method     = lsp.get('method')
    dur_inner  = lsp.get('duration_s')

    dur_str = f'{dur_inner:.1f}s' if isinstance(dur_inner, (int, float)) else 'n/a'
    print(f'\n{dt:.1f}s total (internal {dur_str}), method={method}')
    if short_geos:
        print(f'\n{len(short_geos)} geodesic clusters:')
        for i, g in enumerate(short_geos):
            if 'length_real' in g:
                print(f'  #{i:2d}  L = {g["length_real"]:.10f} + {g["length_imag"]:+.10f}i  '
                      f'(mult {g["multiplicity"]})')
            else:
                print(f'  #{i:2d}  L = {g.get("length", "?")}  (mult {g.get("mult", "?")})')
    else:
        errs = {k: v for k, v in lsp.items() if k.startswith('err')}
        print(f'\nNO geodesics found. errors:')
        for k, v in errs.items():
            print(f'  {k}: {v}')

    h['length_spectrum_short']      = short_geos
    h['length_spectrum_method']     = method
    h['length_spectrum_duration_s'] = dur_inner
    report['raw_subprocess_outputs']['length_spectrum'] = lsp

    out_path = p if args.in_place else p.with_suffix('.lspec.json')
    with out_path.open('w') as f:
        json.dump(report, f, indent=2, default=str)
    print(f'\nWritten: {out_path}')


if __name__ == '__main__':
    main()
