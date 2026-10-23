#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
preflight.py — Analyses an .fseries before the pipeline starts.

Output (--json for wrapper scripts):
  {
    "fseries_path": "...",
    "knot_name": "orbit",
    "n_fourier_terms": 20,
    "offset": [0.0, 0.0, 0.0],
    "nonzero_j_indices": [1, 3, 5, ..., 19],
    "all_zero_j_indices": [2, 4, ..., 20],
    "symmetry_class": "pi_antisymmetric" | "pi_periodic" | "generic",
    "scanner_period_pi_flag": true | false,
    "tracker_fold_pi_antisym_flag": true | false,
    "rationale": "..."
  }

Symmetry classification (Discovery #70):
  - all odd j non-zero, all even j zero               → pi_antisymmetric
  - all even j non-zero, all odd j zero               → pi_periodic
  - otherwise                                          → generic
"""
import argparse
import json
import os
import sys


def parse_fseries(path):
    """As in scanner / kinematic_resolve, with length discriminator (Discovery #69)."""
    coeffs = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('%') or line.startswith('#'):
                continue
            try:
                row = [float(x) for x in line.split()]
            except ValueError:
                continue
            coeffs.append(row)
    if not coeffs:
        raise ValueError(f"No data in {path}")
    if len(coeffs[0]) == 3:
        offset = tuple(coeffs[0])
        harmonics = coeffs[1:]
    else:
        offset = (0.0, 0.0, 0.0)
        harmonics = coeffs
    return offset, harmonics


def classify_symmetry(harmonics, eps=1e-12):
    """Detects pi antisymmetry / pi periodicity from the non-zero j indices."""
    nonzero_j = []
    all_zero_j = []
    for idx, row in enumerate(harmonics):
        j = idx + 1
        if any(abs(v) > eps for v in row):
            nonzero_j.append(j)
        else:
            all_zero_j.append(j)

    odd_nonzero  = [j for j in nonzero_j if j % 2 == 1]
    even_nonzero = [j for j in nonzero_j if j % 2 == 0]

    if nonzero_j and not even_nonzero:
        sym = "pi_antisymmetric"
        rationale = (f"all non-zero j are odd ({odd_nonzero}); "
                     f"all even j zero. → x(s+pi) = -x(s); Discovery #70.")
    elif nonzero_j and not odd_nonzero:
        sym = "pi_periodic"
        rationale = (f"all non-zero j are even ({even_nonzero}); "
                     f"all odd j zero. → x(s+pi) = +x(s); Discovery #70.")
    else:
        sym = "generic"
        rationale = (f"mixed j parity (odd: {odd_nonzero}, "
                     f"even: {even_nonzero}); no pi symmetry.")

    return sym, nonzero_j, all_zero_j, rationale


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('fseries', help='Path to the .fseries file')
    ap.add_argument('--json', action='store_true',
                    help='Machine-readable JSON output (for wrappers)')
    args = ap.parse_args()

    if not os.path.exists(args.fseries):
        print(f"ERROR: file not found: {args.fseries}", file=sys.stderr)
        sys.exit(2)

    offset, harmonics = parse_fseries(args.fseries)
    sym, nonzero_j, all_zero_j, rationale = classify_symmetry(harmonics)

    knot_name = os.path.splitext(os.path.basename(args.fseries))[0]

    period_pi_flag       = (sym == "pi_periodic")
    fold_pi_antisym_flag = (sym == "pi_antisymmetric")

    result = {
        "fseries_path":                args.fseries,
        "knot_name":                   knot_name,
        "n_fourier_terms":             len(harmonics),
        "offset":                      list(offset),
        "nonzero_j_indices":           nonzero_j,
        "all_zero_j_indices":          all_zero_j,
        "symmetry_class":              sym,
        "scanner_period_pi_flag":      period_pi_flag,
        "tracker_fold_pi_antisym_flag": fold_pi_antisym_flag,
        "rationale":                   rationale,
    }

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print(f"[PREFLIGHT] File:             {args.fseries}")
        print(f"[PREFLIGHT] Knot name:        {knot_name}")
        print(f"[PREFLIGHT] Fourier terms:    {len(harmonics)}")
        print(f"[PREFLIGHT] Offset:           {offset}")
        print(f"[PREFLIGHT] non-zero j:       {nonzero_j}")
        print(f"[PREFLIGHT] all-zero j:       {all_zero_j}")
        print(f"[PREFLIGHT] Symmetry class:   {sym}")
        print(f"[PREFLIGHT] Rationale:        {rationale}")
        print(f"[PREFLIGHT] Recommendation scanner --period-pi:        "
              f"{'YES' if period_pi_flag else 'NO'}")
        print(f"[PREFLIGHT] Recommendation tracker --fold-pi-antisym:  "
              f"{'YES' if fold_pi_antisym_flag else 'NO'}")


if __name__ == '__main__':
    main()
