#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
detect_target.py — Reads crossings.parquet, computes the n_cross distribution per
frame and proposes target_count + frame_cap for kinematic_resolve.py.

Strategy (as of Discovery #108):
  - target_count = MIN(n_cross_distribution) — Saturn/orbit experience:
    the minimal n_cross value is usually the number of diagram crossings
    on the optimal frame (even if the mode/median of the distribution
    lies clearly higher).
  - frame_cap = target_count + 2 — tolerates near-optimum frames.
  - n_min, n_max for aggressive_simplify: [target_count, q75]
    (lower bound = true min, upper = 75% quantile = realistic
    reduction search).

Output (--json for wrappers):
  {
    "n_frames": 131072,
    "n_crossings_total": 3051374,
    "per_frame_n_cross": {
      "min": 14, "max": 46, "median": 22, "mean": 23.28,
      "modus": 22, "q05": 16, "q25": 18, "q75": 26, "q95": 34
    },
    "all_even_n_cross": true,         // π-antisym signature
    "all_odd_n_cross":  false,
    "target_count":     14,
    "frame_cap":        16,
    "agg_n_min":        14,
    "agg_n_max":        26,
    "rationale":        "..."
  }
"""
import argparse
import json
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('crossings_parquet')
    ap.add_argument('--json', action='store_true')
    args = ap.parse_args()

    df = pd.read_parquet(args.crossings_parquet)
    if 't_step' not in df.columns:
        print("ERROR: t_step column missing in the parquet", file=sys.stderr)
        sys.exit(2)

    g = df.groupby('t_step').size()

    n_min    = int(g.min())
    n_max    = int(g.max())
    n_median = float(g.median())
    n_mean   = float(g.mean())
    n_modus  = int(g.mode().iloc[0])
    n_q05    = int(g.quantile(0.05))
    n_q25    = int(g.quantile(0.25))
    n_q75    = int(g.quantile(0.75))
    n_q95    = int(g.quantile(0.95))

    n_even = int((g % 2 == 0).sum())
    n_odd  = int((g % 2 == 1).sum())
    all_even = (n_odd == 0)
    all_odd  = (n_even == 0)

    target_count = n_min
    frame_cap    = n_min + 2
    agg_n_min    = n_min
    # Discovery 2026-05-10 (Saturn c=12 finding): q75 is too narrow, because the
    # ~3% of frames that reach a reduction often lie in the right tail
    # (Discovery #92). q95 covers 90% more sweep material and is still
    # clearly below mean+outliers (max), which are only singular-band
    # frames.
    agg_n_max    = max(n_q95, n_min + 4)

    rationale_parts = [
        f"target_count = min(n_cross) = {n_min} "
        f"(Saturn/orbit pattern: the optimal frame gives the minimal diagram).",
        f"frame_cap = target+2 = {frame_cap}.",
        f"aggressive_simplify Range: [{agg_n_min}, {agg_n_max}] "
        f"(min..q95={n_q95}; q75={n_q75} was too narrow — Saturn c=12 finding "
        f"2026-05-10: frames with 3% reduction often lie above q75).",
    ]
    if all_even:
        rationale_parts.append(
            "ALL frame n_cross even ⇒ pi-antisymmetry signature "
            "(every real crossing appears twice via antipodal pairs).")
    elif all_odd:
        rationale_parts.append(
            "ALL frame n_cross odd ⇒ unusual, check by hand.")

    result = {
        "n_frames":          int(g.count()),
        "n_crossings_total": int(g.sum()),
        "per_frame_n_cross": {
            "min": n_min, "max": n_max, "median": n_median, "mean": n_mean,
            "mode_name": n_modus, "q05": n_q05, "q25": n_q25, "q75": n_q75, "q95": n_q95,
        },
        "all_even_n_cross":  all_even,
        "all_odd_n_cross":   all_odd,
        "target_count":      target_count,
        "frame_cap":         frame_cap,
        "agg_n_min":         agg_n_min,
        "agg_n_max":         agg_n_max,
        "rationale":         " ".join(rationale_parts),
    }

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print(f"[DETECT_TARGET] Frames:      {result['n_frames']:,}")
        print(f"[DETECT_TARGET] Crossings:   {result['n_crossings_total']:,}")
        print(f"[DETECT_TARGET] n_cross per frame:")
        print(f"  min={n_min}, modus={n_modus}, median={n_median}, "
              f"mean={n_mean:.1f}, max={n_max}")
        print(f"  Quantiles: q05={n_q05}, q25={n_q25}, q75={n_q75}, q95={n_q95}")
        print(f"[DETECT_TARGET] Parity:      "
              f"{'all even' if all_even else 'all odd' if all_odd else 'mixed'}")
        print(f"[DETECT_TARGET] Recommendations:")
        print(f"  --target-count       = {target_count}")
        print(f"  --frame-cap          = {frame_cap}")
        print(f"  agg --n-min          = {agg_n_min}")
        print(f"  agg --n-max          = {agg_n_max}")


if __name__ == '__main__':
    main()
