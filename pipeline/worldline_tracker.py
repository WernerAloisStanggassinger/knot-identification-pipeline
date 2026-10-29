#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Part of the Knot-Identification Pipeline.
"""worldline_tracker.py — prototype (Werner 2026-05-27).

Diagnosis: the DBSCAN tracker averages over/under across all frames (vote=.mean())
and merges crossings via connected components → at high density the consensus is
≈0.5 and the resolve discards most frame crossings → virtual diagram.

Worldline approach: a real projection frame of a closed 3D knot
ALWAYS yields a classical diagram — as long as ALL its crossings are taken
(over/under from the z of THIS frame, no mean). The temporal accounting serves
(a) the frame CHOICE from the camera wreath and (b) the persistence confirmation
of which crossings are real 3D features (stable (ti,tj) worldline over many
frames) vs. camera-specific tangential artefacts.

NO discard logic via a broken cluster match. Diagram = full frame.

Call: python worldline_tracker.py <crossings.parquet> [--name X] [--tol 0.01]
"""
import sys, json, argparse
import site_config
import numpy as np, pandas as pd

PIPE = site_config.PIPELINE_DIR
sys.path.insert(0, PIPE)
from kinematic_resolve import build_signed_gauss, regina_check  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('crossings')
    ap.add_argument('--name', default='knot')
    ap.add_argument('--tol', type=float, default=0.01,
                    help='(ti,tj) tolerance for the persistence histogram')
    ap.add_argument('--persist-min', type=float, default=0.02,
                    help='min. frame fraction from which a crossing counts as a real '
                         'worldline (persistent)')
    ap.add_argument('--output', default=None)
    args = ap.parse_args()

    df = pd.read_parquet(args.crossings)
    df['absc'] = df['cross2d'].abs()
    nframes = df['t_step'].nunique()
    print(f'[WL] {args.name}: {len(df):,} crossings, {nframes:,} frames')

    # ── (a) Frame choice: min crossing count, then max median|cross2d| (transverse) ──
    g = df.groupby('t_step')
    stats = pd.DataFrame({'n': g.size(), 'medabs': g['absc'].median()}).reset_index()
    nmin = int(stats['n'].min())
    cand = stats[stats['n'] == nmin].sort_values('medabs', ascending=False)
    fstar = int(cand.iloc[0]['t_step'])
    frame = df[df['t_step'] == fstar].copy()
    print(f'[WL] min n_cross/frame = {nmin}  → chosen frame t_step={fstar} '
          f'(med|cross2d|={cand.iloc[0]["medabs"]:.3f})')

    # ── (b) Persistence histogram (worldline confirmation) ──
    G = max(2, int(round(1.0 / args.tol)))
    ti_all = df['t_norm_i'].to_numpy(); tj_all = df['t_norm_j'].to_numpy()
    hist = np.zeros((G, G), dtype=np.int64)
    np.add.at(hist,
              (np.clip((ti_all * G).astype(int), 0, G - 1),
               np.clip((tj_all * G).astype(int), 0, G - 1)), 1)

    frame_obs, persist = [], []
    for _, c in frame.iterrows():
        cx = min(int(c['t_norm_i'] * G), G - 1)
        cy = min(int(c['t_norm_j'] * G), G - 1)
        # sum the 3x3 neighbourhood (a worldline may drift slightly)
        x0, x1 = max(0, cx-1), min(G, cx+2); y0, y1 = max(0, cy-1), min(G, cy+2)
        p = hist[x0:x1, y0:y1].sum() / nframes
        persist.append(p)
        frame_obs.append({'cluster_id': len(frame_obs) + 1,
                          'ti': float(c['t_norm_i']), 'tj': float(c['t_norm_j']),
                          'z_i': float(c['z_i']), 'z_j': float(c['z_j']),
                          'cross2d': float(c['cross2d'])})
    persist = np.array(persist)
    n_persist = int((persist >= args.persist_min).sum())
    print(f'[WL] frame crossings: {len(frame_obs)}  | persistent (≥{args.persist_min} '
          f'of the frames): {n_persist}  | persistence min/med/max = '
          f'{persist.min():.3f}/{np.median(persist):.3f}/{persist.max():.3f}')

    # ── Diagram from the FULL frame (no discarding) ──
    sg, tokens, _ = build_signed_gauss(frame_obs)
    print(f'[WL] signed Gauss ({len(frame_obs)} crossings, full frame):')
    print('     ' + (sg[:140] + ('…' if len(sg) > 140 else '')))

    r = regina_check(sg)
    rep = {k: r[k] for k in r if k != 'signed_gauss'}
    print('[WL] regina:', json.dumps(rep))

    out = {'name': args.name, 'frame': fstar, 'n_cross_frame': len(frame_obs),
           'n_persistent': n_persist, 'persist_min': args.persist_min,
           'signed_gauss': sg, 'regina': rep}
    if args.output:
        json.dump(out, open(args.output, 'w'), indent=1)
        print(f'[WL] → {args.output}')


if __name__ == '__main__':
    main()
