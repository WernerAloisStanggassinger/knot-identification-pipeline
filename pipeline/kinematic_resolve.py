#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
kinematic_resolve.py
====================
crossings.parquet + tracks_dbscan.parquet  →  diagramm.json (+ Regina verification)

Pipeline B closure: from the DBSCAN cluster trajectories (camera-independently
identified topological crossing classes) and a frame choice in which
ALL clusters are transversally visible, the signed Gauss code of the
diagram is built and handed to Regina.

Concept
-------
1. DBSCAN identifies the topological classes in (ti, tj) space
   (camera-independent). Per class: one representative (ti_repr, tj_repr).
2. A consensus vote on over/under and sign is automatically 50/50 FOR
   π-symmetric knots (see the Discovery #69-#72 block — antipodal hemispheres
   see the same topological crossing with mirrored Z order).
   → the vote consensus does NOT give the topological annotation.
3. Instead choose a representative FRAME: a camera angle in
   which all clusters are transversally visible (and few additional
   pseudo-crossings appear). In this frame there is a consistent
   "from above" view and hence a unique over/under annotation per
   visible cluster class.
4. The signed Gauss code is built from the frame observations by
   sorting the tokens along the arc parameter.

Frame choice
------------
Heuristic:
  • n_crossings ≈ number of clusters (few pseudo-crossings)
  • min sin_angle high (all crossings transversal)
Ranking: first by n_crossings (ascending), on a tie
by min sin_angle (descending). --frame-cap caps n_crossings from above.

Antisymmetry
------------
With --fold-pi-antisym: the (folded) cluster set has already resolved the
antisymmetry; in the frame the N crossings correspond to the N (folded)
clusters via lookup on frame (ti, tj) modulo the antisymmetry action.
Without fold: the unfolded cluster list contains antipode pairs;
in the frame the N crossings are distributed over 2N folded representatives.
"""

import argparse
import site_config
import json
import math
import os
import subprocess
import sys
import time
import numpy as np
import pandas as pd

sys.path.insert(0, site_config.FREEZE_DIR)
from pipeline_seal import seal_all
seal_all(verbose=False)


def parse_fseries(path):
    coeffs = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('%'):
                continue
            vals = [float(v) for v in line.split()]
            coeffs.append(vals)
    if len(coeffs[0]) == 3:
        offset = np.array(coeffs[0])
        harmonics = coeffs[1:]
    else:
        offset = np.zeros(3)
        harmonics = coeffs
    a = np.array(harmonics)
    return offset, a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4], a[:, 5]


def reconstruct_curve(offset, ax, bx, ay, by, az, bz, N, period_pi=False):
    period = math.pi if period_pi else 2 * math.pi
    s = np.linspace(0, period, N, endpoint=False)
    js = np.arange(1, len(ax) + 1)
    angles = s[:, None] * js[None, :]
    ca = np.cos(angles)
    sa = np.sin(angles)
    x = ca @ ax + sa @ bx + offset[0]
    y = ca @ ay + sa @ by + offset[1]
    z = ca @ az + sa @ bz + offset[2]
    return np.stack([x, y, z], axis=1).astype(np.float32)


def fold_pi_antisym_pt(ti, tj):
    """As in the tracker: ti < 0.5 as fundamental domain, tj free in [0, 1)."""
    if ti >= 0.5:
        return (ti - 0.5, (tj + 0.5) % 1.0)
    return (ti, tj)


def torus_dist(a, b):
    """Distance on the [0,1) torus."""
    d = abs(a - b)
    return min(d, 1.0 - d)


def match_cluster(ti_obs, tj_obs, cluster_reps, fold_pi):
    """Find nearest cluster_id for an observed (ti_obs, tj_obs).
    cluster_reps: DataFrame with cluster_id, ti_repr, tj_repr.
    fold_pi: if True, compare modulo the antisymmetry action.
    Returns: (best_cluster_id, distance, used_antipode_form)."""
    if fold_pi:
        ti_obs_f, tj_obs_f = fold_pi_antisym_pt(ti_obs, tj_obs)
        d_ti = np.array([torus_dist(ti_obs_f, ti_r) for ti_r in cluster_reps['ti_repr']])
        d_tj = np.array([torus_dist(tj_obs_f, tj_r) for tj_r in cluster_reps['tj_repr']])
    else:
        d_ti_o = np.array([torus_dist(ti_obs, ti_r) for ti_r in cluster_reps['ti_repr']])
        d_tj_o = np.array([torus_dist(tj_obs, tj_r) for tj_r in cluster_reps['tj_repr']])
        ti_a = (ti_obs + 0.5) % 1.0
        tj_a = (tj_obs + 0.5) % 1.0
        d_ti_a = np.array([torus_dist(ti_a, ti_r) for ti_r in cluster_reps['ti_repr']])
        d_tj_a = np.array([torus_dist(tj_a, tj_r) for tj_r in cluster_reps['tj_repr']])
        d_ti = np.minimum(d_ti_o, d_ti_a)
        d_tj = np.minimum(d_tj_o, d_tj_a)
    d = np.maximum(d_ti, d_tj)  # Chebyshev distance on the torus
    best = int(np.argmin(d))
    return int(cluster_reps.iloc[best]['cluster_id']), float(d[best])


def select_frame(df_cross, sin_angle, target_count, frame_cap):
    """Choose a frame with a crossing count close to target_count and max min(sin_angle).
    Returns: list of t_steps sorted by suitability (best first)."""
    df = df_cross.copy()
    df['sin_angle'] = sin_angle
    grp = df.groupby('t_step')
    stats = pd.DataFrame({
        'n_cross': grp.size(),
        'min_sin': grp['sin_angle'].min(),
        'med_sin': grp['sin_angle'].median(),
    }).reset_index()

    if frame_cap is not None:
        stats = stats[stats['n_cross'] <= frame_cap].copy()

    # Distance to the target, then max min_sin
    stats['n_dist'] = (stats['n_cross'] - target_count).abs()
    stats = stats.sort_values(['n_dist', 'min_sin'],
                              ascending=[True, False]).reset_index(drop=True)
    return stats


def build_signed_gauss(frame_obs):
    """Build signed Gauss code from frame observations.
    frame_obs: list of dicts with keys cluster_id, ti, tj, z_i, z_j, cross2d.
    Every frame crossing is its own diagram crossing (even if several
    belong to the same cluster class — mega-clusters absorb several genuine
    3D crossings that are distinct in the frame). cluster_id is carried along as an
    annotation per crossing, NOT as the diagram index.

    Tokens per crossing: 'O' (over) for the strand with higher z, 'U' for the
    other. Sign: + if cross2d > 0, − otherwise.
    Returns: signed Gauss string, list of tokens, mapping diag_idx → cluster_id."""
    tokens = []
    diag_to_cluster = {}
    for diag_idx, c in enumerate(frame_obs, start=1):
        diag_to_cluster[diag_idx] = c['cluster_id']
        i_over = c['z_i'] > c['z_j']
        # Crossing sign (standard knot-theory convention):
        #   + ↔ (T_O × T_U) · n > 0 (right-handed; viewer on the +n side,
        #   over = larger z). FIXED 2026-06-05 (mirror finding, error log
        #   ~08:45): the old formula here was (T_U × T_O) — exactly reversed,
        #   proven convention-free via the Gauss linking integral
        #   (Hopf link Lk=-1: det[T_O,T_U]→-1 ✓, det[T_U,T_O]→+1 ✗) plus
        #   regina calibration (trefoilRight = all '+'). Consequence of the old formula:
        #   all chain diagrams were globally mirrored relative to the true
        #   .fseries curve (corpus healed by the mirror sweep 2026-06-05).
        # Scanner cross2d = (T_i × T_j) · n. If i_over, T_O=T_i,
        # T_U=T_j → (T_O × T_U) · n = +cross2d. If j_over, T_O=T_j,
        # T_U=T_i → (T_O × T_U) · n = -cross2d.
        signed_cross = c['cross2d'] if i_over else -c['cross2d']
        sign = '+' if signed_cross > 0 else '-'
        tokens.append({
            'pos':  c['ti'],
            'crossing_idx': diag_idx,
            'overpass': i_over,
            'sign': sign,
            'cluster_id': c['cluster_id'],
            'role': 'i',
        })
        tokens.append({
            'pos':  c['tj'],
            'crossing_idx': diag_idx,
            'overpass': not i_over,
            'sign': sign,
            'cluster_id': c['cluster_id'],
            'role': 'j',
        })
    tokens.sort(key=lambda t: t['pos'])
    parts = []
    for t in tokens:
        ch = 'O' if t['overpass'] else 'U'
        parts.append(f"{ch}{t['crossing_idx']}{t['sign']}")
    return ''.join(parts), tokens, diag_to_cluster


# --- High-crossing scaling (Werner 2026-05-27; double-6 ~116 crossings) ---
# HOMFLY/Jones are exponential in the crossing number -> only attempt them on SMALL
# simplified diagrams; timeout high so that simplify() on a
# 100+-crossing knot has room. Reversible: .bak-pre-highcross-*.
REGINA_TIMEOUT = 3600   # was 600 s
REGINA_POLY_MAX = 35    # HOMFLY/Jones only if simplified size <= this


def regina_check(signed_gauss):
    """Spawn a base-env subprocess that imports regina and parses the signed Gauss
    code. Returns dict with keys: parsed, classical, virtualGenus, simplified_size,
    homfly, jones."""
    code = (
        "import json, sys\n"
        "sys.path.insert(0, '@FREEZE_DIR@')\n"
        "from pipeline_seal import seal_all\n"
        "seal_all(verbose=False)\n"
        "import regina\n"
        "sg = ".replace('@FREEZE_DIR@', site_config.FREEZE_DIR) + repr(signed_gauss) + "\n"
        "out = {'signed_gauss': sg}\n"
        "try:\n"
        "    L = regina.Link.fromSignedGauss(sg)\n"
        "    out['parsed'] = True\n"
        "    out['classical'] = bool(L.isClassical())\n"
        "    out['virtualGenus'] = int(L.virtualGenus())\n"
        "    if out['classical']:\n"
        "        L.simplify()\n"
        "        out['simplified_size'] = int(L.size())\n"
        f"        if L.size() <= {REGINA_POLY_MAX}:\n"
        "            try: out['homfly'] = str(L.homfly(regina.Algorithm.Treewidth))\n"
        "            except Exception as e: out['homfly_err'] = str(e)\n"
        "            try: out['jones']  = str(L.jones(regina.Algorithm.Treewidth))\n"
        "            except Exception as e: out['jones_err']  = str(e)\n"
        f"        else: out['poly_skipped'] = 'simplified_size=' + str(L.size()) + ' > {REGINA_POLY_MAX}: HOMFLY/Jones skipped (exponential)'\n"
        "except Exception as e:\n"
        "    out['parsed'] = False\n"
        "    out['error']  = str(e)\n"
        "print(json.dumps(out))\n"
    )
    res = subprocess.run(
        ['conda', 'run', '-n', 'base', '--no-capture-output', 'python', '-c', code],
        capture_output=True, text=True, timeout=REGINA_TIMEOUT,
    )
    if res.returncode != 0:
        return {'parsed': False, 'error': f'subprocess: rc={res.returncode}',
                'stderr': res.stderr}
    # Parse last JSON line
    for line in reversed(res.stdout.strip().splitlines()):
        line = line.strip()
        if line.startswith('{') and line.endswith('}'):
            return json.loads(line)
    return {'parsed': False, 'error': 'no json output', 'stdout': res.stdout,
            'stderr': res.stderr}


def main():
    parser = argparse.ArgumentParser(
        description='Pipeline B closure: build the signed Gauss code from DBSCAN clusters + '
                    'frame choice, and cross-check with Regina.')
    parser.add_argument('crossings', help='crossings.parquet (from knot_scanner.py)')
    parser.add_argument('--tracks',  required=True,
                        help='tracks_dbscan.parquet (from knot_tracker_dbscan.py)')
    parser.add_argument('--fseries', required=True,
                        help='Path to the .fseries (for sin angle + curve reconstruction)')
    parser.add_argument('--camera',  required=True,
                        help='camera_*.npz (for theta/phi lookup)')
    parser.add_argument('--N',       type=int, default=1024,
                        help='Curve sampling density (must match the scanner)')
    parser.add_argument('--period-pi', action='store_true',
                        help='Sample the curve on [0,π) (π-PERIODIC knots only)')
    parser.add_argument('--fold-pi-antisym', action='store_true',
                        help='Clusters are in the fold-pi domain. Frame crossings '
                             'are matched modulo the antisymmetry action.')
    parser.add_argument('--target-count', type=int, default=None,
                        help='Expected crossing count in the ideal frame (default: '
                             'number of clusters in tracks). The frame choice optimises '
                             'around this size.')
    parser.add_argument('--frame-cap', type=int, default=None,
                        help='Hard cap on n_crossings in the frame (default: '
                             '2 × target-count).')
    parser.add_argument('--top-frames', type=int, default=10,
                        help='How many top frames to try (default 10). '
                             'The first with a classical diagram wins.')
    parser.add_argument('--match-tolerance', type=float, default=0.02,
                        help='Maximum (ti, tj) distance for a cluster match '
                             '(Chebyshev, default 0.02). Crossings without a match '
                             '= pseudo-crossings, filtered.')
    parser.add_argument('--output', default='kinematic_result.json',
                        help='Output JSON path')
    parser.add_argument('--no-regina', action='store_true',
                        help='Skip the Regina verification (Gauss code construction only)')
    args = parser.parse_args()

    t_total = time.time()

    print(f"[KINEMATIC] Loading tracks: {args.tracks}")
    tracks = pd.read_parquet(args.tracks)
    n_clusters = len(tracks)
    print(f"[KINEMATIC]   {n_clusters} clusters")

    print(f"[KINEMATIC] Loading crossings: {args.crossings}")
    df = pd.read_parquet(args.crossings)
    print(f"[KINEMATIC]   {len(df):,} crossings, {df['t_step'].nunique():,} frames")

    print(f"[KINEMATIC] Loading fseries + curve (N={args.N})")
    coeffs = parse_fseries(args.fseries)
    curve = reconstruct_curve(*coeffs, N=args.N, period_pi=args.period_pi)
    d = np.roll(curve, -1, axis=0) - curve
    dmag = np.linalg.norm(d, axis=1).astype(np.float32)

    seg_i = df['seg_i'].to_numpy()
    seg_j = df['seg_j'].to_numpy()
    c2d   = df['cross2d'].to_numpy()
    sin_a = np.abs(c2d) / (dmag[seg_i].astype(np.float64) *
                            dmag[seg_j].astype(np.float64) + 1e-30)

    target = args.target_count if args.target_count else n_clusters
    cap    = args.frame_cap    if args.frame_cap    else 2 * target
    print(f"[KINEMATIC] Frame choice: target n_cross ≈ {target}, hard-cap {cap}")

    stats = select_frame(df, sin_a, target_count=target, frame_cap=cap)
    print(f"[KINEMATIC]   Top-{min(args.top_frames, len(stats))} frame candidates:")
    print(stats.head(args.top_frames).to_string(index=False))

    cam = np.load(args.camera)
    theta_all = cam['theta']
    phi_all   = cam['phi']

    cluster_reps = tracks[['cluster_id', 'ti_repr', 'tj_repr']].copy()

    chosen = None
    chosen_result = None
    for rank, row in stats.head(args.top_frames).iterrows():
        t_step = int(row['t_step'])
        n_cross = int(row['n_cross'])
        min_sin = float(row['min_sin'])

        frame_df = df[df['t_step'] == t_step].copy()
        frame_df['sin_angle'] = sin_a[df['t_step'].to_numpy() == t_step]

        # Match each frame crossing to the nearest cluster
        matches = []
        for _, fc in frame_df.iterrows():
            cid, dist = match_cluster(float(fc['t_norm_i']), float(fc['t_norm_j']),
                                      cluster_reps, args.fold_pi_antisym)
            matches.append({
                'cluster_id': cid,
                'distance':   dist,
                'ti':         float(fc['t_norm_i']),
                'tj':         float(fc['t_norm_j']),
                'z_i':        float(fc['z_i']),
                'z_j':        float(fc['z_j']),
                'cross2d':    float(fc['cross2d']),
                'sin_angle':  float(fc['sin_angle']),
                'seg_i':      int(fc['seg_i']),
                'seg_j':      int(fc['seg_j']),
            })

        kept = [m for m in matches if m['distance'] <= args.match_tolerance]
        dropped = len(matches) - len(kept)

        # NO cluster dedup: every frame crossing is a distinct diagram crossing.
        # cluster_id serves only for phantom filtering (distance > tolerance).
        signed_gauss, tokens, diag_to_cluster = build_signed_gauss(kept)
        n_unique_clusters = len({m['cluster_id'] for m in kept})

        print(f"\n[KINEMATIC] [{rank+1}] t_step={t_step}  "
              f"theta={theta_all[t_step]:.4f}  phi={phi_all[t_step]:.4f}  "
              f"n={n_cross} min_sin={min_sin:.3f}")
        print(f"[KINEMATIC]   Matches: {len(kept)}/{len(matches)} "
              f"({dropped} pseudo-crossings discarded)")
        print(f"[KINEMATIC]   Frame crossings in the diagram: {len(kept)}, "
              f"of which {n_unique_clusters} distinct clusters (of {n_clusters})")
        print(f"[KINEMATIC]   Signed Gauss ({len(kept)} crossings):")
        print(f"     {signed_gauss}")

        if args.no_regina:
            chosen = {
                't_step': t_step, 'theta': float(theta_all[t_step]),
                'phi': float(phi_all[t_step]),
                'rank': rank+1, 'n_cross_in_frame': n_cross,
                'matches_kept': len(kept), 'matches_dropped': dropped,
                'distinct_clusters_in_frame': n_unique_clusters,
                'n_clusters_total': int(n_clusters),
                'signed_gauss': signed_gauss,
                'diagram_idx_to_cluster': {str(k): v for k, v in diag_to_cluster.items()},
            }
            chosen_result = {'regina_skipped': True}
            break

        print(f"[KINEMATIC]   Regina verification …")
        regina_out = regina_check(signed_gauss)
        if regina_out.get('parsed') and regina_out.get('classical'):
            print(f"[KINEMATIC]   ✓ CLASSICAL, simplified_size="
                  f"{regina_out.get('simplified_size')}")
            print(f"[KINEMATIC]     HOMFLY = {regina_out.get('homfly', '?')}")
            print(f"[KINEMATIC]     Jones  = {regina_out.get('jones', '?')}")
            chosen = {
                't_step': t_step, 'theta': float(theta_all[t_step]),
                'phi': float(phi_all[t_step]),
                'rank': rank+1, 'n_cross_in_frame': n_cross,
                'matches_kept': len(kept), 'matches_dropped': dropped,
                'distinct_clusters_in_frame': n_unique_clusters,
                'n_clusters_total': int(n_clusters),
                'signed_gauss': signed_gauss,
                'diagram_idx_to_cluster': {str(k): v for k, v in diag_to_cluster.items()},
            }
            chosen_result = regina_out
            break
        else:
            if regina_out.get('parsed'):
                print(f"[KINEMATIC]   ✗ virtual  virtualGenus="
                      f"{regina_out.get('virtualGenus')}")
            else:
                print(f"[KINEMATIC]   ✗ parser failed: "
                      f"{regina_out.get('error', 'unknown')}")

    if chosen is None:
        print(f"\n[KINEMATIC] NO top-{args.top_frames} frame yields a classical diagram. "
              f"Output contains only the last attempt.", file=sys.stderr)
        chosen = {
            't_step': int(stats.iloc[args.top_frames-1]['t_step']),
            'rank': args.top_frames, 'no_classical_found': True,
            'signed_gauss': signed_gauss,
            'diagram_idx_to_cluster': {str(k): v for k, v in diag_to_cluster.items()},
        }
        chosen_result = regina_out if not args.no_regina else {'regina_skipped': True}

    output = {
        'crossings_path':   args.crossings,
        'tracks_path':      args.tracks,
        'fseries_path':     args.fseries,
        'n_clusters_dbscan': int(n_clusters),
        'fold_pi_antisym':  bool(args.fold_pi_antisym),
        'period_pi':        bool(args.period_pi),
        'chosen': chosen,
        'regina': chosen_result,
        'duration_seconds': time.time() - t_total,
    }
    with open(args.output, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\n[KINEMATIC] Writing {args.output}")
    print(f"[KINEMATIC] Done in {time.time()-t_total:.1f}s.")


if __name__ == '__main__':
    main()
