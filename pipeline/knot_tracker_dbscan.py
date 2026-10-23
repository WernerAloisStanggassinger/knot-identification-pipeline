#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
knot_tracker_dbscan.py
======================
crossings.parquet  →  tracks_dbscan.parquet

Cell-based connected-component clustering on the camera-independent
arc-parameter axes (t_norm_i, t_norm_j) — the Discovery #33 strategy as a
clean pipeline B tracker.

Core idea
---------
A genuine 3D crossing configuration (local minimum of the curve's self-distance)
has a fixed (ti, tj) point in arc-parameter space. It is visible as a 2D
crossing in one hemisphere of the Fibonacci cameras — the
individual observations across the cameras form a compact cluster
in (ti, tj) space.

Tangencies / camera-specific pseudo-crossings (two segments that are not
close in 3D but happen to overlap from this one perspective)
are spread over the (ti, tj) plane — they form no compact
cluster and are weeded out by the min-cell-count filter.

Anti-pattern avoided here
-------------------------
- 5D clusters over (sx, sy, sz, ti, tj): Discovery #74 — mega-clusters or
  over-fragmentation depending on radii, due to (sx,sy,sz) trajectory drift.
- ti/tj averaging over cluster frames: Discovery #38 — virtualises the
  diagram.

Here (ti, tj) is the cluster identifier (camera-independent), and the
consensus per cluster is formed over all frames, not by averaging.

Options
-------
--fold-pi-antisym : Before clustering, fold (ti, tj) onto the fundamental domain
                    of the π-antisymmetry (treat (s, t) and (s+π, t+π) as
                    equivalent). Clean for Saturn-like knots,
                    assumes antisymmetry.

Output (tracks_dbscan.parquet, one row per cluster, sorted by n_obs):
  cluster_id           int
  ti_repr, tj_repr     float (representative crossing from the max-|cross2d| frame)
  ti_centroid,         float (mean of ti/tj in the cluster — diagnostic only,
   tj_centroid                NOT for Gauss code construction!)
  n_obs                int   (number of observations in the cluster)
  n_unique_t_steps     int
  vote_over            float (fraction of frames with z_i > z_j)
  vote_sign            float (fraction of frames with cross2d > 0)
  rep_t_step           int   (frame with max |cross2d|)
  rep_seg_i, rep_seg_j int
  rep_z_i, rep_z_j     float
  rep_cross2d          float
  fold_partner         int   (cluster_id of the antipodal partner; -1 if not
                              found / not applied in fold mode)
"""

import sys
import site_config
import os
import math
import time
import argparse
import numpy as np
import pandas as pd
import cupy as cp

sys.path.insert(0, site_config.FREEZE_DIR)
from pipeline_seal import seal_all
seal_all(verbose=False)


CELL_COUNT_KERNEL = r"""
extern "C" __global__
void cell_assign(
    const float* __restrict__ ti,
    const float* __restrict__ tj,
    int*         cell_x_out,
    int*         cell_y_out,
    int*         cell_id_out,
    int*         counts,
    const int    n_grid,
    const float  cell_size,
    const int    N
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;
    float t1 = ti[idx];
    float t2 = tj[idx];
    int cx = (int)floorf(t1 / cell_size);
    int cy = (int)floorf(t2 / cell_size);
    if (cx < 0) cx = 0; if (cx >= n_grid) cx = n_grid - 1;
    if (cy < 0) cy = 0; if (cy >= n_grid) cy = n_grid - 1;
    int cid = cx * n_grid + cy;
    cell_x_out[idx] = cx;
    cell_y_out[idx] = cy;
    cell_id_out[idx] = cid;
    atomicAdd(&counts[cid], 1);
}
"""


def fold_pi_antisym(ti, tj):
    """Map (ti, tj) into the fundamental domain of the π-antisymmetry
    (X(s+π) = -X(s)). Each topological crossing has an antipode at
    ((ti+0.5) mod 1, (tj+0.5) mod 1).

    Fundamental domain: {(ti, tj) : 0 ≤ ti < 0.5}.
    In the scanner output ti < tj holds (seg_i < seg_j); this order
    is preserved under the fold — hence no swap needed."""
    flip = ti >= 0.5
    ti_a = np.where(flip, ti - 0.5, ti)
    tj_a = np.where(flip, (tj + 0.5) % 1.0, tj)
    return ti_a.astype(np.float32), tj_a.astype(np.float32)


def cell_connected_components(cells_active, n_grid, periodic=True):
    """8-neighbour connected components on an n_grid × n_grid cell image.
    cells_active: bool array, shape (n_grid, n_grid).
    periodic: wraparound on both axes ([0,1) arc-parameter periodicity).
    Returns: (label_grid: int32 array shape (n_grid, n_grid), n_components).
    Cells with cells_active=False get label = -1."""
    labels = np.full((n_grid, n_grid), -1, dtype=np.int32)
    next_label = 0
    stack = []
    if periodic:
        nbrs = [(-1, -1), (-1, 0), (-1, 1),
                ( 0, -1),          ( 0, 1),
                ( 1, -1), ( 1, 0), ( 1, 1)]
    else:
        nbrs = [(-1, -1), (-1, 0), (-1, 1),
                ( 0, -1),          ( 0, 1),
                ( 1, -1), ( 1, 0), ( 1, 1)]
    for x in range(n_grid):
        for y in range(n_grid):
            if not cells_active[x, y] or labels[x, y] != -1:
                continue
            stack.append((x, y))
            labels[x, y] = next_label
            while stack:
                cx, cy = stack.pop()
                for dx, dy in nbrs:
                    nx = cx + dx
                    ny = cy + dy
                    if periodic:
                        nx %= n_grid
                        ny %= n_grid
                    elif nx < 0 or nx >= n_grid or ny < 0 or ny >= n_grid:
                        continue
                    if cells_active[nx, ny] and labels[nx, ny] == -1:
                        labels[nx, ny] = next_label
                        stack.append((nx, ny))
            next_label += 1
    return labels, next_label


def find_pi_antipodal_partners(reps, eps_repr=0.01):
    """For each cluster, search the antipodal partner ((ti+0.5)%1, (tj+0.5)%1).
    Match tolerance eps_repr on both arc-parameter axes.
    Returns: int array with the fold_partner cluster_id per cluster (-1 if none)."""
    n = len(reps)
    ti = reps['ti_repr'].to_numpy()
    tj = reps['tj_repr'].to_numpy()
    ti_anti = (ti + 0.5) % 1.0
    tj_anti = (tj + 0.5) % 1.0
    swap = ti_anti > tj_anti
    ti_anti_s = np.where(swap, tj_anti, ti_anti)
    tj_anti_s = np.where(swap, ti_anti, tj_anti)
    swap_self = ti > tj
    ti_self = np.where(swap_self, tj, ti)
    tj_self = np.where(swap_self, ti, tj)
    partners = np.full(n, -1, dtype=np.int32)
    for i in range(n):
        d_ti = np.abs(ti_self - ti_anti_s[i])
        d_tj = np.abs(tj_self - tj_anti_s[i])
        d_ti = np.minimum(d_ti, 1.0 - d_ti)
        d_tj = np.minimum(d_tj, 1.0 - d_tj)
        match = (d_ti < eps_repr) & (d_tj < eps_repr)
        cand = np.where(match)[0]
        cand = cand[cand != i]
        if len(cand) > 0:
            partners[i] = int(cand[0])
    return partners


def track(input_parquet, output_path, eps, min_cell_count, fold_pi,
          consensus_threshold, merge_pi_antipodal):

    t_total = time.time()

    print(f"[DBSCAN-TRACKER] Loading {input_parquet}")
    df = pd.read_parquet(input_parquet)
    N = len(df)
    n_frames_total = int(df['t_step'].max()) + 1
    print(f"[DBSCAN-TRACKER]   {N:,} crossings, {n_frames_total:,} frames in the input")

    ti = df['t_norm_i'].to_numpy().astype(np.float32)
    tj = df['t_norm_j'].to_numpy().astype(np.float32)

    if fold_pi:
        print(f"[DBSCAN-TRACKER] --fold-pi-antisym active: folding (ti,tj) onto "
              f"fundamental domain [0, 0.5) × [0, 1) (ti < 0.5)")
        ti_c, tj_c = fold_pi_antisym(ti, tj)
    else:
        ti_c = ti.copy()
        tj_c = tj.copy()

    n_grid = int(math.ceil(1.0 / eps))
    cell_size = 1.0 / n_grid
    print(f"[DBSCAN-TRACKER] Grid: {n_grid}×{n_grid} = {n_grid*n_grid:,} cells, "
          f"cell_size={cell_size:.5f}, eps={eps:.5f}")

    print(f"[DBSCAN-TRACKER] GPU cell assign + count …")
    free_b, total_b = cp.cuda.runtime.memGetInfo()
    print(f"[DBSCAN-TRACKER]   VRAM free: {free_b/1e9:.2f} GB / {total_b/1e9:.2f} GB")

    ti_gpu = cp.asarray(ti_c)
    tj_gpu = cp.asarray(tj_c)
    cell_x = cp.zeros(N, dtype=cp.int32)
    cell_y = cp.zeros(N, dtype=cp.int32)
    cell_id = cp.zeros(N, dtype=cp.int32)
    counts = cp.zeros(n_grid * n_grid, dtype=cp.int32)

    mod = cp.RawModule(code=CELL_COUNT_KERNEL)
    fn = mod.get_function('cell_assign')
    block = 256
    grid = (N + block - 1) // block
    t0 = time.time()
    fn((grid,), (block,),
       (ti_gpu, tj_gpu, cell_x, cell_y, cell_id, counts,
        np.int32(n_grid), np.float32(cell_size), np.int32(N)))
    cp.cuda.Stream.null.synchronize()
    print(f"[DBSCAN-TRACKER]   Cell assign in {time.time()-t0:.2f}s")

    counts_cpu = counts.get().reshape(n_grid, n_grid)
    cells_active = counts_cpu >= min_cell_count
    n_active = int(cells_active.sum())
    print(f"[DBSCAN-TRACKER]   Active cells (≥{min_cell_count} points): "
          f"{n_active} of {n_grid*n_grid}")

    if n_active == 0:
        print(f"[DBSCAN-TRACKER] No active cells. min_cell_count {min_cell_count} "
              f"too high or eps {eps} too small.", file=sys.stderr)
        sys.exit(1)

    print(f"[DBSCAN-TRACKER] Cell connected components (periodic) …")
    t0 = time.time()
    cell_labels, n_components = cell_connected_components(
        cells_active, n_grid, periodic=True)
    print(f"[DBSCAN-TRACKER]   {n_components} components in {time.time()-t0:.2f}s")

    print(f"[DBSCAN-TRACKER] Assigning points → clusters …")
    cx_cpu = cell_x.get()
    cy_cpu = cell_y.get()
    point_cluster = cell_labels[cx_cpu, cy_cpu]

    keep = point_cluster >= 0
    df_kept = df.loc[keep].copy()
    df_kept['cluster'] = point_cluster[keep]
    df_kept['ti_folded'] = ti_c[keep]
    df_kept['tj_folded'] = tj_c[keep]
    n_kept = len(df_kept)
    print(f"[DBSCAN-TRACKER]   {n_kept:,} of {N:,} points in a cluster "
          f"({100*n_kept/N:.1f}%)")

    return _aggregate_and_write(df_kept, output_path, eps, N,
                                    consensus_threshold, merge_pi_antipodal,
                                    fold_pi, t_total)



def _aggregate_and_write(df_kept, output_path, eps, N,
                             consensus_threshold, merge_pi_antipodal,
                             fold_pi, t_total):
    """Shared aggregation and write part for track() and
    track_streaming() (factored out 2026-08-24, Disc #703)."""
    # Per-cluster aggregation
    print(f"[DBSCAN-TRACKER] Per-cluster consensus voting …")
    t0 = time.time()

    df_kept['_over_i'] = (df_kept['z_i'] > df_kept['z_j']).astype(np.int32)
    df_kept['_pos_sign'] = (df_kept['cross2d'] > 0).astype(np.int32)
    df_kept['_abs_c2d'] = np.abs(df_kept['cross2d'].to_numpy())

    gb = df_kept.groupby('cluster')
    rep_idx = gb['_abs_c2d'].idxmax()
    reps = df_kept.loc[rep_idx, ['cluster', 't_step', 'seg_i', 'seg_j',
                                 't_norm_i', 't_norm_j',
                                 'z_i', 'z_j', 'cross2d',
                                 'ti_folded', 'tj_folded']].copy()
    reps = reps.rename(columns={
        't_step': 'rep_t_step',
        'seg_i':  'rep_seg_i',
        'seg_j':  'rep_seg_j',
        't_norm_i': 'ti_repr',
        't_norm_j': 'tj_repr',
        'z_i': 'rep_z_i',
        'z_j': 'rep_z_j',
        'cross2d': 'rep_cross2d',
    }).set_index('cluster')

    stats = pd.DataFrame({
        'n_obs':            gb.size(),
        'n_unique_t_steps': gb['t_step'].nunique(),
        'vote_over':        gb['_over_i'].mean(),
        'vote_sign':        gb['_pos_sign'].mean(),
        'ti_centroid':      gb['ti_folded'].mean(),
        'tj_centroid':      gb['tj_folded'].mean(),
    })

    return _write_from_aggregate(reps, stats, output_path, eps, N,
                                 consensus_threshold, merge_pi_antipodal,
                                 fold_pi, t_total)


def _write_from_aggregate(reps, stats, output_path, eps, N,
                           consensus_threshold, merge_pi_antipodal,
                           fold_pi, t_total):
    """Shared write part (factored out 2026-08-24, Disc #703)."""
    import time as _t
    t0 = _t.time()
    out = reps.join(stats).reset_index().rename(columns={'cluster': 'cluster_raw'})
    out = out.sort_values('n_obs', ascending=False).reset_index(drop=True)
    out['cluster_id'] = np.arange(len(out), dtype=np.int32)
    print(f"[DBSCAN-TRACKER]   Aggregation in {time.time()-t0:.2f}s")

    # ─── π-antipodal merge (post-DBSCAN) ───
    if merge_pi_antipodal:
        if fold_pi:
            print(f"[DBSCAN-TRACKER] --merge-pi-antipodal combined with --fold-pi-antisym: "
                  f"--fold has already resolved the antisymmetry, no merge needed.")
            out['fold_partner'] = -1
        else:
            print(f"[DBSCAN-TRACKER] Searching π-antipodal cluster pairs …")
            partners = find_pi_antipodal_partners(out, eps_repr=max(eps*4, 0.01))
            out['fold_partner'] = partners
            n_pair = int((partners >= 0).sum())
            n_lone = int((partners < 0).sum())
            print(f"[DBSCAN-TRACKER]   {n_pair} clusters with antipodes, {n_lone} without")
    else:
        out['fold_partner'] = -1

    # Reorder columns
    out_cols = ['cluster_id',
                'ti_repr', 'tj_repr', 'ti_centroid', 'tj_centroid',
                'n_obs', 'n_unique_t_steps',
                'vote_over', 'vote_sign',
                'rep_t_step', 'rep_seg_i', 'rep_seg_j',
                'rep_z_i', 'rep_z_j', 'rep_cross2d',
                'fold_partner']
    out_df = out[out_cols].copy()

    print(f"[DBSCAN-TRACKER] Writing {output_path}")
    if output_path.endswith('.csv'):
        out_df.to_csv(output_path, index=False)
    else:
        out_df.to_parquet(output_path, compression='snappy', index=False)

    sz = os.path.getsize(output_path) if os.path.exists(output_path) else 0
    print(f"[DBSCAN-TRACKER] Done in {time.time()-t_total:.1f}s.  "
          f"Output: {output_path} ({sz/1e6:.2f} MB)")
    print()
    print(f"[DBSCAN-TRACKER] Consensus statistics (≥{consensus_threshold:.2f} = stable):")
    n = len(out_df)
    over_stable = ((out_df['vote_over'] >= consensus_threshold) |
                   (out_df['vote_over'] <= 1 - consensus_threshold)).sum()
    sign_stable = ((out_df['vote_sign'] >= consensus_threshold) |
                   (out_df['vote_sign'] <= 1 - consensus_threshold)).sum()
    print(f"  Over/under stable: {over_stable}/{n}")
    print(f"  Sign stable      : {sign_stable}/{n}")

    print()
    print(f"[DBSCAN-TRACKER] Top-{min(40, len(out_df))} Cluster:")
    pd.set_option('display.width', 220)
    pd.set_option('display.max_columns', 20)
    pd.set_option('display.float_format', lambda v: f'{v:.4f}')
    cols_show = ['cluster_id', 'ti_repr', 'tj_repr',
                 'n_obs', 'n_unique_t_steps',
                 'vote_over', 'vote_sign',
                 'fold_partner']
    print(out_df[cols_show].head(40).to_string(index=False))

    return out_df



# ──────────────────────────────────────────────────────────────────────
# STREAMING VARIANT (2026-08-24, Disc #703)
# Two passes over the file instead of one load. Memory depends
# on the batch size and on the number of KEPT points, no longer on N.
# ──────────────────────────────────────────────────────────────────────

def _cells(ti, tj, cell_size, n_grid):
    """Exactly the formula of the CUDA kernel cell_assign: floor(t/cell_size),
    then clamp. float32 so that the rounding is the same."""
    cx = np.floor(ti.astype(np.float32) / np.float32(cell_size)).astype(np.int32)
    cy = np.floor(tj.astype(np.float32) / np.float32(cell_size)).astype(np.int32)
    np.clip(cx, 0, n_grid - 1, out=cx)
    np.clip(cy, 0, n_grid - 1, out=cy)
    return cx, cy


def track_streaming(input_parquet, output_path, eps, min_cell_count, fold_pi,
                    consensus_threshold, merge_pi_antipodal, batch_rows=4_000_000):
    import pyarrow.parquet as _pq

    t_total = time.time()
    pf = _pq.ParquetFile(input_parquet)
    N = pf.metadata.num_rows
    n_grid = int(math.ceil(1.0 / eps))
    cell_size = 1.0 / n_grid
    print(f"[DBSCAN-TRACKER] STREAMING over {input_parquet}")
    print(f"[DBSCAN-TRACKER]   {N:,} crossings, {pf.num_row_groups} row groups, "
          f"batch {batch_rows:,}")
    print(f"[DBSCAN-TRACKER] Grid: {n_grid}x{n_grid} = {n_grid*n_grid:,} cells, "
          f"cell_size={cell_size:.5f}, eps={eps:.5f}")

    # ── Pass 1: count only. Memory = batch + grid. ──
    t0 = time.time()
    counts = np.zeros(n_grid * n_grid, dtype=np.int64)
    max_step = 0
    for b in pf.iter_batches(batch_size=batch_rows,
                             columns=['t_step', 't_norm_i', 't_norm_j']):
        ts = b.column('t_step').to_numpy(zero_copy_only=False)
        if len(ts):
            max_step = max(max_step, int(ts.max()))
        ti = b.column('t_norm_i').to_numpy(zero_copy_only=False).astype(np.float32)
        tj = b.column('t_norm_j').to_numpy(zero_copy_only=False).astype(np.float32)
        if fold_pi:
            ti, tj = fold_pi_antisym(ti, tj)
        cx, cy = _cells(ti, tj, cell_size, n_grid)
        counts += np.bincount(cx.astype(np.int64) * n_grid + cy,
                              minlength=n_grid * n_grid)
    n_frames_total = max_step + 1
    print(f"[DBSCAN-TRACKER]   Pass 1 (count) in {time.time()-t0:.2f}s, "
          f"{n_frames_total:,} frames")

    counts_cpu = counts.reshape(n_grid, n_grid)
    cells_active = counts_cpu >= min_cell_count
    n_active = int(cells_active.sum())
    print(f"[DBSCAN-TRACKER]   Active cells (>={min_cell_count} points): "
          f"{n_active} of {n_grid*n_grid}")
    if n_active == 0:
        print(f"[DBSCAN-TRACKER] No active cells. min_cell_count {min_cell_count} "
              f"too high or eps {eps} too small.", file=sys.stderr)
        sys.exit(1)

    print(f"[DBSCAN-TRACKER] Cell connected components (periodic) ...")
    t0 = time.time()
    cell_labels, n_components = cell_connected_components(
        cells_active, n_grid, periodic=True)
    print(f"[DBSCAN-TRACKER]   {n_components} components in {time.time()-t0:.2f}s")

    # ── Pass 2: AGGREGATE, do not collect (Disc #703). ──
    # From here on memory depends on the number of clusters, not on that of points.
    t0 = time.time()
    K = int(n_components)
    n_worte = (n_frames_total + 63) // 64
    accum = {
        'n': np.zeros(K, dtype=np.int64),
        'over': np.zeros(K, dtype=np.int64),
        'sign': np.zeros(K, dtype=np.int64),
        'sum_ti': np.zeros(K, dtype=np.float64),
        'sum_tj': np.zeros(K, dtype=np.float64),
        'maxabs': np.full(K, -1.0, dtype=np.float64),
        'frames': np.zeros((K, n_worte), dtype=np.uint64),
    }
    rep_fields = ['t_step', 'seg_i', 'seg_j', 't_norm_i', 't_norm_j',
                  'z_i', 'z_j', 'cross2d', 'ti_folded', 'tj_folded']
    rep = {f: [None] * K for f in rep_fields}
    n_kept = 0

    for b in pf.iter_batches(batch_size=batch_rows):
        d = b.to_pandas()
        ti = d['t_norm_i'].to_numpy().astype(np.float32)
        tj = d['t_norm_j'].to_numpy().astype(np.float32)
        if fold_pi:
            ti_f, tj_f = fold_pi_antisym(ti, tj)
        else:
            ti_f, tj_f = ti, tj
        cx, cy = _cells(ti_f, tj_f, cell_size, n_grid)
        cl = cell_labels[cx, cy]
        keep = cl >= 0
        if not keep.any():
            continue
        c = cl[keep].astype(np.int64)
        n_kept += int(c.size)
        zi = d['z_i'].to_numpy()[keep]
        zj = d['z_j'].to_numpy()[keep]
        c2 = d['cross2d'].to_numpy()[keep]
        ts = d['t_step'].to_numpy()[keep]
        tif = ti_f[keep]
        tjf = tj_f[keep]

        np.add.at(accum['n'], c, 1)
        np.add.at(accum['over'], c, (zi > zj).astype(np.int64))
        np.add.at(accum['sign'], c, (c2 > 0).astype(np.int64))
        np.add.at(accum['sum_ti'], c, tif.astype(np.float64))
        np.add.at(accum['sum_tj'], c, tjf.astype(np.float64))

        # Frames: first deduplicate per batch, then set bits — otherwise 42 million
        # scattered accesses instead of a few tens of thousands.
        key_s = np.unique(c * np.int64(n_frames_total) + ts.astype(np.int64))
        kc = (key_s // n_frames_total).astype(np.int64)
        kt = (key_s % n_frames_total).astype(np.int64)
        np.bitwise_or.at(accum['frames'], (kc, kt >> 6),
                         (np.uint64(1) << (kt & 63).astype(np.uint64)))

        # Representative: running maximum, STRICTLY greater (pandas idxmax takes the
        # first maximum — with >= one would get the last).
        aabs = np.abs(c2)
        order = np.argsort(c, kind='stable')
        cs, as_ = c[order], aabs[order]
        limits = np.flatnonzero(np.r_[True, cs[1:] != cs[:-1]])
        for g0, g1 in zip(limits, np.r_[limits[1:], cs.size]):
            k = int(cs[g0])
            loc_i = int(np.argmax(as_[g0:g1]))          # first maximum in the batch
            if as_[g0 + loc_i] > accum['maxabs'][k]:
                accum['maxabs'][k] = as_[g0 + loc_i]
                i = int(order[g0 + loc_i])
                rep['t_step'][k] = ts[i]
                rep['seg_i'][k] = d['seg_i'].to_numpy()[keep][i]
                rep['seg_j'][k] = d['seg_j'].to_numpy()[keep][i]
                rep['t_norm_i'][k] = d['t_norm_i'].to_numpy()[keep][i]
                rep['t_norm_j'][k] = d['t_norm_j'].to_numpy()[keep][i]
                rep['z_i'][k] = zi[i]
                rep['z_j'][k] = zj[i]
                rep['cross2d'][k] = c2[i]
                rep['ti_folded'][k] = tif[i]
                rep['tj_folded'][k] = tjf[i]

    print(f"[DBSCAN-TRACKER]   Pass 2 (aggregate) in {time.time()-t0:.2f}s")
    print(f"[DBSCAN-TRACKER]   {n_kept:,} of {N:,} points in a cluster "
          f"({100*n_kept/N:.1f}%)")
    if n_kept == 0:
        print("[DBSCAN-TRACKER] No point in an active cell.", file=sys.stderr)
        sys.exit(1)

    da = accum['n'] > 0
    idx = np.flatnonzero(da)
    bits = np.zeros(K, dtype=np.int64)
    for k in idx:
        w = accum['frames'][k]
        bits[k] = int(sum(int(x).bit_count() for x in w))
    reps = pd.DataFrame({f: [rep[f][k] for k in idx] for f in rep_fields})
    reps = reps.rename(columns={'t_step': 'rep_t_step', 'seg_i': 'rep_seg_i',
                                'seg_j': 'rep_seg_j', 't_norm_i': 'ti_repr',
                                't_norm_j': 'tj_repr', 'z_i': 'rep_z_i',
                                'z_j': 'rep_z_j', 'cross2d': 'rep_cross2d'})
    reps.index = pd.Index(idx, name='cluster')
    stats = pd.DataFrame({
        'n_obs':            accum['n'][idx],
        'n_unique_t_steps': bits[idx],
        'vote_over':        accum['over'][idx] / accum['n'][idx],
        'vote_sign':        accum['sign'][idx] / accum['n'][idx],
        # float32 BECAUSE pandas' mean() on a float32 column returns float32.
        # ⚠️ DELIBERATELY NOT BIT-IDENTICAL to the classic route (Disc #703): these
        # two columns differ by up to ONE float32 ULP (measured 5.96e-08),
        # because here the sum is taken per batch and pandas sums the column in one go.
        # Werner's decision 2026-08-24 "leave it and document it", admissible
        # BECAUSE checked: kinematic_resolve.py:338 reads only cluster_id/ti_repr/
        # tj_repr, and those are exactly equal. The centroids are reporting quantities.
        # Computing in float64 is more accurate, but then the output is no longer
        # character-identical with the old route — and the counter-case rightly rejects it.
        'ti_centroid':      (accum['sum_ti'][idx] / accum['n'][idx]).astype(np.float32),
        'tj_centroid':      (accum['sum_tj'][idx] / accum['n'][idx]).astype(np.float32),
    }, index=reps.index)

    return _write_from_aggregate(reps, stats, output_path, eps, N,
                                  consensus_threshold, merge_pi_antipodal,
                                  fold_pi, t_total)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='DBSCAN-like cell connected-component tracker on (t_norm_i, t_norm_j). '
                    'Discovery #33 strategy as the pipeline B tracker.')
    parser.add_argument('input',
                        help='Path to crossings.parquet (from knot_scanner.py)')
    parser.add_argument('--output', default='tracks_dbscan.parquet',
                        help='.parquet (default) or .csv')
    parser.add_argument('--eps', type=float, default=0.005,
                        help='Cell size / cluster radius in arc-parameter space '
                             '(default 0.005). Should be clearly smaller than '
                             'the minimal distance between two crossings.')
    parser.add_argument('--min-cell-count', type=int, default=100,
                        help='Cell activation threshold: number of points in a cell '
                             'for it to belong to a cluster (default 100). '
                             'Filters tangential / camera-specific pseudo-crossings.')
    parser.add_argument('--fold-pi-antisym', action='store_true',
                        help='Before clustering: fold (ti,tj) onto [0,0.5)² '
                             '(π-antisymmetry as an assumption). Clean for Saturn-like '
                             'knots; without the flag antisymmetry is verified as a finding.')
    parser.add_argument('--no-merge-pi-antipodal', action='store_true',
                        help='Without --fold-pi-antisym: by default, after DBSCAN '
                             'the antipodal partner is searched for every cluster. '
                             '--no-merge-pi-antipodal disables this diagnostic.')
    parser.add_argument('--consensus-threshold', type=float, default=0.95,
                        help='Threshold (default 0.95) from which a consensus vote counts as '
                             '"stable" (reporting only).')
    parser.add_argument('--klassisch', action='store_true',
                        help='old route: load everything at once (Disc #703). '
                             'Kept only for the counter-case — at 137 million '
                             'crossings it needs over 20 GB and gets killed.')
    parser.add_argument('--batch-rows', type=int, default=4_000_000,
                        help='Batch size for the streaming route (default 4 million).')
    args = parser.parse_args()

    if args.klassisch:
        track(input_parquet=args.input,
              output_path=args.output,
              eps=args.eps,
              min_cell_count=args.min_cell_count,
              fold_pi=args.fold_pi_antisym,
              consensus_threshold=args.consensus_threshold,
              merge_pi_antipodal=not args.no_merge_pi_antipodal)
    else:
        track_streaming(input_parquet=args.input,
                        output_path=args.output,
                        eps=args.eps,
                        min_cell_count=args.min_cell_count,
                        fold_pi=args.fold_pi_antisym,
                        consensus_threshold=args.consensus_threshold,
                        merge_pi_antipodal=not args.no_merge_pi_antipodal,
                        batch_rows=args.batch_rows)
