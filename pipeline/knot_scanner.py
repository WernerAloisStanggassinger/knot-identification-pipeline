#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
knot_scanner.py
===============
.fseries + camera_*.npz -> crossings.parquet  (or .npz / .csv)

Stages 1+2 of the pipeline: pure data acquisition with performance-tuned output.
No interpretation: no over/under, no sign, no clustering. Only raw observations
per camera point.

The output format is chosen by the file extension:
  *.parquet -> pyarrow.parquet, snappy-compressed (default, fast and small)
  *.npz     -> np.savez_compressed with one 1D array per column
  *.csv     -> human-readable, for small smoke tests (slow)

Architecture
------------
The camera flight is a separate artifact (see make_camera_trajectory.py) and is
not computed inline. For every camera point t the scanner applies the rotation
R(t) to the curve (p -> R(t) p), projects onto the xy-plane (i.e. it views the
fixed knot along the direction R(t)^T e_z) and collects all 2D self-intersections
of the projected polygon.

Output columns (sorted by t_step, then by seg_i, seg_j):
  t_step       int32   camera point number along the trajectory (= time stamp)
  seg_i, seg_j int32   segment indices of the sampled curve (0..N_curve-1)
  t_norm_i,    float32 sub-segment position of the intersection point:
   t_norm_j             = (i + t) / N_curve
  z_i, z_j     float32 z-buffer: depth of both strands at the intersection
                       (linearly interpolated on each segment)
  cross2d      float32 determinant of the two 2D segment directions at i and j
                       (prepares the sign; the sign itself is computed later)

NOT in the output (belongs to later pipeline stages):
  - over/under (follows from z_i > z_j)
  - sign +/-  (follows from cross2d and over/under)
  - cluster_id, frame representative, Gauss code, ... (later evaluation)

theta/phi/cam_x/y/z are NOT repeated in the crossings output; they are obtained
from camera_*.npz/csv through t_step (lookup join).

Performance
-----------
CSV writing was the dominant bottleneck in profiling (~82 % of the total time).
NPZ/Parquet output removes the per-row Python loop and stores the columns as
compact 1D arrays.

In the authors' pipeline this script runs under a "blind-test seal" that blocks
knot-table lookups during identification; the scanner itself never uses knot tables,
so the seal is omitted in this standalone release.
"""

import sys
import os
import math
import csv
import time
import argparse
import numpy as np
import cupy as cp


# ──────────────────────────────────────────────
# 1. FSERIES + CURVE
# ──────────────────────────────────────────────

def parse_fseries(path):
    coeffs = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('%'):
                continue
            vals = [float(v) for v in line.split()]
            coeffs.append(vals)
    # The first row is recognised by its length:
    #   3 values -> explicit 3D offset row, harmonics start at row 1
    #   6 values -> no offset row, every row is a Fourier harmonic j = 1, 2, ...
    if len(coeffs[0]) == 3:
        offset = np.array(coeffs[0])
        harmonics = coeffs[1:]
    else:
        offset = np.zeros(3)
        harmonics = coeffs
    ax = np.array([r[0] for r in harmonics])
    bx = np.array([r[1] for r in harmonics])
    ay = np.array([r[2] for r in harmonics])
    by = np.array([r[3] for r in harmonics])
    az = np.array([r[4] for r in harmonics])
    bz = np.array([r[5] for r in harmonics])
    return offset, ax, bx, ay, by, az, bz


def reconstruct_curve(offset, ax, bx, ay, by, az, bz, N, period_pi=False):
    period = math.pi if period_pi else 2 * math.pi
    s = np.linspace(0, period, N, endpoint=False)
    js = np.arange(1, len(ax) + 1)
    angles = s[:, None] * js[None, :]
    cos_a = np.cos(angles)
    sin_a = np.sin(angles)
    x = cos_a @ ax + sin_a @ bx + offset[0]
    y = cos_a @ ay + sin_a @ by + offset[1]
    z = cos_a @ az + sin_a @ bz + offset[2]
    return np.stack([x, y, z], axis=1).astype(np.float32)


# ──────────────────────────────────────────────
# 2. CUDA KERNEL — crossings + z-buffer + cross2d
# ──────────────────────────────────────────────

CROSSING_KERNEL = r"""
extern "C" __global__
void find_crossings(
    const float* __restrict__ proj,    // [N, 3] rotated curve in camera frame
    const int N,
    int*   out_i,
    int*   out_j,
    float* out_t,
    float* out_u,
    float* out_zi,
    float* out_zj,
    float* out_c2d,
    int*   count,
    const int max_cross
) {
    long long k = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    long long total = (long long)N * (N - 1) / 2;
    if (k >= total) return;

    long long lo = 0, hi = N - 1;
    while (lo < hi) {
        long long mid = (lo + hi) / 2;
        long long cnt = mid * ((long long)N - 1) - mid * (mid - 1) / 2;
        if (cnt <= k) lo = mid + 1;
        else hi = mid;
    }
    long long i = lo - 1;
    long long row_start = i * ((long long)N - 1) - i * (i - 1) / 2;
    long long j = i + 1 + (k - row_start);

    if (j == i + 1) return;
    if (i == 0 && j == (long long)(N - 1)) return;

    int ip1 = (int)((i + 1) % N);
    int jp1 = (int)((j + 1) % N);

    float p0x = proj[3*i      ];
    float p0y = proj[3*i    +1];
    float p0z = proj[3*i    +2];
    float p1x = proj[3*ip1    ];
    float p1y = proj[3*ip1  +1];
    float p1z = proj[3*ip1  +2];
    float p2x = proj[3*j      ];
    float p2y = proj[3*j    +1];
    float p2z = proj[3*j    +2];
    float p3x = proj[3*jp1    ];
    float p3y = proj[3*jp1  +1];
    float p3z = proj[3*jp1  +2];

    float dx1 = p1x - p0x, dy1 = p1y - p0y;
    float dx2 = p3x - p2x, dy2 = p3y - p2y;
    float denom = dx1 * dy2 - dy1 * dx2;       // = cross2d
    if (fabsf(denom) < 1e-10f) return;

    float dx0 = p2x - p0x, dy0 = p2y - p0y;
    float t = (dx0 * dy2 - dy0 * dx2) / denom;
    float u = (dx0 * dy1 - dy0 * dx1) / denom;

    if (t > 0.0f && t < 1.0f && u > 0.0f && u < 1.0f) {
        int idx = atomicAdd(count, 1);
        if (idx < max_cross) {
            float zi = p0z * (1.0f - t) + p1z * t;
            float zj = p2z * (1.0f - u) + p3z * u;
            out_i  [idx] = (int)i;
            out_j  [idx] = (int)j;
            out_t  [idx] = t;
            out_u  [idx] = u;
            out_zi [idx] = zi;
            out_zj [idx] = zj;
            out_c2d[idx] = denom;
        }
    }
}
"""


# ──────────────────────────────────────────────
# 3. SCAN
# ──────────────────────────────────────────────

def scan(fseries_path, camera_npz, output_path,
         N_curve=1024, max_cross_per_frame=512, period_pi=False):

    op = output_path.lower()
    if op.endswith('.parquet'):
        fmt = 'parquet'
    elif op.endswith('.npz'):
        fmt = 'npz'
    else:
        fmt = 'csv'

    print(f"[SCANNER] Loading .fseries: {fseries_path}")
    coeffs = parse_fseries(fseries_path)
    curve  = reconstruct_curve(*coeffs, N=N_curve, period_pi=period_pi)
    period_str = "[0,π)" if period_pi else "[0,2π)"
    print(f"[SCANNER]   Curve: {N_curve} points on {period_str}, "
          f"{len(coeffs[1])} Fourier terms")

    print(f"[SCANNER] Loading camera flight: {camera_npz}")
    cam_data = np.load(camera_npz)
    rotmats  = cam_data['rotmat']
    n_frames = len(rotmats)
    print(f"[SCANNER]   {n_frames} camera points")

    print(f"[SCANNER] Output: {output_path}  (Format: {fmt.upper()})")

    mod    = cp.RawModule(code=CROSSING_KERNEL)
    kernel = mod.get_function('find_crossings')

    curve_gpu = cp.asarray(curve)                  # [N, 3]
    inv_N     = np.float32(1.0 / float(N_curve))
    total_pairs = N_curve * (N_curve - 1) // 2
    block = 256
    grid  = (total_pairs + block - 1) // block
    print(f"[SCANNER] GPU: {total_pairs} segment pairs/frame, "
          f"Grid={grid}, Block={block}")

    # Per-frame buffers (reused)
    out_i  = cp.zeros(max_cross_per_frame, dtype=cp.int32)
    out_j  = cp.zeros(max_cross_per_frame, dtype=cp.int32)
    out_t  = cp.zeros(max_cross_per_frame, dtype=cp.float32)
    out_u  = cp.zeros(max_cross_per_frame, dtype=cp.float32)
    out_zi = cp.zeros(max_cross_per_frame, dtype=cp.float32)
    out_zj = cp.zeros(max_cross_per_frame, dtype=cp.float32)
    out_c2d= cp.zeros(max_cross_per_frame, dtype=cp.float32)
    count  = cp.zeros(1, dtype=cp.int32)

    rotmats_gpu = cp.asarray(rotmats)  # [n_frames, 3, 3] resident

    # Per-format initialisation
    if fmt in ('parquet', 'npz'):
        # Collect lists, concatenate at the end and write binary
        all_t_step = []
        all_seg_i  = []
        all_seg_j  = []
        all_t_i    = []
        all_t_j    = []
        all_z_i    = []
        all_z_j    = []
        all_c2d    = []
        csvfile = None
        writer  = None
    else:
        csvfile = open(output_path, 'w', newline='')
        writer  = csv.writer(csvfile)
        writer.writerow(['t_step', 'seg_i', 'seg_j',
                         't_norm_i', 't_norm_j',
                         'z_i', 'z_j', 'cross2d'])

    total_crossings = 0
    overflow_frames = 0
    t0 = time.time()
    report_every = max(1, n_frames // 50)

    try:
        for t_step in range(n_frames):
            R = rotmats_gpu[t_step]                       # [3, 3]
            rotated = curve_gpu @ R.T                     # [N, 3]
            proj = cp.ascontiguousarray(rotated.ravel()) # [N*3]
            count.fill(0)
            kernel(
                (grid,), (block,),
                (proj, np.int32(N_curve),
                 out_i, out_j, out_t, out_u,
                 out_zi, out_zj, out_c2d,
                 count, np.int32(max_cross_per_frame))
            )
            n = int(count[0])
            if n > max_cross_per_frame:
                overflow_frames += 1
                n = max_cross_per_frame

            if n > 0:
                ii = out_i [:n].get()
                jj = out_j [:n].get()
                tt = out_t [:n].get()
                uu = out_u [:n].get()
                zi = out_zi[:n].get()
                zj = out_zj[:n].get()
                c2 = out_c2d[:n].get()
                # Deterministic order: by (seg_i, seg_j)
                order = np.lexsort((jj, ii))
                ii = ii[order]; jj = jj[order]
                tt = tt[order]; uu = uu[order]
                zi = zi[order]; zj = zj[order]
                c2 = c2[order]
                t_norm_i = (ii.astype(np.float32) + tt) * inv_N
                t_norm_j = (jj.astype(np.float32) + uu) * inv_N

                if fmt in ('parquet', 'npz'):
                    all_t_step.append(np.full(n, t_step, dtype=np.int32))
                    all_seg_i.append(ii)
                    all_seg_j.append(jj)
                    all_t_i.append(t_norm_i)
                    all_t_j.append(t_norm_j)
                    all_z_i.append(zi)
                    all_z_j.append(zj)
                    all_c2d.append(c2)
                else:
                    for k in range(n):
                        writer.writerow([
                            t_step,
                            int(ii[k]), int(jj[k]),
                            f"{t_norm_i[k]:.8f}", f"{t_norm_j[k]:.8f}",
                            f"{zi[k]:.6f}",      f"{zj[k]:.6f}",
                            f"{c2[k]:.6e}",
                        ])
                total_crossings += n

            if (t_step + 1) % report_every == 0 or t_step + 1 == n_frames:
                elapsed = time.time() - t0
                pct = (t_step + 1) / n_frames * 100
                eta = elapsed * (n_frames - (t_step + 1)) / (t_step + 1)
                print(f"  [{pct:5.1f}%] t_step={t_step+1}/{n_frames} | "
                      f"Crossings: {total_crossings:,} | "
                      f"{elapsed:.1f}s elapsed, ~{eta:.1f}s ETA",
                      flush=True)

        # ────────── write output ──────────
        if fmt in ('parquet', 'npz'):
            print(f"[SCANNER] Concatenating {total_crossings:,} crossings ...")
            tcat = time.time()
            if all_t_step:
                t_step_arr = np.concatenate(all_t_step)
                seg_i_arr  = np.concatenate(all_seg_i)
                seg_j_arr  = np.concatenate(all_seg_j)
                t_i_arr    = np.concatenate(all_t_i)
                t_j_arr    = np.concatenate(all_t_j)
                z_i_arr    = np.concatenate(all_z_i)
                z_j_arr    = np.concatenate(all_z_j)
                c2d_arr    = np.concatenate(all_c2d)
            else:
                t_step_arr = np.zeros(0, dtype=np.int32)
                seg_i_arr  = np.zeros(0, dtype=np.int32)
                seg_j_arr  = np.zeros(0, dtype=np.int32)
                t_i_arr    = np.zeros(0, dtype=np.float32)
                t_j_arr    = np.zeros(0, dtype=np.float32)
                z_i_arr    = np.zeros(0, dtype=np.float32)
                z_j_arr    = np.zeros(0, dtype=np.float32)
                c2d_arr    = np.zeros(0, dtype=np.float32)
            tcat_done = time.time() - tcat
            twrite = time.time()
            if fmt == 'parquet':
                import pyarrow as pa
                import pyarrow.parquet as pq
                table = pa.Table.from_arrays(
                    [pa.array(t_step_arr), pa.array(seg_i_arr), pa.array(seg_j_arr),
                     pa.array(t_i_arr),    pa.array(t_j_arr),
                     pa.array(z_i_arr),    pa.array(z_j_arr),
                     pa.array(c2d_arr)],
                    names=['t_step', 'seg_i', 'seg_j',
                           't_norm_i', 't_norm_j',
                           'z_i', 'z_j', 'cross2d'],
                )
                print(f"[SCANNER] Writing Parquet ({tcat_done:.1f}s concat) ...")
                pq.write_table(table, output_path, compression='snappy')
                print(f"[SCANNER]   Parquet write {time.time()-twrite:.1f}s")
            else:
                print(f"[SCANNER] Writing NPZ ({tcat_done:.1f}s concat) ...")
                np.savez_compressed(
                    output_path,
                    t_step=t_step_arr,
                    seg_i=seg_i_arr,
                    seg_j=seg_j_arr,
                    t_norm_i=t_i_arr,
                    t_norm_j=t_j_arr,
                    z_i=z_i_arr,
                    z_j=z_j_arr,
                    cross2d=c2d_arr,
                )
                print(f"[SCANNER]   NPZ write {time.time()-twrite:.1f}s")
    finally:
        if csvfile is not None:
            csvfile.close()

    total_time = time.time() - t0
    sz = os.path.getsize(output_path) if os.path.exists(output_path) else 0
    print(f"[SCANNER] Done.")
    print(f"  Total Crossings : {total_crossings:,}")
    print(f"  Overflow frames : {overflow_frames}  "
          f"(use --max-cross-per-frame {max_cross_per_frame*2} if > 0)")
    print(f"  Run time        : {total_time:.1f}s "
          f"({n_frames/total_time:.0f} t_step/s)")
    print(f"  Output          : {output_path}  ({sz/1e6:.1f} MB)")


# ──────────────────────────────────────────────
# 4. ENTRY
# ──────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('fseries',  help='path to the .fseries file')
    parser.add_argument('--camera', required=True,
                        help='path to camera_*.npz (see make_camera_trajectory.py)')
    parser.add_argument('--output', default='crossings.parquet',
                        help='output path. Extension .parquet (default, fast and small), '
                             '.npz (binary, no extra package) or '
                             '.csv (human-readable, slow).')
    parser.add_argument('--N',      type=int, default=1024,
                        help='curve sampling density (default 1024). '
                             'Consider N >= 4096 for dense knots.')
    parser.add_argument('--max-cross-per-frame', type=int, default=512,
                        help='output buffer per frame; increase on overflow.')
    parser.add_argument('--period-pi', action='store_true',
                        help='sample the curve on [0,pi) instead of [0,2pi). '
                             'Required for pi-periodic curves (only even '
                             'harmonics j in the .fseries). Otherwise the doubly '
                             'covered curve is scanned, every crossing is found 4 times '
                             'and triple-point phantoms appear.')
    args = parser.parse_args()

    scan(
        fseries_path=args.fseries,
        camera_npz=args.camera,
        output_path=args.output,
        N_curve=args.N,
        max_cross_per_frame=args.max_cross_per_frame,
        period_pi=args.period_pi,
    )
