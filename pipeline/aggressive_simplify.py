#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
aggressive_simplify.py
======================
Over all frames with n_cross in the chosen range:
  1. Build all signed Gauss codes (in parallel via multiprocessing)
  2. Deduplicate (symmetries often produce identical codes)
  3. Batch subprocess in base (regina) and snappy_env (spherogram, all modes)
     with the loop INSIDE — cuts subprocess overhead drastically
  4. Track the minimal crossing result over all frames + all simplify modes

Output: aggressive_simplify_result.json with:
  - best_signed_gauss        (smallest classical diagram found)
  - best_simplified_size     (true crossing-count lower bound)
  - best_origin_frame
  - best_homfly, best_jones
  - distribution             (histogram of the crossing counts found)
"""

import argparse
import site_config
import json
import math
import multiprocessing as mp
import os
import subprocess
import concurrent.futures as cf
import sys
import time
import numpy as np
import pandas as pd

sys.path.insert(0, site_config.FREEZE_DIR)
sys.path.insert(0, site_config.PIPELINE_DIR)
# Disc #222/#586: subprocess.run + an intermediate layer (conda run, wrapper script) -> the
# timeout kill only hits the intermediate layer, the grandchildren keep running.
from pipeline_subprocess import run_capture_with_group_kill as _run_with_pgkill
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


def reconstruct_curve(offset, ax, bx, ay, by, az, bz, N):
    s = np.linspace(0, 2 * math.pi, N, endpoint=False)
    js = np.arange(1, len(ax) + 1)
    angles = s[:, None] * js[None, :]
    ca = np.cos(angles)
    sa = np.sin(angles)
    return np.stack([
        ca @ ax + sa @ bx + offset[0],
        ca @ ay + sa @ by + offset[1],
        ca @ az + sa @ bz + offset[2],
    ], axis=1).astype(np.float32)


def build_signed_gauss_for_frame(frame_data):
    """frame_data: (t_step, list-of-(seg_i, seg_j, ti, tj, z_i, z_j, cross2d))
    Returns (t_step, signed_gauss_string, n_cross)."""
    t_step, crossings = frame_data
    tokens = []
    for diag_idx, (seg_i, seg_j, ti, tj, z_i, z_j, c2d) in enumerate(crossings, start=1):
        i_over = z_i > z_j
        # FIXED 2026-06-05 (mirror finding, error log ~08:45): + <-> det[T_O,T_U]>0
        # — the old formula was the wrong way round (proven via the Gauss integral), same
        # correction as kinematic_resolve.build_signed_gauss.
        signed_cross = c2d if i_over else -c2d
        sign = '+' if signed_cross > 0 else '-'
        tokens.append((ti, diag_idx, i_over, sign))
        tokens.append((tj, diag_idx, not i_over, sign))
    tokens.sort(key=lambda t: t[0])
    parts = []
    for pos, idx, overpass, sign in tokens:
        ch = 'O' if overpass else 'U'
        parts.append(f"{ch}{idx}{sign}")
    return (t_step, ''.join(parts), len(crossings))


SUBPROCESS_REGINA = """
import sys, json
sys.path.insert(0, '@FREEZE_DIR@')
from pipeline_seal import seal_all
seal_all(verbose=False)
import regina
import os
# R3_SKIP_HOMFLY: homfly/jones are ONLY a by-product here (best_homfly/best_jones, NOT
# the frame selection criterion — that is size_simp_regina, line ~481). On densified
# high-c diagrams regina.homfly() explodes exponentially (15GB/50min/frame -> OOM/orphan).
# The adjudicator path (det-only) sets R3_SKIP_HOMFLY=1 -> identical selection, no hog.
_SKIP_HJ = bool(os.environ.get('R3_SKIP_HOMFLY'))

for line in sys.stdin:
    line = line.strip()
    if not line:
        continue
    try:
        item = json.loads(line)
    except Exception:
        continue
    sg = item['sg']
    out = {'t_step': item.get('t_step'), 'sg_input_len': len(sg)//3}
    try:
        L = regina.Link.fromSignedGauss(sg)
        out['classical'] = bool(L.isClassical())
        out['vGenus']    = int(L.virtualGenus())
        out['size_in']   = int(L.size())
        L.simplify()
        out['size_simp_regina'] = int(L.size())
        if out['classical']:
            try: out['signed_gauss_simp'] = str(L.signedGauss())
            except Exception as e: out['sg_simp_err'] = str(e)
            if not _SKIP_HJ:
                try: out['homfly'] = str(L.homfly(regina.Algorithm.Treewidth))
                except Exception as e: out['homfly_err'] = str(e)
                try: out['jones']  = str(L.jones(regina.Algorithm.Treewidth))
                except Exception as e: out['jones_err'] = str(e)
            try: out['pd'] = str(L.pd())
            except Exception as e: out['pd_err'] = str(e)
    except Exception as e:
        out['parse_err'] = str(e)
    print(json.dumps(out), flush=True)
""".replace('@FREEZE_DIR@', site_config.FREEZE_DIR)

SUBPROCESS_SPHEROGRAM = """
import sys, json, re
sys.path.insert(0, '@FREEZE_DIR@')
from pipeline_seal import seal_all
seal_all(verbose=False)
import snappy, spherogram

def parse_pd(pd_str):
    crossings = re.findall(r'X[\\(\\[](\\-?\\d+(?:,\\s*\\-?\\d+){3})[\\)\\]]', pd_str)
    if not crossings:
        return None
    return [tuple(int(v.strip()) for v in cr.split(',')) for cr in crossings]

for line in sys.stdin:
    line = line.strip()
    if not line:
        continue
    try:
        item = json.loads(line)
    except Exception:
        continue
    pd_str = item['pd']
    out = {'t_step': item.get('t_step')}
    try:
        pd_list = parse_pd(pd_str)
        if pd_list is None:
            out['err'] = 'pd_parse'
            print(json.dumps(out), flush=True); continue
        L0 = spherogram.Link(pd_list)
        out['n_in'] = len(L0.crossings)
        # Disc #562 (2026-07-25): return NOT just the number, but the RECEIPT.
        # Before, this worker delivered only len(Lm.crossings) -> the
        # smaller diagram was counted and thrown away, and Regina's larger code
        # was passed on. Result: crossing_number_min lay above the PROVEN truth on 18 of 40
        # calibration-ladder rungs (by up to +9).
        best_n, best_mode, best_pd = out['n_in'], None, None
        for mode in ('basic', 'level', 'pickup', 'global'):
            try:
                Lm = L0.copy()
                Lm.simplify(mode=mode)
                n = len(Lm.crossings)
                out[f'n_{mode}'] = n
                if n < best_n and n > 0:
                    try:
                        best_pd = str(Lm.PD_code(KnotTheory=True, min_strand_index=1))
                        best_n, best_mode = n, mode
                    except Exception as e:
                        out['pd_best_err'] = str(e)
            except Exception as e:
                out[f'err_{mode}'] = str(e)
        if best_pd:
            out['pd_best']      = best_pd
            out['pd_best_mode'] = best_mode
            out['pd_best_n']    = best_n
    except Exception as e:
        out['err'] = str(e)
    print(json.dumps(out), flush=True)
""".replace('@FREEZE_DIR@', site_config.FREEZE_DIR)

# Disc #562: third mini pass — spherogram-reduced PD back to regina, in order to
# get a SIGNED GAUSS code on the smaller diagram (stage 5 speaks
# signed Gauss, spherogram speaks PD). Round trip verified 2026-07-25 on
# trainstation: fromPD(spherogram KT PD) -> 18 crossings, classical, Jones
# bit-identical to the stored record.
# GUARD against convention/conversion errors (Frankenstein doctrine Disc #256):
# det (= |Delta(-1)|) of the converted diagram MUST equal the det of the Regina diagram
# of the same frame. Unequal -> the handoff is discarded, the old path stays.
SUBPROCESS_REGINA_FROMPD = """
import sys, json
sys.path.insert(0, '@FREEZE_DIR@')
from pipeline_seal import seal_all
seal_all(verbose=False)
import regina

def det_of(L):
    a = L.alexander()
    # regina.Integer breaks int() directly -> go via str() (footgun class #358/#373)
    return abs(sum(int(str(a[i])) * (-1) ** i for i in range(a.degree() + 1)))

for line in sys.stdin:
    line = line.strip()
    if not line:
        continue
    try:
        item = json.loads(line)
    except Exception:
        continue
    out = {'t_step': item.get('t_step')}
    try:
        L = regina.Link.fromPD(item['pd'])
        out['classical'] = bool(L.isClassical())
        out['size_in']   = int(L.size())
        L.simplify()
        out['size']      = int(L.size())
        if out['classical']:
            out['signed_gauss'] = str(L.signedGauss())
            out['det']          = det_of(L)
    except Exception as e:
        out['err'] = str(e)
    try:
        if item.get('sg_ref'):
            R = regina.Link.fromSignedGauss(item['sg_ref'])
            R.simplify()
            out['det_ref']  = det_of(R)
            out['size_ref'] = int(R.size())
    except Exception as e:
        out['det_ref_err'] = str(e)
    print(json.dumps(out), flush=True)
""".replace('@FREEZE_DIR@', site_config.FREEZE_DIR)


def run_batch_subprocess(env, code, items, timeout=900):
    """Send every item dict as a JSON line via stdin to the subprocess; collect
    JSON lines from stdout. Raises subprocess.TimeoutExpired when the timeout is exceeded."""
    cmd = ['conda', 'run', '-n', env, '--no-capture-output', 'python', '-c', code]
    inp = '\n'.join(json.dumps(it) for it in items) + '\n'
    res = _run_with_pgkill(cmd, timeout=timeout, input=inp)
    if res.returncode != 0:
        print(f"[Subprocess {env}] rc={res.returncode}, stderr tail:\n{res.stderr[-500:]}",
              file=sys.stderr)
    out = []
    for line in res.stdout.splitlines():
        line = line.strip()
        if line.startswith('{') and line.endswith('}'):
            try:
                out.append(json.loads(line))
            except Exception:
                pass
    return out


def run_minibatch_parallel(env, code, items, mini_batch_size=20,
                           mini_batch_timeout=60, max_workers=12,
                           global_timeout=600, label='Subprocess'):
    """Discovery #149 (variant 1, stage 4 refactor): parallel mini-batches.

    Instead of one large subprocess call (the old run_batch_subprocess), which kills
    the whole batch on one hanging code and thereby loses all the healthy
    codes too, small mini-batches of `mini_batch_size`
    codes are sent in parallel over `max_workers` threads, each to its own subprocess.

    Each mini-batch has its own timeout (`mini_batch_timeout`). Hanging
    mini-batches are killed, other mini-batches keep running. In addition
    a global wall-clock cap (`global_timeout`): when it expires, open
    futures are cancelled and partial results are returned.

    Returns: (results_list, status_dict)
    """
    items = list(items)
    n_total = len(items)
    if n_total == 0:
        return [], {'status': 'empty', 'n_total_codes': 0,
                    'n_codes_processed': 0, 'n_codes_timeout': 0,
                    'partial': False, 'duration_s': 0.0}
    mini_batches = [items[i:i + mini_batch_size]
                    for i in range(0, n_total, mini_batch_size)]
    n_batches = len(mini_batches)

    results = []
    n_batches_done = 0
    n_codes_done = 0
    n_codes_timeout = 0

    def _worker(batch):
        try:
            return ('OK', run_batch_subprocess(env, code, batch,
                                               timeout=mini_batch_timeout))
        except subprocess.TimeoutExpired:
            return ('TIMEOUT', len(batch))
        except Exception as e:
            return ('ERROR', repr(e), len(batch))

    t_start = time.time()
    print(f"[{label}] Mini-batches: {n_batches} x {mini_batch_size} codes, "
          f"workers={max_workers}, per-batch-timeout={mini_batch_timeout}s, "
          f"global-timeout={global_timeout}s")

    with cf.ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = {ex.submit(_worker, b): b for b in mini_batches}
        global_timed_out = False
        try:
            for fut in cf.as_completed(futures, timeout=global_timeout):
                batch = futures[fut]
                try:
                    tag, *rest = fut.result()
                except Exception as e:
                    print(f"[{label}] worker exception: {e!r}", file=sys.stderr)
                    n_codes_timeout += len(batch)
                    continue
                if tag == 'OK':
                    batch_results = rest[0]
                    results.extend(batch_results)
                    n_codes_done += len(batch)
                elif tag == 'TIMEOUT':
                    n_codes_timeout += rest[0]
                else:  # ERROR
                    n_codes_timeout += rest[1]
                n_batches_done += 1
                elapsed = time.time() - t_start
                if n_batches_done % max(1, n_batches // 20) == 0 or \
                        n_batches_done == n_batches:
                    eta = elapsed * (n_batches - n_batches_done) / max(n_batches_done, 1)
                    print(f"[{label}]   {n_batches_done}/{n_batches} batches "
                          f"({n_codes_done}/{n_total} codes done, "
                          f"{n_codes_timeout} timeout/err), "
                          f"{elapsed:.1f}s elapsed, ETA {eta:.1f}s", flush=True)
        except cf.TimeoutError:
            global_timed_out = True
            elapsed = time.time() - t_start
            print(f"[{label}] GLOBAL TIMEOUT after {elapsed:.1f}s "
                  f"(limit {global_timeout}s). Cancelling "
                  f"{n_batches - n_batches_done} open mini-batches.",
                  flush=True)
            for fut, batch in futures.items():
                if not fut.done():
                    fut.cancel()
                    n_codes_timeout += len(batch)

    status = {
        'n_total_codes': n_total,
        'n_mini_batches': n_batches,
        'mini_batch_size': mini_batch_size,
        'mini_batch_timeout_s': mini_batch_timeout,
        'global_timeout_s': global_timeout,
        'n_codes_processed': n_codes_done,
        'n_codes_timeout': n_codes_timeout,
        'duration_s': time.time() - t_start,
        'partial': n_codes_done < n_total,
        'global_timed_out': global_timed_out,
    }
    return results, status


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('crossings', help='crossings.parquet from knot_scanner.py')
    parser.add_argument('--fseries', required=True)
    parser.add_argument('--N', type=int, default=4096,
                        help='Curve sampling density (must match the scanner).')
    parser.add_argument('--n-min', type=int, default=14,
                        help='Consider frames with n_cross >= n_min.')
    parser.add_argument('--n-max', type=int, default=20,
                        help='Consider frames with n_cross <= n_max.')
    parser.add_argument('--max-frames', type=int, default=20000,
                        help='Hard cap on the number of frames (sample randomly if there are more).')
    parser.add_argument('--processes', type=int, default=12,
                        help='multiprocessing pool size for code generation.')
    parser.add_argument('--regina-batch-size', type=int, default=0,
                        help='Number of codes per Regina mini-batch. '
                             '0 (default) = adaptive: c<=30 → 200, 30<c<=50 → 50, c>50 → 10.')
    parser.add_argument('--regina-timeout', type=int, default=0,
                        help='Timeout per Regina mini-batch (seconds). '
                             '0 (default) = adaptive: c<=30 → 60s, 30<c<=50 → 120s, c>50 → 240s.')
    parser.add_argument('--global-timeout', type=int, default=600,
                        help='Global wall-clock timeout per stage phase '
                             '(Regina, Spherogram). Discovery #149.')
    parser.add_argument('--max-workers', type=int, default=12,
                        help='Number of parallel mini-batch workers (ThreadPool).')
    parser.add_argument('--max-c-multiplier', type=float, default=0.0,
                        help='Pre-filter: codes with c > max_c_multiplier * n_min '
                             'are discarded. 0 (default) = off. '
                             'Recommended: 2.0 — cuts absurdly dense diagrams.')
    parser.add_argument('--spherogram-top-n', type=int, default=200,
                        help='Pass the top-N diagrams (smallest Regina reduction) on to Spherogram.')
    parser.add_argument('--output', default='aggressive_simplify_result.json')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    t_total = time.time()

    print(f"[AGG] Loading {args.crossings}")
    df = pd.read_parquet(args.crossings)
    print(f"[AGG]   {len(df):,} Crossings, {df['t_step'].nunique():,} Frames")

    # Filter frames with n_cross ∈ [n_min, n_max]
    counts = df.groupby('t_step').size()
    frames = counts[(counts >= args.n_min) & (counts <= args.n_max)].index.tolist()
    print(f"[AGG] Frames with n_cross ∈ [{args.n_min}, {args.n_max}]: {len(frames):,}")

    if len(frames) > args.max_frames:
        rng = np.random.default_rng(args.seed)
        frames = sorted(rng.choice(frames, size=args.max_frames, replace=False).tolist())
        print(f"[AGG]   Sampled down to {len(frames):,} frames")

    df_sub = df[df['t_step'].isin(set(frames))].copy()

    # Pre-filter: only transversal crossings (sin_angle > 0.05)
    coeffs = parse_fseries(args.fseries)
    curve = reconstruct_curve(*coeffs, N=args.N)
    d = np.roll(curve, -1, axis=0) - curve
    dmag = np.linalg.norm(d, axis=1).astype(np.float32)
    seg_i = df_sub['seg_i'].to_numpy()
    seg_j = df_sub['seg_j'].to_numpy()
    c2d_  = df_sub['cross2d'].to_numpy()
    sin_a = np.abs(c2d_) / (dmag[seg_i].astype(np.float64) *
                             dmag[seg_j].astype(np.float64) + 1e-30)
    df_sub['sin_angle'] = sin_a

    # Frame-by-frame code generation
    print(f"[AGG] Building {len(frames):,} signed Gauss codes (multiprocessing.Pool({args.processes}))")
    t0 = time.time()
    frame_data_list = []
    grouped = df_sub.groupby('t_step')
    for t_step in frames:
        grp = grouped.get_group(t_step)
        crs = list(zip(grp['seg_i'].to_numpy().tolist(),
                       grp['seg_j'].to_numpy().tolist(),
                       grp['t_norm_i'].to_numpy().tolist(),
                       grp['t_norm_j'].to_numpy().tolist(),
                       grp['z_i'].to_numpy().tolist(),
                       grp['z_j'].to_numpy().tolist(),
                       grp['cross2d'].to_numpy().tolist()))
        frame_data_list.append((int(t_step), crs))

    with mp.Pool(args.processes) as pool:
        results = pool.map(build_signed_gauss_for_frame, frame_data_list,
                           chunksize=200)
    sg_records = [{'t_step': t, 'sg': sg, 'n': n} for (t, sg, n) in results]
    print(f"[AGG]   {len(sg_records):,} Codes in {time.time()-t0:.1f}s")

    # Deduplicate identical signed_gauss strings
    seen = {}
    for r in sg_records:
        seen.setdefault(r['sg'], r)  # first frame wins
    sg_unique = list(seen.values())
    print(f"[AGG]   Unique codes: {len(sg_unique):,} (of {len(sg_records):,})")

    # Optional pre-filter: drop codes with c > max_c_multiplier * n_min
    # (absurdly dense diagrams rarely reduce to the minimum anyway)
    if args.max_c_multiplier > 0 and sg_unique:
        c_threshold = int(args.max_c_multiplier * args.n_min)
        before = len(sg_unique)
        sg_unique = [r for r in sg_unique if r['n'] <= c_threshold]
        print(f"[AGG]   Pre-filter: c > {c_threshold} discarded "
              f"({before} → {len(sg_unique)})")

    # Adaptive batch size + timeout by mean crossing count
    if sg_unique:
        c_mean = sum(r['n'] for r in sg_unique) / len(sg_unique)
        c_median = sorted(r['n'] for r in sg_unique)[len(sg_unique)//2]
    else:
        c_mean = 0
        c_median = 0
    c_ref = max(c_mean, c_median)  # conservative: the larger of the two

    if args.regina_batch_size > 0:
        regina_batch_size = args.regina_batch_size
        batch_src = 'user'
    else:
        # Discovery #149: smaller mini-batches, run in parallel via a ThreadPool instead.
        if c_ref <= 30:
            regina_batch_size = 200
        elif c_ref <= 50:
            regina_batch_size = 50
        else:
            regina_batch_size = 10
        batch_src = f'adaptiv (c_ref={c_ref:.1f})'

    if args.regina_timeout > 0:
        regina_timeout = args.regina_timeout
        timeout_src = 'user'
    else:
        # Per-mini-batch timeout, NOT a total timeout. Global cap via
        # --global-timeout.
        if c_ref <= 30:
            regina_timeout = 60
        elif c_ref <= 50:
            regina_timeout = 120
        else:
            regina_timeout = 240
        timeout_src = f'adaptiv (c_ref={c_ref:.1f})'

    print(f"[AGG]   c_mean={c_mean:.1f}, c_median={c_median} "
          f"→ mini-batch-size={regina_batch_size} [{batch_src}], "
          f"per-batch-timeout={regina_timeout}s [{timeout_src}], "
          f"global-timeout={args.global_timeout}s, workers={args.max_workers}")

    # Regina mini-batches in parallel (Discovery #149)
    print(f"[AGG] Regina simplify on {len(sg_unique):,} unique codes")
    regina_items = [{'t_step': r['t_step'], 'sg': r['sg']} for r in sg_unique]
    regina_results, regina_status = run_minibatch_parallel(
        env='base', code=SUBPROCESS_REGINA, items=regina_items,
        mini_batch_size=regina_batch_size,
        mini_batch_timeout=regina_timeout,
        max_workers=args.max_workers,
        global_timeout=args.global_timeout,
        label='AGG/Regina',
    )
    print(f"[AGG]   Regina done: {regina_status['n_codes_processed']}/"
          f"{regina_status['n_total_codes']} codes processed "
          f"({regina_status['n_codes_timeout']} timeout/err) "
          f"in {regina_status['duration_s']:.1f}s"
          f"{' [PARTIAL/global timeout]' if regina_status['global_timed_out'] else ''}")
    classical = [r for r in regina_results if r.get('classical')]
    print(f"[AGG]   Classical diagrams: {len(classical)}/{len(regina_results)}")

    if not classical:
        print("[AGG] NO classical diagram found in the multi-frame sweep "
              f"(run {'partial' if regina_status['partial'] else 'complete'}).",
              file=sys.stderr)
        report = {
            'no_classical_found': True,
            'regina_results_count': len(regina_results),
            'regina_status': regina_status,
            'stage4_partial': regina_status['partial'],
            'duration_seconds': time.time() - t_total,
        }
        with open(args.output, 'w') as f:
            json.dump(report, f, indent=2, default=str)
        # Discovery #149: no hard rc=3 any more. The pipeline wrapper lets stage 5
        # run on the kinematic_resolve diagram (the algebra is diagram-independent).
        # If REALLY NO classical diagram was found, that is a
        # topological finding (single-frame stagnation), not a pipeline bug.
        print("[AGG] -> no hard exit, the pipeline can continue on the kinematic diagram.")
        sys.exit(0)

    # Sort by Regina simplified size
    classical.sort(key=lambda r: r.get('size_simp_regina', 999))
    sizes = [r['size_simp_regina'] for r in classical]
    from collections import Counter
    size_dist = Counter(sizes)
    print(f"[AGG]   Regina simplified-size distribution: {dict(sorted(size_dist.items()))}")
    min_regina = min(sizes)
    print(f"[AGG]   MIN Regina-simplified-size: {min_regina}")

    # Pass the N smallest on to Spherogram (parallel mini-batches, Discovery #149)
    top = classical[:args.spherogram_top_n]
    print(f"\n[AGG] Spherogram multi-mode simplify on the top {len(top)} diagrams")
    sphe_items = [{'t_step': r['t_step'], 'pd': r['pd']} for r in top
                  if r.get('pd')]
    # adaptive: same scaling as Regina, somewhat more conservative
    sphe_batch_size = max(5, regina_batch_size // 2)
    sphe_timeout = max(60, regina_timeout * 2)  # spherogram global is more expensive
    sphe_results, sphe_status = run_minibatch_parallel(
        env='snappy_env', code=SUBPROCESS_SPHEROGRAM, items=sphe_items,
        mini_batch_size=sphe_batch_size,
        mini_batch_timeout=sphe_timeout,
        max_workers=args.max_workers,
        global_timeout=args.global_timeout,
        label='AGG/Spherogram',
    )
    print(f"[AGG]   Spherogram done: {sphe_status['n_codes_processed']}/"
          f"{sphe_status['n_total_codes']} codes in "
          f"{sphe_status['duration_s']:.1f}s"
          f"{' [PARTIAL]' if sphe_status['global_timed_out'] else ''}")

    # Merge Spherogram results back into classical (per t_step)
    sphe_by_t = {r['t_step']: r for r in sphe_results}
    for r in top:
        sp = sphe_by_t.get(r['t_step'])
        if sp:
            for k, v in sp.items():
                if k != 't_step':
                    r[f'sphe_{k}'] = v

    # Best diagram: minimum over Regina + Spherogram modes
    def best_size(r):
        candidates = [r.get('size_simp_regina', 999)]
        for mode in ('basic', 'level', 'pickup', 'global'):
            v = r.get(f'sphe_n_{mode}')
            if v is not None:
                candidates.append(v)
        return min(candidates)

    classical.sort(key=best_size)
    best = classical[0]
    best_min = best_size(best)
    # R3_DUMP_FRAMES: write out ALL classical frames (t_step/size/signed_gauss) —
    # for the per-frame det distribution diagnosis (several distinct knots at min-cross?).
    if os.environ.get('R3_DUMP_FRAMES'):
        dump = [{'t_step': r.get('t_step'), 'size_simp_regina': r.get('size_simp_regina'),
                 'best_size': best_size(r), 'signed_gauss_simp': r.get('signed_gauss_simp')}
                for r in classical if r.get('signed_gauss_simp')]
        dpath = args.output.rsplit('.json', 1)[0] + '_FRAMES.json'
        with open(dpath, 'w') as f:
            json.dump(dump, f)
        print(f"[AGG] R3_DUMP_FRAMES -> {len(dump)} classical frames to {dpath}")
    print(f"\n[AGG] BEST diagram: t_step={best['t_step']}, "
          f"min_size={best_min}")
    print(f"  Regina simplified: {best.get('size_simp_regina')}")
    for mode in ('basic', 'level', 'pickup', 'global'):
        v = best.get(f'sphe_n_{mode}')
        if v is not None:
            print(f"  Sphe {mode:>7s}: {v}")
    if best.get('signed_gauss_simp'):
        print(f"  Signed Gauss (Regina-reduced): {best['signed_gauss_simp']}")
    if best.get('homfly'):
        print(f"  HOMFLY: {best['homfly']}")
    if best.get('jones'):
        print(f"  Jones : {best['jones']}")

    # Global distribution over ALL crossing counts
    all_min_sizes = [best_size(r) for r in classical]
    all_dist = Counter(all_min_sizes)
    print(f"\n[AGG] Global min-size distribution: {dict(sorted(all_dist.items()))}")

    # ── Disc #562: HANDOFF WITH RECEIPT ───────────────────────────────────────
    # Until 2026-07-25 merge_aggressive always passed on best.signed_gauss_simp (REGINA),
    # although spherogram measurably got smaller in the same run -> the number was known,
    # the code was missing. Exactly this gap is closed here: the smallest
    # spherogram diagram goes back through regina (fromPD), the result is a
    # signed Gauss code on the smaller diagram. The det gate decides on acceptance.
    handoff = {
        'size':         best.get('size_simp_regina'),
        'signed_gauss': best.get('signed_gauss_simp'),
        'source':       'regina',
        't_step':       best.get('t_step'),
    }
    handoff_log = []
    cands = [r for r in classical
             if r.get('sphe_pd_best') and r.get('signed_gauss_simp')
             and r.get('sphe_pd_best_n') is not None
             and r['sphe_pd_best_n'] < r.get('size_simp_regina', 999)]
    cands.sort(key=lambda r: r['sphe_pd_best_n'])
    cands = cands[:int(os.environ.get('AGG_HANDOFF_TOP_N', '6'))]
    if cands:
        print(f"\n[AGG] Handoff pass (Disc #562): {len(cands)} candidates with "
              f"spherogram < regina (smallest: {cands[0]['sphe_pd_best_n']} vs "
              f"{cands[0].get('size_simp_regina')})")
        ho_items = [{'t_step': r['t_step'], 'pd': r['sphe_pd_best'],
                     'sg_ref': r['signed_gauss_simp']} for r in cands]
        ho_results, ho_status = run_minibatch_parallel(
            env='base', code=SUBPROCESS_REGINA_FROMPD, items=ho_items,
            mini_batch_size=max(2, regina_batch_size // 4),
            mini_batch_timeout=max(60, regina_timeout),
            max_workers=args.max_workers,
            global_timeout=max(120, args.global_timeout // 4),
            label='AGG/Handoff',
        )
        accepted = []
        for r in ho_results:
            size = r.get('size')
            reason = None
            if r.get('err'):                      reason = f"fromPD error: {r['err']}"
            elif not r.get('classical'):          reason = 'not classical'
            elif not r.get('signed_gauss'):       reason = 'no signedGauss'
            elif r.get('det') is None or r.get('det_ref') is None:
                                                  reason = 'det not computable'
            elif r['det'] != r['det_ref']:        reason = (f"det-GATE FAIL: {r['det']} != "
                                                            f"{r['det_ref']} (Referenz)")
            elif size is None or size >= handoff['size']:
                                                  reason = f"not smaller ({size} >= {handoff['size']})"
            if reason:
                handoff_log.append({'t_step': r.get('t_step'), 'rejected': reason})
                print(f"[AGG]   rejected t_step={r.get('t_step')}: {reason}")
            else:
                accepted.append(r)
        if accepted:
            accepted.sort(key=lambda r: r['size'])
            a0 = accepted[0]
            handoff = {
                'size':         a0['size'],
                'signed_gauss': a0['signed_gauss'],
                'source':       'spherogram+regina_fromPD',
                't_step':       a0.get('t_step'),
                'det':          a0.get('det'),
                'det_ref':      a0.get('det_ref'),
                'regina_size':  best.get('size_simp_regina'),
            }
            print(f"[AGG]   ✅ Handoff accepted: {handoff['size']} crossings "
                  f"(regina had {best.get('size_simp_regina')}), det={a0.get('det')} "
                  f"== reference det ✓")
        handoff_log.append({'status': ho_status})
    else:
        print(f"\n[AGG] Handoff pass (Disc #562): no candidate "
              f"(spherogram never smaller than regina) -> the regina code stays.")
    if handoff['size'] is not None and best_min < handoff['size']:
        print(f"[AGG]   ⚠️  measured min_size_global={best_min} lies BELOW the "
              f"receipted handoff={handoff['size']} — the receipt for the optimum is missing "
              f"(not an error, but substance given away).")

    # Discovery #149: stage4_partial marker if Regina or Spherogram
    # timed out globally. The mandatory-invariants stage downstream may evaluate it.
    stage4_partial = regina_status.get('partial') or sphe_status.get('partial')
    report = {
        'n_frames_processed': len(sg_records),
        'n_unique_codes':     len(sg_unique),
        'n_classical':        len(classical),
        'best': best,
        'handoff':            handoff,
        'handoff_log':        handoff_log,
        'min_size_global':    best_min,
        'min_size_distribution': dict(sorted(all_dist.items())),
        'sphe_top_n_used':    len(sphe_items),
        'duration_seconds':   time.time() - t_total,
        'stage4_partial':     stage4_partial,
        'regina_status':      regina_status,
        'sphe_status':        sphe_status,
    }
    with open(args.output, 'w') as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\n[AGG] Writing {args.output}")
    print(f"[AGG] Total {time.time()-t_total:.1f}s")


if __name__ == '__main__':
    main()
