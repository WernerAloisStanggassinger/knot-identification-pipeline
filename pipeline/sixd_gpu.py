#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Knot-Identification Pipeline — Werner Alois Stanggassinger & Claude Code (Anthropic).
# ============================================================================
# sixd_gpu.py — optional GPU fast path for the SONO flow (canonical_flow).
#
# PROFILE FINDING (double-6, n=4096, 2026-06-08, cProfile):
#   overlap_pairs = 81% of the runtime (1117s/1374s, 99055 calls). The
#   original O(n^2) target (pair_dist/gauss_wr_acn/min_nonneighbor_dist)
#   together was <3%. Hence the fast path targets EXCLUSIVELY
#   overlap_pairs — the per-iteration cell-list overlap detection.
#
# METHOD: brute-force O(n^2) pair mask on the GPU via Gram matrix
#   D2_ij = |x_i|^2 + |x_j|^2 - 2 x_i.x_j  (matmul, the GPU crushes that),
#   masked with ring distance > excl and upper triangle half. On CPU
#   brute force was ~2s/pass (470MB temporaries, calibration-1 finding) -> hence
#   there the cell list; on the GPU the (n,n) matrix is cheap (134MB @
#   448GB/s) and the kernel launch overhead amortizes from large n on.
#
# IMPORTANT (doctrine): NOT bit-identical to the CPU cell list — different pair
#   order -> different np.add.at float summation -> ULP drift that over
#   thousands of flow steps diverges chaotically. Validation is therefore
#   via INVARIANT tolerance (RL/Wr/ACN/K_tot), not bit-exact. The
#   swept pairs (dist<D, ring>excl) are as a SET exactly the same as
#   the cell list (exact space partition, not approximate); dist is
#   recomputed exactly for the kept pairs (same float formula as
#   CPU) -> identical semantics, only the order is free.
#
# Hook-up: cf.overlap_pairs = sixd_gpu.make_dispatcher(cf.overlap_pairs).
# Threshold GPU_MIN_N: below n=2048 the CPU cell list stays (the GPU overhead
# only pays off at large n). Set via env SIXD_GPU_MIN_N.
# ============================================================================
import os
import numpy as np

GPU_MIN_N = int(os.environ.get("SIXD_GPU_MIN_N", "2048"))

_cp = None
_OK = None
_RING_CACHE = {}     # n -> (n,n) int32 ring distance on the GPU


def _try_cupy():
    global _cp, _OK
    if _OK is not None:
        return _OK
    try:
        import cupy as cp
        # active smoke test (doctrine: probe the installed tool, do not guess)
        a = cp.arange(8.0)
        float((a * a).sum())
        _cp = cp
        _OK = True
    except Exception:
        _cp = None
        _OK = False
    return _OK


def available():
    return _try_cupy()


def device_name():
    if not _try_cupy():
        return None
    return _cp.cuda.runtime.getDeviceProperties(0)["name"].decode()


def _ring_gpu(n):
    cp = _cp
    rc = _RING_CACHE.get(n)
    if rc is None:
        idx = cp.arange(n, dtype=cp.int64)
        b = cp.abs(idx[:, None] - idx[None, :])
        rc = cp.minimum(b, n - b).astype(cp.int32)
        _RING_CACHE[n] = rc
    return rc


def overlap_pairs_gpu(x, n, excl, D):
    """Drop-in for canonical_flow.overlap_pairs — same return semantics
    (ii, jj, dist) for unordered pairs with dist<D and ring distance>excl,
    each pair once (i<j). brute-force O(n^2) Gram mask on the GPU."""
    cp = _cp
    xg = cp.asarray(x, dtype=cp.float64)              # (n,3)
    sq = (xg * xg).sum(1)                             # |x_i|^2
    G = xg @ xg.T                                     # (n,n) Gram
    D2 = sq[:, None] + sq[None, :] - 2.0 * G          # ~ pair distance^2
    rd = _ring_gpu(n)
    # upper triangle half: i<j ; ring>excl ; D2<D^2
    iu = cp.arange(n, dtype=cp.int64)
    mask = (D2 < D * D) & (rd > excl) & (iu[:, None] < iu[None, :])
    ii, jj = cp.where(mask)
    if ii.size == 0:
        return (np.empty(0, np.int64), np.empty(0, np.int64),
                np.empty(0, np.float64))
    # recompute dist exactly (same formula as CPU, no Gram cancellation)
    d = xg[ii] - xg[jj]
    dist = cp.sqrt((d * d).sum(1))
    keep = dist < D                                   # exact final filter
    ii, jj, dist = ii[keep], jj[keep], dist[keep]
    return (cp.asnumpy(ii).astype(np.int64),
            cp.asnumpy(jj).astype(np.int64),
            cp.asnumpy(dist))


def make_dispatcher(cpu_fn, min_n=None):
    """Returns an overlap_pairs function that uses the GPU from min_n and
    otherwise (or without cupy) falls back to the original CPU cell list."""
    thr = GPU_MIN_N if min_n is None else min_n

    def dispatch(x, n, excl, D):
        if n >= thr and _try_cupy():
            return overlap_pairs_gpu(x, n, excl, D)
        return cpu_fn(x, n, excl, D)

    return dispatch
