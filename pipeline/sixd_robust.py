#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Knot-Identification Pipeline — Werner Alois Stanggassinger & Claude Code (Anthropic).
# ============================================================================
# sixd_robust.py — order-ROBUST SONO flow (Werner's route C, Disc #263).
#
# PROBLEM (error log 2026-06-08): the flow is chaotic (positive Lyapunov);
#   np.add.at accumulates the overlap pushes in DISCOVERY order ->
#   float associativity -> ULP -> over thousands of iterations 5–40% drift +
#   branch flips. CPU (cell list) and GPU (Gram mask) deliver
#   the same pair SET, but a different ORDER -> a different canonical form.
#
# FIX (route C):
#   (1) The GPU finds ONLY the pair SET (the O(n^2) part). x and ALL
#       geometry (dist, push, smoothing, ironing) stay numpy -> backend-
#       independently bit-stable. dist is ALWAYS recomputed in numpy from x
#       (never the GPU dist return value -> no cupy sqrt/sum ULP).
#   (2) det_accumulate: scatter-add in FULLY CANONICAL order
#       (lexsort by (target, partner) + cumsum segment difference) instead of
#       np.add.at. Depends only on the pair SET, not on the discovery
#       order -> every backend (CPU/GPU) computes identically.
#
# CONSEQUENCE: CPU-robust == GPU-robust (bit-identical). The robust form differs
#   from the old add.at corpus (a different fixed order) -> recalibration needed,
#   but the form is now WELL-DEFINED + backend-portable + GPU-capable.
#   OPEN (deeper): whether the flow finds the TRUE ideal knot minimum (not
#   just a reproducible plateau representative) — tighter convergence
#   is a separate question (Disc #263 addendum).
#
# Hook-up: sixd_robust.install(backend='gpu'|'cpu'); then cf.run_record.
# ============================================================================
import sys
import site_config

import numpy as np

SCANNER_DIR = site_config.SIXD_DIR
sys.path.insert(0, SCANNER_DIR)
import canonical_flow as cf      # noqa: E402


def det_accumulate(target, partner, vals, n):
    """Order-INDEPENDENT scatter-add: out[target] += vals, summed in
    fully canonical order (lexsort target, then partner). Pure
    numpy -> bit-stable independent of the backend. vals: (m,3)."""
    out = np.zeros((n, 3))
    if target.size == 0:
        return out
    order = np.lexsort((partner, target))     # primary target, secondary partner
    t = target[order]
    v = vals[order]
    csum = np.cumsum(v, axis=0)
    starts = np.flatnonzero(np.concatenate(([True], t[1:] != t[:-1])))
    ends = np.concatenate((starts[1:], [t.size]))
    seg = csum[ends - 1].copy()
    seg[1:] -= csum[starts[1:] - 1]
    out[t[starts]] = seg
    return out


def wall_block_robust(x, n, h, D, circle=False):
    """cf.wall_block, but backend-independent: overlap_pairs delivers only the
    SET, dist/push/accumulation in numpy, det_accumulate instead of np.add.at."""
    excl = min(int(1.5 * D / h) + 2, cf.excl_cap(n, circle))
    cap = min(cf.DISP_CAP * D, 1.0 * h)
    ok, worked = False, False
    for _pass in range(cf.MAX_PASSES + 1):
        ii, jj, _ = cf.overlap_pairs(x, n, excl, D)   # only the pair SET
        if len(ii) == 0:
            ok = True
            break
        d = x[ii] - x[jj]                              # dist ALWAYS numpy from x
        dist = np.sqrt((d * d).sum(1))
        keep = dist < D                               # exact numpy filter
        if not keep.all():
            ii, jj, d, dist = ii[keep], jj[keep], d[keep], dist[keep]
        if ii.size == 0:
            ok = True
            break
        if _pass == cf.MAX_PASSES:
            ok = float(D - dist.min()) <= cf.TOL_DEPTH * D
            break
        worked = True
        u = d / np.maximum(dist, 1e-12)[:, None]
        push = 0.55 * (D - dist)[:, None] * u
        target = np.concatenate([ii, jj])
        partner = np.concatenate([jj, ii])
        vals = np.concatenate([push, -push])
        disp = det_accumulate(target, partner, vals, n)
        disp = cf.spread_disp(disp)
        dn = np.linalg.norm(disp, axis=1)
        scale = np.minimum(1.0, cap / np.maximum(dn, 1e-12))
        x = x + disp * scale[:, None]
        x, _ = cf.curvature_passes(x, n, h, D, max_passes=cf.CURV_MINI)
    x, curv_resid = cf.curvature_passes(x, n, h, D, max_passes=cf.CURV_CLEAN)
    return x, (ok and curv_resid <= cf.CURV_RESID_OK), worked, curv_resid


def install(backend="gpu"):
    """Installs the robust wall_block + selects the overlap backend.
    backend='gpu': GPU Gram mask (min_n=0, always GPU); 'cpu': original cell
    list. det_accumulate makes both routes bit-identical."""
    if not hasattr(cf, "_orig_overlap"):
        cf._orig_overlap = cf.overlap_pairs
    if backend == "gpu":
        import sixd_gpu
        if not sixd_gpu.available():
            raise SystemExit("sixd_robust: backend=gpu, but cupy not available")
        cf.overlap_pairs = sixd_gpu.make_dispatcher(cf._orig_overlap, min_n=0)
    else:
        cf.overlap_pairs = cf._orig_overlap
    cf.wall_block = wall_block_robust
    return backend
