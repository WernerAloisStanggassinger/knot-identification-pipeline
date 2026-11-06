#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""embedding.py — real embedding test: SEGMENT-to-SEGMENT, not point-to-point.

WHY a dedicated tool: a point-sample distance is NOT an embedding proof.
If the found distance is below the sample spacing, a strand can pass between two
samples — the curve looks clean and is not
(Nyquist). Here every segment is computed exactly against every non-adjacent segment,
candidates presorted via a grid.

--selbsttest builds two counter-cases: a curve that SURELY intersects, and
one that surely does not. A test that does not bite is no test.
"""
from __future__ import annotations
import site_config
import sys
import numpy as np
from scipy.spatial import cKDTree


def _seg_seg(p1, q1, p2, q2):
    """Exact distance of two segments (Ericson, Real-Time Collision Detection)."""
    d1 = q1 - p1; d2 = q2 - p2; r = p1 - p2
    a = np.einsum('ij,ij->i', d1, d1)
    e = np.einsum('ij,ij->i', d2, d2)
    f = np.einsum('ij,ij->i', d2, r)
    b = np.einsum('ij,ij->i', d1, d2)
    c = np.einsum('ij,ij->i', d1, r)
    den = a * e - b * b
    s = np.where(den > 1e-30, np.clip((b * f - c * e) / np.where(den > 1e-30, den, 1), 0, 1), 0.0)
    t = (b * s + f) / np.where(e > 1e-30, e, 1)
    t_c = np.clip(t, 0, 1)
    s = np.clip((b * t_c - c) / np.where(a > 1e-30, a, 1), 0, 1)
    diff = (p1 + s[:, None] * d1) - (p2 + t_c[:, None] * d2)
    return np.linalg.norm(diff, axis=1)


def min_segmentabstand(P, aus_bogen_frac=0.02, M=6000):
    """Smallest distance of non-adjacent segments of the closed polyline P.

    EXHAUSTIVE and BOUNDED: the curve is thinned to M segments, then
    every segment is computed exactly, block by block, against every other. No search radius,
    no neighbor count — both have proven to be traps:
      * growing radius -> query_pairs returns O(n^2) pairs -> hang at load ~0
        (20 workers, 2026-08-22, all .fseries already written);
      * fixed k-nearest neighbors -> at fine sampling ALL k are arc neighbors
        -> result "inf" that looks like a proof.
    The neighbor exclusion counts ARC LENGTH, not indices: with index exclusion
    one finds the same strand a few segments further, and the number halves
    obediently with every doubling of N (measured, ratio constant 5.54).

    AUTHENTICITY CHECK: the number must be STABLE across N. If it is not,
    one measures the sampling, not the curve.
    """
    n = len(P)
    schritt_idx = max(1, n // M)
    Pd = P[::schritt_idx]
    m = len(Pd)
    A = Pd
    B = np.roll(Pd, -1, axis=0)
    laengen = np.linalg.norm(B - A, axis=1)
    radius = float(laengen.max())
    L = float(laengen.sum())
    aus_bogen = aus_bogen_frac * L
    bogen = np.concatenate([[0.0], np.cumsum(laengen)])[:m]

    best, wo = np.inf, None
    blk = 512
    for a0 in range(0, m, blk):
        a1 = min(m, a0 + blk)
        ii = np.arange(a0, a1)
        for b0 in range(a0, m, blk):
            b1 = min(m, b0 + blk)
            jj = np.arange(b0, b1)
            I, J = np.meshgrid(ii, jj, indexing="ij")
            I = I.ravel(); J = J.ravel()
            k = I < J
            I, J = I[k], J[k]
            if len(I) == 0:
                continue
            db = np.abs(bogen[I] - bogen[J])
            db = np.minimum(db, L - db)
            k = db > aus_bogen
            I, J = I[k], J[k]
            if len(I) == 0:
                continue
            d = _seg_seg(A[I], B[I], A[J], B[J])
            q = int(np.argmin(d))
            if d[q] < best:
                best, wo = float(d[q]), (int(I[q]), int(J[q]))
    return best, wo, radius


def selbsttest():
    t = np.linspace(0, 2 * np.pi, 4000, endpoint=False)
    kreis = np.stack([np.cos(t), np.sin(t), 0.02 * np.sin(3 * t)], axis=1)
    d1, wo1, h1 = min_segmentabstand(kreis)
    ok1 = d1 > 1e-3
    # Counter-case: a figure eight in the plane intersects itself at the origin
    acht = np.stack([np.sin(t), np.sin(t) * np.cos(t), np.zeros_like(t)], axis=1)
    d2, wo2, h2 = min_segmentabstand(acht)
    ok2 = d2 < 1e-6
    print(f"  Probe 1  round circle          d={d1:.6f}  (expected > 0)      {'✅' if ok1 else '🔴'}")
    print(f"  Probe 2  planar eight (crossing) d={d2:.3e}  (expected ~0)       {'✅' if ok2 else '🔴'}")
    # Counter-case 3: the OLD point test would consider the figure eight clean
    # if sampled coarsely enough — that is the reason for this tool.
    grob = acht[::37]
    dd = np.inf
    for a in range(len(grob)):
        for b in range(a + 3, len(grob) - 2):
            dd = min(dd, np.linalg.norm(grob[a] - grob[b]))
    print(f"  Probe 3  point test on coarse eight d={dd:.5f}  (>0 ⟹ it LIES) "
          f"{'✅ trap proven' if dd > 1e-6 else '🔴'}")
    return 0 if (ok1 and ok2 and dd > 1e-6) else 1


if __name__ == "__main__":
    if "--selbsttest" in sys.argv:
        sys.exit(selbsttest())
    sys.path.insert(0, (site_config.PIPELINE_DIR))
    from aggressive_simplify import parse_fseries, reconstruct_curve
    for pfad in [a for a in sys.argv[1:] if not a.startswith("-")]:
        off, ax, bx, ay, by, az, bz = parse_fseries(pfad)
        for N in (60000, 120000, 240000):
            P = reconstruct_curve(off, ax, bx, ay, by, az, bz, N).astype(np.float64)
            d, wo, h = min_segmentabstand(P)
            L = float(np.linalg.norm(np.diff(np.vstack([P, P[:1]]), axis=0), axis=1).sum())
            print(f"{pfad}  N={N:>7}  segment_length≈{h:.2e}  "
                  f"min_segment_distance={d:.6f}  d/L={d/L:.3e}  "
                  f"{'EMBEDDED' if d > 0 else 'PENETRATED'}")
