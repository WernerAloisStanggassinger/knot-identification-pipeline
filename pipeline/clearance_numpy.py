#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""clearance_numpy.py — segment-exact self-approach distance, numpy ONLY.

Why this exists twice: `clearance_segment.py` needs `scipy.spatial.cKDTree`,
and scipy is not in the base env that `publish_site.sh` uses to run the site
generator. This version replaces the KD tree by a spatial grid (dict of
cells) and needs only numpy.

It is NOT the reference — the reference is `clearance_segment.py`, whose four
counter-cases hit analytically known distances. This version is calibrated
AGAINST that one (`--calibrate`), not against itself (Disc #745: a test that only
measures reproducibility never catches a construction error).

Why segments instead of points (Disc #746/#747): the distance between sample
POINTS overestimates, because the closest spot of the curve almost never lies on a
sample point, and it converges only linearly. The distance between the SEGMENTS is
a lower bound and converges quadratically.
"""
import argparse, sys
import site_config
import numpy as np


def _seg_seg_dist(p1, q1, p2, q2):
    """Shortest distance between segments p1->q1 and p2->q2, vectorised.
       Method after Ericson, Real-Time Collision Detection, ch. 5.1.9.
       Taken over bit-identically from clearance_segment.py."""
    d1 = q1 - p1
    d2 = q2 - p2
    r = p1 - p2
    a = np.einsum('ij,ij->i', d1, d1)
    e = np.einsum('ij,ij->i', d2, d2)
    f = np.einsum('ij,ij->i', d2, r)
    c = np.einsum('ij,ij->i', d1, r)
    b = np.einsum('ij,ij->i', d1, d2)
    denom = a * e - b * b
    s = np.where(denom > 1e-300, np.clip((b * f - c * e) / np.where(denom > 1e-300, denom, 1), 0, 1), 0.0)
    t = (b * s + f) / np.where(e > 1e-300, e, 1)
    t_c = np.clip(t, 0, 1)
    s = np.clip((b * t_c - c) / np.where(a > 1e-300, a, 1), 0, 1)
    t = np.clip((b * s + f) / np.where(e > 1e-300, e, 1), 0, 1)
    diff = (p1 + s[:, None] * d1) - (p2 + t[:, None] * d2)
    return np.linalg.norm(diff, axis=1)


def _point_mca_chunked(points, ring_skip, chunk=2048):
    """Point distance as a COARSE upper bound — chunked, so no n^2 memory.
       (The old version built the full matrix and died at n=32768 on 24 GB.)"""
    n = len(points)
    idx = np.arange(n)
    best = np.inf
    for i0 in range(0, n, chunk):
        i1 = min(i0 + chunk, n)
        d = np.linalg.norm(points[i0:i1, None, :] - points[None, :, :], axis=2)
        rd = np.abs(idx[i0:i1, None] - idx[None, :])
        rd = np.minimum(rd, n - rd)
        d[rd < ring_skip] = np.inf
        m = float(d.min())
        if m < best:
            best = m
    return best


def clearance(points, ring_skip, radius=None):
    """Minimal self-approach distance of the closed polygon, segment-exact.

    points     — (n,3) vertices of the closed polygon
    ring_skip  — how many neighbours along the curve count as 'adjacent'
    radius     — search radius; None => coarse point bound, taken generously
    """
    n = len(points)
    a = points.astype(float)
    b = np.roll(a, -1, axis=0)
    mitten = 0.5 * (a + b)
    lengths = np.linalg.norm(b - a, axis=1)

    if radius is None:
        radius = _point_mca_chunked(a, ring_skip)
    # Two segments can only come closer than `radius` if their
    # midpoints lie at most radius + (half length + half length) apart.
    # The longest edge serves as the bound.
    cell = float(radius + lengths.max())
    if not np.isfinite(cell) or cell <= 0:
        cell = float(np.linalg.norm(a.max(axis=0) - a.min(axis=0))) or 1.0

    # spatial grid: cell -> segment indices
    key_name = np.floor(mitten / cell).astype(np.int64)
    grid = {}
    for i, k in enumerate(map(tuple, key_name)):
        grid.setdefault(k, []).append(i)

    idx = np.arange(n)
    best = np.inf
    neighbours = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)]
    pairs_i, pairs_j = [], []
    for k, ilist in grid.items():
        for d in neighbours:
            nk = (k[0] + d[0], k[1] + d[1], k[2] + d[2])
            jlist = grid.get(nk)
            if not jlist:
                continue
            ii = np.repeat(np.asarray(ilist), len(jlist))
            jj = np.tile(np.asarray(jlist), len(ilist))
            m = ii < jj                      # each pair exactly once
            if not m.any():
                continue
            pairs_i.append(ii[m]); pairs_j.append(jj[m])
    if not pairs_i:
        return float('inf'), 0
    I = np.concatenate(pairs_i); J = np.concatenate(pairs_j)
    # discard neighbours ALONG the curve (the same rule as for the points)
    rd = np.abs(I - J)
    rd = np.minimum(rd, n - rd)
    m = rd >= ring_skip
    I, J = I[m], J[m]
    if len(I) == 0:
        return float('inf'), 0
    # compute in blocks, so that even millions of pairs do not blow the memory
    BL = 2_000_000
    for s in range(0, len(I), BL):
        i, j = I[s:s + BL], J[s:s + BL]
        d = _seg_seg_dist(a[i], b[i], a[j], b[j])
        mm = float(d.min())
        if mm < best:
            best = mm
    return best, len(I)


def _circle(n, r=1.0, dz=0.0, phase=0.0):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False) + phase
    return np.stack([r * np.cos(t), r * np.sin(t), np.full(n, dz)], axis=1)


def _brute(points, ring_skip):
    """All pairs, no grid — the truth the grid is checked against."""
    n = len(points)
    a = points.astype(float); b = np.roll(a, -1, axis=0)
    I, J = np.triu_indices(n, k=1)
    rd = np.minimum(np.abs(I - J), n - np.abs(I - J))
    m = rd >= ring_skip
    I, J = I[m], J[m]
    return float(_seg_seg_dist(a[I], b[I], a[J], b[J]).min())


def calibrate():
    """Against analytics, against brute force and against the scipy reference."""
    print('CALIBRATION clearance_numpy — five counter-cases\n')
    failures = 0

    print('1) THE SEGMENT FUNCTION AGAINST ANALYTICS — two offset circles, distance = dz')
    print('   (as ONE closed curve this does not work: the jump edges between the')
    print('    rings would come arbitrarily close to themselves. Hence the segment function directly.)')
    for dz in (0.5, 0.05, 0.005):
        k1 = _circle(64, 1.0, 0.0)
        k2 = k1 + np.array([0, 0, dz])
        d = float(_seg_seg_dist(k1, np.roll(k1, -1, 0), k2, np.roll(k2, -1, 0)).min())
        ok = abs(d - dz) < 1e-9
        print(f'   {"OK " if ok else "ERROR"} dz={dz:<8} measured {d:.12f}')
        failures += (not ok)

    print('\n2) THE SHARP ONE — segment beats point at coarse, twisted sampling')
    k1 = _circle(16, 1.0, 0.0)
    k2 = _circle(16, 1.0, 0.01, phase=np.pi / 16)
    d_seg = float(_seg_seg_dist(k1, np.roll(k1, -1, 0), k2, np.roll(k2, -1, 0)).min())
    d_pkt = float(np.linalg.norm(k1[:, None, :] - k2[None, :, :], axis=2).min())
    ok = d_seg < d_pkt and abs(d_seg - 0.01) < 1e-3
    print(f'   {"OK " if ok else "ERROR"} segment {d_seg:.6f} · point {d_pkt:.6f} · true 0.010000')
    failures += (not ok)

    print('\n3) THE GRID MUST NOT LOSE A PAIR — against brute force over ALL pairs')
    print('   (the real counter-case of this version: a cell that is too small would miss the')
    print('    closest pair and report a distance that is TOO LARGE — i.e. look safe.)')
    rng = np.random.default_rng(20260904)
    worst_k = 0
    for k in range(12):
        p = rng.normal(size=(400, 3))
        g, _ = clearance(p, ring_skip=3)
        bf = _brute(p, ring_skip=3)
        if abs(g - bf) > 1e-12:
            worst_k += 1
            print(f'   ERROR polygon {k}: grid {g:.12f} vs brute {bf:.12f}')
    print(f'   {"OK " if worst_k == 0 else "ERROR"} {12 - worst_k}/12 polygons bit-identical to brute force')
    failures += (worst_k != 0)

    print('\n4) LOWER BOUND — segment must never be LARGER than the point distance')
    worst_k = 0
    for _ in range(20):
        p = rng.normal(size=(300, 3))
        ds, _ = clearance(p, ring_skip=3)
        dp = _point_mca_chunked(p, ring_skip=3)
        if ds > dp + 1e-12:
            worst_k += 1
    print(f'   {"OK " if worst_k == 0 else "ERROR"} {20 - worst_k}/20 hold the bound')
    failures += (worst_k != 0)

    print('\n5) AGAINST THE scipy REFERENCE — the same number to 1e-12')
    try:
        sys.path.insert(0, site_config.TOOLS_DIR)
        from clearance_segment import clearance as ref_clearance
        worst_k = 0
        for k in range(5):
            p = rng.normal(size=(1500, 3))
            mine, _ = clearance(p, ring_skip=4)
            ref = ref_clearance(p, ring_skip=4)
            ref = ref[0] if isinstance(ref, tuple) else ref
            if abs(mine - ref) > 1e-12:
                worst_k += 1
                print(f'   ERROR polygon {k}: numpy {mine:.12f} · scipy {ref:.12f}')
        print(f'   {"OK " if worst_k == 0 else "ERROR"} {5 - worst_k}/5 polygons bit-identical to the scipy reference')
        failures += (worst_k != 0)
    except ImportError as e:
        print(f'   SKIPPED — scipy reference not importable ({e})')
        print('     (run in gnottheorie_env, scipy is there)')

    print('\n' + ('ALL COUNTER-CASES PASSED' if failures == 0 else f'{failures} ERRORS'))
    return failures


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--calibrate', action='store_true')
    args = ap.parse_args()
    if args.calibrate:
        sys.exit(1 if calibrate() else 0)
    ap.print_help()
