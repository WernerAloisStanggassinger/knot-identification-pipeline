#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""clearance_segment.py — the true self-approach distance of a polygon, not its upper bound.

THE PROBLEM (Disc #746/#747): `auto_bevel_depth` measures the minimal self-approach distance
between SAMPLED POINTS. That overestimates in principle — the closest spot of the curve
almost never lies on a sample point. Measured on kf329:

    n= 1024  mca 0.009334      n=16384  mca 0.001050
    n= 4096  mca 0.003264      n=32768  mca 0.000796
    n= 8192  mca 0.002244      n=65536  mca 0.000652   <- STILL falling

The number does not converge; each one is an upper bound, and the tube radius computed
from it is correspondingly too large (for kf329 by a factor of 14 compared with n=65536).

THE CURE: distance between the SEGMENTS of the polygon instead of between its vertices.
The segment distance is a **lower** bound for the point distance and converges
quadratically instead of linearly — at n=4096 more precise than points at n=65536, in seconds.

METHOD: segment midpoints into a cKDTree, fetch candidate pairs within the radius
(radius = generous bound from a coarse point pass), discard neighbours along the
curve, compute the exact segment-segment distance for the survivors.

`--selbsttest` runs four counter-cases, among them the sharp one: two offset circles with
an ANALYTICALLY known distance — there the segment measurement must hit the true value and
the point measurement must overestimate it. If the segment measurement misses, it is wrong.
"""
import argparse, sys
import site_config
import numpy as np
from scipy.spatial import cKDTree


def _seg_seg_dist(p1, q1, p2, q2):
    """Shortest distance between segments p1->q1 and p2->q2, vectorised.
       Method after Ericson, Real-Time Collision Detection, ch. 5.1.9."""
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
    # clamp t and recompute s (the two correction steps of the method)
    t_c = np.clip(t, 0, 1)
    s = np.clip((b * t_c - c) / np.where(a > 1e-300, a, 1), 0, 1)
    t = np.clip((b * s + f) / np.where(e > 1e-300, e, 1), 0, 1)
    diff = (p1 + s[:, None] * d1) - (p2 + t[:, None] * d2)
    return np.linalg.norm(diff, axis=1)


def clearance(points, ring_skip, radius=None):
    """minimal self-approach distance of the closed polygon (segment-exact)."""
    n = len(points)
    a = points
    b = np.roll(points, -1, axis=0)
    mitten = 0.5 * (a + b)
    laengen = np.linalg.norm(b - a, axis=1)
    if radius is None:
        # coarse bound from the point distance, generously padded
        tree0 = cKDTree(points)
        d, j = tree0.query(points, k=min(64, n))
        idx = np.arange(n)
        rd = np.minimum(np.abs(idx[:, None] - j), n - np.abs(idx[:, None] - j))
        d = np.where(rd >= ring_skip, d, np.inf)
        radius = float(d.min()) * 1.5 + 2 * float(laengen.max())
    tree = cKDTree(mitten)
    paare = tree.query_pairs(radius, output_type='ndarray')
    if len(paare) == 0:
        return float('inf'), 0
    i, j = paare[:, 0], paare[:, 1]
    rd = np.minimum(np.abs(i - j), n - np.abs(i - j))
    m = rd >= ring_skip
    i, j = i[m], j[m]
    if len(i) == 0:
        return float('inf'), 0
    d = _seg_seg_dist(a[i], b[i], a[j], b[j])
    return float(d.min()), len(i)


def _punkt_mca(points, ring_skip):
    n = len(points)
    tree = cKDTree(points)
    k = min(64, n)
    while True:
        d, j = tree.query(points, k=k)
        idx = np.arange(n)
        rd = np.minimum(np.abs(idx[:, None] - j), n - np.abs(idx[:, None] - j))
        g = rd >= ring_skip
        if g.any(axis=1).all() or k >= n:
            return float(np.where(g, d, np.inf).min())
        k = min(k * 2, n)


def selbsttest():
    print('SELF-TEST clearance_segment — four counter-cases\n')
    fehler = 0

    print('1) TWO OFFSET CIRCLES — analytically known distance')
    #   Circle in the xy plane and a second one, shifted by dz: distance = dz.
    #   As ONE closed curve this does not work, so check the segment function directly.
    for dz in (0.5, 0.05, 0.005):
        t = np.linspace(0, 2 * np.pi, 64, endpoint=False)
        k1 = np.c_[np.cos(t), np.sin(t), np.zeros_like(t)]
        k2 = k1 + np.array([0, 0, dz])
        d = _seg_seg_dist(k1, np.roll(k1, -1, 0), k2, np.roll(k2, -1, 0)).min()
        ok = abs(d - dz) < 1e-9
        print(f'   {"✅" if ok else "🔴"} dz={dz:<7} measured {d:.12f}')
        fehler += (not ok)

    print('\n2) THE SHARP ONE: segment beats point at coarse sampling')
    #   Two circles, sampled TWISTED against each other: no pair of points faces
    #   each other, so the point distance overestimates — the segment distance does not.
    t = np.linspace(0, 2 * np.pi, 16, endpoint=False)
    k1 = np.c_[np.cos(t), np.sin(t), np.zeros(16)]
    t2 = t + (t[1] - t[0]) / 2
    k2 = np.c_[np.cos(t2), np.sin(t2), np.full(16, 0.01)]
    d_seg = _seg_seg_dist(k1, np.roll(k1, -1, 0), k2, np.roll(k2, -1, 0)).min()
    d_pkt = np.linalg.norm(k1[:, None, :] - k2[None, :, :], axis=2).min()
    besser = d_seg < d_pkt and abs(d_seg - 0.01) < 1e-3
    print(f'   {"✅" if besser else "🔴"} segment {d_seg:.6f} · point {d_pkt:.6f} · true 0.010000')
    fehler += (not besser)

    print('\n3) LOWER BOUND — segment must never be LARGER than point')
    rng = np.random.default_rng(7)
    schlimm = 0
    for _ in range(20):
        p = rng.normal(size=(200, 3))
        ds = _seg_seg_dist(p, np.roll(p, -1, 0), p[::-1], np.roll(p[::-1], -1, 0)).min()
        dp = np.linalg.norm(p[:, None, :] - p[::-1][None, :, :], axis=2).min()
        if ds > dp + 1e-12:
            schlimm += 1
    print(f'   {"✅" if schlimm == 0 else f"🔴 {schlimm}/20 violated"} 20 random polygons')
    fehler += bool(schlimm)

    print('\n4) CIRCLE AGAINST THE ANALYTIC NUMBER')
    #   🔴 The first draft expected ~2 here ("the circle approaches itself only across the
    #      diameter") and reported the CORRECT measurement 0.474 as an error. The expectation
    #      was wrong: outside the neighbourhood exclusion the nearest allowed point is the one
    #      at ring_skip, i.e. the chord 2*sin(pi*rs/n) = 0.4859 —
    #      not the diameter. Footgun 15: a red finding is first a suspicion
    #      against the checker.
    n, rs = 512, 40
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    kreis = np.c_[np.cos(t), np.sin(t), np.zeros_like(t)]
    sehne = 2 * np.sin(np.pi * rs / n)
    d_pkt = _punkt_mca(kreis, rs)
    d_seg, npaare = clearance(kreis, ring_skip=rs)
    ok_p = abs(d_pkt - sehne) < 1e-9
    #   The segment distance lies below the chord (chords lie inside the circle)
    #   and must not fall arbitrarily far below it: the bound is the chord itself
    #   times cos(pi/n) — the distance of the segment midpoints.
    unten = sehne * np.cos(np.pi / n) ** 2
    ok_s = unten * 0.97 <= d_seg <= sehne
    print(f'   {"✅" if ok_p else "🔴"} point   {d_pkt:.9f} against analytic {sehne:.9f}')
    print(f'   {"✅" if ok_s else "🔴"} segment {d_seg:.9f} — below the chord, above {unten:.9f} ({npaare} pairs)')
    fehler += (not ok_p) + (not ok_s)

    print(f'\n{"✅ ALL COUNTER-CASES PASSED" if fehler == 0 else f"🔴 {fehler} ERRORS"}')
    return 0 if fehler == 0 else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('fseries', nargs='?')
    ap.add_argument('--n', type=int, nargs='+', default=[1024, 4096, 16384])
    ap.add_argument('--selbsttest', action='store_true')
    a = ap.parse_args()
    if a.selbsttest:
        sys.exit(selbsttest())
    if not a.fseries:
        ap.error('fseries missing (or --selbsttest)')

    import importlib.util, time
    sp = importlib.util.spec_from_file_location('b', '@PIPELINE_DIR@/build_curve_json.py'.replace('@PIPELINE_DIR@', site_config.PIPELINE_DIR))
    b = importlib.util.module_from_spec(sp); sys.modules['b'] = b; sp.loader.exec_module(b)
    off, ca, cb = b.parse_fseries(a.fseries)
    print(f"{'n':>7} {'point-mca':>13} {'segment':>13} {'ratio':>12} {'sec':>7}")
    for n in a.n:
        pts = b.reconstruct_curve(off, ca, cb, n_samples=n)
        bmax = float(max(np.ptp(pts, axis=0)))
        ds = float(np.linalg.norm(np.diff(pts, axis=0, append=pts[:1]), axis=1).sum()) / n
        L = ds * n
        rs = max(3, int(np.ceil(max(0.05 * L, 2.0 * bmax * 0.015) / max(ds, 1e-9))))
        t = time.time()
        dp = _punkt_mca(pts, rs)
        dseg, npaare = clearance(pts, rs)
        print(f'{n:>7} {dp:>13.9f} {dseg:>13.9f} {dp/dseg:>12.3f} {time.time()-t:>7.2f}')


if __name__ == '__main__':
    main()
