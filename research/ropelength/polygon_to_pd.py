#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""polygon_to_pd.py — PD code (KnotTheory convention X[i,j,k,l]) of a closed one-component space polygon.

Projection along a random direction (seeded), all crossings of non-adjacent edges, events sorted along the
curve, edge k = arc from event k to event k+1. At a crossing: i = incoming under edge, k = outgoing under edge,
j, l = over edges counter-clockwise after i.  Used to check the knot type of constructed and tightened curves.

Usage (in the sage_stable env):  python -I polygon_to_pd.py file.vect [seed]   -> prints crossings + Alexander polynomial
"""
import sys
import numpy as np


def rotation(seed):
    rng = np.random.default_rng(seed)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    return q


def pd_code(P, seed=1):
    P = P @ rotation(seed).T
    n = len(P)
    A, B = P, np.roll(P, -1, axis=0)
    ev = []                                   # (position along curve, crossing id, is_over)
    dirs = {}
    cid = 0
    a2, d2 = A[:, :2], (B - A)[:, :2]
    for i in range(n):
        j = np.arange(i + 2, n)
        if i == 0:
            j = j[j != n - 1]
        if len(j) == 0:
            continue
        r = a2[j] - a2[i]
        den = d2[i, 0] * d2[j, 1] - d2[i, 1] * d2[j, 0]
        ok = np.abs(den) > 1e-14
        s = np.where(ok, (r[:, 0] * d2[j, 1] - r[:, 1] * d2[j, 0]) / np.where(ok, den, 1), -1)
        t = np.where(ok, (r[:, 0] * d2[i, 1] - r[:, 1] * d2[i, 0]) / np.where(ok, den, 1), -1)
        hit = ok & (s > 0) & (s < 1) & (t > 0) & (t < 1)
        for jj, ss, tt in zip(j[hit], s[hit], t[hit]):
            zi = A[i, 2] + ss * (B[i, 2] - A[i, 2])
            zj = A[jj, 2] + tt * (B[jj, 2] - A[jj, 2])
            assert abs(zi - zj) > 1e-9, 'non-generic projection'
            over_i = zi > zj
            ev.append((i + ss, cid, over_i))
            ev.append((jj + tt, cid, not over_i))
            dirs[cid] = (d2[i], d2[jj], over_i)
            cid += 1
    ev.sort()
    N = len(ev)
    pos = {}                                  # (cid, is_over) -> event index
    for k, (_, c, ov) in enumerate(ev):
        pos[(c, ov)] = k
    pd = []
    for c in range(cid):
        di, dj, over_i = dirs[c]
        du, do = (dj, di) if over_i else (di, dj)
        a, b = pos[(c, False)], pos[(c, True)]
        i_in, k_out = (a - 1) % N, a
        b_in, b_out = (b - 1) % N, b
        cross = du[0] * do[1] - du[1] * do[0]   # > 0: over strand goes right -> left of the under strand
        j, l = (b_in, b_out) if cross > 0 else (b_out, b_in)
        pd.append((i_in, j, k_out, l))
    return pd


if __name__ == '__main__':
    import os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from ropelength_check import read_vect
    comps = read_vect(sys.argv[1])
    assert len(comps) == 1
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    pd = pd_code(comps[0], seed)
    import spherogram
    L = spherogram.Link(pd)
    print(sys.argv[1], 'seed', seed, 'projected crossings', len(pd))
    L.simplify('global')
    print('after simplify', len(L.crossings))
    print('alexander', L.alexander_polynomial())
