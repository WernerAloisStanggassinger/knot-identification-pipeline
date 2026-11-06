# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 126 — widening_the_loom, LOOM-6 (six-lane loom).

Genuine broadening of the Loom-4 geometry to 6 strands: the same
formation law, two additional lanes + two return rails.
  LANES[i] = ((N-1)/2 - i) * S_LANE     (symmetric about 0, step S_LANE)
  RAILS[i] = 1.35 + 0.45*i              (staggered return rails)
Loom-4 is the special case N=4 (bit-compatible: same constants).

Word-agnostic like Loom-4 (build_station takes any braid word as
[(generator_index, sign), ...], generators σ_1..σ_5). build_curve checks
that the closure is a 6-cycle (= a knot, not a link).

PURPOSE 2026-06-09: make the fourth named θ-blind deficit pair 11n_73/74
(6 strands, 25-letter braids) buildable — Loom-4 was not enough.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
YARD = HERE.parent
sys.path.insert(0, str(YARD))

# --- Loom-N geometry (117 scale, generalised) -------------------------
N_STR  = 6
S_LANE = 0.55
H_Z    = 0.22
X_STN  = 2.4
X_RUN  = 3.0
BULGE  = 0.6
NW     = 25

LANES = tuple(((N_STR - 1) / 2.0 - i) * S_LANE for i in range(N_STR))
RAILS = tuple(1.35 + 0.45 * i for i in range(N_STR))


def word_perm(word, n_str=N_STR):
    p = list(range(n_str))
    for g, _s in word:
        p[g - 1], p[g] = p[g], p[g - 1]
    return p


def is_homogeneous(word):
    signs = {}
    for g, s in word:
        if g in signs and signs[g] != s:
            return False
        signs[g] = s
    return True


def build_station(word, x_stn=X_STN, x_run=X_RUN):
    """N strand paths through the fixed station (117 logic, occupant
    generalised: uninvolved lanes = loop)."""
    n_c = len(word)
    w = 2 * x_stn / n_c
    occupant = list(range(N_STR))
    paths = {s: [] for s in range(N_STR)}
    for s in range(N_STR):
        for x in np.linspace(-x_run, -x_stn, 8, endpoint=False):
            paths[s].append((x, LANES[s], 0.0))
    x0 = -x_stn
    ts = np.linspace(0.0, 1.0, NW, endpoint=False)
    for gi, sign in word:
        a, b = gi - 1, gi
        others = [g for g in range(N_STR) if g not in (a, b)]
        sa, sb = occupant[a], occupant[b]
        for t in ts:
            x = x0 + w * t
            ramp = (1 - np.cos(np.pi * t)) / 2
            ya = LANES[a] + (LANES[b] - LANES[a]) * ramp
            yb = LANES[b] + (LANES[a] - LANES[b]) * ramp
            z = H_Z * np.sin(np.pi * t)
            paths[sa].append((x, ya, +z if sign > 0 else -z))
            paths[sb].append((x, yb, -z if sign > 0 else +z))
            for g in others:
                paths[occupant[g]].append((x, LANES[g], 0.0))
        occupant[a], occupant[b] = sb, sa
        x0 += w
    for lane in range(N_STR):
        s = occupant[lane]
        for x in np.linspace(x_stn, x_run, 8, endpoint=False):
            paths[s].append((x, LANES[lane], 0.0))
    end_lane = {occupant[lane]: lane for lane in range(N_STR)}
    return paths, end_lane


def connector(lane: int, x_run=X_RUN):
    """Return rail lane `lane` (117 logic, N rails)."""
    y, R = LANES[lane], RAILS[lane]
    pts = []
    for t in np.linspace(0, 1, 24, endpoint=False)[1:]:
        pts.append((x_run + BULGE * np.sin(np.pi * t), y,
                    R * (1 - np.cos(np.pi * t)) / 2))
    for x in np.linspace(x_run, -x_run, 48, endpoint=False)[1:]:
        pts.append((x, y, R))
    for t in np.linspace(0, 1, 24, endpoint=False)[1:]:
        pts.append((-x_run - BULGE * np.sin(np.pi * t), y,
                    R * (1 + np.cos(np.pi * t)) / 2))
    return pts


def station_crossings(P: np.ndarray, x_stn=X_STN) -> int:
    """Crossings of the z projection, station only (|x|≤x_stn, |z|≤0.5)."""
    n = len(P)
    nxt = (np.arange(n) + 1) % n
    inside = (np.abs(P[:, 0]) <= x_stn + 1e-9) & (np.abs(P[:, 2]) <= 0.5)
    segs = np.where(inside & inside[nxt])[0]
    A, B = P[segs][:, :2], P[nxt[segs]][:, :2]
    cnt = 0
    for ii in range(len(segs)):
        jj = np.arange(ii + 1, len(segs))
        adj = np.abs(segs[jj] - segs[ii]) <= 1
        p, r = A[ii], B[ii] - A[ii]
        q, s = A[jj], B[jj] - A[jj]
        denom = r[0] * s[:, 1] - r[1] * s[:, 0]
        dq = q - p
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (dq[:, 0] * s[:, 1] - dq[:, 1] * s[:, 0]) / denom
            u = (dq[:, 0] * r[1] - dq[:, 1] * r[0]) / denom
        hit = (~adj) & (np.abs(denom) > 1e-12) \
            & (t > 1e-9) & (t < 1 - 1e-9) & (u > 1e-9) & (u < 1 - 1e-9)
        cnt += int(hit.sum())
    return cnt


def outside_station(P: np.ndarray, x_stn=X_STN) -> np.ndarray:
    m = (np.abs(P[:, 0]) > x_stn + 1e-9) | (np.abs(P[:, 2]) > 0.5)
    return P[m]


def build_curve(word):
    """word = [(gi, sign), ...]; closure must be an N_STR cycle."""
    paths, end_lane = build_station(word)
    pts, cur, visited = [], 0, []
    for _ in range(N_STR):
        visited.append(cur)
        pts += paths[cur]
        lane = end_lane[cur]
        pts += connector(lane)
        cur = lane
    if cur != 0 or sorted(visited) != list(range(N_STR)):
        raise RuntimeError("permutation not a 6-cycle — closure not a knot")
    return np.asarray(pts, float)
