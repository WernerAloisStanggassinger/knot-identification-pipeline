#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
polygon_resolve.py — RECIPE #3 (error log 2026-06-04 ~22:40, Werner's go 2026-06-05)
=================================================================================
.fseries → polygon sampling → own projection → signed Gauss → stage-5 input.

Background (tr14/tr15e finding): the pipeline route via the raw scanner
diagram balloons to >28G on polyc07 curves (stage 4) or kills stage-5 regina
(rc=137 on a ~139-crossing Gauss code). The exact-layer measurement on the SAME
curve (polygon → signed Gauss → regina, 4k/8k points) ran with little RAM and a
bit-confirmed det. It is not the knot that is hard — the raw scanner-diagram route
is. With POLYGON_PD=1 this script replaces stages 1–4a and writes
a kinematic_result.json that stage 5+ consumes unchanged.

Method (proven in claude_setup_v2/tools/t67_harmonic_ladder.py, 2026-06-04):
  - sample the curve with ALL harmonics at N points
  - project along the FLATTEST axis first (shortest span axis =
    minimal-diagram doctrine, memory feedback-pd-view-along-flattest-axis),
    fall back to the remaining axes; right-handed frames (chirality!)
  - exact signed Gauss (segment intersections, chunked numpy, abort on tangency)
  - NYQUIST GATE (memory feedback-resolution-nyquist-triangulation): two
    consecutive resolutions N and 2N must be bit-identical in (det, Alexander),
    otherwise keep doubling. No single-resolution claim.
  - export: regina-simplified signed Gauss (stage-5 friendly) + raw Gauss
    as provenance.

Runs in the base env (regina + numpy). Pipeline blind-test compliant: seal_all()
active, no table lookups — everything from the .fseries itself.
"""

import argparse
import site_config
import json
import os
import sys
import time

sys.path.insert(0, site_config.FREEZE_DIR)
from pipeline_seal import seal_all  # noqa: E402
seal_all(verbose=False)

sys.path.insert(0, site_config.PIPELINE_DIR)
import numpy as np  # noqa: E402
import regina  # noqa: E402
from build_curve_json import parse_fseries  # noqa: E402

# Projecting out axis k = permuting coordinates so that k becomes the new
# z axis, RIGHT-HANDED (otherwise mirror drift in the chirality):
#   x out → (y, z, x) · y out → (z, x, y) · z out → (x, y, z)
AXIS_PERM = {0: (1, 2, 0), 1: (2, 0, 1), 2: (0, 1, 2)}
AXIS_NAME = {0: 'x', 1: 'y', 2: 'z'}


def sample_curve(off, A, B, n, phase=0.0):
    """Full curve (all harmonics) at n points; phase = deterministic
    jitter against exactly symmetric tangency degenerations."""
    t = np.linspace(0, 2 * np.pi, n, endpoint=False) + phase
    P = np.tile(np.asarray(off, float), (n, 1))
    for k in range(A.shape[0]):
        P += np.outer(np.cos((k + 1) * t), A[k]) + np.outer(np.sin((k + 1) * t), B[k])
    return P


def signed_gauss(P):
    """Projection along z; segment intersections (chunked numpy); O/U via z,
    sign via crossing handedness. Returns: regina signed-Gauss string,
    or None on suspected tangency (intersection too close to a segment end).
    Identical to the exact layer t67_harmonic_ladder.py (measured 2026-06-04)."""
    n = len(P)
    S = P[:, :2]; E = np.roll(S, -1, axis=0)            # 2D segments i: S[i]->E[i]
    D = E - S
    passages = []                                        # (param, cid, over?, sign)
    cid = 0
    CH = 512
    for i0 in range(0, n, CH):
        i = np.arange(i0, min(i0 + CH, n))
        for j0 in range(i0, n, CH):
            j = np.arange(j0, min(j0 + CH, n))
            I, J = np.meshgrid(i, j, indexing="ij")
            mask = (J > I + 1) & ~((I == 0) & (J == n - 1))
            if not mask.any(): continue
            Ii, Jj = I[mask], J[mask]
            d = (D[Ii, 0] * D[Jj, 1] - D[Ii, 1] * D[Jj, 0])
            ok = np.abs(d) > 1e-15
            Ii, Jj, d = Ii[ok], Jj[ok], d[ok]
            R = S[Jj] - S[Ii]
            t_ = (R[:, 0] * D[Jj, 1] - R[:, 1] * D[Jj, 0]) / d
            u_ = (R[:, 0] * D[Ii, 1] - R[:, 1] * D[Ii, 0]) / d
            hit = (t_ > 0) & (t_ < 1) & (u_ > 0) & (u_ < 1)
            for ii, jj, tt, uu in zip(Ii[hit], Jj[hit], t_[hit], u_[hit]):
                if min(tt, 1 - tt, uu, 1 - uu) < 1e-9:
                    return None                          # tangency/edge case -> increase resolution
                z1 = P[ii, 2] + tt * (P[(ii + 1) % n, 2] - P[ii, 2])
                z2 = P[jj, 2] + uu * (P[(jj + 1) % n, 2] - P[jj, 2])
                if abs(z1 - z2) < 1e-12:
                    return None
                od, ud = (ii, jj) if z1 > z2 else (jj, ii)
                s = 1 if (D[od, 0] * D[ud, 1] - D[od, 1] * D[ud, 0]) > 0 else -1
                cid += 1
                passages.append((ii + tt, cid, z1 > z2, s))
                passages.append((jj + uu, cid, z2 > z1, s))
    passages.sort()
    toks = [f"{'O' if over else 'U'}{c}{'+' if s > 0 else '-'}" for _, c, over, s in passages]
    return "".join(toks)


def measure(off, A, B, axis, n, phase=0.0):
    """One (axis, resolution) measurement. Returns a dict or None (tangency /
    non-classical / empty diagram)."""
    P = sample_curve(off, A, B, n, phase)[:, AXIS_PERM[axis]]
    sg = signed_gauss(P)
    if not sg:                                           # None = tangency, "" = 0 crossings
        return None
    L = regina.Link.fromSignedGauss(sg)
    if not L.isClassical():
        return None
    n_raw = int(L.size())
    L.simplify()
    p = L.alexander()
    det = abs(sum(int(str(p[k])) * (-1) ** k for k in range(p.degree() + 1)))
    return {
        'signed_gauss_raw': sg,
        'n_cross_raw': n_raw,
        'signed_gauss_simplified': str(L.signedGauss()),
        'n_cross_simplified': int(L.size()),
        'det': det,
        'alexander': str(p),
    }


def main():
    ap = argparse.ArgumentParser(description='Recipe #3: polygon diagram straight from the .fseries for stage 5.')
    ap.add_argument('fseries')
    ap.add_argument('--name', default='knot')
    ap.add_argument('--output', required=True, help='kinematic_result.json-compatible output')
    ap.add_argument('--axis', default='auto', choices=['auto', 'x', 'y', 'z'],
                    help='Projection axis; auto = flattest first, the rest as fallback')
    ap.add_argument('--n-start', type=int, default=4096)
    ap.add_argument('--n-max', type=int, default=32768)
    args = ap.parse_args()

    t_total = time.time()
    off, A, B = parse_fseries(args.fseries)
    A = np.asarray(A, float); B = np.asarray(B, float)
    print(f"[POLYGON] Name: {args.name}")
    print(f"[POLYGON] Source: {args.fseries}  ({A.shape[0]} harmonics)")

    # Axis order: flattest (shortest span) first (PD doctrine)
    probe = sample_curve(off, A, B, 2048)
    spans = probe.max(axis=0) - probe.min(axis=0)
    if args.axis == 'auto':
        axes = list(np.argsort(spans))
    else:
        a0 = 'xyz'.index(args.axis)
        axes = [a0] + [a for a in np.argsort(spans) if a != a0]
    print(f"[POLYGON] Spans (x,y,z): {np.round(spans, 4).tolist()} → axis order: "
          f"{[AXIS_NAME[a] for a in axes]}")

    accepted = None
    attempts = []
    for axis in axes:
        prev = None
        n = args.n_start
        attempt_phase = 0
        while n <= args.n_max:
            t0 = time.time()
            r = measure(off, A, B, axis, n, phase=attempt_phase * 1e-4)
            dt = time.time() - t0
            if r is None:
                print(f"[POLYGON]   {AXIS_NAME[axis]}@{n}: tangent/degenerate ({dt:.1f}s) → doubling")
                attempts.append({'axis': AXIS_NAME[axis], 'n': n, 'result': 'tangent/degenerate'})
                prev = None
                n *= 2
                attempt_phase += 1
                continue
            print(f"[POLYGON]   {AXIS_NAME[axis]}@{n}: raw={r['n_cross_raw']} cr, "
                  f"simplified={r['n_cross_simplified']} cr, det={r['det']} ({dt:.1f}s)")
            attempts.append({'axis': AXIS_NAME[axis], 'n': n,
                             'det': r['det'], 'n_cross_raw': r['n_cross_raw'],
                             'n_cross_simplified': r['n_cross_simplified']})
            if prev is not None and prev['det'] == r['det'] and prev['alexander'] == r['alexander']:
                accepted = (axis, prev['n'], n, r)
                print(f"[POLYGON]   ✓ NYQUIST GATE: (det, Alexander) bit-identical at "
                      f"{prev['n']} and {n} → accepted")
                break
            if prev is not None:
                print(f"[POLYGON]   ✗ Nyquist mismatch {prev['n']}↔{n} "
                      f"(det {prev['det']}≠{r['det']} or alex) → doubling")
            prev = dict(r, n=n)
            n *= 2
        if accepted:
            break
        print(f"[POLYGON]   axis {AXIS_NAME[axis]} up to n-max={args.n_max} without "
              f"Nyquist confirmation → next axis")

    if not accepted:
        print(f"[POLYGON] ERROR: no axis yields a Nyquist-stable classical "
              f"diagram up to n-max={args.n_max}. Exit 4.", file=sys.stderr)
        sys.exit(4)

    axis, n1, n2, r = accepted
    output = {
        'fseries_path': args.fseries,
        'polygon_pd': True,                              # recipe #3 marker (provenance)
        'method': ('Rezept #3: Polygon-Sampling -> eigene Projektion (flachste Achse, '
                   'rechtshändiger Frame) -> signed Gauss exakt -> regina; '
                   'Nyquist-Gate (det+Alexander bit-gleich bei zwei Auflösungen)'),
        'chosen': {
            'signed_gauss': r['signed_gauss_simplified'],
            'signed_gauss_raw': r['signed_gauss_raw'],
            'projection_axis': AXIS_NAME[axis],
            'n_samples': n2,
            'n_cross_in_frame': r['n_cross_raw'],
            'n_cross_simplified': r['n_cross_simplified'],
        },
        'nyquist': {'n_low': n1, 'n_high': n2, 'det': r['det'], 'alexander': r['alexander']},
        'attempts': attempts,
        'duration_seconds': time.time() - t_total,
    }
    with open(args.output, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"[POLYGON] det={r['det']}  alexander={r['alexander'][:80]}")
    print(f"[POLYGON] Writing {args.output}")
    print(f"[POLYGON] Done in {time.time()-t_total:.1f}s.")


if __name__ == '__main__':
    main()
