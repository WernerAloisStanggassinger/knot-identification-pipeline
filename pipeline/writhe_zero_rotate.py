#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# writhe_zero_rotate.py — prove that the zero line is realisation, not knot:
# take irrational (T(3,7), Wr=+17.25) and coil a sub-arc with opposite handedness until Wr=0.
# Coiling = local isotopy (twisting a rope) -> knot type stays T(3,7), writhe wanders.
import site_config
import sys, json, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from helicity_spectrum import parse_fseries, reconstruct, writhe

FS = (site_config.ROOT + '/own_knots/irrational/irrational.fseries')


def frame_parallel(Q):
    """Parallel-transport frame (u1,u2) perpendicular to the tangent along Q."""
    n = len(Q)
    T = np.gradient(Q, axis=0); T /= (np.linalg.norm(T, axis=1, keepdims=True) + 1e-12)
    u1 = np.zeros_like(Q)
    # Start normal
    a = np.array([1.0, 0, 0])
    if abs(np.dot(a, T[0])) > 0.9:
        a = np.array([0, 1.0, 0])
    u1[0] = a - np.dot(a, T[0]) * T[0]; u1[0] /= np.linalg.norm(u1[0])
    for i in range(1, n):
        v = u1[i-1] - np.dot(u1[i-1], T[i]) * T[i]
        nv = np.linalg.norm(v)
        u1[i] = v/nv if nv > 1e-9 else u1[i-1]
    u2 = np.cross(T, u1)
    return u1, u2


def coil_subarc(C, i0, i1, n_turns, amp, sign):
    """Coil the arc C[i0:i1] into n_turns helix turns, amplitude amp (windowed to 0 at the ends)."""
    C = C.copy()
    Q = C[i0:i1].copy()
    m = len(Q)
    u1, u2 = frame_parallel(Q)
    k = np.arange(m)
    env = np.sin(np.pi * k / (m - 1)) ** 2            # 0 at both ends
    phi = sign * 2 * np.pi * n_turns * (k / (m - 1))
    off = (amp * env)[:, None] * (np.cos(phi)[:, None] * u1 + np.sin(phi)[:, None] * u2)
    C[i0:i1] = Q + off
    return C


def pigtail(C, j, n_turns, sign, m_out=1600, m_ret=800):
    """Insert at C[j] a pigtail into the empty space: radially out as a helix
    (n_turns turns, opposite handedness sign), straight back along the axis.
    In empty space -> cannot hit the knot body -> T(3,7) stays; writhe ~ sign*n_turns."""
    ext = np.abs(C).max()
    p = C[j].copy()
    u = p / (np.linalg.norm(p) + 1e-12)              # radially outward
    e1 = np.cross(u, [0, 0, 1.0]);  e1 /= (np.linalg.norm(e1) + 1e-12)
    e2 = np.cross(u, e1)
    D = 1.6 * ext                                    # far out, beyond the knot
    R = 0.35 * ext                                   # coil radius, large in the empty space
    s = np.linspace(0, 1, m_out)
    env = np.sin(np.pi * s)                          # R=0 at p and at the turning point
    phi = sign * 2 * np.pi * n_turns * s
    axis = p[None, :] + (D * s)[:, None] * u[None, :]
    helix = axis + (R * env)[:, None] * (np.cos(phi)[:, None] * e1[None, :] + np.sin(phi)[:, None] * e2[None, :])
    far = p + D * u
    ret = np.linspace(far, p, m_ret, endpoint=False)  # straight back along the axis
    return np.vstack([C[:j+1], helix, ret, C[j+1:]])


def main():
    import sys as _s
    _s.path.insert(0, (site_config.MACHINE_YARD + '/the_knotfamily_machine'))
    from fseries_io import polyline_to_fourier, write_fseries

    a, b = parse_fseries(FS)
    C0 = reconstruct(a, b, 6000)
    Wr0 = writhe(C0)
    j = int(np.argmax(np.linalg.norm(C0, axis=1)))
    print("Wr0 (irrational, served) = %+.4f   Extent=%.3f" % (Wr0, np.abs(C0).max()), flush=True)

    # Bisection: sign=+1 subtracts ~1.1/turn -> Wr(N)=0 near N~16
    f = lambda N: writhe(pigtail(C0, j, N, sign=+1))
    lo, hi = 5.0, 30.0
    flo, fhi = f(lo), f(hi)
    print("Bisection  f(%.0f)=%+.3f  f(%.0f)=%+.3f" % (lo, flo, hi, fhi), flush=True)
    for _ in range(28):
        mid = 0.5*(lo+hi); fm = f(mid)
        if (fm > 0) == (flo > 0): lo, flo = mid, fm
        else: hi, fhi = mid, fm
    Nstar = 0.5*(lo+hi)
    C1 = pigtail(C0, j, Nstar, sign=+1)
    Wr1 = writhe(C1)
    print("N* = %.4f turns  ->  Wr1 (Polyline) = %+.5f" % (Nstar, Wr1), flush=True)

    # fseries round trip: enough harmonics so that the coil (and thus Wr=0) survives
    print("\nfseries round trip (Wr after reconstruction from N terms):", flush=True)
    best = None
    for nt in [200, 400, 600, 900]:
        aa, bb = polyline_to_fourier(C1, n_terms=nt, n_samples=max(65536, 40*nt), center=True)
        Cr = reconstruct(aa, bb, 20000)
        wr = writhe(Cr)
        print("  n_terms=%-4d  Wr=%+8.4f" % (nt, wr), flush=True)
        if best is None or abs(wr) < abs(best[2]): best = (nt, (aa, bb), wr)
    nt, (aa, bb), wr = best
    out = (site_config.ROOT + '/own_knots/irrational/irrational_wr0.fseries')
    write_fseries(out, aa, bb, knot_name="irrational_wr0",
                  source="irrational (T(3,7)) turned to Wr=0 with a pigtail, N*=%.3f turns, n_terms=%d" % (Nstar, nt))
    print("\nwritten: %s  (n_terms=%d, Wr=%+.4f)" % (out, nt, wr), flush=True)
    print("COMPARISON:  irrational served Wr=%+.3f   ->   irrational_wr0 Wr=%+.3f   (same knot T(3,7)?)" % (Wr0, wr), flush=True)


if __name__ == "__main__":
    main()
