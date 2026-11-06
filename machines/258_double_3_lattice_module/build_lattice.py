#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""build_lattice.py — Werner's Double-3 lattice module, closed at last.

Construction (Werner 18.08.2026): Double-3 is a connecting piece with docking
turning points. At k of them a MIRRORED and ENLARGED Double-3 is
spliced in each. Result is ONE closed loop.

Doctrine that counts here:
  * Source is own_knots/double-3/double-3.fseries — NOT the
    double3_module_dense*.npy from the phase-C era (Disc #251/#252: there the
    algebra of a dense module object was stapled onto a foreign curve).
  * NO number travels from the construction into a record. Only a .fseries
    is output; it is identified anew by the pipeline.
  * "no self intersection" is Werner's CLAIM. Here it is checked —
    and on the curve reconstructed from the .fseries, not on the
    build polyline (the record object is the served curve).
"""
from __future__ import annotations
import site_config
import argparse, json, math, sys
import numpy as np

sys.path.insert(0, (site_config.PIPELINE_DIR))
sys.path.insert(0, (site_config.MACHINE_YARD + '/the_knotfamily_machine'))
from aggressive_simplify import parse_fseries, reconstruct_curve
from fseries_io import polyline_to_fourier, write_fseries, resample_uniform

FS = (site_config.ROOT + '/own_knots/double-3/double-3.fseries')


# ---------------------------------------------------------------- Load module
def modul(n=40000):
    off, ax, bx, ay, by, az, bz = parse_fseries(FS)
    P = reconstruct_curve(off, ax, bx, ay, by, az, bz, n).astype(np.float64)
    C = P - P.mean(0)
    T = C.T @ C / len(C)
    w, V = np.linalg.eigh(T)
    V = V[:, np.argsort(w)[::-1]]
    e3 = V[:, 2] / np.linalg.norm(V[:, 2])
    e1 = V[:, 0] - (V[:, 0] @ e3) * e3
    e1 /= np.linalg.norm(e1)
    Q = np.stack([e1, np.cross(e3, e1), e3], axis=1)
    return C @ Q                              # curve in its own frame, centred


def aussenecken(M):
    """The six outer turning points: local radius maxima, the top six."""
    r = np.hypot(M[:, 0], M[:, 1])
    n = len(r)
    kand = [i for i in range(n) if r[i] >= r[i-1] and r[i] > r[(i+1) % n]]
    # merge plateaus
    grp, cur = [], [kand[0]]
    for i in kand[1:]:
        (cur if i - cur[-1] <= n // 400 else grp.append(cur) or (cur := [i]))
        if i - cur[-1] > n // 400:
            pass
    grp2, cur = [], [kand[0]]
    for i in kand[1:]:
        if i - cur[-1] <= n // 400:
            cur.append(i)
        else:
            grp2.append(cur); cur = [i]
    grp2.append(cur)
    if len(grp2) > 1 and (grp2[0][0] + n - grp2[-1][-1]) <= n // 400:
        grp2[0] = grp2[-1] + grp2[0]; grp2.pop()
    mitten = [g[int(np.argmax(r[g]))] for g in grp2]
    mitten.sort(key=lambda i: -r[i])
    sechs = sorted(mitten[:6])
    return sechs


def rot_z(w):
    c, s = math.cos(w), math.sin(w)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])


# ---------------------------------------------------------------- Build

def modul_lage(k=3, skala=1.35, spalt=0.10, nmod=40000, richtung="aussen",
               phase=0.0, kipp=0.0, hub=0.0):
    """Centre + k attached modules as CLOSED curves, BEFORE splicing.
    Exactly the placement that baue() then cuts open — so that the linking number measures
    the same object that is built, and not a reconstruction of it."""
    M = modul(nmod)
    ecken = aussenecken(M)
    schritt = max(1, len(ecken) // k)
    gewaehlt = sorted(ecken[(i * schritt) % len(ecken)] for i in range(k))
    spiegel = np.diag([1.0, 1.0, -1.0])
    anhaenge = []
    for i in gewaehlt:
        A = (M @ spiegel) * skala
        p_i = M[i]
        u = np.array([p_i[0], p_i[1], 0.0]); u /= np.linalg.norm(u)
        # DIRECTION: where does the body of the attached module fall?
        #   "aussen" — radial of the attachment points inward, body grows
        #      OUTWARD. When turning up it moves away from the central body: at k=3
        #      measured up to scale 14 NEVER interlocked with the center (2026-08-22).
        #   "ueber"  — radial points outward, body falls BACK OVER the
        #      center. That is Werner's picture: "so that the INNER bodies
        #      are penetrated" — the larger module is turned over the
        #      smaller one. With "aussen" that is geometrically impossible.
        # CORNER (formerly wrongly "phase"): WHICH of the six outer corners of the
        # attached module docks. A rotation angle about z was a NO-OP here —
        # the subsequent alignment to the radial of the docking corner undoes it
        # exactly. The search grid gave it away: four phase values,
        # four times bit-identical numbers (2026-08-22). The control is discrete.
        ecken_A = aussenecken(A)
        c = ecken_A[int(phase) % len(ecken_A)]
        q = A[c]
        v = np.array([q[0], q[1], 0.0]); v /= np.linalg.norm(v)
        ziel_radial = -u if richtung == "aussen" else u
        w = math.atan2(ziel_radial[1], ziel_radial[0]) - math.atan2(v[1], v[0])
        A = A @ rot_z(w).T
        A = A + ((p_i + spalt * u) - A[c])
        # TILT: about the radial axis THROUGH the docking site. Breaks the
        # coplanarity that lays the strands against each other.
        if kipp:
            achse = u
            K = np.cos(kipp) * np.eye(3) + np.sin(kipp) * np.array(
                [[0, -achse[2], achse[1]], [achse[2], 0, -achse[0]], [-achse[1], achse[0], 0]]
            ) + (1 - np.cos(kipp)) * np.outer(achse, achse)
            A = (A - A[c]) @ K.T + A[c]
        # HUB: offset perpendicularly. The cheapest way to open a narrow passage.
        if hub:
            A = A + np.array([0.0, 0.0, hub])
        anhaenge.append(A)
    return M, anhaenge, gewaehlt


def baue(k=3, skala=1.35, spalt=0.10, fenster=0.012, nmod=40000, richtung="aussen",
         phase=0.0, kipp=0.0, hub=0.0):
    """Splice in k modules. skala>1 = attached ones larger (Werner's control).

    richtung: "aussen" | "ueber"  — see modul_lage(). The placement lives
    ONLY there; baue() does not rebuild it. (It once stood here twice, and
    the direction correction of 22.08. therefore at first acted only half.)
    """
    M, anhaenge, gewaehlt = modul_lage(k, skala, spalt, nmod, richtung, phase, kipp, hub)
    n = len(M)
    w = max(2, int(fenster * n))

    schnitte = [(i - w, i + w) for i in gewaehlt]
    boegen = []
    for j in range(k):
        start = (schnitte[j][1]) % n
        ende  = (schnitte[(j + 1) % k][0]) % n
        idx = np.arange(start, ende + n if ende < start else ende) % n
        boegen.append(M[idx])

    anhang = []
    for A in anhaenge:
        c = aussenecken(A)[0]
        wA = max(2, int(fenster * len(A)))
        s2, e2 = (c + wA) % len(A), (c - wA) % len(A)
        idx = np.arange(s2, e2 + len(A) if e2 < s2 else e2) % len(A)
        anhang.append(A[idx])

    # --- Bridges (Hermite) + assembly ---------------------------------------
    def bruecke(p0, t0, p1, t1, m=140):
        d = np.linalg.norm(p1 - p0)
        s = np.linspace(0, 1, m)[:, None]
        h00 = 2*s**3 - 3*s**2 + 1; h10 = s**3 - 2*s**2 + s
        h01 = -2*s**3 + 3*s**2;    h11 = s**3 - s**2
        return h00*p0 + h10*(d*t0) + h01*p1 + h11*(d*t1)

    teile = []
    for j in range(k):
        B = boegen[j]; A = anhang[j]
        teile.append(B)
        teile.append(bruecke(B[-1], (B[-1]-B[-2])/np.linalg.norm(B[-1]-B[-2]),
                             A[0],  (A[1]-A[0])/np.linalg.norm(A[1]-A[0])))
        teile.append(A)
        Bn = boegen[(j + 1) % k]
        teile.append(bruecke(A[-1], (A[-1]-A[-2])/np.linalg.norm(A[-1]-A[-2]),
                             Bn[0], (Bn[1]-Bn[0])/np.linalg.norm(Bn[1]-Bn[0])))
    return np.vstack(teile), gewaehlt


# ---------------------------------------------------------------- Embedding
def min_abstand(P, ausschluss_frac=0.01):
    """Smallest distance between NON-adjacent segments (sample)."""
    n = len(P)
    m = max(1, n // 4000)
    Q = P[::m]
    nq = len(Q)
    aus = max(2, int(ausschluss_frac * nq))
    best, wo = np.inf, None
    for i in range(nq):
        d = np.linalg.norm(Q - Q[i], axis=1)
        for o in range(-aus, aus + 1):
            d[(i + o) % nq] = np.inf
        j = int(np.argmin(d))
        if d[j] < best:
            best, wo = float(d[j]), (i, j)
    return best, wo


def laenge(P):
    return float(np.linalg.norm(np.diff(np.vstack([P, P[:1]]), axis=0), axis=1).sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--skala", type=float, default=1.35)
    ap.add_argument("--spalt", type=float, default=0.10)
    ap.add_argument("--terme", type=int, default=1400)
    ap.add_argument("--name", default=None)
    ap.add_argument("--out", default=".")
    a = ap.parse_args()

    name = a.name or f"gitter_d3_k{a.k}_s{a.skala:g}".replace(".", "p")
    P, ecken = baue(a.k, a.skala, a.spalt)
    print(f"[{name}] build polyline: {len(P)} points, length {laenge(P):.3f}, "
          f"corners used {ecken}")


    ns = max(8192, 4 * a.terme)
    fa, fb = polyline_to_fourier(P, n_terms=a.terme, n_samples=ns)
    pfad = f"{a.out}/{name}.fseries"
    write_fseries(pfad, fa, fb, knot_name=name, source="build_lattice.py (Maschine 258)")

    off, ax, bx, ay, by, az, bz = parse_fseries(pfad)
    R = reconstruct_curve(off, ax, bx, ay, by, az, bz, 60000).astype(np.float64)
    Pu = resample_uniform(P, 60000); Pu -= Pu.mean(0)
    fehler = float(np.abs(np.linalg.norm(R, axis=1) - np.linalg.norm(Pu, axis=1)).max())
    from embedding import min_segmentabstand
    d_rek, wo_r, _h = min_segmentabstand(R)
    print(f"   .fseries: {a.terme} terms -> reconstructed 60k points")
    print(f"   smallest distance RECONSTRUCTED = {d_rek:.5f}  at {wo_r}")
    print(f"   ratio distance/length           = {d_rek/laenge(R):.3e}")
    urteil = "EINGEBETTET" if d_rek > 0 else "DURCHDRUNGEN"
    print(f"   ⟹ {urteil}")

    json.dump(dict(name=name, k=a.k, skala=a.skala, spalt=a.spalt, terme=a.terme,
                   punkte_bau=int(len(P)), ecken=[int(e) for e in ecken],
                    min_abstand_rekonstruiert=d_rek,
                   laenge=laenge(R), urteil=urteil,
                   quelle=FS, doktrin="Record-Objekt = servierte .fseries; keine Zahl aus der Konstruktion"),
              open(f"{a.out}/{name}.bau.json", "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print(f"   → {pfad}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
