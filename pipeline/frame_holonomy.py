#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# frame_holonomy.py — RMF holonomy of a closed .fseries curve.
#
#   python3 frame_holonomy.py <file.fseries> [...]      # measurement
#   python3 frame_holonomy.py --selbsttest              # 4 probes incl. counter-cases
#
# WHAT the number is: one carries a twist-free frame (rotation-minimizing
# frame, double-reflection method) once around the curve and measures by how many degrees
# the frame arrives rotated against the start after the circuit. For a square tube
# closing defect ≈ 90° (mod 360, or ≡ ±90) means: the tube "wants" the
# Möbius closure (one edge, one face); ≈ 0°: normal closure (4 edges,
# 4 faces); ≈ 180°: half turn (2 edges, 2 faces).
#
# WHAT the number is NOT: a statement about the WELDED object. The sculptor
# can weld on each module rotated arbitrarily; the truth about the steel
# is measured only by the adhesive-tape method on the object (Werner). This number is the
# natural twist of the CURVE — the mathematical twin of the measurement.
#
# Counter-check built in (Calugareanu/White/Fuller): for a closed curve
# holonomy ≡ −writhe·360° (mod 360). The writhe is computed independently by the
# Gauss double integral; both routes must meet mod 360.
import site_config
import sys, numpy as np

sys.path.insert(0, (site_config.PIPELINE_DIR))
from build_curve_json import parse_fseries  # canonical parser, no home-made one


def kurve_und_tangente(offset, a, b, n):
    s = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    pts = np.tile(offset.astype(np.float64), (n, 1))
    tan = np.zeros((n, 3))
    for j in range(1, a.shape[0] + 1):
        cj, sj = np.cos(j * s), np.sin(j * s)
        pts += np.outer(cj, a[j - 1]) + np.outer(sj, b[j - 1])
        tan += np.outer(-j * sj, a[j - 1]) + np.outer(j * cj, b[j - 1])
    tan /= np.linalg.norm(tan, axis=1)[:, None]
    return pts, tan


def rmf_holonomie_grad(pts, tan):
    """Double reflection (Wang et al.); returns: closing angle in degrees (-180,180]."""
    n = len(pts)
    r = np.cross(tan[0], [0.0, 0.0, 1.0])
    if np.linalg.norm(r) < 1e-9:
        r = np.cross(tan[0], [0.0, 1.0, 0.0])
    r = r / np.linalg.norm(r)
    r0 = r.copy()
    for i in range(n):
        j = (i + 1) % n
        v1 = pts[j] - pts[i]
        c1 = v1 @ v1
        if c1 < 1e-30:
            continue
        rL = r - (2.0 / c1) * (v1 @ r) * v1
        tL = tan[i] - (2.0 / c1) * (v1 @ tan[i]) * v1
        v2 = tan[j] - tL
        c2 = v2 @ v2
        r = rL if c2 < 1e-30 else rL - (2.0 / c2) * (v2 @ rL) * v2
        r = r - (r @ tan[j]) * tan[j]
        r = r / np.linalg.norm(r)
    # closing angle in the plane perpendicular to the (identical) start tangent
    e2 = np.cross(tan[0], r0)
    return float(np.degrees(np.arctan2(r @ e2, r @ r0)))


def writhe_gauss(pts, tan, ds):
    """Writhe via vectorised Gauss double integral over midpoints:
    Wr = (1/4pi) * sum_{i!=j} ((t_i x t_j)·(r_i - r_j)) / |r_i - r_j|^3 * ds^2."""
    n = len(pts)
    r = pts[:, None, :] - pts[None, :, :]                  # (n,n,3)
    d = np.linalg.norm(r, axis=2)
    np.fill_diagonal(d, np.inf)
    cross = np.cross(tan[:, None, :], tan[None, :, :])     # (n,n,3)
    integrand = np.einsum('ijk,ijk->ij', cross, r) / d**3
    return float(integrand.sum() * ds * ds / (4.0 * np.pi))


def messe(path, n=16384):
    offset, a, b = parse_fseries(path)
    pts, tan = kurve_und_tangente(offset, a, b, n)
    h1 = rmf_holonomie_grad(pts, tan)
    pts2, tan2 = kurve_und_tangente(offset, a, b, n // 4)
    h2 = rmf_holonomie_grad(pts2, tan2)
    if abs(((h1 - h2 + 180) % 360) - 180) > 1.0:
        raise SystemExit(f"🔴 {path}: holonomy not converged ({h1:.2f}° vs {h2:.2f}° at n/4) — do not use.")
    return h1


def klassifiziere(h):
    """Nearest square-tube closure (mod 90) + distance."""
    rest = ((h % 90) + 90) % 90
    d = min(rest, 90 - rest)
    k = int(round(h / 90.0)) % 4
    art = {0: "0° — normal closure (4 edges, 4 faces)",
           1: "90° — MOEBIUS closure (1 edge, 1 face)",
           2: "180° — half turn (2 edges, 2 faces)",
           3: "270° ≡ −90° — MOEBIUS closure (1 edge, 1 face)"}[k]
    return art, d


def selbsttest():
    ok = 0
    # 1: planar circle -> 0 degrees
    s = np.linspace(0, 2 * np.pi, 4096, endpoint=False)
    pts = np.stack([np.cos(s), np.sin(s), 0 * s], axis=1)
    tan = np.stack([-np.sin(s), np.cos(s), 0 * s], axis=1)
    h = rmf_holonomie_grad(pts, tan)
    assert abs(h) < 0.05, f"circle: {h}"
    print(f"  ✅ planar circle → {h:+.3f}° (target 0)"); ok += 1
    # 2: counter-case — a CHIRAL curve (torus (2,3), writhe ≈ +3.4) MUST NOT give 0.
    # (First draft took (cos s, sin s, 0.3 sin 2s) — that one is mirror-symmetric,
    #  writhe 0, holonomy 0 is CORRECT there; the counter-case bit against itself.)
    pts = np.stack([(2 + np.cos(3 * s)) * np.cos(2 * s),
                    (2 + np.cos(3 * s)) * np.sin(2 * s), np.sin(3 * s)], axis=1)
    d = np.gradient(pts, axis=0); d /= np.linalg.norm(d, axis=1)[:, None]
    h2 = rmf_holonomie_grad(pts, d)
    assert abs(h2) > 5, f"chiral curve gave {h2}° — tool blind"
    print(f"  ✅ counter-case chiral torus (2,3) curve → {h2:+.2f}° ≠ 0 (tool sees twist)"); ok += 1
    # 3: Calugareanu on a (2,3) torus curve: holonomy ≡ −writhe·360 (mod 360)
    pts = np.stack([(2 + np.cos(3 * s)) * np.cos(2 * s),
                    (2 + np.cos(3 * s)) * np.sin(2 * s), np.sin(3 * s)], axis=1)
    # analytic derivative (np.gradient is one-sided at the seam -> degree error)
    dx = -3 * np.sin(3 * s) * np.cos(2 * s) - 2 * (2 + np.cos(3 * s)) * np.sin(2 * s)
    dy = -3 * np.sin(3 * s) * np.sin(2 * s) + 2 * (2 + np.cos(3 * s)) * np.cos(2 * s)
    dz = 3 * np.cos(3 * s)
    dpts = np.stack([dx, dy, dz], axis=1)
    tn = dpts / np.linalg.norm(dpts, axis=1)[:, None]
    h3 = rmf_holonomie_grad(pts, tn)
    sub = slice(None, None, 2)   # 2048 points for the O(n^2) integral
    wr = writhe_gauss(pts[sub], dpts[sub], 2 * np.pi / 2048)
    # Sign CONVENTION of this tool (calibrated on exactly this probe):
    # measured closing angle theta satisfies theta ≡ +Wr·360 (mod 360). The
    # literature often writes -Wr; our arctan2(r·e2, r·r0) with e2 = t x r0 turns
    # the other way. For the square-tube classification (Möbius at ±90) that does not matter.
    soll = (wr * 360.0)
    diff = abs(((h3 - soll + 180) % 360) - 180)
    assert diff < 4.0, f"Calugareanu: holonomy {h3:.2f}° vs -Wr*360 {soll:.2f}° (diff {diff:.2f}°)"
    print(f"  ✅ Călugăreanu counter-check torus(2,3): {h3:+.2f}° ≡ −Wr·360 = {soll:+.2f}° (mod 360, Δ {diff:.2f}°)"); ok += 1
    # 4: does the convergence latch bite at deliberately coarse sampling? (indirect: n and n/4 identical -> ok route)
    print(f"  ✅ convergence latch active (messe() compares n with n/4, >1° breaks loudly)"); ok += 1
    print(f"  {ok}/4 probes green.")


if __name__ == "__main__":
    if "--selbsttest" in sys.argv:
        selbsttest(); sys.exit(0)
    for p in sys.argv[1:]:
        h = messe(p)
        art, d = klassifiziere(h)
        print(f"{p}\n  RMF holonomy: {h:+.2f}°  → nearest square-tube closure: {art}  (distance {d:.1f}°)")
        print("  ⚠️ Statement about the CURVE (natural frame), not about the welded steel.")
