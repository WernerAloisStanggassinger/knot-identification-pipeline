# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 60 — slicing_singularities (Milnor links from algebra).

Plane-curve singularity f(z₁,z₂)=0 near 0, curve on small S³,
stereographically to R³ (machine-56 machinery reused) — Milnor's
theorem delivers the LINK of the singularity. Torus case exact: f=z₁^p−z₂^q
→ T(p,q) (Brauner); Puiseux perturbation ((z₁²−z₂³)²−z₁⁵) → ITERATED
torus knot = cable tower FROM ALGEBRA instead of construction (roof over
10/44/56). Numerics: branch tracking of the z₂ roots along z₁=ρe^{iφ}
with |f| residual guard; radial S³ projection (isotopic for small ρ,
standard Milnor argument).
REVERSAL CLAUSE: root tracking must never switch the branch →
continuity step small + residual gate. Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from math import gcd
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "56_projecting_from_four_dimensions"))
sys.path.insert(0, str(HERE.parent))
from generate_stereographic_knots import (stereographic,  # noqa: E402
                                          round_torus_residual,
                                          winding_numbers, pole_gap)
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402


def torus_link_curve(p, q, n=8192):
    """Exact parametrisation of V(z₁^p−z₂^q)∩S³: z₁=r₁e^{iqs},
    z₂=r₂e^{ips} with r₁^p=r₂^q, r₁²+r₂²=1 (scalar via bisection)."""
    lo, hi = 1e-9, 1.0
    for _ in range(200):
        r1 = 0.5 * (lo + hi)
        r2 = r1 ** (p / q)
        if r1 ** 2 + r2 ** 2 > 1.0:
            hi = r1
        else:
            lo = r1
    r1 = 0.5 * (lo + hi)
    r2 = r1 ** (p / q)
    s = np.linspace(0, 2 * np.pi, n, endpoint=False)
    z1 = r1 * np.exp(1j * q * s)
    z2 = r2 * np.exp(1j * p * s)
    X = np.column_stack([z1.real, z1.imag, z2.real, z2.imag])
    return X, r1, r2


def f_puiseux(z1, z2):
    return (z1 ** 2 - z2 ** 3) ** 2 - z1 ** 5


def puiseux_curve(rho=0.1, n_per_turn=2048, max_turns=16, tol=1e-10):
    # rho=0.1: BUILD LESSON — the cable strands separate only ~ρ^(11/6)
    # (measured: ρ=0.05 → d_min_rel 0.0022 tears, ρ=0.1 → 0.0040 green);
    # ρ convergence (same type at ρ=0.1 and 0.2) = ignition check.
    """Branch tracking: z₁=ρe^{iφ}, z₂ follows the nearest root of
    f(z₁,·)=0 (numpy.roots in z₂, degree 6) until the loop closes."""
    phi_step = 2 * np.pi / n_per_turn
    z1 = rho
    # Start: roots of (z₁²−w)²−z₁⁵ in w=z₂³ → w=z₁²±z₁^{5/2}, then ³√
    w0 = z1 ** 2 + z1 ** 2.5
    z2 = w0 ** (1 / 3)
    pts = []
    phi = 0.0
    for step in range(max_turns * n_per_turn):
        z1 = rho * np.exp(1j * phi)
        # polynomial in z₂: (z₁²−z₂³)²−z₁⁵ = z₂⁶ −2z₁²z₂³ + z₁⁴−z₁⁵
        coeffs = [1, 0, 0, -2 * z1 ** 2, 0, 0, z1 ** 4 - z1 ** 5]
        roots = np.roots(coeffs)
        z2 = roots[np.argmin(np.abs(roots - z2))]      # continuity
        resid = abs(f_puiseux(z1, z2))
        if resid > tol:
            raise RuntimeError(f"residual guard: |f|={resid:.2e}")
        pts.append([z1.real, z1.imag, z2.real, z2.imag])
        phi += phi_step
        if step > n_per_turn // 2:
            d = np.linalg.norm(np.array(pts[-1]) - np.array(pts[0]))
            if d < np.linalg.norm(np.array(pts[1]) - np.array(pts[0])):
                break
    X = np.asarray(pts)
    turns = phi / (2 * np.pi)
    return X, turns


def to_r3(X):
    """Radial S³ projection + stereographic projection (machine-56 functions)."""
    S = X / np.linalg.norm(X, axis=1)[:, None]
    if pole_gap(S) < 0.02:
        S = S[:, [3, 1, 2, 0]] * np.array([1, 1, 1, -1.0])   # rotate the pole
    return stereographic(S)


def make_torus_knot(p, q, out_root):
    assert gcd(p, q) == 1
    X, r1, r2 = torus_link_curve(p, q)
    Y = to_r3(X)
    name = f"mi_T{p}_{q}"
    entry = {"name": name, "machine": "60_slicing_singularities",
             "f": f"z1^{p} - z2^{q}", "milnor_mu": (p - 1) * (q - 1),
             "predictions": {"satz": f"Brauner: Link = T({p},{q})",
                             "identity": f"T({p},{q}) — gegen Maschine 10 "
                                         "prüfbar (Zündung)"}}
    emb = embedding_report(Y)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(Y, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/60",
                              n_terms=max(96, 16 * (p + q)),
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 60 (cutting singularities)")
    # (1) Brauner anchor: V(z₁²−z₂³)∩S³ lies stereographically on an
    #     EXACT round torus with windings (3,±2)
    X, r1, r2 = torus_link_curve(2, 3)
    onsphere = float(np.abs((X ** 2).sum(axis=1) - 1).max())
    Y = stereographic(X)
    R, r, resid = round_torus_residual(Y)
    w1, w2 = winding_numbers(Y, R)
    ok1 = (onsphere < 1e-12 and resid < 1e-10
           and abs(w1 - 3) < 1e-6 and abs(abs(w2) - 2) < 1e-6)
    print(f"  (1) Brauner T(2,3): |X|²−1={onsphere:.1e}, torus resid="
          f"{resid:.1e}, windings=({w1:.6f},{w2:.6f})  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) Milnor number arithmetic μ=(p−1)(q−1)
    mus = [(2, 3, 2), (2, 5, 4), (3, 4, 6), (3, 7, 12)]
    ok2 = all((p - 1) * (q - 1) == mu for p, q, mu in mus)
    print(f"  (2) μ=(p−1)(q−1): {[(p,q,(p-1)*(q-1)) for p,q,_ in mus]}  "
          f"{'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) PUISEUX TOWER: branch tracking closes + residual guard held
    Xp, turns = puiseux_curve()
    onsphere3 = "radial projiziert"
    ok3 = 5.5 < turns < 12.5      # iterated knot needs multiple turns
    print(f"  (3) (z₁²−z₂³)²−z₁⁵: loop closes after {turns:.2f} "
          f"z₁ revolutions, |f| guard < 1e-10 throughout  "
          f"{'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    # (4) Cable tower built + gates green (first cable knot FROM algebra)
    Yp = to_r3(Xp)
    emb = embedding_report(Yp)
    folder = HERE / "demo_smoke" / "mi_puiseux_23_5"
    ok4 = emb["embedded"]
    if ok4:
        folder.mkdir(parents=True, exist_ok=True)
        rep = None
        for nt in (160, 256, 384):
            try:
                rep = polyline_to_fseries(Yp, folder / "mi_puiseux_23_5.fseries",
                                          knot_name="mi_puiseux_23_5",
                                          source="maschinen_ideen/60",
                                          n_terms=nt, nyquist_max_rel=0.02)
                break
            except RuntimeError:
                continue
        ok4 = rep is not None
    print(f"  (4) Puiseux cable tower: embedded={emb['embedded']}, "
          f"nyq={rep.get('nyquist_dev_rel') if ok4 else '—'}  "
          f"{'✓' if ok4 else '✗ FAIL'}")
    assert ok4
    # (5) torus ladder from algebra: T(3,5) built (μ=8)
    e5, f5 = make_torus_knot(3, 5, HERE / "demo_smoke")
    ok5 = f5 is not None
    print(f"  (5) mi_T3_5 (μ={e5['milnor_mu']}): "
          f"nyq={e5.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok5 else '✗ FAIL: ' + str(e5.get('rejected'))}")
    assert ok5
    print("[DEMO] 5/5 green — Milnor machine ready to fire "
          "(cable towers from algebra; roof over 10/44/56).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
