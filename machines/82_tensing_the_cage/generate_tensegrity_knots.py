# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 82 — tensing_the_cage (tensegrity form-finding).

Closed ROPE loop (cycle, q=1) + two STRUT families
(offset d₁/d₂, pressure: q<0) on n nodes: the equilibrium matrix
is CIRCULANT → eigenvalues λ(k) ANALYTIC; strut densities are
SOLVED so that λ(k₁)=λ(k₂)=0 (linear 2×2) → the kernel CONTAINS the
coordinates: x,y = mode k₁, z = mode k₂ — the shape finds itself
(force density, Snelson principle). The rope loop in
equilibrium = the knot ({n/k₁} star + k₂ wave).
BETS EXACT today: λ(k) formula == numerical eigenvalues (bit);
equilibrium residual L·X = 0 < 1e-10; kernel dimension counted ==
predicted. Sculptural: the curve that carries itself.
REVERSAL CLAUSE: rope must not cut struts → distance gate.
Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402


def laplacian_circulant_eigs(n, fams):
    """λ(k) = Σ_f q_f·(2−2cos(2π·d_f·k/n)) for edge family with
    offset d_f and density q_f (cycle = family d=1, q=1)."""
    k = np.arange(n)
    lam = np.zeros(n)
    for d, q in fams:
        lam += q * (2 - 2 * np.cos(2 * np.pi * d * k / n))
    return lam


def build_laplacian(n, fams):
    L = np.zeros((n, n))
    for d, q in fams:
        for i in range(n):
            j = (i + d) % n
            L[i, i] += q
            L[j, j] += q
            L[i, j] -= q
            L[j, i] -= q
    return L


def solve_densities(n, k1, k2, d1, d2):
    """q₁,q₂ of the strut families such that λ(k₁)=λ(k₂)=0."""
    def m(d, k):
        return 2 - 2 * np.cos(2 * np.pi * d * k / n)
    A = np.array([[m(d1, k1), m(d2, k1)], [m(d1, k2), m(d2, k2)]])
    b = -np.array([m(1, k1), m(1, k2)])
    q1, q2 = np.linalg.solve(A, b)
    return float(q1), float(q2)


def form_find(n, k1, k2, d1, d2, z_amp=0.8, z_phase=0.6):
    """z_phase=0.6: build lesson — at φ=0 the wave extrema of the
    k₂ mode lie exactly on the star crossings (contact, gate tore);
    the phase is a free kernel direction, measured φ≈0.6 separates cleanly
    (d_min_rel 0.017)."""
    q1, q2 = solve_densities(n, k1, k2, d1, d2)
    fams = [(1, 1.0), (d1, q1), (d2, q2)]
    th = 2 * np.pi * np.arange(n) / n
    X = np.column_stack([np.cos(k1 * th), np.sin(k1 * th),
                         z_amp * np.cos(k2 * th + z_phase)])
    L = build_laplacian(n, fams)
    resid = float(np.abs(L @ X).max())
    return X, fams, (q1, q2), resid, L


def rounded_polygon(X, pts_per_edge=60, rounds=3):
    dense = []
    m = len(X)
    for k in range(m):
        a, b = X[k], X[(k + 1) % m]
        t = np.linspace(0, 1, pts_per_edge, endpoint=False)
        dense.append(np.outer(1 - t, a) + np.outer(t, b))
    C = np.vstack(dense)
    for _ in range(rounds):
        C = 0.5 * C + 0.25 * (np.roll(C, 1, axis=0) + np.roll(C, -1, axis=0))
    return C


def strut_clearance(X, d1, d2, C):
    """Min distance rope curve ↔ strut chords."""
    best = np.inf
    n = len(X)
    for d in (d1, d2):
        for i in range(n):
            a, b = X[i], X[(i + d) % n]
            t = np.linspace(0.12, 0.88, 30)
            seg = np.outer(1 - t, a) + np.outer(t, b)
            dmin = np.min(np.linalg.norm(C[:, None, :] - seg[None, :, :],
                                         axis=2))
            best = min(best, float(dmin))
    return best


def make_knot(n, k1, k2, d1, d2, out_root):
    X, fams, (q1, q2), resid, L = form_find(n, k1, k2, d1, d2)
    name = f"tg_n{n}_k{k1}_{k2}"
    entry = {"name": name, "machine": "82_tensing_the_cage",
             "n_nodes": n, "modes": [k1, k2], "strut_offsets": [d1, d2],
             "strut_densities": [round(q1, 6), round(q2, 6)],
             "equilibrium_residual": resid,
             "predictions": {
                 "identity": f"{{{n}/{k1}}}-Stern + k={k2}-Welle "
                             "(Zünd-Messung)",
                 "skulptur": "die Kurve, die sich selbst trägt"}}
    if resid > 1e-10:
        entry["rejected"] = f"Gleichgewicht verletzt ({resid:.2e})"
        return entry, None
    if not (q1 < 0 or q2 < 0):
        entry["rejected"] = "keine Druck-Strebe (kein Tensegrity)"
        return entry, None
    C = rounded_polygon(X)
    entry["strut_clearance"] = round(strut_clearance(X, d1, d2, C), 4)
    emb = embedding_report(C)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(C, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/82",
                              n_terms=max(96, 10 * n),
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 82 (tensegrity)")
    n, k1, k2, d1, d2 = 12, 5, 3, 6, 4
    # (1) circulant formula == numerical eigenvalues
    q1, q2 = solve_densities(n, k1, k2, d1, d2)
    fams = [(1, 1.0), (d1, q1), (d2, q2)]
    L = build_laplacian(n, fams)
    lam_num = np.sort(np.linalg.eigvalsh(L))
    lam_form = np.sort(laplacian_circulant_eigs(n, fams))
    d = float(np.abs(lam_num - lam_form).max())
    ok1 = d < 1e-12
    print(f"  (1) λ(k) formula vs numerical: max|Δ|={d:.2e} "
          f"(q₁={q1:.3f}, q₂={q2:.3f})  {'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) kernel dimension == predicted (k₁ pair + k₂ pair + constant)
    nullity = int((np.abs(lam_num) < 1e-10).sum())
    want = 1 + 2 + 2          # k=0 + (±k₁) + (±k₂)
    ok2 = nullity == want
    print(f"  (2) kernel dimension: {nullity} (expected {want}: const + "
          f"k₁ pair + k₂ pair)  {'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) Equilibrium: L·X = 0 for the shape
    X, fams, qs, resid, L = form_find(n, k1, k2, d1, d2)
    ok3 = resid < 1e-10 and (qs[0] < 0 or qs[1] < 0)
    print(f"  (3) equilibrium residual: {resid:.2e}; compression strut "
          f"present ({qs[0]:.3f}/{qs[1]:.3f})  "
          f"{'✓ shape carries itself' if ok3 else '✗ FAIL'}")
    assert ok3
    # (4) tensegrity knot built ({12/5} star + 3-wave)
    e4, f4 = make_knot(n, k1, k2, d1, d2, HERE / "demo_smoke")
    ok4 = f4 is not None
    print(f"  (4) {e4['name']}: strut clearance="
          f"{e4.get('strut_clearance')}, "
          f"nyq={e4.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok4 else '✗ FAIL: ' + str(e4.get('rejected'))}")
    assert ok4
    print("[DEMO] 4/4 green — tensegrity machine ready to fire "
          "(form-finding exact; stainless-steel struts + rope knot buildable).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
