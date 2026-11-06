# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 81 — clasping_the_double (Whitehead doubles).

The thread runs THERE (+d) and BACK in the Bishop tube of a companion
(−d); one end closes plainly, the other as CLASP (cap winds
once AROUND the double band — the hook). ALEXANDER-BLIND BY CONSTRUCTION:
Δ(W(K)) = 1, det = 1 for EVERY companion — the series factory for
θ rail (does θ separate the doubles of different companions? Dror wire).
BETS EXACT today: Bishop frame is TWIST-FREE (∮⟨e₁',e₂⟩ds ≈ 0,
measured — the machine-44 holonomy lesson as a number); closure exact;
clasp region survives the embedding gate. det=1 = ignition bet.
REVERSAL CLAUSE: clasp = densest point → own radius budget.
Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "10_climbing_the_torus_ladder"))
sys.path.insert(0, str(HERE.parent))
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402
from the_knotfamily_machine.fseries_io import reconstruct_from_fourier  # noqa: E402


def bishop_frame(X):
    """Parallel-transport frame (twist-free) along a closed curve."""
    n = len(X)
    T = np.roll(X, -1, axis=0) - np.roll(X, 1, axis=0)
    T /= np.linalg.norm(T, axis=1)[:, None]
    e1 = np.zeros_like(X)
    v = np.array([0.0, 0.0, 1.0])
    e1[0] = v - (v @ T[0]) * T[0]
    e1[0] /= np.linalg.norm(e1[0])
    for i in range(1, n):
        # apply minimal rotation T[i-1]→T[i] to e1
        a, b = T[i - 1], T[i]
        c = np.cross(a, b)
        s = np.linalg.norm(c)
        if s < 1e-14:
            e1[i] = e1[i - 1]
        else:
            c = c / s
            cosang = np.clip(a @ b, -1, 1)
            ang = np.arctan2(s, cosang)
            K = np.array([[0, -c[2], c[1]], [c[2], 0, -c[0]],
                          [-c[1], c[0], 0]])
            R = (np.eye(3) + np.sin(ang) * K
                 + (1 - np.cos(ang)) * (K @ K))
            e1[i] = R @ e1[i - 1]
        e1[i] -= (e1[i] @ T[i]) * T[i]
        e1[i] /= np.linalg.norm(e1[i])
    e2 = np.cross(T, e1)
    return T, e1, e2


def frame_twist(T, e1, e2):
    """∮⟨de₁, e₂⟩ — Bishop must give ≈0 (twist-free)."""
    de1 = np.roll(e1, -1, axis=0) - e1
    return float(np.abs((de1 * e2).sum(axis=1)[:-1].sum()))


def companion(kind="trefoil", n=2048):
    if kind == "trefoil":
        from generate_torus_ladder import torus_coeffs
        a, b = torus_coeffs(2, 3)
        return reconstruct_from_fourier(a, b, n_samples=n)
    if kind == "t25":
        from generate_torus_ladder import torus_coeffs
        a, b = torus_coeffs(2, 5)
        return reconstruct_from_fourier(a, b, n_samples=n)
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.column_stack([3 * np.cos(t), 3 * np.sin(t), 0 * t])


def whitehead_double(X, rho=0.05, n_sigma=4096, wdir=+1):
    """THIRD ARCHITECTURE (lessons 1+2 in the footer): continuous
    σ parametrisation σ∈[0,2): s(σ) = L·(1−cos πσ)/2 (smooth U-turn,
    no hairpin kink), cross-section angle ψ = πσ (the punctures
    ±ρe₁ are DISTINCT — lesson 2: fixed ±ρe₁ loops end twice at the
    same point). CLASP: in the turn window around σ=1 the thread winds in
    the cross-section disc once around the opposite puncture (+ρe₁) —
    the hook. Winding number 0 + clasp = Whitehead."""
    T, e1, e2 = bishop_frame(X)
    n = len(X)
    sig = np.linspace(0, 2, n_sigma, endpoint=False)
    s_idx = (1 - np.cos(np.pi * sig)) / 2 * n          # s(σ) as index
    i0 = np.floor(s_idx).astype(int) % n
    i1 = (i0 + 1) % n
    fr = (s_idx - np.floor(s_idx))[:, None]
    Xs = X[i0] * (1 - fr) + X[i1] * fr
    E1s = e1[i0] * (1 - fr) + e1[i1] * fr
    E2s = e2[i0] * (1 - fr) + e2[i1] * fr
    psi = np.pi * sig
    base = Xs + rho * (np.cos(psi)[:, None] * E1s
                       + np.sin(psi)[:, None] * E2s)
    # CLASP via the OUTER CORRIDOR (lesson 4: the crowd at the
    # turn is not entered at all): from the turn point −ρe₁ radially
    # out to 2.6ρ (nobody parks there), outside over the +e₂ side to the
    # +e₁ outer point, wind ONCE around the closing strand pair
    # (radius 0.9ρ around +ρe₁), outside back over the −e₂ side.
    # Outgoing and return arcs separated in T by ±0.5ρ.
    T0, E1, E2 = T[0], e1[0], e2[0]
    X0 = X[0]
    # Lesson 6 (root of the 0.05ρ collision): at the turn ds/dσ = 0 —
    # insertion point and follow-up point lie 0.0015ρ apart, the
    # feet of outgoing and return corridor are forced onto each other.
    # Hence cut out a SMALL window w'=0.08: the docking points
    # then lie at ψ = π(1∓w′) → ±0.26ρ separated in e₂.
    wprime = 0.08
    i1 = int((1 - wprime) / 2 * n_sigma)
    i2 = int((1 + wprime) / 2 * n_sigma)
    A_pt, B_pt = base[i1], base[i2]
    R_out = 2.6 * rho
    rc = 0.9 * rho
    C_wrap = X0 + rho * E1

    def disk_pt(radius, ang, t_off):
        return (X0 + radius * (np.cos(ang) * E1 + np.sin(ang) * E2)
                + t_off * T0)

    def arc(radius, a0, a1, t_off, m=60):
        aa = np.linspace(a0, a1, m)
        return np.array([disk_pt(radius, a, t_off) for a in aa])

    def chord(p, q, m=16):
        tt = np.linspace(0, 1, m + 2)[1:-1]
        return np.outer(1 - tt, p) + np.outer(tt, q)

    # Lesson 5: BOTH σ branches approach the turn from the −T side (s(σ) is
    # maximal at σ=1) — the free lane lies on +T. First along +T
    # out of the crowd, then radially; outgoing and return path separated in T (1.1ρ/1.8ρ)
    # AND in angle (π∓0.25), otherwise the chords cross.
    a_hin, a_rueck = np.pi - 0.25, np.pi + 0.25
    t_hin, t_rueck = 1.1 * rho, 1.8 * rho
    lift = chord(A_pt, A_pt + (t_hin - (A_pt - X0) @ T0) * T0, m=12)
    out1 = chord(lift[-1], disk_pt(R_out, a_hin, t_hin))
    swing1 = arc(R_out, a_hin, 0.0, t_hin)
    # Winding around the closing strand pair at +ρe₁ (T centered around t_hin)
    # wdir=−1 reverses the winding direction (clasp sign flip;
    # docking points φ=0/±2π identical — machine 113 uses this)
    C_w = C_wrap + t_hin * T0
    phi = wdir * np.linspace(0, 2 * np.pi, 121)
    wrap = (C_w[None, :]
            + rc * np.outer(np.cos(phi), E1)
            + rc * np.outer(np.sin(phi), E2)
            + np.outer(np.linspace(-0.15, 0.15, len(phi)) * rho, T0))
    in1 = chord(disk_pt(R_out, 0.0, t_hin), wrap[0])
    out2 = chord(wrap[-1], disk_pt(R_out, 0.0, t_rueck))
    swing2 = arc(R_out, 2 * np.pi, a_rueck, t_rueck)
    in2 = chord(disk_pt(R_out, a_rueck, t_rueck),
                B_pt + (t_rueck * 0.35) * T0, m=20)
    drop = chord(B_pt + (t_rueck * 0.35) * T0, B_pt, m=10)
    detour = np.vstack([lift, out1, swing1, in1, wrap, out2, swing2,
                        in2, drop])
    D = np.vstack([base[:i1 + 1], detour, base[i2:]])
    # Smooth transitions
    for _ in range(3):
        D = 0.5 * D + 0.25 * (np.roll(D, 1, axis=0)
                              + np.roll(D, -1, axis=0))
    return D


def roll_to_max_clearance(X, n_probe=1024):
    """Rotate the companion so that X[0] lies at the point of MAXIMAL self-
    distance (fix building block 2: the turn gets the freest spot)."""
    from the_knotfamily_machine.fseries_io import resample_uniform
    Q = resample_uniform(X, n_probe)
    n = len(Q)
    idx = np.arange(n)
    best_i, best_d = 0, -1.0
    for i in range(0, n, 4):
        d = np.linalg.norm(Q - Q[i], axis=1)
        ring = np.abs(idx - i)
        ring = np.minimum(ring, n - ring)
        d[ring <= 16] = np.inf
        dm = float(d.min())
        if dm > best_d:
            best_d, best_i = dm, i
    # transfer to X index
    j = int(round(best_i / n * len(X))) % len(X)
    return np.roll(X, -j, axis=0), best_d


def make_knot(kind, out_root, rho=None):
    X = companion(kind)
    X, local_clear = roll_to_max_clearance(X)
    name = f"wh_{kind}"
    entry = {"name": name, "machine": "81_clasping_the_double",
             "companion": kind,
             "predictions": {
                 "det_wette": "det(W(K)) = 1, Δ = 1 — für JEDEN "
                              "Begleiter (Zündung)",
                 "theta_frage": "trennt θ Doubles verschiedener "
                                "Begleiter? (Dror-Schiene)"}}
    # FIX 2026-06-08 (Disc #264 + addendum): TWO leaks at once —
    # (1) rho=0.3*local_clear was TOO THIN, the clasp does not hook -> raw curve
    #     is the UNKNOT (embedding gate let it pass: embedded != knotted).
    # (2) Even with a knotting raw curve the fseries EXPORT KILLS the knot:
    #     the clasp is a sub-resolution feature -> Nyquist-faithful (0.3%) is
    #     NOT topology-faithful (measured: clasp needs ~H=768; H response non-
    #     monotonic, 1024 wrong again). Hence search over (rho, n_terms) where the
    #     EXPORTED fseries really knots -> topo gate on the FSERIES (not only
    #     the raw curve). Fatter rho first (better representable). None ->
    #     reject honestly (clasp redesign needed).
    from machine_topo_gate import reduced_crossings, fseries_curve
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    fpath = folder / f"{name}.fseries"
    rhos = ([rho] if rho is not None else
            [round(f * local_clear, 4) for f in
             (0.90, 0.80, 0.70, 0.65, 0.60, 0.55, 0.52, 0.50, 0.45)])  # fat->thin
    chosen = None
    for rc_rho in rhos:
        Dc = whitehead_double(X, rc_rho)
        if not embedding_report(Dc)["embedded"]:
            continue
        if reduced_crossings(Dc) == 0:           # raw curve must knot first
            continue
        for nt in (768, 640, 896, 512, 1024, 384):   # H search: clasp ~768
            try:
                rep = polyline_to_fseries(Dc, fpath, knot_name=name,
                                          source="maschinen_ideen/81",
                                          n_terms=nt, nyquist_max_rel=0.05)
            except RuntimeError:
                continue
            if reduced_crossings(fseries_curve(str(fpath))) > 0:  # GATE: does the fseries knot?
                chosen = (rc_rho, Dc, nt, rep)
                break
        if chosen:
            break
    if chosen is None:
        entry["rejected"] = "topo gate (keine (rho,H)-Kombi liefert knotende fseries)"
        return entry, None
    rho, D, nt, rep = chosen
    entry["rho"] = rho
    entry["n_terms"] = nt
    entry["embedding"] = embedding_report(D)
    entry["reduced_crossings_fseries"] = reduced_crossings(fseries_curve(str(fpath)))
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 81 (Whitehead doubles)")
    X = companion("trefoil")
    T, e1, e2 = bishop_frame(X)
    # (1) Frame quality: orthonormal + twist CONVERGENCE (the residual
    #     twist is discretization error and must fall with n —
    #     law instead of wishful number, machine-72 lesson reused)
    ortho = max(float(np.abs((T * e1).sum(axis=1)).max()),
                float(np.abs((e1 * e2).sum(axis=1)).max()))
    tws = []
    for nn in (2048, 8192):
        Xn = companion("trefoil", n=nn)
        Tn, e1n, e2n = bishop_frame(Xn)
        tws.append(frame_twist(Tn, e1n, e2n))
    ok1 = ortho < 1e-12 and tws[1] < tws[0] / 3 and tws[1] < 1e-3
    print(f"  (1) Bishop frame: ortho {ortho:.1e}, twist "
          f"{tws[0]:.2e}→{tws[1]:.2e} (n×4 ⇒ falls — twist-free in "
          f"the limit)  {'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) Tube fidelity: base part (far from the hook) keeps distance ρ
    #     to the companion curve (σ architecture: position rotates in the tube)
    X2, _ = roll_to_max_clearance(X)
    rho = 0.16
    D = whitehead_double(X2, rho)
    probe = D[200:1500:50]                   # σ ≈ 0.1..0.73, without hook
    dmin = np.array([np.linalg.norm(X2 - q, axis=1).min() for q in probe])
    dev = float(np.abs(dmin - rho).max() / rho)
    ok2 = dev < 0.05
    print(f"  (2) tube fidelity (base): |distance−ρ|/ρ ≤ {dev:.4f}  "
          f"{'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) W(trefoil) built + gates (clasp survives)
    e3, f3 = make_knot("trefoil", HERE / "demo_smoke")
    ok3 = f3 is not None
    print(f"  (3) {e3['name']}: "
          f"nyq={e3.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓ clasp embedded' if ok3 else '✗ FAIL: ' + str(e3.get('rejected'))}")
    assert ok3
    # (4) THE SERIES: second companion (T(2,5)) — the blind family grows
    e4, f4 = make_knot("t25", HERE / "demo_smoke")
    ok4 = f4 is not None
    print(f"  (4) {e4['name']}: "
          f"nyq={e4.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓ series factory runs' if ok4 else '✗ FAIL: ' + str(e4.get('rejected'))}")
    assert ok4
    print("[DEMO] 4/4 green — Whitehead machine ready to fire "
          "(Alexander-blind series; det=1 bet goes to ignition).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
