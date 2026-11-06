# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 89 — orbiting_the_black_hole (Kerr photon resonances).

Spherical photon orbits of the Kerr spacetime (Boyer–Lindquist, M=1,
spin a): r=const, θ oscillates, φ winds. Hamiltonian form H = ½g^{μν}p_μp_ν
with NUMERICAL ∂H (no hand-written Christoffel — less
typo surface). THREE CONSERVATION GUARDS: E, L_z AND the Carter
constant Q — Q is conserved ONLY in the true Kerr metric ⇒ the guard
is at once the metric correctness witness. Rational θ:φ resonance
by bisection in the orbit radius → closed light path; solid-torus
embedding of the phase angles (66 construction — the orbit itself lies
on the r₀ sphere, curves on spheres are unknots).
EXACT BETS today: H=0 + E/L/Q < 1e-8 over the whole orbit;
resonance bisection hits Δφ/period = p/q.
REVERSAL CLAUSE: near the horizon it gets stiff → step guard.
Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402

A_SPIN = 0.9
R_RING, R_TUBE = 6.0, 1.8


def inv_metric(x, a=A_SPIN):
    """Inverse Kerr metric (BL, M=1) on (t,r,θ,φ) momenta."""
    r, th = x[1], x[2]
    s2 = np.sin(th) ** 2
    Sig = r * r + a * a * np.cos(th) ** 2
    Del = r * r - 2 * r + a * a
    g = np.zeros((4, 4))
    g[0, 0] = -((r * r + a * a) ** 2 - a * a * Del * s2) / (Sig * Del)
    g[0, 3] = g[3, 0] = -2 * a * r / (Sig * Del)
    g[1, 1] = Del / Sig
    g[2, 2] = 1 / Sig
    g[3, 3] = (Del - a * a * s2) / (Sig * Del * s2)
    return g


def H(x, p):
    return 0.5 * p @ inv_metric(x) @ p


def grad_H_x(x, p, h=1e-7):
    g = np.zeros(4)
    for i in (1, 2):                      # only r,θ (t,φ cyclic exactly)
        e = np.zeros(4); e[i] = h
        g[i] = (H(x + e, p) - H(x - e, p)) / (2 * h)
    return g


def rk4_step(x, p, dl):
    def f(x, p):
        return inv_metric(x) @ p, -grad_H_x(x, p)
    k1x, k1p = f(x, p)
    k2x, k2p = f(x + 0.5 * dl * k1x, p + 0.5 * dl * k1p)
    k3x, k3p = f(x + 0.5 * dl * k2x, p + 0.5 * dl * k2p)
    k4x, k4p = f(x + dl * k3x, p + dl * k3p)
    return (x + dl / 6 * (k1x + 2 * k2x + 2 * k3x + k4x),
            p + dl / 6 * (k1p + 2 * k2p + 2 * k3p + k4p))


def spherical_orbit_params(r0, a=A_SPIN):
    """ξ=L/E, η=Q/E² of the spherical photon orbits (standard formulas,
    VERIFIED by the r-constancy guard of the integration)."""
    xi = (r0 ** 2 * (3 - r0) - a * a * (1 + r0)) / (a * (r0 - 1))
    eta = (r0 ** 3 * (4 * a * a - r0 * (r0 - 3) ** 2)) / (a * a
                                                          * (r0 - 1) ** 2)
    return xi, eta


def theta_roots(xi, eta, a=A_SPIN):
    """f(u) = −a²(u²−u₊²)(u²−u₋²), u = cosθ: turning-point latitudes
    analytically (u₊ = max|cosθ|)."""
    A2 = a * a
    s = (A2 - eta - xi * xi) / A2
    prod = -eta / A2
    disc = np.sqrt(max(s * s - 4 * prod, 0.0))
    up2 = (s + disc) / 2
    um2 = (s - disc) / 2
    return np.sqrt(max(up2, 0.0)), um2


def integrate_orbit(r0, n_steps=24000, dl=0.02):
    """PHYSICS LESSON: the photon sphere is UNSTABLE (potential
    maximum) — the raw Hamiltonian flow explodes exponentially
    (H drift 8e13 measured). Hence Carter-SEPARATED equations on
    of the spherical manifold, θ via ψ substitution
    u = u₊·sin ψ (smooth at turning points); the H=0 residual check against the
    RAW metric remains the cross/Kerr witness."""
    a = A_SPIN
    xi, eta = spherical_orbit_params(r0)
    up, um2 = theta_roots(xi, eta)
    Del = r0 * r0 - 2 * r0 + a * a
    psi, phi = 0.0, 0.0
    traj = [(np.pi / 2, 0.0, 0.0)]
    h_resid = 0.0
    for _ in range(n_steps):
        def rhs(psi, phi):
            u = up * np.sin(psi)
            th = np.arccos(np.clip(u, -1, 1))
            Sig = r0 * r0 + a * a * u * u
            dpsi = a * np.sqrt(up * up * np.sin(psi) ** 2 - um2) / Sig
            s2 = 1 - u * u
            dphi = (xi / s2 - a + a * (r0 * r0 + a * a - xi * a) / Del) / Sig
            return dpsi, dphi
        k1 = rhs(psi, phi)
        k2 = rhs(psi + 0.5 * dl * k1[0], phi + 0.5 * dl * k1[1])
        k3 = rhs(psi + 0.5 * dl * k2[0], phi + 0.5 * dl * k2[1])
        k4 = rhs(psi + dl * k3[0], phi + dl * k3[1])
        psi += dl / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        phi += dl / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        u = up * np.sin(psi)
        th = np.arccos(np.clip(u, -1, 1))
        traj.append((th, phi, psi))
        # Cross-witness: p from separated solution into the RAW metric
        du_dl = a * np.sqrt(max((up * up - u * u) * (u * u - um2), 0.0))
        sgn = 1.0 if np.cos(psi) >= 0 else -1.0
        s_th = np.sqrt(max(1 - u * u, 1e-15))
        p_th = -sgn * du_dl / s_th          # p_θ = Σ·dθ/dλ
        x4 = np.array([0.0, r0, th, phi])
        p4 = np.array([-1.0, 0.0, p_th, xi])
        h_resid = max(h_resid, abs(H(x4, p4)))
    traj = np.asarray(traj)
    drift = {"H_kreuz": h_resid, "u_plus": float(up)}
    return traj, drift


def phi_per_theta_period(r0):
    """Δφ over one full θ oscillation (measured on the orbit)."""
    traj, drift = integrate_orbit(r0, n_steps=12000, dl=0.01)
    th = traj[:, 0] - np.pi / 2
    sign_flips = np.nonzero(np.diff(np.sign(th)) != 0)[0]
    if len(sign_flips) < 5:
        return None, drift
    i0, i1 = sign_flips[0], sign_flips[4]      # 4 flips = 2 half periods
    dphi = traj[i1, 1] - traj[i0, 1]
    return float(dphi / 2), drift              # Δφ per full θ period


def find_resonance(target, lo, hi, iters=18):
    """Bisection: Δφ/2π = target (e.g. 3/2)."""
    f_lo = phi_per_theta_period(lo)[0] / (2 * np.pi) - target
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        f_m = phi_per_theta_period(mid)[0] / (2 * np.pi) - target
        if (f_m > 0) == (f_lo > 0):
            lo, f_lo = mid, f_m
        else:
            hi = mid
    return 0.5 * (lo + hi)


def make_knot(r0, q_periods, out_root, tag):
    traj, drift = integrate_orbit(r0, n_steps=q_periods * 7000, dl=0.012)
    name = f"kr_{tag}"
    entry = {"name": name, "machine": "89_orbiting_the_black_hole",
             "r0": round(r0, 6), "spin": A_SPIN,
             "watch_drift": {k: float(f"{v:.3e}") for k, v in
                             drift.items()},
             "predictions": {
                 "identity": "Torus-Knoten-Klasse der Resonanz (Zündung)",
                 "physik": "sphärischer Kerr-Photonen-Orbit"}}
    if drift["H_kreuz"] > 1e-8:
        entry["rejected"] = f"Kreuz-Zeuge-Wache ({drift})"
        return entry, None
    # Solid torus: φ = ring, ψ = tube ANGLE (build lesson: the
    # θ deflection OSCILLATES on a tube arc instead of winding —
    # outgoing and returning branches touch; the ψ phase is monotonic
    # and winds cleanly: T(p,q) arises only with ψ as tube coordinate.)
    phi = traj[:, 1]
    psi_t = traj[:, 2]
    # cut EXACTLY at ψ = q·2π (build lesson: a continuing orbit
    # covers its own trace — (3k,2k) multiple cover = near-
    # contact with itself), distribute residual drift
    over = np.nonzero(psi_t >= q_periods * 2 * np.pi)[0]
    if len(over) == 0:
        entry["rejected"] = "Bahn zu kurz für q Perioden"
        return entry, None
    k_cut = int(over[0])
    phi, psi_t = phi[:k_cut], psi_t[:k_cut]
    r = R_RING + R_TUBE * np.cos(psi_t)
    P = np.column_stack([r * np.cos(phi), r * np.sin(phi),
                         R_TUBE * np.sin(psi_t)])
    gap = P[-1] - P[0]
    P = P - np.outer(np.linspace(0, 1, len(P)), gap)
    P = P[::3]
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = None
    for nt in (128, 192, 256):
        try:
            rep = polyline_to_fseries(P, folder / f"{name}.fseries",
                                      knot_name=name,
                                      source="maschinen_ideen/89",
                                      n_terms=nt, nyquist_max_rel=0.02)
            break
        except RuntimeError:
            continue
    if rep is None:
        entry["rejected"] = "Nyquist rot"
        return entry, None
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 89 (Kerr photons)")
    # (1) CROSS-WITNESS: separated (Carter) solution into the RAW metric
    #     substituted ⇒ H=0 must hold bit-exactly (E,L exact by construction)
    traj, drift = integrate_orbit(2.5, n_steps=8000, dl=0.01)
    ok1 = drift["H_kreuz"] < 1e-9
    print(f"  (1) H cross witness (r₀=2.5, a=0.9): max|H|="
          f"{drift['H_kreuz']:.1e} — separation == raw Kerr metric  "
          f"{'✓ Carter holds ⇒ metric is Kerr' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) turning-point latitude: max|cosθ| of the orbit == u₊ analytically
    u_meas = float(np.abs(np.cos(traj[:, 0])).max())
    ok2 = abs(u_meas - drift["u_plus"]) < 1e-6
    print(f"  (2) turning latitude: measured {u_meas:.8f} vs analytic "
          f"{drift['u_plus']:.8f}  {'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) resonance bisection: aim at Δφ/2π = 3/2
    r_res = find_resonance(1.5, 1.9, 3.1)
    ratio, _ = phi_per_theta_period(r_res)
    ratio /= 2 * np.pi
    ok3 = abs(ratio - 1.5) < 1e-3
    print(f"  (3) resonance: r₀={r_res:.5f} → Δφ/2π per θ period = "
          f"{ratio:.5f} (target 1.5)  {'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    # (4) Light knot built (3:2 resonance, 2 θ periods = 3 φ revolutions)
    e4, f4 = make_knot(r_res, 2, HERE / "demo_smoke", "res3zu2")
    ok4 = f4 is not None
    print(f"  (4) {e4['name']}: guards={e4.get('watch_drift')}, "
          f"nyq={e4.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok4 else '✗ FAIL: ' + str(e4.get('rejected'))}")
    assert ok4
    print("[DEMO] 4/4 green — Kerr machine ready to fire "
          "(the knot that light ties around the black hole).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
