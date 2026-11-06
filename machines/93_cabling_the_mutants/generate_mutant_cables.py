# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 93 — cabling_the_mutants (satellite separator).

Classic: cables/satellites SEPARATE mutants where the polynomials of the
pair itself collide (2-cable HOMFLY). Here: (2,·) Möbius band
boundary cable (σ∈[0,2), s=σ·L linear, ψ=πσ — ONE strand, Bishop tube
from 81) onto BOTH curves of a 90er mutant pair.
BETS EXACT today: tube fidelity (distance ρ); closure (ψ=2π meets
ψ=0); both cables embedded as a pair. IGNITION: do the cable
invariants separate the pair? (If yes: mutants proven distinct + θ-on-
satellites ammunition for Dror.) REVERSAL CLAUSE: ρ ≤ 0.25·clearance
(81er lesson), mutant curves are tight. Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "81_clasping_the_double"))
sys.path.insert(0, str(HERE.parent / "90_rotating_the_tangle_ball"))
sys.path.insert(0, str(HERE.parent))
from generate_whitehead_doubles import bishop_frame  # noqa: E402
from generate_mutant_pairs import mutant_pair  # noqa: E402
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402
from the_knotfamily_machine.fseries_io import resample_uniform  # noqa: E402


def moebius_cable(X, rho, n_sigma=8192):
    """(2,·) band boundary: s = σ·L linear, ψ = πσ — closes after σ=2."""
    X = resample_uniform(np.asarray(X, float), 2048)
    T, e1, e2 = bishop_frame(X)
    n = len(X)
    sig = np.linspace(0, 2, n_sigma, endpoint=False)
    s_idx = (sig * n / 2 * 2) % n              # σ·L as index (L↔n)
    i0 = np.floor(s_idx).astype(int) % n
    i1 = (i0 + 1) % n
    fr = (s_idx - np.floor(s_idx))[:, None]
    Xs = X[i0] * (1 - fr) + X[i1] * fr
    E1s = e1[i0] * (1 - fr) + e1[i1] * fr
    E2s = e2[i0] * (1 - fr) + e2[i1] * fr
    psi = np.pi * sig
    return Xs + rho * (np.cos(psi)[:, None] * E1s
                       + np.sin(psi)[:, None] * E2s)


def make_cable_pair(seed, out_root):
    K1, K2, _, _ = mutant_pair(seed)
    name = f"mc_pair_s{seed}"
    entry = {"name": name, "machine": "93_cabling_the_mutants",
             "seed": seed,
             "predictions": {
                 "zuendung": "Kabel-Invarianten trennen das Paar? "
                             "(Mutanten-Detektor, Dror-Klassiker)",
                 "typ": "(2,·)-Band-Rand-Satellit"}}
    folder = Path(out_root) / name
    built = 0
    for tag, K in (("K", K1), ("Kmut", K2)):
        clear = embedding_report(resample_uniform(K, 2048),
                                 n_samples=2048)["d_min"]
        rho = 0.22 * clear
        Cb = moebius_cable(K, rho)
        entry[f"rho_{tag}"] = round(float(rho), 5)
        # measure tube fidelity
        Xr = resample_uniform(K, 2048)
        probe = Cb[500:4000:250]
        dmin = np.array([np.linalg.norm(Xr - q, axis=1).min()
                         for q in probe])
        entry[f"tube_dev_{tag}"] = round(
            float(np.abs(dmin - rho).max() / rho), 4)
        emb = embedding_report(Cb)
        if not emb["embedded"]:
            entry["rejected"] = f"embedding gate ({tag})"
            return entry, None
        folder.mkdir(parents=True, exist_ok=True)
        rep = None
        for nt in (192, 256, 384):
            try:
                rep = polyline_to_fseries(
                    Cb, folder / f"{name}_{tag}.fseries",
                    knot_name=f"{name}_{tag}",
                    source="maschinen_ideen/93", n_terms=nt,
                    nyquist_max_rel=0.02)
                break
            except RuntimeError:
                continue
        if rep is None:
            entry["rejected"] = f"Nyquist rot ({tag})"
            return entry, None
        entry[f"fseries_{tag}"] = rep
        built += 1
    if built == 2:
        (folder / f"{name}_entry.json").write_text(
            json.dumps(entry, indent=1))
        return entry, folder
    return entry, None


def cmd_demo(args):
    print("[DEMO] Smoke test machine 93 (mutant cables)")
    # (1) closure geometry: cable start meets cable end (ψ closure)
    K1, _, _, _ = mutant_pair(0)
    Cb = moebius_cable(K1, 0.02)
    gap = float(np.linalg.norm(Cb[0] - Cb[-1]))
    step = float(np.linalg.norm(Cb[1] - Cb[0]))
    ok1 = gap < 3 * step
    print(f"  (1) Möbius closure: residual gap {gap:.2e} ≈ step size "
          f"{step:.2e}  {'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) double cover: every companion point is passed TWICE
    #     (σ and σ+1 — band-boundary signature)
    Xr = resample_uniform(K1, 512)
    n_near = []
    ring = np.arange(512)
    for i in range(0, 512, 64):
        q = Xr[i]
        # only samples at places with free self-distance (else
        # the 0.03 ball also sees the cable passage of the NEIGHBOR strand)
        ds = np.linalg.norm(Xr - q, axis=1)
        rd = np.minimum(np.abs(ring - i), 512 - np.abs(ring - i))
        ds[rd <= 12] = np.inf
        if ds.min() < 0.09:
            continue
        d = np.linalg.norm(Cb - q, axis=1)
        idx = np.nonzero(d < 0.03)[0]
        if len(idx):
            n_cl = 1 + int((np.diff(idx) > 50).sum())
            # Wrap seam: cluster across the array end otherwise counts twice
            if idx[0] == 0 and idx[-1] == len(Cb) - 1 and n_cl > 1:
                n_cl -= 1
            n_near.append(n_cl)
    ok2 = len(n_near) >= 4 and all(c == 2 for c in n_near)
    print(f"  (2) band-boundary signature: {n_near} passages per sample "
          f"(target 2)  {'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) cable PAIR built on mutants (scan seeds)
    e3 = f3 = None
    for s in range(20):
        e3, f3 = make_cable_pair(s, HERE / "demo_smoke")
        if f3 is not None:
            break
    ok3 = f3 is not None
    print(f"  (3) {e3['name']}: ρ={e3.get('rho_K')}/{e3.get('rho_Kmut')}, "
          f"tube fidelity {e3.get('tube_dev_K')}/{e3.get('tube_dev_Kmut')}, "
          f"nyq={e3.get('fseries_K', {}).get('nyquist_dev_rel')}/"
          f"{e3.get('fseries_Kmut', {}).get('nyquist_dev_rel')}  "
          f"{'✓ satellite pair in place' if ok3 else '✗ FAIL: ' + str(e3.get('rejected'))}")
    assert ok3
    print("[DEMO] 3/3 green — cable machine ready to fire "
          "(the mutant detector waits for ignition).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
