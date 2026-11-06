# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 115 — forging_the_rotoreflection (S_2k forge).

Generalisation of the 91 point symmetry (S_2) to ROTOREFLECTION
S_2k: σ = Rz(π/k)∘(z→−z), σ-invariant curve ⟹ AMPHICHIRAL BY
GEOMETRY (σ is an orientation-reversing isometry that maps the knot
onto itself — the certificate is the construction itself).
SPECTRALLY EXACT: xy as c(t)=Σ aₙe^{int} with n ≡ 1 mod 2k;
z with n ≡ k mod 2k. σ² = Rz(2π/k) ⟹ Z/k rotation symmetry included.
THIS IS THE double-6 CLASS (amphichiral + cyclic, k=3) AS
SERIES FACTORY — Stoimenow track (amphi specialist) + bridge
sculpture↔machine. Invertibility NOT imposed → candidates
for amphichiral-non-invertible.
BETS EXACT today: σ invariance bit-exact on a 2k-divisible grid;
Spectral support verified; series k=2,3,4 embedding-gated
written. FIRING BET: sig=0 for all; θ amphi-consistent;
6D fingerprint sees the S_2k class. Doctrine: do NOT fire.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402

N = 4320      # divisible by 2k for k=2,3,4 (and 5,6)


def forge(k, seed, n_xy=3, n_z=2, decay=0.55):
    """σ-invariant curve: xy harmonics n≡1 mod 2k, z: n≡k mod 2k."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0, 2 * np.pi, N, endpoint=False)
    # xy: complex harmonics, both directions 1±2k·j
    ns_xy = [1]
    for j in range(1, n_xy + 1):
        ns_xy += [1 + 2 * k * j, 1 - 2 * k * j]
    c = np.zeros(N, complex)
    for i, n in enumerate(ns_xy):
        amp = (1.0 if n == 1 else decay ** abs((n - 1) // (2 * k)))
        a = amp * (rng.normal() + 1j * rng.normal()) * (0.45 if n != 1 else 1.2)
        c += a * np.exp(1j * n * t)
    # z: reelle Harmonische n ≡ k mod 2k
    z = np.zeros(N)
    for j in range(n_z):
        n = k + 2 * k * j
        z += (decay ** j) * 0.7 * (rng.normal() * np.cos(n * t)
                                   + rng.normal() * np.sin(n * t))
    return np.column_stack([c.real, c.imag, z]), ns_xy


def sigma_image(P, k):
    """σ = Rz(π/k) ∘ (z→−z), applied to the point row."""
    a = np.pi / k
    R = np.array([[np.cos(a), -np.sin(a), 0],
                  [np.sin(a), np.cos(a), 0], [0, 0, -1.0]])
    return P @ R.T


def make_knot(k, seed, out_root):
    P, ns = forge(k, seed)
    name = f"rr_S{2*k}_s{seed}"
    entry = {"name": name, "machine": "115_forging_the_rotoreflection",
             "k": k, "group": f"S_{2*k}", "harmonics_xy": ns,
             "predictions": {
                 "amphichiral": "JA per Konstruktion (σ orientierungs-"
                                "umkehrend, Kurve σ-invariant)",
                 "rotation": f"Z/{k} via σ²",
                 "zuendung": "sig=0 + θ amphi-konsistent + 6D sieht "
                             "S-Klasse; Invertierbarkeit OFFEN "
                             "(double-6-Klasse, Stoimenow)"}}
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = None
    for nt in (max(32, 1 + 2 * k * 3), 128, 256):
        try:
            rep = polyline_to_fseries(P, folder / f"{name}.fseries",
                                      knot_name=name,
                                      source="maschinen_ideen/115",
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
    print("[DEMO] smoke test machine 115 (S_2k rotoreflection forge)")
    # (1) σ invariance BIT: σ(P(t)) == P(t + π/k) on the grid
    worst = 0.0
    for k in (2, 3, 4):
        P, _ = forge(k, seed=7)
        shift = N // (2 * k)
        d = float(np.abs(sigma_image(P, k)
                         - np.roll(P, -shift, axis=0)).max())
        worst = max(worst, d)
    ok1 = worst < 1e-12
    print(f"  (1) σ(P(t)) == P(t+π/k) for k=2,3,4: max|Δ| = {worst:.1e}"
          f"  {'✓ Amphi-Zertifikat = Konstruktion' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) Spectral support: FFT only on allowed frequencies
    k = 3
    P, ns = forge(k, seed=7)
    spec_xy = np.fft.fft(P[:, 0] + 1j * P[:, 1]) / N
    live = {int(f) for f in np.nonzero(np.abs(spec_xy) > 1e-10)[0]}
    live = {f if f <= N // 2 else f - N for f in live}
    allowed = set(ns)
    ok2 = live <= allowed
    spec_z = np.fft.fft(P[:, 2]) / N
    live_z = {int(f) for f in np.nonzero(np.abs(spec_z) > 1e-10)[0]}
    live_z = {abs(f if f <= N // 2 else f - N) for f in live_z}
    ok2 = ok2 and all(f % (2 * k) == k for f in live_z)
    print(f"  (2) Spectral support k=3: xy on {sorted(live)} "
          f"(⊆ n≡1 mod 6), z on {sorted(live_z)} (≡3 mod 6)  "
          f"{'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) series k=2,3,4 written, embedding-gated
    built = []
    for k in (2, 3, 4):
        for seed in (7, 8, 9, 10):
            e, f = make_knot(k, seed, HERE / "demo_smoke")
            if f is not None:
                built.append(e["name"])
                break
    ok3 = len(built) == 3
    print(f"  (3) series: {built}  "
          f"{'✓ S_4/S_6/S_8 amphis in place' if ok3 else '✗ FAIL'}")
    assert ok3
    print("[DEMO] 3/3 green — rotoreflection forge ready to fire "
          "(certified amphichirality in series, double-6 class).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
