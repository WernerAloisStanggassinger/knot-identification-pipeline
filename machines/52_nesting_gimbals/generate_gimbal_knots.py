# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 52 — nesting_gimbals (nested rotations, skew axes).

P(t) = R₁(ω₁t)·(a₁ + R₂(ω₂t)·(a₂ + R₃(ω₃t)·a₃)) — clockwork in 3D.
FINITELY HARMONIC by theorem (products of trig = sums): the spectral
support lies in the ω-sum lattice {|±ω₁±ω₂±ω₃|, sub-chains} — measured
to 1e-12. Axis GEOMETRY = new parameter space beyond the 01 world
(there everything is axis-parallel). BET: coplanar/parallel axes collapse
onto a known stratum (z-support classes measurable via the 31 support),
skew axes = new stratum. Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from itertools import product
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402

N_SAMPLES = 8192


def rotmat(axis, ang):
    """Rodrigues, vectorized over ang (…,)→(…,3,3)."""
    axis = np.asarray(axis, float)
    axis = axis / np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    ang = np.asarray(ang, float)[..., None, None]
    I = np.eye(3)
    return I + np.sin(ang) * K + (1 - np.cos(ang)) * (K @ K)


def gimbal_curve(stages, n=N_SAMPLES):
    """stages: list of (axis(3,), omega, arm(3,)) from outside to inside.
    P(t) = R₁(ω₁t)(a₁ + R₂(ω₂t)(a₂ + …))."""
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    P = np.zeros((n, 3))
    for (axis, om, arm) in reversed(stages):
        R = rotmat(axis, om * t)                       # (n,3,3)
        P = np.einsum("nij,nj->ni", R, np.asarray(arm, float) + P)
    return P


def freq_lattice(omegas):
    """Allowed frequencies: |±ω₁ ±ω₂ … (sub-chains from outside)| — stage k
    uses ω₁..ω_k. Constant (0) belongs to it (arm offsets)."""
    allowed = {0}
    for k in range(1, len(omegas) + 1):
        for signs in product([-1, 0, 1], repeat=k):
            if signs and signs[-1] == 0:
                # innermost stage of the sub-chain must rotate, else duplicate
                continue
            allowed.add(abs(int(sum(s * w for s, w in zip(signs, omegas[:k])))))
    return allowed


def spectral_leak(P, allowed, tol_rel=1e-12):
    """Energy outside the allowed lattice, relative to total energy."""
    n = len(P)
    spec = np.fft.rfft(P - P.mean(axis=0), axis=0) / n
    E = (np.abs(spec) ** 2).sum()
    leak = sum(float((np.abs(spec[f]) ** 2).sum())
               for f in range(spec.shape[0]) if f not in allowed)
    return leak / E, leak / E < tol_rel


def z_support(P, eps_rel=1e-9):
    """Frequency support of the z axis (31 class language)."""
    n = len(P)
    spec = np.abs(np.fft.rfft(P[:, 2] - P[:, 2].mean())) / n
    m = spec.max()
    return sorted(int(f) for f in np.nonzero(spec > eps_rel * max(m, 1e-30))[0])


def make_knot(name, stages, out_root):
    P = gimbal_curve(stages)
    omegas = [int(s[1]) for s in stages]
    allowed = freq_lattice(omegas)
    leak, finite = spectral_leak(P, allowed)
    entry = {"name": name, "machine": "52_nesting_gimbals",
             "omegas": omegas,
             "axes": [[round(float(x), 6) for x in s[0]] for s in stages],
             "arms": [[round(float(x), 6) for x in s[2]] for s in stages],
             "spectral_leak_rel": float(leak),
             "finite_harmonic": bool(finite),
             "z_support": z_support(P),
             "predictions": {
                 "satz": "Support ⊂ ω-Summen-Gitter (gemessen)",
                 "stratum": "schiefe Achsen = neu; koplanar kollabiert",
                 "identity": "Zünd-Messung"}}
    if not finite:
        entry["rejected"] = f"Harmonik-Gate (leak={leak:.2e})"
        return entry, None
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    n_terms = max(64, 2 * max(allowed) + 8)
    rep = polyline_to_fseries(P, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/52", n_terms=n_terms,
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def skew_axes(seed=0):
    """Three deliberately skew (non-coplanar, non-orthogonal) axes."""
    rng = np.random.default_rng(seed)
    A = rng.normal(size=(3, 3))
    return [a / np.linalg.norm(a) for a in A]


def cmd_demo(args):
    print("[DEMO] Smoke test machine 52 (gimbals)")
    z = np.array([0.0, 0.0, 1.0])
    # (1) finite-harmonics theorem: skew axes, leak < 1e-12
    ax = skew_axes(seed=3)
    stages = [(ax[0], 2, np.array([1.0, 0, 0])),
              (ax[1], 5, np.array([0.45, 0, 0])),
              (ax[2], 11, np.array([0.18, 0, 0]))]
    P = gimbal_curve(stages)
    leak, finite = spectral_leak(P, freq_lattice([2, 5, 11]))
    print(f"  (1) finite harmonics (skew, ω=2/5/11): leak={leak:.2e}  "
          f"{'✓' if finite else '✗ FAIL'}")
    assert finite
    # (2) Lattice is SHARP: without the sum frequencies the gate tears
    too_small = {0, 2, 5, 11}
    leak2, fin2 = spectral_leak(P, too_small)
    ok2 = not fin2
    print(f"  (2) lattice sharp (only fundamental frequencies allowed): "
          f"leak={leak2:.2e} → gate tears  {'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) Axis-parallel collapse: all axes z ⇒ z-support only {0}
    stages_z = [(z, 2, np.array([1.0, 0, 0.3])),
                (z, 5, np.array([0.45, 0, 0.2])),
                (z, 11, np.array([0.18, 0, 0.1]))]
    zs = z_support(gimbal_curve(stages_z))
    ok3 = zs == [0]
    print(f"  (3) coplanar collapse: z-support={zs} (arm z constant — "
          f"01 world class)  {'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    # (4) Skew stratum: z-support carries mixed frequencies
    zs4 = z_support(P)
    ok4 = len(zs4) > 3
    print(f"  (4) skew stratum: z-support={zs4[:8]}…  "
          f"({len(zs4)} frequencies)  {'✓' if ok4 else '✗ FAIL'}")
    assert ok4
    # (5) Knot built + gates green
    e5, f5 = make_knot("gb_skew_2_5_11_s3", stages, HERE / "demo_smoke")
    ok5 = f5 is not None
    print(f"  (5) gimbal knot: embedded="
          f"{e5.get('embedding', {}).get('embedded')}, "
          f"nyq={e5.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok5 else '✗ FAIL: ' + str(e5.get('rejected'))}")
    assert ok5
    print("[DEMO] 5/5 green — gimbal machine ready to fire "
          "(axis geometry = new parameter space, theorem measured).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
