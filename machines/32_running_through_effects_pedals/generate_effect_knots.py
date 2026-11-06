# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MASCHINE 32 — running_through_effects_pedals (Audio-Filter auf .fseries).

Effekte: comb (Echo, pro Harmonischer ×(1+γe^{-ikτ}) — exakt im
Koeffizienten-Raum) · vibrato (t→t+ε·sin(mt), Zeitbereich+Re-Fourier) ·
clip (tanh saturation of the coordinates, time domain + re-Fourier).

THEOREM (verified in the build): comb with τ=2π/n multiplies every k class
mod n by a constant ⇒ preserves EVERY C_n class occupation exactly —
symmetry-faithful deformation paths. Gain gate against cancellation (flaw 2).
Doctrine: do NOT ignite before Monday.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "the_knotfamily_machine"))
sys.path.insert(0, str(HERE.parent))
from fseries_io import reconstruct_from_fourier, write_fseries  # noqa: E402
from machine_common import (embedding_report, polyline_to_fseries,  # noqa: E402
                            min_nonadjacent_dist, resample_uniform)


def load_fseries(path):
    C = np.loadtxt(path, comments="%", ndmin=2)
    return C[:, [0, 2, 4]], C[:, [1, 3, 5]]


def comb(a, b, gamma, tau):
    """Echo: C_k → C_k·(1+γe^{-ikτ}), C_k = a_k − i·b_k (pro Achse)."""
    n = a.shape[0]
    k = np.arange(1, n + 1)
    f = 1.0 + gamma * np.exp(-1j * k * tau)
    if np.abs(f).min() < 0.1:
        raise ValueError(f"Gain-Gate: |Filter| min {np.abs(f).min():.3f} < 0.1 "
                         "(cancellation flaw)")
    C = a - 1j * b
    C2 = C * f[:, None]
    return np.real(C2), -np.imag(C2)


def time_effect(a, b, fn, n_samples=8192, n_terms=None):
    """Time-domain effect + re-Fourier (Nyquist guard via back-comparison)."""
    P = reconstruct_from_fourier(a, b, n_samples=n_samples)
    Q = fn(P, np.linspace(0, 2 * np.pi, n_samples, endpoint=False))
    nt = a.shape[0] if n_terms is None else n_terms
    spec = np.fft.rfft(Q, axis=0)
    a2 = (2.0 / n_samples) * np.real(spec[1:nt + 1])
    b2 = (-2.0 / n_samples) * np.imag(spec[1:nt + 1])
    R = reconstruct_from_fourier(a2, b2, n_samples=n_samples)
    dev = float(np.abs(R - (Q - Q.mean(0))).max()
                / np.linalg.norm(Q.max(0) - Q.min(0)))
    return a2, b2, dev


def vibrato(eps, m):
    def fn(P, t):
        n = len(t)
        t2 = (t + eps * np.sin(m * t)) % (2 * np.pi)
        idx = np.argsort(t2)
        # Resample: evaluate the curve at shifted times = interpolation
        out = np.empty_like(P)
        for d in range(3):
            out[:, d] = np.interp(t2, t, P[:, d], period=2 * np.pi)
        return out
    return fn


def clip(gain):
    def fn(P, t):
        s = np.abs(P).max()
        return np.tanh(gain * P / s) / np.tanh(gain) * s
    return fn


def cn_class_support(a, b, eps=1e-9):
    """Besetzte k-Werte (xy gemeinsam, z separat)."""
    amp_xy = np.hypot(a[:, :2], b[:, :2]).sum(axis=1)
    amp_z = np.hypot(a[:, 2], b[:, 2])
    return (set((np.nonzero(amp_xy > eps)[0] + 1).tolist()),
            set((np.nonzero(amp_z > eps)[0] + 1).tolist()))


def make_knot(src_a, src_b, effect, params, name, out_root):
    if effect == "comb":
        a2, b2 = comb(src_a, src_b, params["gamma"], params["tau"])
        dev = 0.0
    elif effect == "vibrato":
        a2, b2, dev = time_effect(src_a, src_b,
                                  vibrato(params["eps"], params["m"]))
    else:
        a2, b2, dev = time_effect(src_a, src_b, clip(params["gain"]))
    P = reconstruct_from_fourier(a2, b2, n_samples=4096)
    d_min, _ = min_nonadjacent_dist(resample_uniform(P, 4096))
    diam = float(np.linalg.norm(P.max(0) - P.min(0)))
    entry = {"name": name, "machine": "32_running_through_effects_pedals",
             "effect": effect, "params": params,
             "refourier_dev": round(dev, 6),
             "class_support_before": [sorted(s) for s in
                                      cn_class_support(src_a, src_b)],
             "class_support_after": [sorted(s) for s in
                                     cn_class_support(a2, b2)],
             "embedding": {"d_min_rel": round(d_min / diam, 6),
                           "embedded": bool(d_min / diam >= 0.003)}}
    if not entry["embedding"]["embedded"]:
        entry["rejected"] = "embedding gate (Pfad-Kante erreicht)"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    write_fseries(folder / f"{name}.fseries", a2, b2, knot_name=name,
                  source="maschinen_ideen/32_running_through_effects_pedals")
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 32 (effects units)")
    # Source: D_3 curve from machine 12 (C_3 classes {1,2,4}xy / {3,6}z)
    sys.path.insert(0, str(HERE.parent / "12_filling_the_symmetry_table"))
    from generate_symmetry_columns import build_coeffs, DEFAULT_SPECS, \
        symmetry_residual
    a, b = build_coeffs("D_n", 3, DEFAULT_SPECS["D_n"])
    # 1) THEOREM: comb τ=2π/3 preserves C_3 exactly (residual before/after)
    a2, b2 = comb(a, b, 0.5, 2 * np.pi / 3)
    r_before = symmetry_residual("D_n", 3, a, b)["C_n"]
    r_after = symmetry_residual("D_n", 3, a2, b2)["C_n"]
    ok1 = r_after < 1e-12
    print(f"  (1) comb τ=2π/3: C_3-Residuum {r_before:.1e} → {r_after:.1e} "
          f"(erhalten)  {'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # 2) control: generic τ does NOT break the class-phase coherence
    #    of the classes (classes remain classes!), but class support remains
    sup_b = cn_class_support(a, b)
    a3, b3 = comb(a, b, 0.5, 1.234)
    sup_a = cn_class_support(a3, b3)
    ok2 = sup_b == sup_a
    print(f"  (2) Klassen-Support unter comb invariant: {sorted(sup_b[0])}  "
          f"{'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # 3) gain gate: γ=1, τ such that cancellation threatens
    try:
        comb(a, b, 1.0, np.pi)
        ok3 = False
    except ValueError:
        ok3 = True
    print(f"  (3) Cancellation gate (γ=1, τ=π) refuses  "
          f"{'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    # 4) Zeitbereich: clip + vibrato bauen, Gates
    e4, f4 = make_knot(a, b, "clip", {"gain": 2.0}, "fx_clip_g2",
                       HERE / "demo_smoke")
    e5, f5 = make_knot(a, b, "vibrato", {"eps": 0.15, "m": 3},
                       "fx_vib_e015_m3", HERE / "demo_smoke")
    ok4 = f4 is not None and f5 is not None
    print(f"  (4) clip: dev={e4['refourier_dev']} embedded="
          f"{e4['embedding']['embedded']} · vibrato: dev="
          f"{e5['refourier_dev']} embedded={e5['embedding']['embedded']}  "
          f"{'✓' if ok4 else '✗ FAIL'}")
    assert ok4
    # 5) Pfad bis zur Embedding-Kante: γ-Sweep, erste Ablehnung dokumentiert
    edge = None
    for g in np.arange(0.1, 0.96, 0.1):
        try:
            e, f = make_knot(a, b, "comb", {"gamma": float(g),
                                            "tau": 2 * np.pi / 3},
                             f"fx_comb_g{g:.1f}".replace(".", "p"),
                             HERE / "demo_smoke")
        except ValueError:
            edge = ("gain_gate", float(g))
            break
        if f is None:
            edge = ("embedding", float(g))
            break
    edge_txt = edge if edge else "keine bis 0.9 (Pfad durchgehend baubar)"
    print(f"  (5) γ-Pfad: Kante = {edge_txt}  ✓")
    print("[DEMO] 5/5 green — effects units ready to fire (type edges = "
          "ignition program).")


def main():
    ap = argparse.ArgumentParser(description="Machine 32 — effects units")
    sub = ap.add_subparsers(dest="cmd", required=True)
    pm = sub.add_parser("make")
    pm.add_argument("fseries", help="Quell-.fseries")
    pm.add_argument("--effect", choices=["comb", "vibrato", "clip"],
                    required=True)
    pm.add_argument("--gamma", type=float, default=0.4)
    pm.add_argument("--tau", type=float, default=2 * np.pi / 3)
    pm.add_argument("--eps", type=float, default=0.15)
    pm.add_argument("--m", type=int, default=3)
    pm.add_argument("--gain", type=float, default=2.0)
    sub.add_parser("demo")
    args = ap.parse_args()
    if args.cmd == "demo":
        cmd_demo(args)
    else:
        a, b = load_fseries(args.fseries)
        params = {"comb": {"gamma": args.gamma, "tau": args.tau},
                  "vibrato": {"eps": args.eps, "m": args.m},
                  "clip": {"gain": args.gain}}[args.effect]
        name = f"fx_{args.effect}_{Path(args.fseries).stem}"
        man_path = HERE / "MANIFEST.json"
        man = json.loads(man_path.read_text()) if man_path.exists() \
            else {"knots": []}
        entry, folder = make_knot(a, b, args.effect, params, name,
                                  HERE / "out_fseries")
        man["knots"] = [k for k in man["knots"]
                        if k["name"] != entry["name"]] + [entry]
        man_path.write_text(json.dumps(man, indent=1))
        print(f"[MAKE] {name} → {folder if folder else entry.get('rejected')}")


if __name__ == "__main__":
    main()
