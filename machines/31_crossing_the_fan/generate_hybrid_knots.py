# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 31 — crossing_the_fan (machine hybrids).

Mixer: xy coefficients from source A + z coefficients from source B
(z normalised to the xy span, flaw 2). INTERSECTION PREDICTION (12 style,
exact): the hybrid is C_n-symmetric ⇔ xy classes ⊆ ±1 mod n AND
z classes ⊆ 0 mod n — computed before the build, measured after the build at the
residual. Doctrine: do NOT ignite before Monday.
"""
from __future__ import annotations
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "the_knotfamily_machine"))
sys.path.insert(0, str(HERE.parent))
from fseries_io import reconstruct_from_fourier, write_fseries  # noqa: E402
from machine_common import min_nonadjacent_dist, resample_uniform  # noqa: E402


def load_fseries(path):
    C = np.loadtxt(path, comments="%", ndmin=2)
    return C[:, [0, 2, 4]], C[:, [1, 3, 5]]


def support(a, b, axes, eps=1e-9):
    amp = np.hypot(a[:, axes], b[:, axes]).sum(axis=1)
    return set((np.nonzero(amp > eps)[0] + 1).tolist())


def signed_support(a, b, eps=1e-9):
    """SIGNED xy harmonics j ∈ ±k (build lesson: cos axes
    occupy BOTH signs — |k| classes are not enough; Lissajous
    correctly failed the residual, the predicate was too lax)."""
    out = set()
    for k in range(1, a.shape[0] + 1):
        ax, ay = a[k - 1, 0], a[k - 1, 1]
        bx, by = b[k - 1, 0], b[k - 1, 1]
        c_pos = complex(ax + by, ay - bx) / 2
        c_neg = complex(ax - by, ay + bx) / 2
        if abs(c_pos) > eps:
            out.add(k)
        if abs(c_neg) > eps:
            out.add(-k)
    return out


def surviving_cn(a, b, n_max=12):
    """Exact: C_n ⇔ ∃m, gcd(m,n)=1: all signed xy j ≡ m mod n,
    all z k ≡ 0 mod n."""
    J = signed_support(a, b)
    zz = support(a, b, [2])
    out = []
    for n in range(2, n_max + 1):
        if any(k % n for k in zz):
            continue
        residues = {j % n for j in J}
        if len(residues) == 1 and math.gcd(residues.pop(), n) == 1:
            out.append(n)
    return out


def cn_residual(a, b, n, n_check=512):
    t = np.linspace(0, 2 * np.pi, n_check, endpoint=False)
    nt = a.shape[0]

    def at(tv):
        out = np.zeros((len(tv), 3))
        for k in range(1, nt + 1):
            out += np.outer(np.cos(k * tv), a[k - 1]) \
                 + np.outer(np.sin(k * tv), b[k - 1])
        return out
    al = 2 * np.pi / n
    R = np.array([[np.cos(al), -np.sin(al), 0],
                  [np.sin(al), np.cos(al), 0], [0, 0, 1]])
    P = at(t)
    diam = float(np.linalg.norm(P.max(0) - P.min(0)))
    return float(np.abs(at(t + al) - P @ R.T).max()) / diam


def hybrid(a_xy, b_xy, a_z, b_z):
    n = max(a_xy.shape[0], a_z.shape[0])

    def pad(m, rows):
        out = np.zeros((rows, 3))
        out[:m.shape[0]] = m
        return out
    A, B = pad(a_xy, n), pad(b_xy, n)
    Az, Bz = pad(a_z, n), pad(b_z, n)
    a = np.zeros((n, 3))
    b = np.zeros((n, 3))
    a[:, :2], b[:, :2] = A[:, :2], B[:, :2]
    # normalise z to the xy span (flaw 2)
    P_xy = reconstruct_from_fourier(A, B, n_samples=1024)
    span_xy = float(np.linalg.norm(P_xy[:, :2].max(0) - P_xy[:, :2].min(0)))
    z_amp = np.hypot(Az[:, 2], Bz[:, 2]).sum()
    sc = 0.35 * span_xy / max(z_amp, 1e-12)
    a[:, 2], b[:, 2] = sc * Az[:, 2], sc * Bz[:, 2]
    return a, b


def make_knot(name, a, b, parents, out_root):
    pred = surviving_cn(a, b)
    res = {n: cn_residual(a, b, n) for n in range(2, 8)}
    measured = [n for n, r in res.items() if r < 1e-9]
    entry = {"name": name, "machine": "31_crossing_the_fan",
             "parents": parents,
             "predicted_cn": pred,
             "measured_cn": measured,
             "prediction_holds": bool(set(pred) >= set(measured)
                                      and set(measured) >= {n for n in pred
                                                            if n < 8})}
    P = reconstruct_from_fourier(a, b, n_samples=4096)
    d_min, _ = min_nonadjacent_dist(resample_uniform(P, 4096))
    diam = float(np.linalg.norm(P.max(0) - P.min(0)))
    entry["embedding"] = {"d_min_rel": round(d_min / diam, 6),
                          "embedded": bool(d_min / diam >= 0.003)}
    if not entry["embedding"]["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    write_fseries(folder / f"{name}.fseries", a, b, knot_name=name,
                  source="maschinen_ideen/31_crossing_the_fan")
    return entry, folder


def cmd_demo(args):
    print("[DEMO] smoke test machine 31 (hybrids)")
    sys.path.insert(0, str(HERE.parent / "12_filling_the_symmetry_table"))
    sys.path.insert(0, str(HERE.parent / "01_knotting_oscillations"))
    from generate_symmetry_columns import build_coeffs, DEFAULT_SPECS
    from generate_oscillation_knots import lissajous_coeffs
    # Eltern: A = D_3-Kurve (xy-Klassen {1,2,4}), B = Lissajous (3,4,7)
    aA, bA = build_coeffs("D_n", 3, DEFAULT_SPECS["D_n"])
    aB, bB = lissajous_coeffs((3, 4, 7), (0.7, 0.2, 0.0))
    # 1) intersection prediction: A-xy {1,2,4} + B-z {7}: C_n needs 7≡0 mod n
    #    AND {1,2,4} ⊆ ±1 mod n ⇒ n=7? 1,2,4 mod 7: 1,2,4 — 2,4 ∉ {1,6} ✗
    #    ⇒ prediction: NO C_n survives (empty)
    ah, bh = hybrid(aA, bA, aB, bB)
    pred = surviving_cn(ah, bh)
    ok1 = pred == []
    print(f"  (1) D_3-xy × Lissajous-z(7): predicted intersection = {pred} "
          f"(empty — 7 does not divide, classes collide)  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # 2) constructive case: A-xy {1,2,4} + A's own z {3,6} ⇒ C_3 survives
    ah2, bh2 = hybrid(aA, bA, aA, bA)
    pred2 = surviving_cn(ah2, bh2)
    res3 = cn_residual(ah2, bh2, 3)
    ok2 = 3 in pred2 and res3 < 1e-12
    print(f"  (2) self-hybrid: intersection {pred2}, C_3 residual {res3:.1e}  "
          f"{'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # 3) Lissajous occupies ±k (signed lesson) ⇒ intersection MUST be empty,
    #    although |k| classes would fit — predicate hardening measured
    aC, bC = lissajous_coeffs((5, 7, 6), (0.9, 0.4, 0.0))
    ah3, bh3 = hybrid(aC, bC, aA, bA)
    pred3 = surviving_cn(ah3, bh3)
    ok3 = pred3 == [] and cn_residual(ah3, bh3, 3) > 1e-6
    print(f"  (3) Lissajous(5,7)-xy × D_3-z: Schnitt={pred3} (leer — "
          f"signed ±k collides), residual confirmed  "
          f"{'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    # 4) REAL cross-machine hybrid with survivor: 12×10 —
    #    D_3-xy (rotating pairs, signed ≡ 1 mod 3) + T(2,3)-z ({3})
    sys.path.insert(0, str(HERE.parent / "10_climbing_the_torus_ladder"))
    from generate_torus_ladder import torus_coeffs
    aT, bT = torus_coeffs(2, 3)
    ah4, bh4 = hybrid(aA, bA, aT, bT)
    pred4 = surviving_cn(ah4, bh4)
    e4, f4 = make_knot("hy_d3xy_t23z", ah4, bh4, ["D_3 (12)", "T(2,3)-z (10)"],
                       HERE / "demo_smoke")
    ok4 = 3 in pred4 and 3 in e4["measured_cn"] and f4 is not None
    print(f"  (4) D_3-xy × T(2,3)-z (12×10): Schnitt={pred4}, gemessen="
          f"{e4['measured_cn']}, embedded={e4['embedding']['embedded']}  "
          f"{'✓' if ok4 else '✗ FAIL: ' + str(e4.get('rejected'))}")
    assert ok4
    print("[DEMO] 4/4 green — hybrid machine ready to fire (signed-"
          "intersection theorem predicted+measured, 12×10 hybrid built).")


def main():
    ap = argparse.ArgumentParser(description="Machine 31 — hybrids")
    sub = ap.add_subparsers(dest="cmd", required=True)
    pm = sub.add_parser("make")
    pm.add_argument("fseries_xy")
    pm.add_argument("fseries_z")
    sub.add_parser("demo")
    args = ap.parse_args()
    if args.cmd == "demo":
        cmd_demo(args)
    else:
        aA, bA = load_fseries(args.fseries_xy)
        aB, bB = load_fseries(args.fseries_z)
        a, b = hybrid(aA, bA, aB, bB)
        name = (f"hy_{Path(args.fseries_xy).stem}"
                f"__{Path(args.fseries_z).stem}")
        man_path = HERE / "MANIFEST.json"
        man = json.loads(man_path.read_text()) if man_path.exists() \
            else {"knots": []}
        entry, folder = make_knot(name, a, b,
                                  [args.fseries_xy, args.fseries_z],
                                  HERE / "out_fseries")
        man["knots"] = [k for k in man["knots"]
                        if k["name"] != entry["name"]] + [entry]
        man_path.write_text(json.dumps(man, indent=1))
        print(f"[MAKE] {name}: Schnitt {entry['predicted_cn']} → "
              f"{folder if folder else entry.get('rejected')}")


if __name__ == "__main__":
    main()
