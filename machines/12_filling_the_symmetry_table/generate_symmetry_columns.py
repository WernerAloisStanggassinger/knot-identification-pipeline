# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 12 — filling_the_symmetry_table (new axis-family columns).

Selection rules DERIVED (founding txt, 2026-06-06), build in 01 style
(exact coefficients, no sampling):

  D_n  : x = Σ p_k cos(kt) · y = Σ s_k p_k sin(kt) (s_k = ±1 by class
         k ≡ ±1 mod n) · z = Σ r_k sin(kt), k ≡ 0 mod n
         ⇒ C_2 axis meets the curve at t=0,π ⇒ STRONGLY INVERTIBLE exactly.
  C_nv : like D_n, z = Σ r_k cos(kt) ⇒ spatial mirror ⇒ AMPHICHIRAL exactly.
  C_nh : ONLY n odd (theorem in the txt: reversal variant impossible for n>=3);
         xy rotating pairs k ≡ ±1 mod n AND even, phases free;
         z: k ≡ 0 mod n AND odd ⇒ σ_h via half-circuit ⇒ AMPHICHIRAL.

Every build verifies its symmetry NUMERICALLY (isometry + parameter map,
machine epsilon — 03 doctrine: measure, do not believe). D_nd/D_nh =
Expansion stage. Doctrine: do NOT fire before Monday.
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
from fseries_io import write_fseries, reconstruct_from_fourier  # noqa: E402
from machine_common import min_nonadjacent_dist, resample_uniform  # noqa: E402

N_TERMS = 24
N_CHECK = 1024


def build_coeffs(column, n, spec):
    """spec: list (k, amp[, phase]) for xy + list (k, amp) for z.

    Returns (a, b) in fseries format. Raises on rule violations.
    """
    xy, zz = spec["xy"], spec["z"]
    a = np.zeros((N_TERMS, 3))
    b = np.zeros((N_TERMS, 3))
    for item in xy:
        k, amp = item[0], item[1]
        cls = k % n
        if cls not in (1, n - 1):
            raise ValueError(f"xy harmonic {k}: class {cls} ∉ ±1 mod {n}")
        if column == "C_nh":
            if k % 2 != 0:
                raise ValueError(f"C_nh: xy-k={k} must be even")
            phase = item[2] if len(item) > 2 else 0.0
            s = +1.0 if cls == 1 else -1.0
            # rotating pair with free phase
            a[k - 1, 0] += amp * np.cos(phase)
            b[k - 1, 0] += amp * np.sin(phase)
            a[k - 1, 1] += -s * amp * np.sin(phase)
            b[k - 1, 1] += s * amp * np.cos(phase)
        else:                       # D_n / C_nv: pure cos/sin, class sign
            s = +1.0 if cls == 1 else -1.0
            a[k - 1, 0] += amp          # x: cos
            b[k - 1, 1] += s * amp      # y: ±sin
    for (k, amp) in zz:
        if k % n != 0:
            raise ValueError(f"z harmonic {k} ≢ 0 mod {n}")
        if column == "C_nh":
            if k % 2 == 0:
                raise ValueError(f"C_nh: z-k={k} must be odd")
            b[k - 1, 2] += amp          # phase free would be ok; sin chosen
        elif column == "D_n":
            b[k - 1, 2] += amp          # sin (odd in t)
        else:                           # C_nv
            a[k - 1, 2] += amp          # cos (even in t)
    return a, b


def symmetry_residual(column, n, a, b):
    """Numerical proof: isometry ∘ curve == curve ∘ parameter map."""
    t = np.linspace(0, 2 * np.pi, N_CHECK, endpoint=False)
    P = reconstruct_from_fourier(a, b, n_samples=N_CHECK)
    diam = float(np.linalg.norm(P.max(0) - P.min(0)))

    def at(tv):
        out = np.zeros((len(tv), 3))
        for k in range(1, N_TERMS + 1):
            out += np.outer(np.cos(k * tv), a[k - 1]) \
                 + np.outer(np.sin(k * tv), b[k - 1])
        return out

    checks = {}
    # C_n part: rotation by 2π/n == shift
    al = 2 * np.pi / n
    R = np.array([[np.cos(al), -np.sin(al), 0],
                  [np.sin(al), np.cos(al), 0], [0, 0, 1]])
    checks["C_n"] = float(np.abs(at(t + al) - at(t) @ R.T).max()) / diam
    if column == "D_n":     # (x,−y,−z) == t→−t
        M = np.diag([1.0, -1.0, -1.0])
        checks["C2x"] = float(np.abs(at(-t) - at(t) @ M).max()) / diam
    if column == "C_nv":    # (x,−y,z) == t→−t
        M = np.diag([1.0, -1.0, 1.0])
        checks["sigma_v"] = float(np.abs(at(-t) - at(t) @ M).max()) / diam
    if column == "C_nh":    # (x,y,−z) == t→t+π
        M = np.diag([1.0, 1.0, -1.0])
        checks["sigma_h"] = float(np.abs(at(t + np.pi) - at(t) @ M).max()) / diam
    return checks


DEFAULT_SPECS = {
    "D_n": {"xy": [(1, 1.0), (4, 0.35), (2, 0.3)],
            "z": [(3, 0.5), (6, 0.2)]},
    "C_nv": {"xy": [(1, 1.0), (4, 0.35), (2, 0.3)],
             "z": [(3, 0.5), (6, 0.2)]},
    "C_nh": {"xy": [(2, 1.0, 0.0), (4, 0.4, 1.1), (8, 0.25, 0.4)],
             "z": [(3, 0.5), (9, 0.2)]},
}


def make_knot(column, n, out_root, spec=None):
    if column in ("C_nh", "D_nh"):
        raise ValueError(
            "THEOREM (2026-06-06): for n>=3, C_nh/D_nh are IMPOSSIBLE as the symmetry "
            "of an embedded curve — the σ_h half-circuit forces a "
            "self-intersection at every zero of z (see founding txt)")
    spec = DEFAULT_SPECS[column] if spec is None else spec
    if column in ("D_n", "C_nv"):       # adapt classes to n (default n=3)
        spec = {"xy": [(k, A) for (k, A, *r) in
                       [tuple(x) for x in spec["xy"]]
                       if k % n in (1, n - 1)],
                "z": [(k, A) for (k, A) in spec["z"] if k % n == 0]}
    a, b = build_coeffs(column, n, spec)
    res = symmetry_residual(column, n, a, b)
    name = f"sym_{column.replace('_n', str(n))}"
    bets = {"D_n": {"strongly_invertible": True,
                    "basis": "C2-Achse trifft Kurve t=0,π — exakt"},
            "C_nv": {"amphichiral_unoriented": True,
                     "basis": "Raumspiegel in der Gruppe — exakt"},
            "C_nh": {"amphichiral_unoriented": True,
                     "basis": "σ_h via Halbumlauf (Satz: nur n ungerade)"}}
    entry = {"name": name, "machine": "12_filling_the_symmetry_table",
             "column": column, "n": n,
             "symmetry_residuals": {k: float(v) for k, v in res.items()},
             "symmetry_exact": bool(max(res.values()) < 1e-12),
             "predictions": bets[column]}
    if not entry["symmetry_exact"]:
        entry["rejected"] = f"Symmetrie-Residuum {max(res.values()):.2e}"
        return entry, None
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
                  source="maschinen_ideen/12_filling_the_symmetry_table")
    return entry, folder


def cmd_make(args):
    man_path = HERE / "MANIFEST.json"
    man = json.loads(man_path.read_text()) if man_path.exists() else {"knots": []}
    entry, folder = make_knot(args.column, args.n, HERE / "out_fseries")
    man["knots"] = [k for k in man["knots"] if k["name"] != entry["name"]] + [entry]
    man_path.write_text(json.dumps(man, indent=1))
    print(f"[MAKE] {entry['name']}: residuals={entry['symmetry_residuals']} → "
          f"{folder if folder else entry.get('rejected')}")


def cmd_demo(args):
    print("[DEMO] Smoke test machine 12 (symmetry columns)")
    # 1) D_3: symmetry exact + strongly-invertible bet
    e1, f1 = make_knot("D_n", 3, HERE / "demo_smoke")
    ok1 = f1 is not None and e1["symmetry_exact"]
    print(f"  (1) D_3: residuals={e1['symmetry_residuals']} embedded="
          f"{e1.get('embedding', {}).get('embedded')}  "
          f"{'✓' if ok1 else '✗ FAIL: ' + str(e1.get('rejected'))}")
    assert ok1
    # 2) C_3v: mirror exact (amphi factory)
    e2, f2 = make_knot("C_nv", 3, HERE / "demo_smoke")
    ok2 = f2 is not None and e2["symmetry_exact"]
    print(f"  (2) C_3v: residuals={e2['symmetry_residuals']}  "
          f"{'✓' if ok2 else '✗ FAIL: ' + str(e2.get('rejected'))}")
    assert ok2
    # 3) MEASURE IMPOSSIBILITY THEOREM: C_3h half-circuit curve has d_min == 0
    #    exactly (self-intersection at zeros of z) — 200/200 in the search run
    a3, b3 = build_coeffs("C_nh", 3, DEFAULT_SPECS["C_nh"])
    P3 = reconstruct_from_fourier(a3, b3, n_samples=4096)
    d3, _ = min_nonadjacent_dist(resample_uniform(P3, 4096))
    ok3 = d3 < 1e-9
    print(f"  (3) THEOREM measured: C_3h candidate d_min = {d3:.2e} (exactly 0 — "
          f"self-intersection forced)  {'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    # 4) Theorem gate: C_nh/D_nh build refuses with justification
    try:
        make_knot("C_nh", 3, HERE / "demo_smoke")
        ok4 = False
    except ValueError:
        ok4 = True
    print(f"  (4) C_nh build refused (THEOREM)  {'✓' if ok4 else '✗ FAIL'}")
    assert ok4
    # 5) Rule gate: wrong class is thrown out
    try:
        build_coeffs("D_n", 3, {"xy": [(3, 1.0)], "z": [(3, 0.5)]})
        ok5 = False
    except ValueError:
        ok5 = True
    print(f"  (5) xy-k=3 in D_3 (class 0) refused  "
          f"{'✓' if ok5 else '✗ FAIL'}")
    assert ok5
    print("[DEMO] 5/5 green — TWO new columns ready to fire, ONE "
          "impossibility theorem measured (D_nd = last open column).")


def main():
    ap = argparse.ArgumentParser(description="Machine 12 — symmetry columns")
    sub = ap.add_subparsers(dest="cmd", required=True)
    pm = sub.add_parser("make")
    pm.add_argument("column", choices=["D_n", "C_nv"])
    pm.add_argument("--n", type=int, default=3)
    sub.add_parser("demo")
    args = ap.parse_args()
    {"make": cmd_make, "demo": cmd_demo}[args.cmd](args)


if __name__ == "__main__":
    main()
