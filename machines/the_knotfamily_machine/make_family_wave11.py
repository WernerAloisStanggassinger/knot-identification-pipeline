# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave11.py — Wave 11: fine shear grid at twist=0 (basin transition).

Discovery #180 TODO e: Wave 9/10 showed that shear is the primary topology
selector. At tw=0 we already know:
  sh-1.5, sh-1 → S_24_d49_Z3_chiral (kf103/108/127 all salvageable)
  sh+1         → K12n_837
  sh+1.5       → K12a_1019
Wave 11 sweeps the fine range between ±1 and ±1.5, plus the transition
around sh=0 (S_24 baseline basin), to map the basin boundaries.

Axis design (Discovery #180 TODO e — exact):
  twist  = 0  (fixed)
  shear ∈ {-1.5, -1.25, -1, -0.75, -0.5,
           +0.5, +0.75, +1, +1.25, +1.5}  (10 values)

= 10 variants kf147..kf156.

Expected findings:
  - sh ∈ {-1.5, -1.25, -1, -0.75}: S_24_d49 or K12n_706 (basin expectation)
  - sh ∈ {-0.5, +0.5}: S_24 baseline basin (D6 amphi)
  - sh ∈ {+0.75, +1}: K12n_837
  - sh ∈ {+1.25, +1.5}: K12a_1019 or K12n_839 (bifurcation range)

Failure expectation: shear=-1.5 (kf103) was a failure, sh=-1.0 (kf108, kf127)
also — neg-shear asymmetry from waves 9/10. Wave-11 failures on the neg side
plausible; salvage via skill `/salvage-failures` (stage-5+theta) as standard.

Convention: first shear (z'=z+slope·x), then helix twist=0 — identical to
Wave 8/9/10 (perturb_shear_then_twist).

Output: erzeugte_knoten/kf{147..156}_<descr>/<name>.fseries + wave11_index.json.
Sweep via the separate sweep_wave11.sh.
"""
from __future__ import annotations
import site_config
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from knot_machine import build_inner, rot_z  # noqa: E402
from fseries_io import polyline_to_fourier, write_fseries  # noqa: E402

OUT_BASE = Path((site_config.GENERATED_KNOTS_DIR))
MODULE_NPY = ROOT / "double3_module_extracted.npy"


def perturb_zshear(M: np.ndarray, slope: float) -> np.ndarray:
    Mp = M.copy()
    Mp[1:-1, 2] = M[1:-1, 2] + slope * M[1:-1, 0]
    return Mp


def perturb_helix_twist(M: np.ndarray, twist_total: float) -> np.ndarray:
    Mp = M.copy()
    L = len(M) - 1
    for i in range(1, L):
        ang = twist_total * (i / L)
        R = rot_z(ang)
        Mp[i] = R @ M[i]
    return Mp


def perturb_shear_then_twist(M: np.ndarray, slope: float,
                             twist_total: float) -> np.ndarray:
    Mp = perturb_zshear(M, slope)
    Mp = perturb_helix_twist(Mp, twist_total)
    return Mp


def _fmt(x: float) -> str:
    if x == 0:
        return "0"
    return f"{x:+.2f}".rstrip("0").rstrip(".")


SHEARS = [-1.5, -1.25, -1.0, -0.75, -0.5,
          +0.5, +0.75, +1.0, +1.25, +1.5]
TWIST = 0.0  # fixed


def build_variants():
    out = []
    kf = 147
    for sh in SHEARS:
        name = f"kf{kf:03d}_sh{_fmt(sh)}_tw{_fmt(TWIST)}"
        out.append((name, sh, TWIST))
        kf += 1
    return out


def main():
    if not MODULE_NPY.exists():
        print(f"MISSING: {MODULE_NPY}", file=sys.stderr)
        sys.exit(2)
    M = np.load(MODULE_NPY)
    R120 = rot_z(2 * np.pi / 3)

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    index = []
    variants = build_variants()
    print(f"{'name':<32}  {'sh':>5}  {'tw':>5}  {'pts':>5}  {'bbox':>6}"
          f"  {'odd':>7}  {'even':>7}")

    for name, sh, tw in variants:
        Mp = perturb_shear_then_twist(M, sh, tw)
        ce = float(np.linalg.norm(R120 @ Mp[0] - Mp[-1]))
        P = build_inner(Mp, drop_seam_duplicates=True, closure_tol=ce + 1e-4)
        try:
            a, b = polyline_to_fourier(P, n_terms=20, n_samples=2048,
                                       center=True)
        except Exception as e:
            print(f"  {name}: Fourier error: {e}")
            index.append({"name": name, "error": str(e)})
            continue

        out_dir = OUT_BASE / name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{name}.fseries"
        descr = (f"module_shear_then_twist "
                 f"{{'slope': {sh}, 'twist_total': {tw}}}")
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.v1 wave11 ({descr})")

        amps = np.sqrt((a ** 2 + b ** 2).sum(axis=1))
        odd_sum = float(amps[::2].sum())
        even_sum = float(amps[1::2].sum())
        bbox = float(np.abs(P).max())

        index.append({
            "name": name, "slope": sh, "twist_total": tw, "descr": descr,
            "polyline_points": len(P), "bbox": bbox,
            "odd_sum": odd_sum, "even_sum": even_sum,
            "amps_top5": [float(x) for x in amps[:5]],
            "out": str(out_path),
        })
        print(f"  {name:<32}  {sh:+.2f}  {tw:+.2f}  {len(P):>5}"
              f"  {bbox:>6.3f}  {odd_sum:>7.3f}  {even_sum:>7.3f}")

    idx_path = OUT_BASE / "wave11_index.json"
    idx_path.write_text(json.dumps(index, indent=2))
    print(f"\n→ {idx_path}")
    print(f"→ {len(variants)} fseries files in erzeugte_knoten/kf{{147..156}}_*/")


if __name__ == "__main__":
    main()
