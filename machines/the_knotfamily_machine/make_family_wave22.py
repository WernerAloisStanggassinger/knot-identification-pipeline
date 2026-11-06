# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave22.py — Wave 22: basin-switch diagonal fine grid.

Wave 17 (Disc #198) has a **basin-switch diagonal** at tw+2.5
localized between sh+1.75 and sh+1.875:
  - kf219 (sh+1.75,  tw+2.5) → S_21_d1369_Z3_chiral (post-Salvage Disc #201)
  - kf222 (sh+1.875, tw+2.5) → 12a_1019 basin

Wave 22 refines the boundary in both axes with step 0.03125 in sh
and step 0.25 in tw, 3x3 grid:
  sh ∈ {+1.78125, +1.8125, +1.84375}   (step 0.03125, between 1.75/1.875)
  tw ∈ {+2.25, +2.5, +2.75}            (step 0.25, around tw+2.5)

9 knots kf259..kf267.

Hypotheses:
  - The basin boundary is sharp (3 sh points land either entirely
    in S_21_d1369 or entirely in 12a_1019, with a jump in between).
  - OR the boundary is a true bifurcation curve (tw-dependent
    shifted) — then one sees a diagonal pattern in the 3x3 grid.
  - OR between the basins lie stage-6 failures (bifurcation
    points are stage-6 killers, see Disc #200/#201).
  - **Bonus finding possible:** a third cluster between S_21_d1369
    and 12a_1019 (e.g. new Z3 family with different det) as a bridging stone.
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


def perturb_zshear(M, slope):
    Mp = M.copy()
    Mp[1:-1, 2] = M[1:-1, 2] + slope * M[1:-1, 0]
    return Mp


def perturb_helix_twist(M, twist_total):
    Mp = M.copy()
    L = len(M) - 1
    for i in range(1, L):
        ang = twist_total * (i / L)
        R = rot_z(ang)
        Mp[i] = R @ M[i]
    return Mp


def perturb_shear_then_twist(M, slope, twist_total):
    return perturb_helix_twist(perturb_zshear(M, slope), twist_total)


def _fmt(x):
    if x == 0:
        return "0"
    return f"{x:+.3f}".rstrip("0").rstrip(".")


SHEARS = [+1.78125, +1.8125, +1.84375]
TWISTS = [+2.25, +2.5, +2.75]


def build_variants():
    out = []
    kf = 259
    for sh in SHEARS:
        for tw in TWISTS:
            name = f"kf{kf:03d}_sh{_fmt(sh)}_tw{_fmt(tw)}"
            out.append((name, sh, tw))
            kf += 1
    return out


def main():
    if not MODULE_NPY.exists():
        print(f"FEHLT: {MODULE_NPY}", file=sys.stderr)
        sys.exit(2)
    M = np.load(MODULE_NPY)
    R120 = rot_z(2 * np.pi / 3)

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    index = []
    variants = build_variants()
    print(f"{'name':<32}  {'sh':>7}  {'tw':>5}  {'pts':>5}  {'bbox':>6}"
          f"  {'odd':>7}  {'even':>7}")

    for name, sh, tw in variants:
        Mp = perturb_shear_then_twist(M, sh, tw)
        ce = float(np.linalg.norm(R120 @ Mp[0] - Mp[-1]))
        P = build_inner(Mp, drop_seam_duplicates=True, closure_tol=ce + 1e-4)
        try:
            a, b = polyline_to_fourier(P, n_terms=20, n_samples=2048, center=True)
        except Exception as e:
            print(f"  {name}: Fourier error: {e}")
            index.append({"name": name, "error": str(e)})
            continue
        out_dir = OUT_BASE / name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{name}.fseries"
        descr = f"module_shear_then_twist {{'slope': {sh}, 'twist_total': {tw}}}"
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.v1 wave22 ({descr})")
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
        print(f"  {name:<32}  {sh:+.5f}  {tw:+.2f}  {len(P):>5}"
              f"  {bbox:>6.3f}  {odd_sum:>7.3f}  {even_sum:>7.3f}")

    idx_path = OUT_BASE / "wave22_index.json"
    idx_path.write_text(json.dumps(index, indent=2))
    print(f"\n→ {idx_path}")
    print(f"→ {len(variants)} fseries-Files in erzeugte_knoten/kf{{259..267}}_*/")


if __name__ == "__main__":
    main()
