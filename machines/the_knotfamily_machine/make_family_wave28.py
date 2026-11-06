# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave28.py — Wave 28: D7 cluster extension + crash retry.

Cycle 4, wave 3.

Wave 25 has established S_21_d169_D7_chiral with only 2 members (kf287/kf289),
both at tw=-1.0. Hypothesis:
  - negative tw is a stable basin: grazing more finely tw ∈ {-1.5, -1.25}
    should produce further D7 cluster members.
  - tw between 0 and +1 is a crash column: graze more finely tw ∈ {-0.75, -0.25, +0.25}
    shows whether there is a narrow stable strip before the crash.
  - the kf284 unknot degeneration (sh=-0.5/tw=-1.0) stays out:
    ONLY sh ∈ {0, +0.5} (no more sh=-0.5).

Operator: identical to wave 25 — module shear_then_helix-twist, then
snap_closure and build_polygon_symmetric(n=7).

  shear in {0.0, +0.5}                              (Step 0.5)
  twist in {-1.5, -1.25, -0.75, -0.25, +0.25}       (mixed fine)

Centre (wave-25 knots) left out -> 10 new knots kf311..kf320.
"""
from __future__ import annotations
import site_config
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from knot_machine import (  # noqa: E402
    build_polygon_symmetric, polygon_generator, rot_z,
)
from fseries_io import polyline_to_fourier, write_fseries  # noqa: E402
from make_family_polygon import snap_closure  # noqa: E402

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


SHEARS = [0.0, +0.5]
TWISTS = [-1.50, -1.25, -0.75, -0.25, +0.25]
# Center is conceptually „the wave-25 sweep" (sh ∈ {-0.5, 0, +0.5}, tw ∈ {-1, 0, +1});
# all wave-28 points are outside the wave-25 grid -> no SKIP.


def _fmt(x):
    if x == 0:
        return "0"
    return f"{x:+.3f}".rstrip("0").rstrip(".")


def build_variants():
    out = []
    kf = 311
    for sh in SHEARS:
        for tw in TWISTS:
            name = f"kf{kf:03d}_polyc07_sh{_fmt(sh)}_tw{_fmt(tw)}"
            out.append((name, sh, tw))
            kf += 1
    return out


def main():
    if not MODULE_NPY.exists():
        print(f"MISSING: {MODULE_NPY}", file=sys.stderr)
        sys.exit(2)
    M_base = np.load(MODULE_NPY)
    G7 = polygon_generator(7, False)

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    index = []
    variants = build_variants()
    print(f"{'name':<40}  {'sh':>5}  {'tw':>6}  {'pts':>5}  {'bbox':>6}"
          f"  {'odd':>7}  {'even':>7}")

    for name, sh, tw in variants:
        Mp = perturb_shear_then_twist(M_base, sh, tw)
        Ms = snap_closure(Mp, G7)
        ce = float(np.linalg.norm(G7 @ Ms[0] - Ms[-1]))
        try:
            P = build_polygon_symmetric(
                Ms, n=7, mirror=False,
                drop_seam_duplicates=True,
                closure_tol=max(ce + 1e-4, 1e-4),
            )
        except Exception as e:
            print(f"  {name}: build error: {e}")
            index.append({"name": name, "error": str(e)})
            continue
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
        descr = (f"poly_c07_shear_then_twist "
                 f"{{'slope': {sh}, 'twist_total': {tw}}}")
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.v1 wave28 ({descr})")
        amps = np.sqrt((a ** 2 + b ** 2).sum(axis=1))
        odd_sum = float(amps[::2].sum())
        even_sum = float(amps[1::2].sum())
        bbox = float(np.abs(P).max())
        index.append({
            "name": name, "slope": sh, "twist_total": tw, "n": 7,
            "mirror": False, "descr": descr,
            "closure_err_after_snap": ce,
            "polyline_points": len(P), "bbox": bbox,
            "odd_sum": odd_sum, "even_sum": even_sum,
            "amps_top5": [float(x) for x in amps[:5]],
            "out": str(out_path),
        })
        print(f"  {name:<40}  {sh:>+5.2f}  {tw:>+6.3f}  {len(P):>5}"
              f"  {bbox:>6.3f}  {odd_sum:>7.3f}  {even_sum:>7.3f}")

    idx_path = OUT_BASE / "wave28_index.json"
    idx_path.write_text(json.dumps(index, indent=2))
    print(f"\n-> {idx_path}")
    print(f"-> {len(variants)} fseries files in erzeugte_knoten/kf{{311..320}}_*/")


if __name__ == "__main__":
    main()
