# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave19.py — Wave 19: sh+1.375 column gap fill.

Disc #197 g) + wave-16 finding: at sh+1.375 the basins lie
S_18_d1_D3_chiral (negative tw) and S_24_d961_D6_amphi (positive tw)
very close to each other, kf201-FAIL (sh+1.375, tw+1.6875) marks the
basin boundary. Wave 19 fills the gaps along this column in order to
localize the sharp basin boundary.

Present at sh+1.375 (post wave 14/15/16):
  tw = -1.0625 (kf209 → S_18_d1)
  tw = -1      (kf196 → S_18_d1, Wave 15 Foundation)
  tw = -0.9375 (kf210 → S_18_d1)
  tw = +1.6875 (kf201 → FAIL, basin boundary)
  tw = +1.75   (kf191 → S_24_d961, Wave 14 Foundation)
  tw = +1.8125 (kf202 → S_24_d961)

Wave 19 fills 8 gaps (step 0.25 coarse, densified at the kink):
  tw = -0.75, -0.5, -0.25, 0, +0.5, +1.0, +1.25, +1.5

8 knots kf234..kf241, all at sh=+1.375.

Expectation: basin boundary between S_18_d1 and S_24_d961 lies
somewhere between tw=+1.0 and tw=+1.6875 (kf201-FAIL). Wave 19
should localize exactly that — presumably with further stage-6
timeouts at the basin boundary (bifurcation points can kill stage 6,
see wave-13/14 failure theory).
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


SH_FIXED = +1.375
TWISTS = [-0.75, -0.5, -0.25, 0.0, +0.5, +1.0, +1.25, +1.5]


def build_variants():
    out = []
    kf = 234
    for tw in TWISTS:
        name = f"kf{kf:03d}_sh{_fmt(SH_FIXED)}_tw{_fmt(tw)}"
        out.append((name, SH_FIXED, tw))
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
    print(f"{'name':<32}  {'sh':>5}  {'tw':>5}  {'pts':>5}  {'bbox':>6}"
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
                      source=f"knotfamily_machine.v1 wave19 ({descr})")
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

    idx_path = OUT_BASE / "wave19_index.json"
    idx_path.write_text(json.dumps(index, indent=2))
    print(f"\n→ {idx_path}")
    print(f"→ {len(variants)} fseries-Files in erzeugte_knoten/kf{{234..241}}_*/")


if __name__ == "__main__":
    main()
