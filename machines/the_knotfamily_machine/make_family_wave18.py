# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave18.py — Wave 18: Mid-Field Inter-Cluster-Border-Probe.

Disc #195 + #197 g): point 5 in the mail draft Bar-Natan/vdV is factually
wrong, because A vs B differ in polynomials+vol+theta (not
only theta). A real algebra collision (same c+det+sig+sym+Alex+Jones+
HOMFLY+vol but bit-distinct theta) would be a significant finding.

Strategy: sweep at basin borders where different clusters lie close
together. Expectation: at borders, knots with almost
identical polynomial profiles but distinct theta values may appear.

Hot spot: between the sh=0 D6 plateau (S_24_d81_D6_amphi 28 members)
and the sh+0.5 region (12n_837, S_18_d243 etc.). Plus a neg-shear side probe.

8 knots kf226..kf233:
  pos-shear focus (6 knots step 0.125):
    sh∈{+0.125, +0.25, +0.375} × tw∈{0, +1}
  neg-shear side probe (2 knots):
    sh∈{-0.125, -0.25} × tw=0

Expected findings:
  - kf226..kf231 (pos-shear): probably S_24_d81_D6 (D6 plateau
    edge) or transition to 12n_837. Real algebra collision possible
    but unlikely.
  - kf232/233 (neg-shear): probably the S_24_d81_D6 plateau edge
    (wave-11 sh=-0.5 was D6-amphi too).

If a wave-18 probe is theta-distinct from all existing
S clusters but otherwise has an identical profile → REAL algebra collision
for mail point 5.
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


VARIANTS_SPEC = [
    ("kf226", +0.125, 0.0),
    ("kf227", +0.125, +1.0),
    ("kf228", +0.25,  0.0),
    ("kf229", +0.25,  +1.0),
    ("kf230", +0.375, 0.0),
    ("kf231", +0.375, +1.0),
    ("kf232", -0.125, 0.0),
    ("kf233", -0.25,  0.0),
]


def build_variants():
    return [(f"{kf}_sh{_fmt(sh)}_tw{_fmt(tw)}", sh, tw)
            for kf, sh, tw in VARIANTS_SPEC]


def main():
    if not MODULE_NPY.exists():
        print(f"MISSING: {MODULE_NPY}", file=sys.stderr)
        sys.exit(2)
    M = np.load(MODULE_NPY)
    R120 = rot_z(2 * np.pi / 3)

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    index = []
    variants = build_variants()
    print(f"{'name':<32}  {'sh':>6}  {'tw':>5}  {'pts':>5}  {'bbox':>6}"
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
                      source=f"knotfamily_machine.v1 wave18 ({descr})")
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
        print(f"  {name:<32}  {sh:+.3f}  {tw:+.2f}  {len(P):>5}"
              f"  {bbox:>6.3f}  {odd_sum:>7.3f}  {even_sum:>7.3f}")

    idx_path = OUT_BASE / "wave18_index.json"
    idx_path.write_text(json.dumps(index, indent=2))
    print(f"\n→ {idx_path}")
    print(f"→ {len(variants)} fseries files in erzeugte_knoten/kf{{226..233}}_*/")


if __name__ == "__main__":
    main()
