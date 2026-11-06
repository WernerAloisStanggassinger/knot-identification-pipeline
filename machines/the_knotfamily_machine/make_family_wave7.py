# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave7.py — Wave 7: Negative + Low + Extreme Helix-Twist.

Wave 6 (Discovery #wave6) had t∈{1.25, 1.5, 1.75, 2.25, 2.5, 2.75, 3.25, 3.5,
3.75} positive. Wave 7 extends in three directions:

  Negative range (mirror hypothesis):
    t ∈ {-3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5}
    → expects mirrors of the positive-twist counterparts

  Low positive (gap 0..1.0):
    t ∈ {0.25, 0.5, 0.75}
    → fills between identity (t=0) and wave 4's t=1.0 (= cluster member K_24 cluster)

  Extreme positive (beyond wave 6):
    t ∈ {4.5, 5.0}
    → tests whether t=4.0 (kf59=K3_1 trefoil collapse) persists or further
      topology appears

= 12 variants kf74-kf85.

Hypotheses:
  - negative twist yields mirror knots (theta sign flipped) of the
    positive-twist counterparts
  - low positive (t < 1.0) probably stays in the K_24 basin
  - extreme (t > 4.0) probably stays trefoil collapse or escalates into a
    new class

Output: erzeugte_knoten/kf{74..85}_<descr>/<name>.fseries + wave7_index.json.
Sweep via the separate sweep_wave7.sh.
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


def perturb_helix_twist(M: np.ndarray, twist_total: float) -> np.ndarray:
    """Identical to wave 4-6. Rotate interior vertex i in xy by angle
    twist_total*(i/L), L=len(M)-1. Negative twist = mirror rotation direction."""
    Mp = M.copy()
    L = len(M) - 1
    for i in range(1, L):
        ang = twist_total * (i / L)
        R = rot_z(ang)
        Mp[i] = R @ M[i]
    return Mp


VARIANTS = [
    # Negative range (mirror hypothesis)
    ("kf74_helix_twist_t-0.5", -0.5),
    ("kf75_helix_twist_t-1.0", -1.0),
    ("kf76_helix_twist_t-1.5", -1.5),
    ("kf77_helix_twist_t-2.0", -2.0),
    ("kf78_helix_twist_t-2.5", -2.5),
    ("kf79_helix_twist_t-3.0", -3.0),
    ("kf80_helix_twist_t-3.5", -3.5),
    # Low positive (gap 0..1.0)
    ("kf81_helix_twist_t0.25", 0.25),
    ("kf82_helix_twist_t0.5",  0.5),
    ("kf83_helix_twist_t0.75", 0.75),
    # Extreme positive
    ("kf84_helix_twist_t4.5",  4.5),
    ("kf85_helix_twist_t5.0",  5.0),
]


def main():
    if not MODULE_NPY.exists():
        print(f"MISSING: {MODULE_NPY}", file=sys.stderr)
        sys.exit(2)
    M = np.load(MODULE_NPY)

    # C_3 rotation for closure tolerance
    R120 = rot_z(2 * np.pi / 3)

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    index = []
    print(f"{'name':<28}  {'t':>5}  {'pts':>5}  {'bbox':>6}  {'odd':>7}  {'even':>7}")

    for name, twist in VARIANTS:
        Mp = perturb_helix_twist(M, twist)
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
        descr = f"module_perturb_then_inner {{'twist_total': {twist}}}"
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.v1 wave7 ({descr})")

        amps = np.sqrt((a ** 2 + b ** 2).sum(axis=1))
        odd_sum = float(amps[::2].sum())
        even_sum = float(amps[1::2].sum())
        bbox = float(np.abs(P).max())

        index.append({
            "name": name, "twist_total": twist, "descr": descr,
            "polyline_points": len(P), "bbox": bbox,
            "odd_sum": odd_sum, "even_sum": even_sum,
            "amps_top5": [float(x) for x in amps[:5]],
            "out": str(out_path),
        })
        print(f"  {name:<28}  {twist:+.2f}  {len(P):>5}  {bbox:>6.3f}  {odd_sum:>7.3f}  {even_sum:>7.3f}")

    idx_path = OUT_BASE / "wave7_index.json"
    idx_path.write_text(json.dumps(index, indent=2))
    print(f"\n→ {idx_path}")
    print(f"→ {len(VARIANTS)} fseries files in erzeugte_knoten/kf{{74..85}}_*/")


if __name__ == "__main__":
    main()
