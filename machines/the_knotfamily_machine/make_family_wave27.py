# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave27.py — Wave 27: kf53 + kf80 dual foundation.

Cycle 4, wave 2.

Two foundation singletons tested together, both with small
Prime-square determinants (5^2 and 3^2 for the first time as foundation probings).

kf53 (c=24, det=25=5^2, Z/3, sig=0, fib=T) — wave-5 operator polyline_post
with global_z_twist=2*pi. Sweep gtw in {1.5, 1.75, 2.25, 2.5} (step pi/4
without centre) -> kf303..kf306.

kf80 (c=32, det=9=3^2, D3, sig=0, fib=T) — wave-7 operator
module_perturb_then_inner with helix_twist=-3.5. Sweep ht in
{-4.0, -3.75, -3.25, -3.0} (step 0.25 without centre) -> kf307..kf310.

Hypotheses:
  - kf53 family: first 5^2 cluster (so far only 3^k and 7^k prime squares)
  - kf80 family: compact Z/2 or D3 cluster at very small det=9

In total 8 knots kf303..kf310.
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


def global_z_twist(P: np.ndarray, total_twist: float) -> np.ndarray:
    """polyline_post (wave 5/23 convention): rotate points about the z axis
    proportional to the z height."""
    P = P.copy()
    z = P[:, 2]
    z_min, z_max = z.min(), z.max()
    if z_max - z_min < 1e-12:
        return P
    angles = total_twist * (z - z_min) / (z_max - z_min)
    cos_a = np.cos(angles)
    sin_a = np.sin(angles)
    x_new = cos_a * P[:, 0] - sin_a * P[:, 1]
    y_new = sin_a * P[:, 0] + cos_a * P[:, 1]
    P[:, 0] = x_new
    P[:, 1] = y_new
    return P


def perturb_helix_twist(M: np.ndarray, twist_total: float) -> np.ndarray:
    """module_perturb_then_inner (wave 4-7 convention): rotate interior
    vertex i by the angle twist_total*(i/L)."""
    Mp = M.copy()
    L = len(M) - 1
    for i in range(1, L):
        ang = twist_total * (i / L)
        R = rot_z(ang)
        Mp[i] = R @ M[i]
    return Mp


# kf53 family: global_z_twist on the finished polyline, total_twist in units of pi
GTW_FACTORS = [1.50, 1.75, 2.25, 2.50]  # kf303..kf306, step pi/4 without 2.0 (kf53)

# kf80 family: helix_twist on the module, twist_total absolute (rad equivalent)
HTW_TOTALS = [-4.00, -3.75, -3.25, -3.00]  # kf307..kf310, step 0.25 without -3.5 (kf80)


def _fmt_gtw(x):
    return f"{x:.3f}".rstrip("0").rstrip(".")


def _fmt_htw(x):
    if x == 0:
        return "0"
    return f"{x:+.3f}".rstrip("0").rstrip(".")


def build_variants():
    out = []
    kf = 303
    for tf in GTW_FACTORS:
        name = f"kf{kf:03d}_gtw_{_fmt_gtw(tf)}pi"
        out.append((name, "gtw", tf))
        kf += 1
    for ht in HTW_TOTALS:
        name = f"kf{kf:03d}_helix_twist_t{_fmt_htw(ht)}"
        out.append((name, "htw", ht))
        kf += 1
    return out


def main():
    if not MODULE_NPY.exists():
        print(f"MISSING: {MODULE_NPY}", file=sys.stderr)
        sys.exit(2)
    M = np.load(MODULE_NPY)
    R120 = rot_z(2 * np.pi / 3)

    # kf53 family needs the C_3 polyline once, then polyline-post
    P0 = build_inner(M, drop_seam_duplicates=True)

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    index = []
    variants = build_variants()
    print(f"{'name':<36}  {'op':>4}  {'param':>7}  {'pts':>5}  {'bbox':>6}"
          f"  {'odd':>7}  {'even':>7}")

    for name, op, param in variants:
        if op == "gtw":
            total_twist = param * np.pi
            P = global_z_twist(P0, total_twist)
            descr = (f"polyline_post_global_z_twist "
                     f"{{'total_twist_pi': {param}}}")
            extra = {"total_twist_pi": param, "total_twist": total_twist}
        else:  # htw
            Mp = perturb_helix_twist(M, param)
            ce = float(np.linalg.norm(R120 @ Mp[0] - Mp[-1]))
            P = build_inner(Mp, drop_seam_duplicates=True,
                            closure_tol=ce + 1e-4)
            descr = f"module_perturb_then_inner {{'twist_total': {param}}}"
            extra = {"twist_total": param, "closure_err_module": ce}

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
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.v1 wave27 ({descr})")
        amps = np.sqrt((a ** 2 + b ** 2).sum(axis=1))
        odd_sum = float(amps[::2].sum())
        even_sum = float(amps[1::2].sum())
        bbox = float(np.abs(P).max())
        rec = {"name": name, "op": op, "descr": descr,
               "polyline_points": len(P), "bbox": bbox,
               "odd_sum": odd_sum, "even_sum": even_sum,
               "amps_top5": [float(x) for x in amps[:5]],
               "out": str(out_path)}
        rec.update(extra)
        index.append(rec)
        print(f"  {name:<36}  {op:>4}  {param:>+7.3f}  {len(P):>5}"
              f"  {bbox:>6.3f}  {odd_sum:>7.3f}  {even_sum:>7.3f}")

    idx_path = OUT_BASE / "wave27_index.json"
    idx_path.write_text(json.dumps(index, indent=2))
    print(f"\n-> {idx_path}")
    print(f"-> {len(variants)} fseries files in erzeugte_knoten/kf{{303..310}}_*/")


if __name__ == "__main__":
    main()
