# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave31.py — Wave 31: kf304-Foundation-Probe (gtw 1D-Sweep).

Cycle 5, wave 3.

kf304_gtw_1.75pi (c=21, det=243=3^5, sym=D3, sig=+2) is a singleton with
the first 3^5 pattern (so far 3^3=27 in S_18_d27, 3^4=81 in S_30_d81,
3^7=2187 still singleton kf89).

PRE-CHECK 2026-05-17: kf36_module_reverse (c=21, det=243, sig=-2, D3) is
NOT mirror of kf304 (theta sum != 0, vol different). Therefore kf36
NOT in W31, but a pure kf304 gtw sweep wave.

Operator: perturb_global_twist(M, theta_total) on inner module, then
build_inner (C_3-symmetric). 1D sweep gtw in 16 steps from 1.0π to
2.5π (step width 0.1π), centre left out -> 15 new knots kf355..kf369.

  gtw in {1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, [1.75 = kf304], 1.8, 1.9, 2.0,
          2.1, 2.2, 2.3, 2.4, 2.5} * pi
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


def perturb_global_twist(M, theta_total):
    """Rotate all interior vertices by theta_total*z_normalized about the z axis.
    z_normalized = (z - z_min) / (z_max - z_min) on [0, 1]."""
    Mp = M.copy()
    interior = M[1:-1]
    if len(interior) == 0:
        return Mp
    z = interior[:, 2]
    zmin, zmax = z.min(), z.max()
    if zmax == zmin:
        return Mp
    z_norm = (z - zmin) / (zmax - zmin)
    for i, t in enumerate(z_norm, start=1):
        ang = theta_total * t
        R = rot_z(ang)
        Mp[i] = R @ M[i]
    return Mp


PI = np.pi
GTWS = [round(x * 0.1, 1) for x in range(10, 26)]  # 1.0 .. 2.5 step 0.1
SKIP = 1.75  # = kf304


def _fmt_pi(g):
    return f"{g:g}pi".replace(".", "_")


def build_variants():
    out = []
    kf = 355
    for g in GTWS:
        if g == SKIP:
            continue
        name = f"kf{kf:03d}_gtw_{_fmt_pi(g)}"
        out.append((name, g))
        kf += 1
    return out


def main():
    if not MODULE_NPY.exists():
        print(f"FEHLT: {MODULE_NPY}", file=sys.stderr)
        sys.exit(2)
    M_base = np.load(MODULE_NPY)

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    index = []
    variants = build_variants()
    print(f"{'name':<40}  {'gtw':>8}  {'pts':>5}")

    for name, g in variants:
        Mp = perturb_global_twist(M_base, g * PI)
        try:
            P = build_inner(Mp, drop_seam_duplicates=True)
        except Exception as e:
            print(f"  {name}: build error: {e}")
            index.append({"name": name, "error": str(e)})
            continue
        try:
            a, b = polyline_to_fourier(P, n_terms=20, n_samples=2048, center=True)
        except Exception as e:
            print(f"  {name}: Fourier error: {e}")
            index.append({"name": name, "error": f"fourier: {e}"})
            continue
        dest_dir = OUT_BASE / name
        dest_dir.mkdir(parents=True, exist_ok=True)
        fpath = dest_dir / f"{name}.fseries"
        write_fseries(fpath, a, b, knot_name=name, source="make_family_wave31.py")
        print(f"  {name:<40}  {g:>5.2f}*pi  {len(P):>5d}")
        index.append({"name": name, "gtw_in_pi": g,
                      "n_points": len(P), "fseries": str(fpath)})

    idx_path = OUT_BASE / "wave31_index.json"
    json.dump(index, open(idx_path, "w"), indent=2)
    print(f"\nWrote {idx_path} ({len(index)} variants)")


if __name__ == "__main__":
    main()
