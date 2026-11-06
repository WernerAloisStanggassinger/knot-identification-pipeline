# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""make_family_wave30.py — Wave 30: kf31_spiral_overlay-Foundation-Probe.

Cycle 5, wave 2.

kf31_spiral_overlay (c=24, det=2401=7^4, sym=D3, sig=0, amphi) is a singleton
with first 7^4 pattern in the corpus (so far 7^2 in S_18_d49 and S_21_d169
seen, 7^3 in S_18_d8281 together with 13^2, 7^4 never isolated).

Operator: perturb_spiral_overlay(M, amp, freq) on inner module, then
build_inner (C_3-symmetric). 4x4 grid (amp, freq) around (0.2, 2) =
kf31 position, centre left out -> 15 new knots kf340..kf354.

  amp  in {0.10, 0.20, 0.30, 0.40}    (Step 0.10)
  freq in {1, 2, 3, 4}                (integer, Step 1)
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


def perturb_spiral_overlay(M, amp=0.2, freq=2):
    """Spiral disturbance on all interior vertices in xy. (Identical to wave3 variant)."""
    Mp = M.copy()
    n_int = len(M) - 2
    for i in range(1, len(M) - 1):
        t = (i - 1) / max(1, n_int - 1)
        ang = t * freq * 2 * np.pi
        Mp[i, 0] += amp * np.cos(ang)
        Mp[i, 1] += amp * np.sin(ang)
    return Mp


AMPS = [0.10, 0.20, 0.30, 0.40]
FREQS = [1, 2, 3, 4]
SKIP = (0.20, 2)  # = kf31


def _fmt_amp(a):
    return f"{a:.2f}".rstrip("0").rstrip(".")


def build_variants():
    out = []
    kf = 340
    for amp in AMPS:
        for freq in FREQS:
            if (amp, freq) == SKIP:
                continue
            name = f"kf{kf:03d}_spiral_amp{_fmt_amp(amp)}_freq{freq}"
            out.append((name, amp, freq))
            kf += 1
    return out


def main():
    if not MODULE_NPY.exists():
        print(f"MISSING: {MODULE_NPY}", file=sys.stderr)
        sys.exit(2)
    M_base = np.load(MODULE_NPY)

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    index = []
    variants = build_variants()
    print(f"{'name':<40}  {'amp':>5}  {'freq':>4}  {'pts':>5}")

    for name, amp, freq in variants:
        Mp = perturb_spiral_overlay(M_base, amp=amp, freq=freq)
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
        write_fseries(fpath, a, b, knot_name=name, source="make_family_wave30.py")
        print(f"  {name:<40}  {amp:>5.2f}  {freq:>4d}  {len(P):>5d}")
        index.append({"name": name, "amp": amp, "freq": freq,
                      "n_points": len(P), "fseries": str(fpath)})

    idx_path = OUT_BASE / "wave30_index.json"
    json.dump(index, open(idx_path, "w"), indent=2)
    print(f"\nWrote {idx_path} ({len(index)} variants)")


if __name__ == "__main__":
    main()
