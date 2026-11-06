# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
Generate knot family from the extracted double-3 module (12-vertex C₃).

Per variant:
  - perturb inner vertices (M[1..L-1]) — boundary M[0]/M[L] stays fixed,
    so that the closure constraint M[-1] = R_{120°}·M[0] is preserved
  - rigid build_inner → closed 3D polyline
  - Fourier(20) → .fseries → $KNOT_PIPELINE_ROOT/erzeugte_knoten/<name>/

NO pipeline run here — only after user approval.
"""
from __future__ import annotations
import site_config
import numpy as np
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from knot_machine import build_inner, rot_z, make_module_endpoint  # noqa
from fseries_io import polyline_to_fourier, write_fseries, resample_uniform  # noqa

OUT_BASE = Path((site_config.GENERATED_KNOTS_DIR))
MODULE_NPY = ROOT / "double3_module_extracted.npy"


def perturbation_zscale(M: np.ndarray, factor: float) -> np.ndarray:
    """Scale z of the inner vertices by factor; boundary unchanged."""
    Mp = M.copy()
    Mp[1:-1, 2] = M[1:-1, 2] * factor
    return Mp


def perturbation_xyscale(M: np.ndarray, factor: float) -> np.ndarray:
    """Scale xy of the inner vertices by factor; boundary unchanged."""
    Mp = M.copy()
    Mp[1:-1, 0:2] = M[1:-1, 0:2] * factor
    return Mp


def perturbation_xyrotate(M: np.ndarray, angle_deg: float) -> np.ndarray:
    """Rotate inner vertices in the xy plane by angle_deg; boundary fixed."""
    Mp = M.copy()
    R = rot_z(np.deg2rad(angle_deg))
    Mp[1:-1] = M[1:-1] @ R.T
    return Mp


def perturbation_zshear_x(M: np.ndarray, slope: float) -> np.ndarray:
    """Linear shear z' = z + slope·x for inner vertices."""
    Mp = M.copy()
    Mp[1:-1, 2] = M[1:-1, 2] + slope * M[1:-1, 0]
    return Mp


def perturbation_radial_pinch(M: np.ndarray, factor: float) -> np.ndarray:
    """Radial scale factor about the z axis — acts on xy with factor depending
    on |z| (narrower or wider at top/bottom)."""
    Mp = M.copy()
    z = M[1:-1, 2]
    # f(z) = factor if |z| max, 1 if z=0
    z_abs_norm = np.abs(z) / max(np.abs(M[:, 2]).max(), 1e-9)
    scale = 1.0 + (factor - 1.0) * z_abs_norm
    Mp[1:-1, 0:2] = M[1:-1, 0:2] * scale[:, None]
    return Mp


VARIANTS = [
    # (name, callable, kwargs, description)
    ("kf01_zsquash",   perturbation_zscale,    dict(factor=0.5),    "z-Stauchung 0.5x"),
    ("kf02_zstretch",  perturbation_zscale,    dict(factor=1.6),    "z-Streckung 1.6x"),
    ("kf03_zflat",     perturbation_zscale,    dict(factor=0.2),    "z-Stauchung 0.2x (flach)"),
    ("kf04_xyshrink",  perturbation_xyscale,   dict(factor=0.6),    "xy-Schrumpfung 0.6x"),
    ("kf05_xygrow",    perturbation_xyscale,   dict(factor=1.4),    "xy-Aufweitung 1.4x"),
    ("kf06_xyrot+30",  perturbation_xyrotate,  dict(angle_deg=30),  "xy-Rotation +30°"),
    ("kf07_xyrot-30",  perturbation_xyrotate,  dict(angle_deg=-30), "xy-Rotation -30°"),
    ("kf08_shear+",    perturbation_zshear_x,  dict(slope=0.5),     "z-shear +0.5·x"),
    ("kf09_shear-",    perturbation_zshear_x,  dict(slope=-0.5),    "z-shear -0.5·x"),
    ("kf10_pinchnarrow", perturbation_radial_pinch, dict(factor=0.5), "radial: oben/unten enger"),
    ("kf11_pinchwide", perturbation_radial_pinch, dict(factor=1.5), "radial: oben/unten weiter"),
    ("kf12_baseline",  perturbation_zscale,    dict(factor=1.0),    "unverändert (Kontrolle)"),
]


def main():
    M = np.load(MODULE_NPY)
    print(f"Loaded module: {M.shape}")
    print(f"  M[0]  = {M[0]}")
    print(f"  M[-1] = {M[-1]}")
    R120 = rot_z(2 * np.pi / 3)
    print(f"  R₁₂₀·M[0] = {R120 @ M[0]}")
    closure_err = float(np.linalg.norm(R120 @ M[0] - M[-1]))
    print(f"  closure error of the original module: {closure_err:.6f}")
    print()

    OUT_BASE.mkdir(parents=True, exist_ok=True)

    summary = []
    for name, fn, kwargs, descr in VARIANTS:
        Mp = fn(M, **kwargs)
        # check closure (should be preserved, since boundary unchanged)
        ce = float(np.linalg.norm(R120 @ Mp[0] - Mp[-1]))
        # Build closed polyline
        P = build_inner(Mp, drop_seam_duplicates=True, closure_tol=closure_err + 1e-6)
        # Fourier
        a, b = polyline_to_fourier(P, n_terms=20, n_samples=2048, center=True)
        # write .fseries to dedicated folder
        out_dir = OUT_BASE / name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{name}.fseries"
        write_fseries(out_path, a, b, knot_name=name,
                      source=f"knotfamily_machine.v1 ({descr})")
        amps = np.sqrt((a ** 2 + b ** 2).sum(axis=1))
        odd_sum = float(amps[::2].sum())   # j=1,3,5,...
        even_sum = float(amps[1::2].sum()) # j=2,4,6,...
        bbox = float(np.abs(P).max())
        # record
        summary.append({
            "name": name, "descr": descr, "closure_err": ce,
            "bbox": bbox,
            "amps_top5": amps[:5].tolist(),
            "odd_sum": odd_sum, "even_sum": even_sum,
            "out": str(out_path),
        })

    # Tabular overview
    print(f"{'name':<18}  {'bbox':>6}  {'odd':>7}  {'even':>7}  "
          f"{'|c1|':>6}  {'|c2|':>6}  {'|c3|':>6}  {'|c4|':>6}  {'|c5|':>6}  descr")
    for s in summary:
        a = s["amps_top5"]
        print(f"{s['name']:<18}  {s['bbox']:6.3f}  "
              f"{s['odd_sum']:7.4f}  {s['even_sum']:7.4f}  "
              f"{a[0]:6.3f}  {a[1]:6.3f}  {a[2]:6.3f}  {a[3]:6.3f}  {a[4]:6.3f}  "
              f"{s['descr']}")

    # save the JSON index
    import json
    with open(OUT_BASE / "family_index.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nFamily-Index: {OUT_BASE / 'family_index.json'}")
    print(f"\n{len(summary)} modules generated under {OUT_BASE}")


if __name__ == "__main__":
    main()
