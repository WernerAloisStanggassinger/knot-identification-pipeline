# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
Polyline → Fourier series → .fseries (DHF/Z6).

Format matching saturn.fseries:
  header (4 lines with %),
  20 lines with 6 floats each: a_x(j) b_x(j) a_y(j) b_y(j) a_z(j) b_z(j).
  Convention: x(s) = Σ_j [a_x(j) cos(j·s) + b_x(j) sin(j·s)]   for s ∈ [0,2π).
  No DC term (polyline is centered beforehand).
"""
from __future__ import annotations
import numpy as np
from pathlib import Path

DEFAULT_TERMS = 20


def resample_uniform(P: np.ndarray, n_samples: int = 1024) -> np.ndarray:
    """Resample closed 3D polyline uniformly in arc length.

    P: (n,3). It is assumed that P[-1] != P[0] (seam open) — the
    transition P[-1] → P[0] is added as a closing segment during resampling.
    Output: (n_samples, 3) — uniform in arc length, start = P[0].
    """
    P = np.asarray(P, dtype=float)
    Pc = np.vstack([P, P[0:1]])  # close explicitly
    seg = np.diff(Pc, axis=0)
    seg_len = np.linalg.norm(seg, axis=1)
    L = np.concatenate([[0.0], np.cumsum(seg_len)])
    L_total = L[-1]
    if L_total <= 0:
        raise ValueError("polyline with length 0")
    s_target = np.linspace(0.0, L_total, n_samples, endpoint=False)
    out = np.empty((n_samples, 3))
    for d in range(3):
        out[:, d] = np.interp(s_target, L, Pc[:, d])
    return out


def polyline_to_fourier(P: np.ndarray, n_terms: int = DEFAULT_TERMS,
                        n_samples: int = 1024,
                        center: bool = True
                        ) -> tuple[np.ndarray, np.ndarray]:
    """Compute real Fourier coefficients a, b of shape (n_terms, 3).

    Returns (a, b) with:
      a[j-1, d] = a_d(j),  b[j-1, d] = b_d(j),  d∈{x,y,z}, j=1..n_terms.
    """
    Q = resample_uniform(P, n_samples)
    if center:
        Q = Q - Q.mean(axis=0)
    N = Q.shape[0]
    spec = np.fft.rfft(Q, axis=0)  # (N//2+1, 3)
    if n_terms + 1 > spec.shape[0]:
        raise ValueError(f"n_samples too small for {n_terms} terms; need n_samples >= {2*(n_terms+1)}")
    # j=1..n_terms
    a = (2.0 / N) * np.real(spec[1:n_terms + 1])
    b = (-2.0 / N) * np.imag(spec[1:n_terms + 1])
    return a, b


def write_fseries(path: str | Path, a: np.ndarray, b: np.ndarray,
                  knot_name: str = "knot",
                  source: str = "knotfamily_machine") -> Path:
    """Write .fseries in saturn format.

    a, b: (n_terms, 3) arrays. Columns = (x,y,z).
    """
    path = Path(path)
    n_terms = a.shape[0]
    if a.shape != b.shape or a.shape[1] != 3:
        raise ValueError("a and b must have the same (n_terms,3) shape")
    lines = [
        f" % {knot_name} (DHF)",
        f" % process with doknotscad.Z6",
        f"% created by {source}",
        f"% lines  a_x(j) b_x(j)  a_y(j)  b_y(j)  a_z(j)  b_z(j)",
    ]
    for j in range(n_terms):
        ax, ay, az = a[j]
        bx, by, bz = b[j]
        lines.append(
            f"   {ax: .6f}    {bx: .6f}    {ay: .6f}   {by: .6f}    {az: .6f}   {bz: .6f}"
        )
    path.write_text("\n".join(lines) + "\n")
    return path


def reconstruct_from_fourier(a: np.ndarray, b: np.ndarray,
                             n_samples: int = 1024) -> np.ndarray:
    """Reconstruct 3D curve from Fourier coefficients (for verification)."""
    n_terms = a.shape[0]
    s = np.linspace(0.0, 2 * np.pi, n_samples, endpoint=False)
    out = np.zeros((n_samples, 3))
    for j in range(1, n_terms + 1):
        out += np.outer(np.cos(j * s), a[j - 1]) + np.outer(np.sin(j * s), b[j - 1])
    return out
