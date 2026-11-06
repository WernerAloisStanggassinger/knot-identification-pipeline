# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
Common foundation of the machine ideas (from machine 02 on).

Built-in corpus lessons:
  - Densify fix (memory 2026-05-31): polyline is BEFORE the Fourier analysis
    densified + uniformly resampled — never transform a thin raw curve.
  - Nyquist guard (memory resolution_nyquist): after writing, the
    reconstruction is measured against the source polyline; too large a deviation
    = under-resolution (double-6 trap) → hard error instead of silent drift.
  - Embedding gate: curves with near-self-intersection are rejected
    (no Frankenstein feed for the pipeline).
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "the_knotfamily_machine"))
from fseries_io import (resample_uniform, polyline_to_fourier,  # noqa: E402
                        write_fseries, reconstruct_from_fourier)

EMBED_MIN_REL = 0.003   # min non-neighbour distance / diameter


def densify(P: np.ndarray, max_seg_rel: float = 0.002) -> np.ndarray:
    """Subdivide segments longer than max_seg_rel·diameter linearly."""
    P = np.asarray(P, float)
    diam = float(np.linalg.norm(P.max(0) - P.min(0)))
    max_len = max_seg_rel * diam
    Pc = np.vstack([P, P[0:1]])
    out = []
    for i in range(len(P)):
        seg = Pc[i + 1] - Pc[i]
        L = float(np.linalg.norm(seg))
        k = max(1, int(np.ceil(L / max_len)))
        for j in range(k):
            out.append(Pc[i] + seg * (j / k))
    return np.asarray(out)


def min_nonadjacent_dist(P: np.ndarray, exclude: int = 8):
    """Minimal distance of non-adjacent vertices (cyclic), chunked."""
    n = len(P)
    d_min, arg = np.inf, (0, 0)
    idx = np.arange(n)
    for s in range(0, n, 256):
        block = P[s:s + 256]
        D = np.linalg.norm(block[:, None, :] - P[None, :, :], axis=2)
        ring = np.abs(idx[None, :] - idx[s:s + 256, None])
        ring = np.minimum(ring, n - ring)
        D[ring <= exclude] = np.inf
        k = np.unravel_index(np.argmin(D), D.shape)
        if D[k] < d_min:
            d_min, arg = float(D[k]), (s + int(k[0]), int(k[1]))
    return d_min, arg


def embedding_report(P: np.ndarray, n_samples: int = 4096) -> dict:
    Q = resample_uniform(P, n_samples)
    diam = float(np.linalg.norm(Q.max(0) - Q.min(0)))
    d_min, (i, j) = min_nonadjacent_dist(Q)
    rel = d_min / diam if diam > 0 else 0.0
    return {"n_samples": n_samples, "diameter": round(diam, 6),
            "d_min": round(d_min, 6), "d_min_rel": round(rel, 6),
            "closest_pair": [int(i), int(j)],
            "embedded": bool(rel >= EMBED_MIN_REL)}


def polyline_to_fseries(P: np.ndarray, path, knot_name: str, source: str,
                        n_terms: int = 64, n_samples: int = 8192,
                        nyquist_max_rel: float = 0.01) -> dict:
    """Polyline → .fseries with densify + Nyquist guard. Returns report.

    nyquist_max_rel: max. permissible deviation reconstruction↔source relative
    to the diameter. If the guard trips, n_terms is too small (double-6 trap).
    """
    Pd = densify(P)
    Q = resample_uniform(Pd, n_samples)
    Qc = Q - Q.mean(axis=0)
    a, b = polyline_to_fourier(Qc, n_terms=n_terms, n_samples=n_samples,
                               center=False)
    R = reconstruct_from_fourier(a, b, n_samples=n_samples)
    diam = float(np.linalg.norm(Qc.max(0) - Qc.min(0)))
    dev = float(np.abs(R - Qc).max()) / diam
    if dev > nyquist_max_rel:
        raise RuntimeError(
            f"{knot_name}: Nyquist guard red — reconstruction deviation "
            f"{dev:.4f} > {nyquist_max_rel} at n_terms={n_terms}. "
            f"More terms needed (under-resolution trap, cf. double-6).")
    write_fseries(path, a, b, knot_name=knot_name, source=source)
    return {"n_terms": n_terms, "n_samples": n_samples,
            "nyquist_dev_rel": round(dev, 6), "diameter": round(diam, 6)}


def count_projection_crossings(P: np.ndarray, axis: int = 1) -> int:
    """Crossings of the projection along `axis` (segment intersections, O(n²)).

    Intended for coarse polylines (≤ ~2000 corners) as geometry sanity check —
    Target: known diagram crossing number of the construction.
    """
    keep = [d for d in range(3) if d != axis]
    Q = np.asarray(P, float)[:, keep]
    n = len(Q)
    A = Q
    B = Q[(np.arange(n) + 1) % n]
    cnt = 0
    for i in range(n):
        # vectorised against all j > i+1 (neighbour segments excluded)
        j = np.arange(i + 2, n if i > 0 else n - 1)
        if len(j) == 0:
            continue
        p, r = A[i], B[i] - A[i]
        q_, s = A[j], B[j] - A[j]
        denom = r[0] * s[:, 1] - r[1] * s[:, 0]
        dq = q_ - p
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (dq[:, 0] * s[:, 1] - dq[:, 1] * s[:, 0]) / denom
            u = (dq[:, 0] * r[1] - dq[:, 1] * r[0]) / denom
        hit = (np.abs(denom) > 1e-12) & (t > 1e-9) & (t < 1 - 1e-9) \
            & (u > 1e-9) & (u < 1 - 1e-9)
        cnt += int(hit.sum())
    return cnt
