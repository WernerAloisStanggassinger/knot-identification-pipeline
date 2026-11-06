# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
knot_machine — modular knot generation following the "How to do knots" construction rules.

Generic generator path from 2026-05-12:

  build_polygon_symmetric(M, n, mirror)  — n rigid module copies g^k(M).
    mirror=False : g = R_{2π/n}(z)         (cyclic C_n symmetry)
    mirror=True  : g = R_{2π/n}(z) ∘ σ_h   (Schoenflies S_n, n even only)

    Closure constraint on the module: M[-1] = g·M[0]. Each segment traverses the
    module in M order (even if g is orientation-reversing).

  Backwards-compatibility wrappers:
    build_outer(M) == build_polygon_symmetric(M, n=6, mirror=True)   (S_6)
    build_inner(M) == build_polygon_symmetric(M, n=3, mirror=False)  (C_3)
"""
from __future__ import annotations
import numpy as np


def rot_z(angle: float) -> np.ndarray:
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0.0],
                     [s,  c, 0.0],
                     [0., 0., 1.]])


def polygon_generator(n: int, mirror: bool) -> np.ndarray:
    """Generator g for n-fold module chaining.

    mirror=False : g = R_{2π/n}(z)         → pure Cₙ (order n)
    mirror=True  : g = R_{2π/n}(z) ∘ σ_h   → Schoenflies Sₙ (order n, n even)

    With mirror=True and n odd the order of g is 2n instead of n —
    that would be a closed chain with 2n copies, not n. We reject
    this combination; the user should then choose mirror=False+2n or odd
    n without mirror.
    """
    if n < 2:
        raise ValueError(f"n must be >= 2 (got {n})")
    if mirror and n % 2 != 0:
        raise ValueError(
            f"mirror=True requires even n (got n={n}). "
            f"With odd n, g = R_{{2π/n}}·σ_h has order 2n instead of n — "
            f"the chain does not close after n copies."
        )
    G = rot_z(2 * np.pi / n)
    if mirror:
        G[2, 2] = -1.0
    return G


def polygon_orbit(n: int, mirror: bool) -> list[np.ndarray]:
    """n Gruppenelemente {e, g, g², …, g^{n-1}}."""
    G = polygon_generator(n, mirror)
    out = [np.eye(3)]
    for _ in range(n - 1):
        out.append(G @ out[-1])
    return out


def s6_generator() -> np.ndarray:
    """Backwards-Kompat: S_6-Generator (n=6, mirror=True)."""
    return polygon_generator(6, mirror=True)


def s6_orbit() -> list[np.ndarray]:
    """Backwards compatibility: all 6 S_6 group elements."""
    return polygon_orbit(6, mirror=True)


def make_module_endpoint(M0: np.ndarray, *, generator: np.ndarray | None = None,
                         n: int | None = None, mirror: bool | None = None
                         ) -> np.ndarray:
    """Returns M[-1] that satisfies the closure M[-1] = g·M[0].

    Call variants:
      make_module_endpoint(M0, generator=G)        — explicit generator
      make_module_endpoint(M0, n=5, mirror=False)  — generic C_n/S_n
      make_module_endpoint(M0)                     — Default S_6 (legacy)
    """
    if generator is None:
        if n is not None:
            if mirror is None:
                raise ValueError("If n is passed, mirror must also be set.")
            generator = polygon_generator(n, mirror)
        else:
            generator = s6_generator()
    M0 = np.asarray(M0, dtype=float).reshape(3)
    return generator @ M0


def _apply(R: np.ndarray, t: np.ndarray, M: np.ndarray) -> np.ndarray:
    return M @ R.T + t


def _check_closure(M: np.ndarray, generator: np.ndarray, *, tol: float = 1e-6,
                   strict: bool = False) -> float:
    err = float(np.linalg.norm(generator @ M[0] - M[-1]))
    if err > tol:
        msg = (f"Closure-Constraint verletzt: |g·M[0] - M[-1]| = {err:.6g} (tol={tol}).")
        if strict:
            raise ValueError(msg)
    return err


def build_polygon_symmetric(M: np.ndarray, n: int, mirror: bool, *,
                            drop_seam_duplicates: bool = True,
                            strict_closure: bool = False,
                            closure_tol: float = 1e-6) -> np.ndarray:
    """Generic module chaining: n rigid copies g^k(M), k=0..n-1.

    mirror=False : g = R_{2π/n}(z)          → Cₙ symmetry (n unique modules)
    mirror=True  : g = R_{2π/n}(z) ∘ σ_h    → Sₙ (n/2 unique × mirror, n even)

    Examples:
        n=3, mirror=False  → C_3 (like build_inner without z-flip pairing)
        n=5, mirror=False  → C_5  (5 unique modules)
        n=6, mirror=True   → S_6  (like build_outer)
        n=8, mirror=True   → S_8  (4 unique + 4 Mirror)
        n=10, mirror=True  → S_10 (5 unique + 5 Mirror)

    Closure condition on the module: M[-1] = g·M[0].
    """
    M = np.asarray(M, dtype=float)
    if M.ndim != 2 or M.shape[1] != 3 or M.shape[0] < 2:
        raise ValueError("M must be (n,3) with n>=2")
    G = polygon_generator(n, mirror)
    _check_closure(M, G, tol=closure_tol, strict=strict_closure)

    chunks = []
    Gk = np.eye(3)
    for k in range(n):
        seg = M @ Gk.T
        if drop_seam_duplicates and chunks:
            seg = seg[1:]
        chunks.append(seg)
        Gk = G @ Gk
    return np.vstack(chunks)


def build_outer(M: np.ndarray, *, drop_seam_duplicates: bool = True,
                strict_closure: bool = False, closure_tol: float = 1e-6
                ) -> np.ndarray:
    """Backwards compatibility: S_6 chaining. Identical to
    build_polygon_symmetric(M, n=6, mirror=True)."""
    return build_polygon_symmetric(M, n=6, mirror=True,
                                   drop_seam_duplicates=drop_seam_duplicates,
                                   strict_closure=strict_closure,
                                   closure_tol=closure_tol)


def build_inner(M: np.ndarray, *, drop_seam_duplicates: bool = True,
                strict_closure: bool = False, closure_tol: float = 1e-6
                ) -> np.ndarray:
    """Backwards compatibility: C_3 chaining. Identical to
    build_polygon_symmetric(M, n=3, mirror=False)."""
    return build_polygon_symmetric(M, n=3, mirror=False,
                                   drop_seam_duplicates=drop_seam_duplicates,
                                   strict_closure=strict_closure,
                                   closure_tol=closure_tol)


def close_seam(P: np.ndarray, mode: str = "wrap") -> np.ndarray:
    """Closes the seam (P[0] vs P[-1]) of the polyline.

    mode="wrap": appends the start point again at the end (for plot/visualization).
    mode="snap": distributes the offset linearly over the whole path
                 (minimally invasive, smooths mini-discontinuity).
    """
    P = np.asarray(P, dtype=float)
    gap = P[0] - P[-1]
    if mode == "wrap":
        return np.vstack([P, P[0:1]])
    if mode == "snap":
        n = len(P)
        # correct linearly so that P[-1] += gap (= P[0])
        weights = np.linspace(0.0, 1.0, n).reshape(-1, 1)
        return P + weights * gap
    raise ValueError(f"unbekannter mode: {mode}")
