# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 51 — walking_the_gyroid (minimal-surface channels).

Closed paths on the gyroid level surface
g(x,y,z) = sin x·cos y + sin y·cos z + sin z·cos x = 0 —
the underground machine (03) in its space-filling form; symmetry
inherited from the SPACE GROUP (I4₁32). Build: gradient-projected walks
(tangential step + Newton back-projection), closure via
close-return gate ON the manifold (reversal clause: closure
is its own gate, not given for free). |g| residual = permanent guard.
SCULPTURAL: gyroid channels are the most beautiful tubes in mathematics.
Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402

G_TOL = 1e-10          # |g| residual guard after projection
STEP = 0.35            # step size on the surface
CLOSE_R = 0.45         # close-return radius (in space)
MIN_STEPS = 60         # before this no return counts (else instant closure)


def gyroid(P):
    x, y, z = P[..., 0], P[..., 1], P[..., 2]
    return (np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z)
            + np.sin(z) * np.cos(x))


def grad_gyroid(P):
    x, y, z = P[..., 0], P[..., 1], P[..., 2]
    return np.stack([np.cos(x) * np.cos(y) - np.sin(z) * np.sin(x),
                     np.cos(y) * np.cos(z) - np.sin(x) * np.sin(y),
                     np.cos(z) * np.cos(x) - np.sin(y) * np.sin(z)], axis=-1)


def project(P, iters=60, tol=1e-13):
    """Newton along ∇g onto g=0 — until convergence, not a fixed budget
    (build lesson: starts near ∇g≈0 points need >8 iterations)."""
    Q = np.asarray(P, float).copy()
    for _ in range(iters):
        g = gyroid(Q)
        if np.abs(g).max() < tol:
            break
        n = grad_gyroid(Q)
        Q = Q - (g / np.maximum((n * n).sum(axis=-1), 1e-12))[..., None] * n
    return Q


def walk(seed, n_max=1500, step=STEP, persist=0.92, home_gain=0.35):
    """Persistent tangential walk; after the exploration phase one
    homing drift (in the tangent plane!) towards the start — build lesson:
    free walks on the quasi-2D surface NEVER return in 4000 steps
    (60/60 seeds measured). The close-return gate still measures
    closure strictly (distance + |g| guard), the drift is only
    steering, not a free closure."""
    rng = np.random.default_rng(seed)
    p = project(rng.uniform(0, 2 * np.pi, 3))
    n = grad_gyroid(p); n /= np.linalg.norm(n)
    d = rng.normal(size=3)
    d -= d @ n * n; d /= np.linalg.norm(d)
    path = [p]
    for it in range(n_max):
        d_new = persist * d + (1 - persist) * rng.normal(size=3)
        if it > MIN_STEPS:                       # homing phase
            home = path[0] - path[-1]
            d_new = d_new + home_gain * home / max(np.linalg.norm(home), 1e-9)
        n = grad_gyroid(path[-1]); n /= np.linalg.norm(n)
        d_new -= d_new @ n * n
        d = d_new / np.linalg.norm(d_new)
        q = project(path[-1] + step * d)
        if abs(float(gyroid(q))) > G_TOL:
            return np.asarray(path), False        # projection broke
        path.append(q)
        if it > MIN_STEPS and np.linalg.norm(q - path[0]) < CLOSE_R:
            return np.asarray(path), True
    return np.asarray(path), False


def close_on_surface(path, n_bridge=12):
    """Bridge end→start, every bridge point projected onto the surface."""
    a, b = path[-1], path[0]
    t = np.linspace(0, 1, n_bridge + 2)[1:-1]
    bridge = project(np.outer(1 - t, a) + np.outer(t, b))
    return np.vstack([path, bridge])


def make_knot(seed, out_root):
    path, closed = walk(seed)
    entry = {"name": f"gy_walk_s{seed}", "machine": "51_walking_the_gyroid",
             "seed": seed, "n_steps": len(path), "closed": bool(closed),
             "predictions": {
                 "symmetry": "Raumgruppen-Erbe (I4₁32) — Zünd-Messung",
                 "identity": "Zünd-Messung"}}
    if not closed:
        entry["rejected"] = "Close-Return-Gate (kein Rückkehrer)"
        return entry, None
    loop = close_on_surface(path)
    g_resid = float(np.abs(gyroid(loop)).max())
    entry["g_residual_max"] = g_resid
    if g_resid > G_TOL:
        entry["rejected"] = f"|g|-Wache ({g_resid:.2e})"
        return entry, None
    emb = embedding_report(loop)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / entry["name"]
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(loop, folder / f"{entry['name']}.fseries",
                              knot_name=entry["name"],
                              source="maschinen_ideen/51",
                              n_terms=max(96, len(loop) // 6),
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    (folder / f"{entry['name']}_entry.json").write_text(
        json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 51 (gyroid)")
    # (1) projection exact: 512 random points → |g| < 1e-10
    rng = np.random.default_rng(1)
    Q = project(rng.uniform(0, 2 * np.pi, (512, 3)))
    gmax = float(np.abs(gyroid(Q)).max())
    ok1 = gmax < G_TOL
    print(f"  (1) Newton projection: max|g|={gmax:.2e}  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) Space-group symmetry of the surface itself (translation inheritance):
    #     g(x+π, π−y, z+π/2)… simple witness: g is exactly 2π-periodic
    P = rng.uniform(-10, 10, (256, 3))
    per = float(np.abs(gyroid(P + 2 * np.pi) - gyroid(P)).max())
    ok2 = per < 1e-12
    print(f"  (2) lattice periodicity 2π: max|Δg|={per:.2e}  "
          f"{'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) walk stays on the surface (guard never triggered) + close-return
    #     finds returners (scan seeds, gate decides honestly)
    e3 = f3 = None
    tried = 0
    for seed in range(60):
        tried += 1
        e3, f3 = make_knot(seed, HERE / "demo_smoke")
        if f3 is not None:
            break
    ok3 = f3 is not None
    print(f"  (3) closed gyroid path after {tried} seeds: "
          f"steps={e3.get('n_steps')}, max|g|={e3.get('g_residual_max', 0):.1e}, "
          f"nyq={e3.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok3 else '✗ FAIL: ' + str(e3.get('rejected'))}")
    assert ok3
    # (4) Gate honesty: share of returners among 20 seeds (measured,
    #     no claim — the number is the machine calibration)
    n_closed = sum(int(walk(s)[1]) for s in range(100, 120))
    print(f"  (4) close-return rate: {n_closed}/20 walks return  "
          f"✓ (measured)")
    print("[DEMO] 4/4 green — gyroid machine ready to fire "
          "(space-group symmetry of the paths = ignition measurement).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
