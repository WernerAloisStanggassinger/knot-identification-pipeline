# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 56 — projecting_from_four_dimensions (S³ stereography).

Harmonic curves on S³ ⊂ R⁴, stereographically projected to R³.
EXACT ANCHOR: (cos pt, sin pt, cos qt, sin qt)/√2 → T(p,q) — the
Torus ladder (machine 10) drops out as a special case; the image of the
flat torus is a ROUND torus of revolution (bit-measurable). Richer
4D spectra = untrodden stratum; Hopf fibration as a symmetry source.
REVERSAL CLAUSE: curves near the projection pole explode → pole-distance
gate (max x₄ ≤ 1−δ) BEFORE the projection. Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402

POLE_GAP_MIN = 0.08    # min distance from x₄=1 (pole) — reversal clause
N_SAMPLES = 8192


def harmonic4d(spec, n=N_SAMPLES):
    """spec: list per coordinate (0..3) of (freq, amp, phase).
    Delivers raw 4D curve (n,4) — NOT yet on S³."""
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    X = np.zeros((n, 4))
    for d, terms in enumerate(spec):
        for (f, A, ph) in terms:
            X[:, d] += A * np.cos(f * t + ph)
    return X


def to_sphere(X):
    """Normalize radially onto S³ (points near 0 = construction unusable)."""
    r = np.linalg.norm(X, axis=1)
    if r.min() < 1e-6:
        raise ValueError("4D curve runs through the origin — change spectrum")
    return X / r[:, None]


def stereographic(S):
    """S³ → R³ from the north pole (0,0,0,1): y_i = x_i/(1−x₄)."""
    return S[:, :3] / (1.0 - S[:, 3])[:, None]


def inverse_stereographic(Y):
    """R³ → S³ (way back, for the roundtrip proof)."""
    s = (Y ** 2).sum(axis=1)
    X4 = (s - 1.0) / (s + 1.0)
    Xi = 2.0 * Y / (s + 1.0)[:, None]
    return np.hstack([Xi, X4[:, None]])


def pole_gap(S):
    return float(1.0 - S[:, 3].max())


def torus_spec(p, q):
    """Flat torus: (cos pt, sin pt, cos qt, sin qt)/√2 → T(p,q)."""
    A = 1.0 / np.sqrt(2.0)
    return [[(p, A, 0.0)], [(p, A, -np.pi / 2)],
            [(q, A, 0.0)], [(q, A, -np.pi / 2)]]


def round_torus_residual(Y):
    """Does Y lie on a torus of revolution (axis z)? Returns (R, r, resid_rel).
    (u−R)² + z² = r² with u=√(x²+y²); R,r by least squares from u,z."""
    u = np.linalg.norm(Y[:, :2], axis=1)
    z = Y[:, 2]
    # circle fit in (u,z): u²+z² = 2Ru + (r²−R²)  → linear in (R, c)
    Acol = np.column_stack([2 * u, np.ones_like(u)])
    rhs = u ** 2 + z ** 2
    (R, c), *_ = np.linalg.lstsq(Acol, rhs, rcond=None)
    r2 = c + R ** 2
    if r2 <= 0:
        return float(R), 0.0, np.inf
    r = float(np.sqrt(r2))
    resid = np.abs(np.sqrt((u - R) ** 2 + z ** 2) - r)
    return float(R), r, float(resid.max() / r)


def winding_numbers(Y, R):
    """(winding about the z axis, winding about the tube core) — exactly integer
    for torus images: unwrap phases, total rotation/2π."""
    th1 = np.unwrap(np.arctan2(Y[:, 1], Y[:, 0]))
    u = np.linalg.norm(Y[:, :2], axis=1)
    th2 = np.unwrap(np.arctan2(Y[:, 2], u - R))
    w1 = (th1[-1] - th1[0] + (th1[1] - th1[0])) / (2 * np.pi)
    w2 = (th2[-1] - th2[0] + (th2[1] - th2[0])) / (2 * np.pi)
    return float(w1), float(w2)


def make_knot(name, spec, out_root):
    X = to_sphere(harmonic4d(spec))
    gap = pole_gap(X)
    entry = {"name": name, "machine": "56_projecting_from_four_dimensions",
             "spec": [[list(map(float, term)) for term in axis]
                      for axis in spec],
             "pole_gap": round(gap, 6),
             "predictions": {"stratum": "4D-Harmonik, stereographisch",
                             "identity": "Zünd-Messung"}}
    if gap < POLE_GAP_MIN:
        entry["rejected"] = f"Pol-Gate (gap={gap:.4f} < {POLE_GAP_MIN})"
        return entry, None
    Y = stereographic(X)
    emb = embedding_report(Y)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(Y, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/56", n_terms=128,
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 56 (S³ stereography)")
    # (1) roundtrip S³ → R³ → S³ exact
    rng = np.random.default_rng(7)
    S = to_sphere(rng.normal(size=(2048, 4)))
    S = S[S[:, 3] < 0.9]
    err = np.abs(inverse_stereographic(stereographic(S)) - S).max()
    ok1 = err < 1e-12
    print(f"  (1) stereography roundtrip: max|Δ|={err:.2e}  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) T(2,3) anchor: image lies EXACTLY on a round torus of revolution
    X = to_sphere(harmonic4d(torus_spec(2, 3)))
    Y = stereographic(X)
    R, r, resid = round_torus_residual(Y)
    ok2 = resid < 1e-10
    print(f"  (2) T(2,3) anchor: round torus R={R:.6f} r={r:.6f} "
          f"resid_rel={resid:.2e}  {'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) winding numbers = (p,q) integer
    w1, w2 = winding_numbers(Y, R)
    ok3 = abs(w1 - 2) < 1e-6 and abs(abs(w2) - 3) < 1e-6
    print(f"  (3) windings: ({w1:.8f}, {w2:.8f}) ≈ (2, ±3)  "
          f"{'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    # (4) pole gate fires: a curve that comes close to the north pole
    bad = [[(1, 0.1, 0.0)], [(1, 0.1, -np.pi / 2)], [(2, 0.05, 0.0)],
           [(0, 1.0, 0.0)]]   # x₄-dominant → runs close to the pole
    e_bad, f_bad = make_knot("st_polegate_probe", bad, HERE / "demo_smoke")
    ok4 = f_bad is None and "Pol-Gate" in str(e_bad.get("rejected", ""))
    print(f"  (4) pole gate: {e_bad.get('rejected')}  "
          f"{'✓' if ok4 else '✗ FAIL'}")
    assert ok4
    # (5) Rich 4D spectrum (virgin-territory stratum) built + gates green
    rich = [[(2, 0.7, 0.0), (5, 0.35, 1.1)],
            [(2, 0.7, -np.pi / 2), (5, 0.35, 2.6)],
            [(3, 0.7, 0.0), (7, 0.25, 0.7)],
            [(3, 0.7, -np.pi / 2), (4, 0.3, 1.9)]]
    e5, f5 = make_knot("st_rich_2v5_3v7v4", rich, HERE / "demo_smoke")
    ok5 = f5 is not None
    print(f"  (5) new-territory knot: gap={e5['pole_gap']}, "
          f"embedded={e5.get('embedding', {}).get('embedded')}, "
          f"nyq={e5.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok5 else '✗ FAIL: ' + str(e5.get('rejected'))}")
    assert ok5
    print("[DEMO] 5/5 green — stereography machine ready to fire "
          "(torus ladder drops out exactly, Hopf door open).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
