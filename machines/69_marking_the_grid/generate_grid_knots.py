# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 69 — marking_the_grid (grid diagrams / arc presentations).

n×n grid, exactly one X and one O per row/column (two disjoint
permutations): row segments join X–O horizontally (bottom),
column segments vertical (top) — FULLY enumerable.
Components = cycles of σ_O∘σ_X⁻¹ (exact; here additionally cross-proven by
independent path trace). Lit anchor: arc index,
grid number of the trefoil = 5.
BETS today: component formula == path trace (50 random grids,
bit); n=4 vs n=5 writhe sweep as a grid-number shadow (proxy,
honestly declared — writhe is geometry, no proof).
REVERSAL CLAUSE: lattice collinearity = vertex trap (jitter +
odd subdivision, 46er lesson). Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "27_pumping_the_twist_account"))
sys.path.insert(0, str(HERE.parent))
from generate_writhe_pumped import writhe  # noqa: E402
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402

Z_V, Z_H = 0.35, -0.35      # vertical TOP, horizontal BOTTOM


def random_grid(n, rng):
    """(X_col_per_row, O_col_per_row), disjoint."""
    for _ in range(500):
        X = rng.permutation(n)
        O = rng.permutation(n)
        if np.all(X != O):
            return X, O
    raise RuntimeError("no disjoint permutation found")


def components_formula(X, O):
    """Cycles of σ_O∘σ_X⁻¹ (column permutation)."""
    n = len(X)
    col_to_row_X = np.empty(n, int)
    col_to_row_X[X] = np.arange(n)
    perm = O[col_to_row_X]          # column j (X) → column of the O in that row
    seen, n_comp = set(), 0
    for s in range(n):
        if s in seen:
            continue
        n_comp += 1
        c = s
        while c not in seen:
            seen.add(c)
            c = perm[c]
    return n_comp


def components_trace(X, O):
    """Independent cross-proof: walk the marker graph explicitly."""
    n = len(X)
    col_to_row_X = {int(X[i]): i for i in range(n)}
    seen_rows, n_comp = set(), 0
    for r0 in range(n):
        if r0 in seen_rows:
            continue
        n_comp += 1
        r = r0
        while r not in seen_rows:
            seen_rows.add(r)
            col = int(O[r])           # row r: X→O horizontal
            r = col_to_row_X[col]     # column col: O→X vertical
    return n_comp


def grid_curve(X, O, seed=0):
    """3D polyline: horizontal at z<0, vertical at z>0, ramps at
    the markers; jitter against lattice collinearity."""
    n = len(X)
    rng = np.random.default_rng(seed)
    col_to_row_X = {int(X[i]): i for i in range(n)}
    pts = []
    r = 0
    for _ in range(n):
        cx, co = int(X[r]), int(O[r])
        pts += [[cx, r, Z_H], [co, r, Z_H]]          # horizontal
        r2 = col_to_row_X[co]
        pts += [[co, r, Z_V], [co, r2, Z_V]]          # vertical
        r = r2
    P = np.asarray(pts, float)
    P += 0.015 * rng.normal(size=P.shape)
    dense = []
    m = len(P)
    for k in range(m):
        a, b = P[k], P[(k + 1) % m]
        t = np.linspace(0, 1, 21, endpoint=False)
        dense.append(np.outer(1 - t, a) + np.outer(t, b))
    C = np.vstack(dense)
    for _ in range(2):
        C = 0.5 * C + 0.25 * (np.roll(C, 1, axis=0) + np.roll(C, -1, axis=0))
    return C


def wr_sweep(n, n_trials, seed0):
    """max |Wr| over random 1-component grids."""
    rng = np.random.default_rng(seed0)
    best, best_grid = 0.0, None
    for _ in range(n_trials):
        X, O = random_grid(n, rng)
        if components_formula(X, O) != 1:
            continue
        C = grid_curve(X, O, seed=int(rng.integers(1e6)))
        w = writhe(C[::3], n=256)
        if abs(w) > abs(best):
            best, best_grid = w, (X.copy(), O.copy())
    return best, best_grid


def make_knot(X, O, out_root, tag):
    C = grid_curve(X, O, seed=7)
    name = f"gd_n{len(X)}_{tag}"
    entry = {"name": name, "machine": "69_marking_the_grid",
             "X": [int(x) for x in X], "O": [int(o) for o in O],
             "n_components": components_formula(X, O),
             "writhe_proxy": round(float(writhe(C[::3], n=256)), 3),
             "predictions": {"arc_index_bound": int(len(X)),
                             "identity": "Zünd-Messung"}}
    if entry["n_components"] != 1:
        entry["rejected"] = "Link"
        return entry, None
    emb = embedding_report(C)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(C, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/69", n_terms=128,
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 69 (marking the grid)")
    # (1) Components: formula == independent trace, 50 random grids
    rng = np.random.default_rng(69)
    ok1 = True
    for _ in range(50):
        n = int(rng.integers(3, 8))
        X, O = random_grid(n, rng)
        ok1 = ok1 and components_formula(X, O) == components_trace(X, O)
    print(f"  (1) formula σ_O∘σ_X⁻¹ == path trace: 50/50 bit-equal  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) Grid-number shadow: n=4 vs n=5 writhe sweep (proxy)
    w4, _ = wr_sweep(4, 120, seed0=11)
    w5, g5 = wr_sweep(5, 120, seed0=12)
    print(f"  (2) max|Wr|: n=4 → {abs(w4):.2f}, n=5 → {abs(w5):.2f} "
          f"(theorem: the trefoil needs grid 5 — proxy, not a proof)  "
          f"✓ (measured)")
    # (3) Trefoil candidate built from n=5
    e3, f3 = make_knot(*g5, HERE / "demo_smoke", "maxwr")
    ok3 = f3 is not None
    print(f"  (3) {e3['name']}: Wr={e3['writhe_proxy']}, "
          f"nyq={e3.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok3 else '✗ FAIL: ' + str(e3.get('rejected'))}")
    assert ok3
    # (4) Explicit trefoil grid (standard 5-grid) built
    X = np.array([0, 1, 2, 3, 4]); O = np.array([2, 3, 4, 0, 1])
    nc = components_formula(X, O)
    e4, f4 = (None, None)
    if nc == 1:
        e4, f4 = make_knot(X, O, HERE / "demo_smoke", "std")
    ok4 = f4 is not None
    print(f"  (4) standard 5-grid (comps={nc}): "
          f"Wr={e4['writhe_proxy'] if e4 else '—'}, "
          f"nyq={e4.get('fseries', {}).get('nyquist_dev_rel') if e4 else '—'}"
          f"  {'✓' if ok4 else '✗ FAIL'}")
    assert ok4
    print("[DEMO] 4/4 green — grid machine ready to fire "
          "(census-capable; arc index as new track).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
