# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 64 — counting_the_sticks (equilateral stick knots).

Knots from n STRAIGHT equal-length sticks — the square-tube-nearest
machine of all. Theorem anchor: stick number s(3₁)=6, n≤5 ⇒ ALWAYS unknot.
Machine: random search in the space of equilateral closed n-gons
(batch closure by projection iteration), cheap witness = writhe
(Gauss) + projection crossings, harvest of knotted specimens.
EXACT BET: n=5 sweep delivers ZERO knot suspects (theorem as
self-test); n=6 delivers trefoil suspects at a measurable rate (the
rate = new number of our own statistics). REVERSAL CLAUSE: equi-
sidedness when closing = own gate (stick-length residual),
otherwise the sticks are a claim. Doctrine: do NOT ignite.
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

EQ_TOL = 1e-12        # stick-length and closure residual


def random_equilateral(n, batch, seed):
    """Batch of closed equilateral n-gons: unit steps,
    iterative (subtract mean + renormalise) until residual < EQ_TOL."""
    rng = np.random.default_rng(seed)
    E = rng.normal(size=(batch, n, 3))
    E /= np.linalg.norm(E, axis=2, keepdims=True)
    for _ in range(2000):
        drift = E.mean(axis=1, keepdims=True)
        if np.abs(drift).max() < EQ_TOL:
            break
        E = E - drift
        E /= np.linalg.norm(E, axis=2, keepdims=True)
    P = np.cumsum(E, axis=1)
    return P, E


def gates_ok(E):
    """Equilaterality + closure, bit-checked (per polygon)."""
    len_res = np.abs(np.linalg.norm(E, axis=2) - 1.0).max(axis=1)
    close_res = np.linalg.norm(E.sum(axis=1), axis=1)
    return (len_res < 1e-9) & (close_res < 1e-9), len_res, close_res


def writhe_batch(P):
    """Polygon writhe via Gauss over segment midpoints (vectorized
    over the batch; n small ⇒ exact enough as suspect filter)."""
    return np.array([writhe(p, n=256) for p in P])


def min_projection_crossings(P, n_dirs=24, seed=0):
    """Minimum of the crossing numbers over random projections — ≤2 ⇒
    unknot (cheap witness)."""
    from machine_common import count_projection_crossings
    rng = np.random.default_rng(seed)
    best = np.inf
    for _ in range(n_dirs):
        R, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        Q = P @ R
        best = min(best, count_projection_crossings(Q, axis=2))
        if best <= 2:
            break
    return int(best)


def harvest(n, n_batches, batch, seed0, wr_thresh=2.5):
    """Sweep: suspects = |Wr| > threshold AND min crossings ≥ 3."""
    suspects, n_total = [], 0
    for b in range(n_batches):
        P, E = random_equilateral(n, batch, seed0 + b)
        ok, _, _ = gates_ok(E)
        P = P[ok]
        n_total += len(P)
        W = writhe_batch(P)
        for i in np.nonzero(np.abs(W) > wr_thresh)[0]:
            mc = min_projection_crossings(P[i])
            if mc >= 3:
                suspects.append((P[i], float(W[i]), mc))
    return suspects, n_total


def make_knot(P, wr, out_root, tag):
    name = f"st_{tag}"
    # Slightly round the corners (3 smoothing passes on dense polyline)
    dense = []
    m = len(P)
    for k in range(m):
        a, b = P[k], P[(k + 1) % m]
        t = np.linspace(0, 1, 80, endpoint=False)
        dense.append(np.outer(1 - t, a) + np.outer(t, b))
    C = np.vstack(dense)
    for _ in range(3):
        C = 0.5 * C + 0.25 * (np.roll(C, 1, axis=0) + np.roll(C, -1, axis=0))
    entry = {"name": name, "machine": "64_counting_the_sticks",
             "n_sticks": int(len(P)), "writhe": round(wr, 4),
             "predictions": {"satz": "s(3₁)=6 — n=6-Beute = Trefoil-"
                                     "Kandidat", "identity": "Zünd-Messung"}}
    emb = embedding_report(C)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(C, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/64", n_terms=96,
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 64 (counting sticks)")
    # (1) Gates exact: equilaterality + closure
    P, E = random_equilateral(6, 256, seed=1)
    ok, lr, cr = gates_ok(E)
    ok1 = ok.all()
    print(f"  (1) gates: max length res={lr.max():.2e}, "
          f"max closure res={cr.max():.2e}, {int(ok.sum())}/256  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) Theorem self-test: n=5 → ZERO suspects
    sus5, tot5 = harvest(5, n_batches=4, batch=512, seed0=100)
    ok2 = len(sus5) == 0
    print(f"  (2) n=5 sweep ({tot5} polygons): {len(sus5)} suspects "
          f"(theorem: 0)  {'✓' if ok2 else '✗ FAIL — THEOREM VIOLATED?!'}")
    assert ok2
    # (3) n=6 harvest: trefoil suspects + rate
    sus6, tot6 = harvest(6, n_batches=12, batch=1024, seed0=2000)
    rate = len(sus6) / tot6
    ok3 = len(sus6) >= 1
    print(f"  (3) n=6 harvest: {len(sus6)}/{tot6} suspects "
          f"(Rate {rate:.2e}; |Wr|: "
          f"{', '.join(f'{w:.2f}' for _, w, _ in sus6[:4])})  "
          f"{'✓' if ok3 else '✗ FAIL (rate too small for batch?)'}")
    assert ok3
    # (4) Prey built + gates green
    Pb, wr, mc = sus6[0]
    e4, f4 = make_knot(Pb, wr, HERE / "demo_smoke", "hex_trefoil_cand")
    ok4 = f4 is not None
    print(f"  (4) hexagon catch: Wr={e4['writhe']}, minCross={mc}, "
          f"nyq={e4.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok4 else '✗ FAIL: ' + str(e4.get('rejected'))}")
    assert ok4
    print("[DEMO] 4/4 green — stick machine ready to fire "
          "(6 sticks, 1 trefoil candidate — the square-tube machine).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
