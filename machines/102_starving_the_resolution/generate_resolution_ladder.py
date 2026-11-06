# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 102 — starving_the_resolution (Nyquist ladder).

The double-6 lesson (Disc #235) as a controlled MEASUREMENT SERIES: a
dense mother knot (triple sum from machine 100), sampled at
n_terms = 8..64 — the ladder shows where the geometry tips over and from
when it is stable. Firing: identity per stage = the resolution
doctrine as a measurement curve of its own (referee/methods material).
BETS EXACT today: reconstruction deviation monotonically falling;
geometry tipping stage located (embedding/writhe per stage).
Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "100_summing_the_knots"))
sys.path.insert(0, str(HERE.parent / "27_pumping_the_twist_account"))
sys.path.insert(0, str(HERE.parent))
from generate_connected_sums import trefoil, connected_sum  # noqa: E402
from generate_writhe_pumped import writhe  # noqa: E402
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402
from the_knotfamily_machine.fseries_io import (polyline_to_fourier,  # noqa: E402
                                               reconstruct_from_fourier)

LADDER = (8, 12, 16, 24, 32, 48, 64)


def mother_curve():
    """Dense mother knot: 3₁#3₁#3₁ (machine 100)."""
    T = trefoil()
    return connected_sum(connected_sum(T, T), T)


def ladder_report(P):
    """Per stage: reconstruction deviation, embedding, writhe."""
    from the_knotfamily_machine.fseries_io import resample_uniform
    P = resample_uniform(np.asarray(P, float), 4096)
    Pc = P - P.mean(axis=0)
    diam = float(np.linalg.norm(Pc.max(0) - Pc.min(0)))
    rows = []
    for nt in LADDER:
        a, b = polyline_to_fourier(Pc, n_terms=nt, n_samples=4096,
                                   center=False)
        R = reconstruct_from_fourier(a, b, n_samples=4096)
        dev = float(np.abs(R - Pc).max() / diam)
        emb = embedding_report(R, n_samples=2048)
        wr = float(writhe(R, n=384))
        rows.append({"n_terms": nt, "deviation_rel": round(dev, 5),
                     "embedded": emb["embedded"],
                     "d_min_rel": emb["d_min_rel"],
                     "writhe": round(wr, 3)})
    return rows, diam


def make_ladder(out_root):
    P = mother_curve()
    rows, diam = ladder_report(P)
    name = "rl_triple31"
    entry = {"name": name, "machine": "102_starving_the_resolution",
             "ladder": rows,
             "predictions": {
                 "zuendung": "Identität je Stufe — wo kippt sie, ab wann "
                             "stabil? (Auflösungs-Doktrin als Messung)"}}
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    Pc = P - P.mean(axis=0)
    built = []
    for nt in LADDER:
        try:
            rep = polyline_to_fseries(
                Pc, folder / f"{name}_nt{nt}.fseries",
                knot_name=f"{name}_nt{nt}",
                source="maschinen_ideen/102", n_terms=nt,
                nyquist_max_rel=1.0)          # guard OFF: the ladder
            built.append(nt)                   # IS SUPPOSED to be under-resolved!
        except RuntimeError:
            pass
    entry["built_stufen"] = built
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 102 (resolution ladder)")
    P = mother_curve()
    rows, diam = ladder_report(P)
    # (1) Monotonicity: deviation falls with n_terms
    devs = [r["deviation_rel"] for r in rows]
    ok1 = all(devs[i] >= devs[i + 1] - 1e-6 for i in range(len(devs) - 1))
    print(f"  (1) deviation ladder {LADDER}: "
          f"{[f'{d:.4f}' for d in devs]}  "
          f"{'✓ monotonically falling' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) geometry tipping located: at the bottom embedding/writhe tip over
    wr = [r["writhe"] for r in rows]
    wr_dev = [abs(w - wr[-1]) for w in wr]
    kipp = next((LADDER[i] for i in range(len(wr_dev))
                 if wr_dev[i] < 0.3), None)
    ok2 = kipp is not None and wr_dev[0] > 0.3
    print(f"  (2) writhe per stage: {[f'{w:.2f}' for w in wr]} — "
          f"coarse up to n_terms<{kipp}, stable from {kipp}  "
          f"{'✓ tipping measured' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) Ladder written (Nyquist guard deliberately OFF — the
    #     under-resolution IS the object of measurement; declared!)
    e3, f3 = make_ladder(HERE / "demo_smoke")
    ok3 = f3 is not None and len(e3["built_stufen"]) == len(LADDER)
    print(f"  (3) {len(e3['built_stufen'])}/{len(LADDER)} stages "
          f"written (guard declared off — double-6 lesson as "
          f"series)  {'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    print("[DEMO] 3/3 green — resolution ladder ready to fire "
          "(the doctrine gets its own measurement curve).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
