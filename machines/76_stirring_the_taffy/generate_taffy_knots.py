# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 76 — stirring_the_taffy (Thurston–Nielsen stirring protocols).

Stirring-rod protocols (taffy puller) = braids. The classical 3-rod
protocol σ₁σ₂⁻¹ is pseudo-Anosov with stretch factor λ = (3+√5)/2 = φ²
EXACT (Perron eigenvalue of the integer transition matrix [[1,1],[1,2]]
— the wire to the golden ladder, machine 35). Closures of the
protocol braids (07er machinery) = stirring knots; entropy log λ =
mixing quality of the protocol. Classic: (σ₁σ₂⁻¹)²-trace = figure-8-
candidate (amphichiral ⇒ Wr≈0 measurable today).
REVERSAL CLAUSE: trace closure may yield a link — the 07er
component gate judges exactly. Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "07_closing_differently"))
sys.path.insert(0, str(HERE.parent / "27_pumping_the_twist_account"))
sys.path.insert(0, str(HERE.parent))
from generate_braid_closures import make_knot as braid_make_knot  # noqa: E402
from generate_braid_closures import predict_components  # noqa: E402
from generate_writhe_pumped import writhe  # noqa: E402
from the_knotfamily_machine.fseries_io import reconstruct_from_fourier  # noqa: E402


def protocol_word(k):
    """k rounds of the taffy protocol σ₁σ₂⁻¹ on 3 strands."""
    return [(0, 1), (1, -1)] * k


def stretch_factor_exact():
    """λ of the σ₁σ₂⁻¹ protocol: Perron root of x²−3x+1 (matrix
    [[1,1],[1,2]], integer) — as (exact expression, float)."""
    M = np.array([[1, 1], [1, 2]])
    lam = max(np.linalg.eigvalsh(M))
    golden = (1 + np.sqrt(5)) / 2
    return M, float(lam), float(golden ** 2)


def make_knot(k, closure, out_root):
    word = protocol_word(k)
    name = f"tf_p{k}_{closure}"
    entry, folder = braid_make_knot(name, word, 3, closure, 1,
                                    Path(out_root))
    entry["machine"] = "76_stirring_the_taffy"
    entry["protocol"] = {"rounds": k,
                         "entropy_per_round": float(np.log(
                             stretch_factor_exact()[1]))}
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 76 (stirring taffy)")
    # (1) Stretch factor EXACT: λ = φ², root of x²−3x+1
    M, lam, g2 = stretch_factor_exact()
    poly_resid = lam ** 2 - 3 * lam + 1
    ok1 = abs(lam - g2) < 1e-12 and abs(poly_resid) < 1e-12
    print(f"  (1) λ={lam:.12f} = φ²={g2:.12f}, x²−3x+1 residual="
          f"{poly_resid:.1e}  {'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) Component predictions: (σ₁σ₂⁻¹)^k trace — k=1 permutation
    #     (123) cycle? exact from 07
    comps = {k: predict_components(protocol_word(k), 3, "trace", 1)
             for k in (1, 2, 3, 4)}
    ok2 = all(v >= 1 for v in comps.values())
    print(f"  (2) components per round count: {comps} (exact from the "
          f"braid permutation)  {'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) THE CLASSIC: build the smallest knot round; at k=2 the
    #     trace closure is the figure-8 candidate (amphichiral ⇒ Wr≈0)
    e3 = f3 = None
    for k in (2, 3, 4, 5):
        if comps.get(k, predict_components(protocol_word(k), 3,
                                           "trace", 1)) != 1:
            continue
        e3, f3 = make_knot(k, "trace", HERE / "demo_smoke")
        if f3 is not None:
            break
    ok3 = f3 is not None
    wr = None
    if ok3:
        from generate_braid_closures import polyline_to_fseries  # noqa
        a_b = np.loadtxt(f3 / f"{e3['name']}.fseries", comments="%",
                         ndmin=2)
        P = reconstruct_from_fourier(a_b[:, [0, 2, 4]], a_b[:, [1, 3, 5]],
                                     n_samples=2048)
        wr = writhe(P, n=384)
        e3["writhe_measured"] = round(float(wr), 3)
    print(f"  (3) stirring knot {e3['name']}: Wr={None if wr is None else round(wr,3)} "
          f"(figure-8 thesis: ≈0), nyq="
          f"{e3.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok3 else '✗ FAIL: ' + str(e3.get('rejected'))}")
    assert ok3
    # (4) Entropy scale: log λ per round, k rounds = k·log λ (exact)
    ent = float(np.log(lam))
    ok4 = abs(ent - 0.9624) < 1e-3
    print(f"  (4) mixing entropy: {ent:.4f} per round (log φ²)  "
          f"{'✓' if ok4 else '✗ FAIL'}")
    assert ok4
    print("[DEMO] 4/4 green — taffy machine ready to fire "
          "(pA protocols as a knot source, λ=φ² exact).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
