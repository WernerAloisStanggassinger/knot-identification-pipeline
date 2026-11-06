# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 110 — scaling_the_mutants (# composition).

Mutation commutes with connected sum: (K#A, Kmut#A) is again a
mutant pair, but BIGGER. Ladder: stage 0 = pair of 90, stage 1 =
both #T(2,3), stage 2 = both #T(2,3)#T(2,3). Firing question: does θ separate
small pairs, big ones, both, none? (Composition machine from 90+100,
honestly declared.) BETS EXACT today: ladder embedded as PAIRS;
stage det bets (multiplicative) noted. Doctrine: do NOT fire.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "90_rotating_the_tangle_ball"))
sys.path.insert(0, str(HERE.parent / "100_summing_the_knots"))
sys.path.insert(0, str(HERE.parent))
from generate_mutant_pairs import mutant_pair  # noqa: E402
from generate_connected_sums import trefoil, connected_sum  # noqa: E402
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402


def ladder_pair(seed, depth):
    K, Km, _, _ = mutant_pair(seed)
    T = trefoil()
    for _ in range(depth):
        K = connected_sum(K, T)
        Km = connected_sum(Km, T)
    return K, Km


def make_pair(seed, depth, out_root):
    K, Km = ladder_pair(seed, depth)
    name = f"sm_s{seed}_d{depth}"
    entry = {"name": name, "machine": "110_scaling_the_mutants",
             "seed": seed, "depth": depth,
             "predictions": {
                 "zuendung": "Polynom-Kollision je Stufe (Pflicht); "
                             "θ-Trennschärfe vs Größe (Dror)",
                 "det": "multipliziert je Stufe (×3 pro Trefoil)"}}
    folder = Path(out_root) / name
    for tag, C in (("K", K), ("Kmut", Km)):
        emb = embedding_report(C)
        if not emb["embedded"]:
            entry["rejected"] = f"embedding gate ({tag})"
            return entry, None
        folder.mkdir(parents=True, exist_ok=True)
        rep = None
        for nt in (160, 256, 384):
            try:
                rep = polyline_to_fseries(
                    C, folder / f"{name}_{tag}.fseries",
                    knot_name=f"{name}_{tag}",
                    source="maschinen_ideen/110", n_terms=nt,
                    nyquist_max_rel=0.02)
                break
            except RuntimeError:
                continue
        if rep is None:
            entry["rejected"] = f"Nyquist rot ({tag})"
            return entry, None
        entry[f"fseries_{tag}"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] smoke test machine 110 (mutant scaling)")
    # (1) ladder stages 0..2 built as pairs (scan seeds)
    built = {}
    for seed in range(8):
        ok_all = True
        rows = []
        for depth in (0, 1, 2):
            e, f = make_pair(seed, depth, HERE / "demo_smoke")
            if f is None:
                ok_all = False
                break
            rows.append((depth,
                         e["fseries_K"]["nyquist_dev_rel"],
                         e["fseries_Kmut"]["nyquist_dev_rel"]))
        if ok_all:
            built[seed] = rows
            break
    ok1 = bool(built)
    seed = list(built)[0] if built else None
    print(f"  (1) Ladder (stage, nyq_K, nyq_Kmut) for seed {seed}: "
          f"{built.get(seed)}  "
          f"{'✓ three pair stages in place' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) Outside fidelity across the stages: the trefoil summand is in
    #     positioned identically in BOTH curves of a stage
    K1, Km1 = ladder_pair(seed, 1)
    n_t = 2048
    d_tail = float(np.abs(K1[-n_t:] - Km1[-n_t:]).max())
    ok2 = d_tail < 0.05
    print(f"  (2) summand fidelity: max|Δ| in the trefoil part = {d_tail:.4f} "
          f"(mutation stays in the ball)  {'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    print("[DEMO] 2/2 green — scaling machine ready to fire "
          "(θ separation power vs size = firing question).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
