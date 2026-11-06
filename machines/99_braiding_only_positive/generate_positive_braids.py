# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 99 — braiding_only_positive (Bennequin self-test).

Positive braid closures: fibered, Seifert genus EXACTLY
g = (c − n + 1)/2 (Bennequin equality) — the machine knows the
genus BEFOREHAND. IGNITION BET = double in-house self-test:
pipeline genus == formula AND fibered == True over the series.
MAIL THREAD: Stoimenow (positivity) + referee hardening.
EXACT BETS today: component gate; c = word length; series built.
Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "07_closing_differently"))
sys.path.insert(0, str(HERE.parent))
from generate_braid_closures import (make_knot as braid_make_knot,  # noqa: E402
                                     predict_components)


def positive_word(seed, n_strands, length):
    rng = np.random.default_rng(seed)
    return [(int(rng.integers(0, n_strands - 1)), 1)
            for _ in range(length)]


def bennequin_genus(word, n_strands):
    c = sum(abs(k) for (_, k) in word)
    return (c - n_strands + 1) // 2, c


def make_knot(word, n_strands, out_root, tag):
    g, c = bennequin_genus(word, n_strands)
    name = f"pb_{tag}"
    entry, folder = braid_make_knot(name, word, n_strands, "trace", 1,
                                    Path(out_root))
    entry["machine"] = "99_braiding_only_positive"
    entry["bennequin"] = {"crossings": c, "n_strands": n_strands,
                          "genus_exact": g}
    entry["predictions"]["selbsttest"] = (
        f"pipeline genus == {g} AND fibered == True (ignition)")
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 99 (positive braids)")
    # (1) Anchor arithmetic: σ1³ (2 strands) = trefoil, g=(3−2+1)/2=1
    g, c = bennequin_genus([(0, 3)], 2)
    ok1 = g == 1 and c == 3
    print(f"  (1) anchor σ₁³: c={c}, Bennequin g={g} (trefoil genus 1)  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) T(p,q) control: (σ1σ2)^4 on 3 strands = T(3,4)? c=8,
    #     g=(8−3+1)/2=3 — torus formula (p−1)(q−1)/2 = 3 ✓ independent
    w34 = [(0, 1), (1, 1)] * 4
    g2, c2 = bennequin_genus(w34, 3)
    g_torus = (3 - 1) * (4 - 1) // 2
    ok2 = g2 == g_torus == 3
    print(f"  (2) (σ₁σ₂)⁴: Bennequin g={g2} == torus formula {g_torus}  "
          f"{'✓ two formulas, one number' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) Series: 4 positive knot words (component scan), built
    built = []
    seed = 0
    while len(built) < 4 and seed < 80:
        n_str = 3 + (seed % 2)
        word = positive_word(seed, n_str, 8 + 2 * (seed % 3))
        seed += 1
        if predict_components(word, n_str, "trace", 1) != 1:
            continue
        e, f = make_knot(word, n_str, HERE / "demo_smoke",
                         f"s{seed - 1}")
        if f is not None:
            built.append((e["name"], e["bennequin"]["genus_exact"],
                          e["fseries"]["nyquist_dev_rel"]))
    ok3 = len(built) == 4
    print(f"  (3) series (name, g_exact, nyq): {built}  "
          f"{'✓ genus-known class is in place' if ok3 else '✗ FAIL'}")
    assert ok3
    print("[DEMO] 3/3 green — positive-braid machine ready to fire "
          "(the pipeline is tested on knots whose genus we "
          "KNOW BEFOREHAND).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
