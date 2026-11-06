# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 26 — listening_to_number_noise (number sequences as twist words).

Sequences: moebius (μ: +1→σ1, −1→σ1⁻¹, 0→σ2 pause; Mertens function =
carried twist account) · primegap (gaps mod 3) · digits (digits
of a number in base 3). Convention 'squeeze': μ=0 skipped (flaw 2).
Noise trio: 20 (maximally even) / 26 (numbers) / 33 (random) — one apparatus.
Doctrine: do NOT ignite before Monday.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "07_closing_differently"))
sys.path.insert(0, str(HERE.parent))
from generate_braid_closures import (predict_components,  # noqa: E402
                                     make_knot as braid_make_knot)


def moebius_upto(n):
    mu = [1] * (n + 1)
    primes = []
    is_comp = [False] * (n + 1)
    for i in range(2, n + 1):
        if not is_comp[i]:
            primes.append(i)
            mu[i] = -1
        for p in primes:
            if i * p > n:
                break
            is_comp[i * p] = True
            if i % p == 0:
                mu[i * p] = 0
                break
            mu[i * p] = -mu[i]
    return mu[1:]


def primes_upto(n):
    sieve = [True] * (n + 1)
    sieve[0:2] = [False, False]
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = [False] * len(sieve[i * i::i])
    return [i for i, s in enumerate(sieve) if s]


def sequence(kind, length):
    if kind == "moebius":
        return moebius_upto(length)
    if kind == "moebius_squeeze":
        mu = moebius_upto(4 * length)
        return [m for m in mu if m != 0][:length]
    if kind == "primegap":
        ps = primes_upto(20 * length)
        gaps = [(b - a) for a, b in zip(ps, ps[1:])][:length]
        return [g % 3 - 1 for g in gaps]            # → {−1, 0, 1}
    raise ValueError(kind)


def seq_to_word(seq):
    """+1→σ1, −1→σ1⁻¹, 0→σ2 (3 strands)."""
    return [(0, 1) if v > 0 else (0, -1) if v < 0 else (1, 1) for v in seq]


def mertens_track(seq):
    acc, out = 0, []
    for v in seq:
        acc += v
        out.append(acc)
    return out


def strands_for(kind):
    # squeeze = pure σ1 word → strand 3 would stay separate (build finding):
    # at home on 2 strands
    return 2 if kind == "moebius_squeeze" else 3


def make_knot(kind, length, out_root):
    seq = sequence(kind, length)
    word = seq_to_word(seq)
    name = f"nn_{kind}_L{length}"
    entry, folder = braid_make_knot(name, word, strands_for(kind), "trace",
                                    1, Path(out_root))
    entry["machine"] = "26_listening_to_number_noise"
    entry["sequence"] = {"kind": kind, "length": length,
                         "values_head": seq[:16],
                         "drall_konto_final": mertens_track(seq)[-1],
                         "zero_share": round(seq.count(0) / len(seq), 3)}
    return entry, folder


def cmd_demo(args):
    print("[DEMO] smoke test machine 26 (number noise)")
    # 1) Möbius anchor: μ(1..10) textbook
    mu = moebius_upto(10)
    soll = [1, -1, -1, 0, -1, 1, -1, 0, 0, 1]
    ok1 = mu == soll
    print(f"  (1) μ(1..10) = {mu} (textbook)  {'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # 2) Mertens account: M(10) = −1 (textbook)
    ok2 = mertens_track(mu)[-1] == -1
    print(f"  (2) Mertens M(10) = {mertens_track(mu)[-1]} (expected −1)  "
          f"{'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # 3) μ=0 share documented + squeeze convention pause-free
    s = sequence("moebius", 60)
    sq = sequence("moebius_squeeze", 60)
    ok3 = 0.2 < s.count(0) / 60 < 0.55 and sq.count(0) == 0
    print(f"  (3) μ=0 share {s.count(0)/60:.2f} (flaw 2 documented), "
          f"squeeze pause-free  {'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    # 4) build: first knot case per sequence (length scan)
    built = []
    for kind in ("moebius", "moebius_squeeze", "primegap"):
        for L in range(8, 30):
            if predict_components(seq_to_word(sequence(kind, L)),
                                  strands_for(kind), "trace") == 1:
                e, f = make_knot(kind, L, HERE / "demo_smoke")
                if f is not None:
                    built.append((kind, L,
                                  e["sequence"]["drall_konto_final"]))
                break
    ok4 = len(built) == 3
    print(f"  (4) built (kind, L, account): {built}  "
          f"{'✓' if ok4 else '✗ FAIL'}")
    assert ok4
    print("[DEMO] 4/4 green — number noise ready to fire "
          "(spectra comparison = ignition program).")


def main():
    ap = argparse.ArgumentParser(description="Machine 26 — number noise")
    sub = ap.add_subparsers(dest="cmd", required=True)
    pm = sub.add_parser("make")
    pm.add_argument("kind", choices=["moebius", "moebius_squeeze", "primegap"])
    pm.add_argument("--length", type=int, default=20)
    sub.add_parser("demo")
    args = ap.parse_args()
    if args.cmd == "demo":
        cmd_demo(args)
    else:
        man_path = HERE / "MANIFEST.json"
        man = json.loads(man_path.read_text()) if man_path.exists() \
            else {"knots": []}
        entry, folder = make_knot(args.kind, args.length, HERE / "out_fseries")
        man["knots"] = [k for k in man["knots"]
                        if k["name"] != entry["name"]] + [entry]
        man_path.write_text(json.dumps(man, indent=1))
        print(f"[MAKE] {entry['name']} → "
              f"{folder if folder else entry.get('rejected')}")


if __name__ == "__main__":
    main()
