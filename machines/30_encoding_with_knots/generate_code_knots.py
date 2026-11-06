# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 30 — encoding_with_knots (cyclic codes as twist words).

Cyclic [7,4] Hamming code, generator g(x)=x³+x+1 over GF(2):
Codewords = multiples of g mod x⁷−1. Bit 1→σ1, 0→σ2 (3 strands,
trace). Markov ⇒ the machine knots NECKLACE classes (flaw 2, honestly).
Hamming balls around codewords = controlled defect families (23 wire).
Doctrine: do NOT ignite before Monday.
"""
from __future__ import annotations
import argparse
import itertools
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "07_closing_differently"))
sys.path.insert(0, str(HERE.parent))
from generate_braid_closures import (predict_components,  # noqa: E402
                                     make_knot as braid_make_knot)

N, G = 7, 0b1011                 # x³+x+1


def polymul_gf2(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    return r


def codewords():
    """All 16 words of the cyclic [7,4] code (multiples of g)."""
    out = set()
    for m in range(16):
        c = polymul_gf2(G, m)
        c = (c | (c >> N)) & ((1 << N) - 1)      # mod x^7−1
        out.add(c)
    return sorted(out)


def bits(c):
    return [(c >> i) & 1 for i in range(N)]


def hamming(a, b):
    return bin(a ^ b).count("1")


def necklace_class(c):
    rots = [((c << r) | (c >> (N - r))) & ((1 << N) - 1) for r in range(N)]
    return min(rots)


def word_of(c):
    return [(0, 1) if b else (1, 1) for b in bits(c)]


TAILS = ([], [(0, 1)], [(1, 1)], [(0, 1), (1, 1)], [(1, 1), (0, 1)])


def best_closure(c):
    """Closure × closing-syllable scan. MEASURED FINDING 2026-06-06: all 15
    codewords give 2 components without tail in ALL closures —
    the code lives, as a permutation, in the wrong coset. The adaptive
    closing syllable (minimal tail → 3-cycle) is convention like 02's
    odd-CF normalization; code bits stay untouched."""
    for tail in TAILS:
        for closure, q in (("trace", 1), ("shift", 1), ("shift", 2)):
            if predict_components(word_of(c) + tail, 3, closure, q) == 1:
                return closure, q, tail
    return None, None, None


def make_knot(c, out_root, tag="cw"):
    closure, q, tail = best_closure(c)
    if closure is None:
        return {"name": f"cd_{tag}_{c:07b}", "rejected":
                "kein Verschluss+Tail gibt Knoten"}, None
    name = f"cd_{tag}_{c:07b}_{closure}{q if closure=='shift' else ''}"         + (f"_t{len(tail)}" if tail else "")
    entry, folder = braid_make_knot(name, word_of(c) + tail, 3, closure, q,
                                    Path(out_root))
    entry["closing_syllable"] = [list(x) for x in tail]
    entry["machine"] = "30_encoding_with_knots"
    entry["code"] = {"bits": f"{c:07b}", "necklace_class":
                     f"{necklace_class(c):07b}",
                     "is_codeword": c in codewords()}
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 30 (code knots)")
    cw = codewords()
    # 1) code properties: 16 words, minimum distance 3, cyclic
    dists = [hamming(a, b) for a, b in itertools.combinations(cw, 2)]
    rot1 = {((c << 1) | (c >> (N - 1))) & ((1 << N) - 1) for c in cw}
    ok1 = len(cw) == 16 and min(dists) == 3 and rot1 == set(cw)
    print(f"  (1) [7,4] Hamming: {len(cw)} words, d_min={min(dists)}, "
          f"cyclic  {'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # 2) necklace classes: knots see classes, not words (flaw 2)
    classes = {necklace_class(c) for c in cw}
    ok2 = 1 < len(classes) < 16
    print(f"  (2) {len(classes)} necklace classes from 16 codewords "
          f"(machine knots CLASSES)  {'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # 3) build: first knot codeword + Hamming-1 defect neighbour
    built = []
    for c in cw:
        if not c:
            continue
        e, f = make_knot(c, HERE / "demo_smoke")
        if f is not None:
            built.append(e["name"])
            d = c ^ 1                          # bit 0 flipped = distance 1
            e2, f2 = make_knot(d, HERE / "demo_smoke", tag="h1")
            if f2 is not None:
                built.append(e2["name"])
            break
    ok3 = len(built) >= 1
    print(f"  (3) built: {built} (codeword + Hamming-1 neighbour if "
          f"knot)  {'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    print("[DEMO] 3/3 green — code machine ready to fire (code↔knot-"
          "distance = open question in the manifest).")


def main():
    ap = argparse.ArgumentParser(description="Machine 30 — code knots")
    sub = ap.add_subparsers(dest="cmd", required=True)
    pm = sub.add_parser("make")
    pm.add_argument("bits", help="7-bit word, e.g. 1011000")
    sub.add_parser("demo")
    args = ap.parse_args()
    if args.cmd == "demo":
        cmd_demo(args)
    else:
        c = int(args.bits, 2)
        man_path = HERE / "MANIFEST.json"
        man = json.loads(man_path.read_text()) if man_path.exists() \
            else {"knots": []}
        entry, folder = make_knot(c, HERE / "out_fseries")
        man["knots"] = [k for k in man["knots"]
                        if k["name"] != entry["name"]] + [entry]
        man_path.write_text(json.dumps(man, indent=1))
        print(f"[MAKE] {entry['name']} → "
              f"{folder if folder else entry.get('rejected')}")


if __name__ == "__main__":
    main()
