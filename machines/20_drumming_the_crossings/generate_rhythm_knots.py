# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 20 — drumming_the_crossings (Euclidean rhythms → knots).

Rhythm E(k,n) (k beats maximally evenly on n positions,
Bresenham form: beat at (i·k) mod n < k) is read as a braid word:
  --mode pos : beat = σ1, rest = σ2      (positive braids only)
  --mode alt : beat = σ1, rest = σ2⁻¹    (mixed signs)
Both conventions on equal footing (flip clause 1) — the convention is
OUR construction, no literature theorem (honest in the manifest).

Bets in two stages (flaw 2): word periodicity chunk^gcd EXACT (checked);
knot periodicity = ignition-time CANDIDATE. Rotation invariance = Markov
consistency gate, not a finding (flaw 3).

Doctrine: do NOT ignite before Monday.
"""
from __future__ import annotations
import argparse
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "07_closing_differently"))
sys.path.insert(0, str(HERE.parent))
from generate_braid_closures import (predict_components,  # noqa: E402
                                     make_knot as braid_make_knot)


def euclidean_rhythm(k, n):
    """Bresenham form of the Euclidean rhythm: 1 = beat, 0 = rest."""
    return [1 if (i * k) % n < k else 0 for i in range(n)]


def rhythm_word(pattern, mode="alt"):
    """Rhythm → braid word (3 strands): beat σ1, rest σ2^±1."""
    pause = (1, 1) if mode == "pos" else (1, -1)
    return [(0, 1) if b else pause for b in pattern]


def word_chunk_period(word):
    """Smallest period p with word == chunk^(len/p) (word level, exact)."""
    n = len(word)
    for p in range(1, n + 1):
        if n % p == 0 and word == word[:p] * (n // p):
            return p
    return n


def make_knot(k, n, mode, out_root):
    pattern = euclidean_rhythm(k, n)
    word = rhythm_word(pattern, mode)
    g = math.gcd(k, n)
    name = f"rh_E{k}_{n}_{mode}"
    entry, folder = braid_make_knot(name, word, 3, "trace", 1, Path(out_root))
    entry["machine"] = "20_drumming_the_crossings"
    entry["rhythm"] = {"k": k, "n": n, "pattern": pattern, "mode": mode}
    entry["word_period_exact"] = word_chunk_period(word)
    entry["predictions"].update({
        "word_periodicity": f"chunk^{g} (gcd) — Wort-Ebene EXAKT geprüft",
        "knot_period_candidate": g if g > 1 else None,
        "convention_note": "Rhythmus→Knoten ist EIGENE Konstruktion, "
                           "kein Literatur-Satz",
    })
    return entry, folder


# ----------------------------------------------------------------------- CLI

def cmd_plan(args):
    rows = []
    for n in range(4, args.max_n + 1):
        for k in range(2, n - 1):
            pattern = euclidean_rhythm(k, n)
            for mode in ("pos", "alt"):
                w = rhythm_word(pattern, mode)
                if predict_components(w, 3, "trace") == 1:
                    rows.append({"k": k, "n": n, "mode": mode,
                                 "gcd": math.gcd(k, n)})
    out = HERE / f"rhythm_plan_maxn{args.max_n}.json"
    out.write_text(json.dumps({"n_knot_rhythms": len(rows), "rows": rows},
                              indent=1))
    by_gcd = {}
    for r in rows:
        by_gcd[r["gcd"]] = by_gcd.get(r["gcd"], 0) + 1
    print(f"[PLAN] {len(rows)} knot rhythms (n≤{args.max_n}, 2 modes), "
          f"by gcd: {by_gcd} → {out.name}")


def cmd_make(args):
    man_path = HERE / "MANIFEST.json"
    man = json.loads(man_path.read_text()) if man_path.exists() else {"knots": []}
    entry, folder = make_knot(args.k, args.n, args.mode, HERE / "out_fseries")
    man["knots"] = [x for x in man["knots"] if x["name"] != entry["name"]] + [entry]
    man_path.write_text(json.dumps(man, indent=1))
    print(f"[MAKE] {entry['name']} → {folder if folder else entry.get('rejected')}")


def cmd_demo(args):
    print("[DEMO] Smoke test machine 20 (rhythms)")
    # 1) tresillo anchor: E(3,8) = [1,0,0,1,0,0,1,0]
    t = euclidean_rhythm(3, 8)
    ok1 = t == [1, 0, 0, 1, 0, 0, 1, 0] and sum(t) == 3
    print(f"  (1) E(3,8) = {t} (Tresillo)  {'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # 2) gcd word periodicity exact: E(4,8) → chunk², E(3,8) → aperiodic
    w48 = rhythm_word(euclidean_rhythm(4, 8))
    w38 = rhythm_word(euclidean_rhythm(3, 8))
    p48, p38 = word_chunk_period(w48), word_chunk_period(w38)
    ok2 = p48 == 2 and p38 == 8
    print(f"  (2) word periods: E(4,8)→{p48} (gcd 4 ⇒ chunk period 2), "
          f"E(3,8)→{p38} (aperiodic)  {'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # 3) rotation consistency (Markov gate, not a finding): components equal
    rot = w38[3:] + w38[:3]
    ok3 = predict_components(w38, 3, "trace") == predict_components(rot, 3,
                                                                    "trace")
    print(f"  (3) necklace rotation: components invariant  "
          f"{'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    # 4) full build: first knot rhythm from the mini plan
    built = None
    for (k, n, mode) in ((3, 8, "alt"), (3, 8, "pos"), (5, 8, "alt"),
                         (3, 7, "alt"), (5, 7, "pos")):
        if predict_components(rhythm_word(euclidean_rhythm(k, n), mode),
                              3, "trace") == 1:
            e, f = make_knot(k, n, mode, HERE / "demo_smoke")
            if f is not None:
                built = e["name"]
                break
    print(f"  (4) built: {built} (embedded="
          f"{e['embedding']['embedded']}, Nyquist="
          f"{e.get('fseries', {}).get('nyquist_dev_rel')})  "
          f"{'✓' if built else '✗ FAIL'}")
    assert built
    print("[DEMO] 4/4 green — rhythm machine ready to fire "
          "(firing after Monday).")


def main():
    ap = argparse.ArgumentParser(description="Machine 20 — rhythm knots")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("plan"); p1.add_argument("--max-n", type=int, default=12)
    p2 = sub.add_parser("make")
    p2.add_argument("k", type=int)
    p2.add_argument("n", type=int)
    p2.add_argument("--mode", choices=["pos", "alt"], default="alt")
    sub.add_parser("demo")
    args = ap.parse_args()
    {"plan": cmd_plan, "make": cmd_make, "demo": cmd_demo}[args.cmd](args)


if __name__ == "__main__":
    main()
