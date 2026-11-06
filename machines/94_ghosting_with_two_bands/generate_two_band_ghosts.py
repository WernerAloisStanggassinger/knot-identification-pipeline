# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 94 — ghosting_with_two_bands (polynomial ghost family).

6-strand plat: fixed entanglement frame + TWO twist boxes
(2p half twists in σ₂, 2q in σ₄). EVEN twists leave the braid
permutation unchanged ⇒ component count is CONSTANT over the whole
(p,q) family (exactly provable from the permutation).
THE IGNITION BET (in-house, blind-test compliant): do the polynomials
depend only on p+q? Measured by the OWN pipeline — the finding
rests on own numbers (Kanenobu reference = footnote after lit check).
MAIL THREAD: Dror (polynomial ghosts in series for the θ track).
REVERSAL CLAUSE: bet may fail — then the family is a
separation-sensitivity test. Doctrine: do NOT ignite.
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

# Fixed frame (found by permutation scan: first choice with a
# 1-component plat — build lesson: the hand-made frame gave a 2-link);
# twist boxes sit in σ₂ (slot A) and σ₄ (slot B)
FRAME_PRE = [(0, 1), (4, -1), (2, 1), (1, -1), (3, 1)]
FRAME_MID = [(4, 1), (1, 1), (3, 1)]
FRAME_POST = [(0, -1), (0, -1), (4, 1), (2, -1), (1, -1)]


def family_word(p, q):
    """Rahmen + 2p-Twists (σ₂) + 2q-Twists (σ₄)."""
    w = list(FRAME_PRE)
    if p:
        w.append((1, 2 * p))
    w += FRAME_MID
    if q:
        w.append((3, 2 * q))
    w += FRAME_POST
    return w


def make_knot(p, q, out_root):
    word = family_word(p, q)
    name = f"gh_p{p}_q{q}"
    entry, folder = braid_make_knot(name, word, 6, "plat", 1,
                                    Path(out_root))
    entry["machine"] = "94_ghosting_with_two_bands"
    entry["family"] = {"p": p, "q": q, "p_plus_q": p + q}
    entry["predictions"]["geister_wette"] = (
        f"Polynomials depend only on p+q={p+q}? (ignition measures in-house)")
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 94 (two-band ghosts)")
    # (1) Even twists ⇒ component count constant over the family
    comps = {}
    for p in range(0, 4):
        for q in range(0, 4):
            comps[(p, q)] = predict_components(family_word(p, q), 6,
                                               "plat", 1)
    vals = set(comps.values())
    ok1 = len(vals) == 1
    print(f"  (1) components over 16 (p,q): {sorted(vals)} — constant  "
          f"{'✓ (even twists do not change the permutation)' if ok1 else '✗ FAIL'}")
    assert ok1
    n_comp = vals.pop()
    # (2) frame yields a KNOT (1 component)
    ok2 = n_comp == 1
    print(f"  (2) frame choice: {n_comp} component(s)  "
          f"{'✓ knot family' if ok2 else '✗ FAIL — swap the frame'}")
    assert ok2
    # (3) ghost family p+q=4 built: (1,3), (2,2), (3,1)
    built = []
    for (p, q) in ((1, 3), (2, 2), (3, 1)):
        e, f = make_knot(p, q, HERE / "demo_smoke")
        if f is None:
            print(f"      ({p},{q}): ✗ {e.get('rejected')}")
        else:
            built.append((p, q, e["fseries"]["nyquist_dev_rel"]))
    ok3 = len(built) == 3
    print(f"  (3) family p+q=4: {built}  "
          f"{'✓ three ghost candidates in place' if ok3 else '✗ FAIL'}")
    assert ok3
    # (4) control point of another class: p+q=5
    e4, f4 = make_knot(2, 3, HERE / "demo_smoke")
    ok4 = f4 is not None
    print(f"  (4) control (2,3) [p+q=5]: "
          f"nyq={e4.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok4 else '✗ FAIL: ' + str(e4.get('rejected'))}")
    assert ok4
    print("[DEMO] 4/4 green — ghost machine ready to fire "
          "(ignition bet: polynomials == f(p+q)? — own numbers, no lookup).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
