# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 73 — unrolling_continued_fractions (2-bridges from constants).

Continued fraction [a₀;a₁,a₂,…] of π/e/√2 → convergent p/q → 2-bridge
knot b(p,q) as 4-plat: Conway correspondence continued fraction = tangle
word, built as braid σ₂^a σ₁^b σ₂^c … on 4 strands with plat
closure (machinery of 07). DET = p BY CONSTRUCTION — the
sharpest det control machine in the fan (11 built det=p², here p).
„The π knot": 355/113 in stainless steel.
BETS EXACT today: convergent arithmetic bit (22/7, 355/113);
p odd ⇒ knot (07 component gate as convention detector).
det=p = ignition bet. REVERSE: plat word convention (even/odd
layers) must match the Conway normal form — component prediction
reveals any convention drift. Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from fractions import Fraction
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "07_closing_differently"))
sys.path.insert(0, str(HERE.parent))
from generate_braid_closures import make_knot as braid_make_knot  # noqa: E402
from generate_braid_closures import predict_components  # noqa: E402


def cf_terms(x, n_terms):
    """Continued-fraction terms of x (float-safe for small n)."""
    out = []
    for _ in range(n_terms):
        a = int(np.floor(x))
        out.append(a)
        frac = x - a
        if frac < 1e-12:
            break
        x = 1.0 / frac
    return out


def convergent(terms):
    """[a₀;a₁,…] → p/q exact (Fraction)."""
    f = Fraction(terms[-1])
    for a in reversed(terms[:-1]):
        f = a + 1 / f
    return f.numerator, f.denominator


def normalize_odd(terms):
    """4-plat normal form needs ODD CF length (word begins and
    ends with σ₂) — build lesson: even length does NOT realize b(p,q),
    the component detector fired exactly there. [.,a]→[.,a−1,1]
    (value-equal), resp. [.,a,1]→[.,a+1]."""
    t = list(terms)
    if len(t) % 2 == 0:
        if t[-1] >= 2:
            t = t[:-1] + [t[-1] - 1, 1]
        else:
            t = t[:-2] + [t[-2] + 1]
    return t


def plat_word(terms):
    """Murasugi-4-Plat: σ₂^{b₁} σ₁^{−b₂} σ₂^{b₃} … σ₂^{b_k}
    (k odd via normalize_odd)."""
    word = []
    for i, a in enumerate(normalize_odd(terms)):
        gen = 1 if i % 2 == 0 else 0
        sgn = 1 if i % 2 == 0 else -1
        word.append((gen, sgn * a))
    return word


def make_knot(label, terms, out_root):
    p, q = convergent(terms)
    word = plat_word(terms)
    name = f"cf_{label}_{p}_{q}"
    entry, folder = braid_make_knot(name, word, 4, "plat", 1,
                                    Path(out_root))
    entry["machine"] = "73_unrolling_continued_fractions"
    entry["cf"] = {"terms": terms, "p": int(p), "q": int(q)}
    entry["predictions"]["det_wette"] = f"det = p = {p} (Zündung misst)"
    entry["predictions"]["familie"] = "2-Brücken/rational (b(p,q))"
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 73 (unrolling continued fractions)")
    # (1) convergent arithmetic exact: π → 22/7, 355/113
    t_pi = cf_terms(np.pi, 4)            # [3,7,15,1]
    ok1 = (convergent([3, 7]) == (22, 7)
           and convergent([3, 7, 15, 1]) == (355, 113)
           and t_pi == [3, 7, 15, 1])
    print(f"  (1) π convergents: {t_pi} → 22/7 ✓, 355/113 ✓  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) Parity: p odd ⇒ knot — component gate as
    #     convention detector over several convergents
    cases = {"pi22_7": [3, 7], "e_8_3": cf_terms(np.e, 3),
             "sqrt2_7_5": [1, 2, 2]}
    rows = []
    for lab, terms in cases.items():
        p, q = convergent(terms)
        n_comp = predict_components(plat_word(terms), 4, "plat", 1)
        rows.append((lab, p, q, p % 2, n_comp))
    ok2 = all((p % 2 == 1) == (nc == 1) for _, p, q, par, nc in rows)
    print(f"  (2) parity rule (p odd ⇔ 1 component): "
          f"{[(l, p, nc) for l, p, q, _, nc in rows]}  "
          f"{'✓' if ok2 else '✗ FAIL — convention drift!'}")
    assert ok2
    # (3) THE π KNOT: 22/7 built (det bet 22? p=22 even → link!
    #     → take 355/113: p=355 odd = knot)
    e3, f3 = make_knot("pi", [3, 7, 15, 1], HERE / "demo_smoke")
    ok3 = f3 is not None
    print(f"  (3) π knot b(355,113): "
          f"nyq={e3.get('fseries', {}).get('nyquist_dev_rel')}, "
          f"det bet={e3['predictions']['det_wette']}  "
          f"{'✓' if ok3 else '✗ FAIL: ' + str(e3.get('rejected'))}")
    assert ok3
    # (4) Control family: √2 convergent (choose p odd)
    e4, f4 = make_knot("sqrt2", [1, 2, 2], HERE / "demo_smoke")
    ok4 = f4 is not None
    print(f"  (4) √2 knot b({e4['cf']['p']},{e4['cf']['q']}): "
          f"nyq={e4.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok4 else '✗ FAIL: ' + str(e4.get('rejected'))}")
    assert ok4
    print("[DEMO] 4/4 green — continued-fraction machine ready to fire "
          "(det=p by construction; the π knot waits for stainless steel).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
