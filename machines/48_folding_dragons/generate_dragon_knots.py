# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 48 — folding_dragons (automatic sequences as structure words).

Thue-Morse, paperfolding sequence (dragon curve!), Rudin-Shapiro: the class
BETWEEN periodic and random, with self-similarity at WORD level
(13 echo). Words run in 38 (torsion backbone) and 07 (braids).
BET: fold sequence as torsion word in 38 — the sequence is in blocks
±-balanced → ring-resonance candidates for free (pitch law of 38).
Thue-Morse is cube-free — do word laws hit the knot level?
(= ignition question). Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "38_twisting_the_backbone"))
sys.path.insert(0, str(HERE.parent / "07_closing_differently"))
sys.path.insert(0, str(HERE.parent))
from generate_backbone_knots import make_knot as backbone_make_knot  # noqa: E402
from generate_braid_closures import make_knot as braid_make_knot  # noqa: E402


# ----------------------------------------------------- the three sequences exactly

def thue_morse(n):
    """t(k) = bit-parity sequence, as ±1."""
    return [1 - 2 * (bin(k).count("1") % 2) for k in range(n)]


def paperfold(n):
    """Regular paperfolding sequence (dragon curve), k=1..n, as ±1:
    f(k) = 1 if (k / 2^v(k)) ≡ 1 mod 4, else −1."""
    out = []
    for k in range(1, n + 1):
        m = k
        while m % 2 == 0:
            m //= 2
        out.append(1 if m % 4 == 1 else -1)
    return out


def rudin_shapiro(n):
    """s(k) = (−1)^(number of '11' pairs in binary(k))."""
    out = []
    for k in range(n):
        b = bin(k)[2:]
        pairs = sum(1 for i in range(len(b) - 1) if b[i:i + 2] == "11")
        out.append(1 - 2 * (pairs % 2))
    return out


def has_cube(seq):
    """Does www (cube) exist as a factor? (Thue-Morse: NEVER — theorem)."""
    s = "".join("a" if x > 0 else "b" for x in seq)
    n = len(s)
    for L in range(1, n // 3 + 1):
        for i in range(n - 3 * L + 1):
            if s[i:i + L] == s[i + L:i + 2 * L] == s[i + 2 * L:i + 3 * L]:
                return True
    return False


# ------------------------------------------------------------ word routes

def torsion_word(seq, tau=0.85):
    return [tau * x for x in seq]


def braid_word(seq):
    """±1 sequence → braid word on 3 strands: +1→σ₁, −1→σ₂⁻¹
    (mixed generators, otherwise strand 3 stays separate — machine-26 lesson)."""
    return [(0, 1) if x > 0 else (1, -1) for x in seq]


def make_backbone(kind, length, theta_deg, out_root, tau=0.85):
    seq = {"thue_morse": thue_morse, "paperfold": paperfold,
           "rudin_shapiro": rudin_shapiro}[kind](length)
    entry, folder = backbone_make_knot(torsion_word(seq, tau), theta_deg,
                                       out_root)
    entry["machine"] = "48_folding_dragons"
    entry["sequence"] = {"kind": kind, "length": length,
                         "balance": int(sum(seq)),
                         "head": seq[:12]}
    if folder is not None:
        new = f"dr_{kind}_L{length}_th{theta_deg}"
        entry["name"] = new
    return entry, folder


def make_braid(kind, length, out_root):
    seq = {"thue_morse": thue_morse, "paperfold": paperfold,
           "rudin_shapiro": rudin_shapiro}[kind](length)
    name = f"dr_{kind}_braid_L{length}"
    entry, folder = braid_make_knot(name, braid_word(seq), 3, "trace", 1,
                                    Path(out_root))
    entry["machine"] = "48_folding_dragons"
    entry["sequence"] = {"kind": kind, "length": length,
                         "balance": int(sum(seq))}
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 48 (folding dragons)")
    # (1) recurrences exact: TM t(2k)=t(k), t(2k+1)=−t(k);
    #     fold sequence f(2k)=f(k), f(4k+1)=+1, f(4k+3)=−1;
    #     RS s(2k)=s(k), s(2k+1)=(−1)^k s(k)
    tm, pf, rs = thue_morse(4096), paperfold(4096), rudin_shapiro(4096)
    ok_tm = all(tm[2 * k] == tm[k] and tm[2 * k + 1] == -tm[k]
                for k in range(2047))
    okp1 = all(pf[4 * k + 1 - 1] == 1 for k in range(1023))
    okp3 = all(pf[4 * k + 3 - 1] == -1 for k in range(1023))
    okp2 = all(pf[2 * k - 1] == pf[k - 1] for k in range(1, 2048))
    ok_rs = all(rs[2 * k] == rs[k] and rs[2 * k + 1] == (-1) ** k * rs[k]
                for k in range(2047))
    ok1 = ok_tm and okp1 and okp3 and okp2 and ok_rs
    print(f"  (1) recurrences (TM/fold/RS, 4096 terms each): "
          f"{ok_tm}/{okp1 and okp2 and okp3}/{ok_rs}  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) Thue-Morse cube-free (theorem, measured on prefix 256)
    ok2 = not has_cube(thue_morse(256))
    print(f"  (2) TM prefix 256 cube-free (www-free)  "
          f"{'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) balance properties: TM prefix 2^k exactly balanced,
    #     fold sequence bounded-twist (|Σ|/n → small)
    bal_tm = sum(thue_morse(1024))
    bal_pf = sum(paperfold(1024))
    ok3 = bal_tm == 0 and abs(bal_pf) <= 8
    print(f"  (3) Balance: TM(1024)={bal_tm} (exactly 0), "
          f"Fold(1024)={bal_pf} (bounded)  {'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    # (4) automatic sequence as torsion word → ring resonance + embedding.
    #     Build lesson (bet half reversed, MEASURED on the 1424-point scan):
    #     fold-sequence resonances exist (L=24/τ=0.5: err_rel 0.028
    #     @θ=59°), but ALL intersect themselves — the block balance
    #     wrestles, but does not embed. Thue-Morse carries BOTH gates:
    #     L=8/τ=0.3/θ=11.5° → k=4, err_rel 0.0044, d_min_rel 0.0049.
    e4 = f4 = None
    for th in (11.5, 6.5, 3.5, 12.0, 11.0):
        e4, f4 = make_backbone("thue_morse", 8, th, HERE / "demo_smoke",
                               tau=0.3)
        if f4 is not None:
            break
    ok4 = f4 is not None
    if ok4:
        print(f"  (4) Dragon backbone (TM): L={e4['sequence']['length']}, "
              f"θ={e4['theta_deg']}°, Ring k={e4['ring']['k']}, "
              f"err_rel={e4['ring']['err_rel']:.4f}  ✓ "
              f"(fold resonances wrestle, but do not embed — measured)")
    else:
        print(f"  (4) ✗ FAIL: no resonance hit ({e4.get('rejected')})")
    assert ok4
    # (5) Thue-Morse as braid → 07 closure (component gate honest)
    e5 = f5 = None
    for L in (15, 16, 17, 18, 13, 14, 19, 20):
        e5, f5 = make_braid("thue_morse", L, HERE / "demo_smoke")
        if f5 is not None:
            break
    ok5 = f5 is not None
    print(f"  (5) TM braid L={e5['sequence']['length']}: "
          f"nyq={e5.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok5 else '✗ FAIL: ' + str(e5.get('rejected'))}")
    assert ok5
    print("[DEMO] 5/5 green — dragon machine ready to fire "
          "(cube-free words at the knot level = ignition question).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
