# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 19 — knotting_square_roots, STAGE 1 (CF → Lorenz LR word).

√d → periodic continued fraction (exact, ints) → hyperbolic matrix
A = ∏ [[a_i,1],[1,0]] over the (possibly doubled) period →
BIT-VERIFIED decomposition A == R^{a1}·L^{a2}··· in the generators
R=[[1,1],[0,1]], L=[[1,0],[1,1]] → cyclic LR word of the Lorenz/
modular knot (convention Birman arXiv:1201.0214 p.24 after Ghys:
L=UV, R=UV⁻¹; a global L↔R swap = mirror knot, registered as a
convention bet).

BETS: 19_quadratwurzel_maschine.txt (W0 bit checks · W1 Birman-
Williams fibred+positive · W2 first non-torus d · W3 c vs
period length · W4 „root of 13"). Stage 2 = template landing gear.
"""
from __future__ import annotations
import argparse, json, math, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

DS = (2, 3, 5, 6, 7, 10, 11, 13, 14, 19, 22, 31, 43)


def cf_sqrt(d: int):
    """Periodic continued fraction of √d: (a0, [period]) — exact in ints
    (standard recursion m,q; the period ends at a = 2·a0)."""
    a0 = math.isqrt(d)
    if a0 * a0 == d:
        raise ValueError(f"{d} is a square")
    m, q, a = 0, 1, a0
    period = []
    while True:
        m = q * a - m
        q = (d - m * m) // q
        a = (a0 + m) // q
        period.append(a)
        if a == 2 * a0:
            return a0, period


def cf_selftest(d: int, a0: int, period: list[int]) -> bool:
    """Pell self-test (exact): the convergent p_{l−1}/q_{l−1} (a0 followed
    from the first l−1 period terms, l = period length) satisfies
    p² − d·q² = (−1)^l. Lesson from the first flight: off-by-one — the fundamental
    solution sits at the convergent l−1, not l."""
    l = len(period)
    p0, p1 = 1, a0
    q0, q1 = 0, 1
    for a in period[:l - 1]:
        p0, p1 = p1, a * p1 + p0
        q0, q1 = q1, a * q1 + q0
    return p1 * p1 - d * q1 * q1 == (-1) ** l


def mat_mul(A, B):
    return ((A[0][0] * B[0][0] + A[0][1] * B[1][0],
             A[0][0] * B[0][1] + A[0][1] * B[1][1]),
            (A[1][0] * B[0][0] + A[1][1] * B[1][0],
             A[1][0] * B[0][1] + A[1][1] * B[1][1]))


R = ((1, 1), (0, 1))
L = ((1, 0), (1, 1))
I2 = ((1, 0), (0, 1))


def hyperbolic_matrix(period: list[int]):
    """A = ∏ [[a_i,1],[1,0]] over the period (even length forced:
    odd period is doubled — geodesic word needs
    R-block/L-block pairs)."""
    per = period if len(period) % 2 == 0 else period * 2
    A = I2
    for a in per:
        A = mat_mul(A, ((a, 1), (1, 0)))
    return A, per


def rl_word_matrix(per: list[int]):
    """W = R^{a1}·L^{a2}·R^{a3}··· as a matrix (bit-check counterpart)."""
    W = I2
    for i, a in enumerate(per):
        G = R if i % 2 == 0 else L
        for _ in range(a):
            W = mat_mul(W, G)
    return W


def lorenz_entry(d: int):
    a0, period = cf_sqrt(d)
    pell_ok = cf_selftest(d, a0, period)
    A, per = hyperbolic_matrix(period)
    W = rl_word_matrix(per)
    bit_ok = (A == W)
    tr = A[0][0] + A[1][1]
    word = "".join(("R" if i % 2 == 0 else "L") * a
                   for i, a in enumerate(per))
    trip = sum(1 for i in range(len(word))
               if word[i] == "L" and word[(i + 1) % len(word)] == "R")
    return {"d": d, "cf": f"[{a0}; {period}]",
            "period_len": len(period), "word_period": per,
            "pell_selftest": pell_ok,
            "matrix_bit_check_A_eq_RLword": bit_ok,
            "trace": tr, "hyperbolic": tr > 2,
            "lorenz_word": word, "word_len": len(word),
            "trip_number": trip,
            "predictions": {
                "W1": "gefasert + positiver Zopf (Birman-Williams)",
                "W2": "Torus-Klasse falls trip<=1; sonst Kandidat",
                "W3": "c waechst mit period_len (zweiseitig)",
                "konvention": "globale L<->R-Spiegelung = Spiegel-"
                              "Knoten (Birman S.24, L=UV R=UV^-1)"}}


def cmd_demo(_args):
    print("[DEMO] Machine 19, stage 1 — CF(√d) → Lorenz LR word")
    rows, ok = [], 0
    s1 = s2 = True
    for d in DS:
        e = lorenz_entry(d)
        rows.append(e)
        s1 &= e["pell_selftest"]
        s2 &= e["matrix_bit_check_A_eq_RLword"] and e["hyperbolic"]
        print(f"  d={d:3d}: CF={e['cf']:>22s} · word={e['lorenz_word'][:26]}"
              f"{'…' if len(e['lorenz_word']) > 26 else ''}"
              f" · trip={e['trip_number']}"
              f" · Pell {'✓' if e['pell_selftest'] else '✗'}"
              f" · A==RL word {'✓' if e['matrix_bit_check_A_eq_RLword'] else '✗'}"
              f" (tr={e['trace']})")
    print(f"  (S1) Pell/convergent self-test all d: "
          f"{'✓' if s1 else '✗'}")
    ok += s1
    print(f"  (S2) matrix BIT CHECK A == R^a1·L^a2·… all d + "
          f"hyperbolic: {'✓' if s2 else '✗'}")
    ok += s2
    s3 = len(rows) == len(DS)
    (HERE / "wurzel_wett_protokoll_20260607.json").write_text(
        json.dumps(rows, indent=1))
    print(f"  (S3) bet protocol {len(rows)} roots → "
          f"wurzel_wett_protokoll_20260607.json {'✓' if s3 else '✗'}")
    ok += s3
    w13 = next(e for e in rows if e["d"] == 13)
    print(f"  ★ √13 (Werner's sentence): word={w13['lorenz_word']} "
          f"(trip={w13['trip_number']}) — waits for the stage-2 landing gear")
    print(f"[DEMO] {ok}/3 green — {'root converter stands.' if ok == 3 else 'RED.'}"
          f" Doctrine: do NOT ignite; stage 2 = template curve.")
    return 0 if ok == 3 else 1


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    a = ap.parse_args()
    if a.cmd == "demo":
        sys.exit(cmd_demo(a))


if __name__ == "__main__":
    main()
