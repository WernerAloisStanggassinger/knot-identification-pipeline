#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Determinant of a knot directly from the PD code — without Sage, without Singular.

WHY OWN: spherogram computes invariants only "inside Sage", and Sage needs Singular here.
A second, independent route is mandatory anyway (triple counter-check) — and this one
is elementary linear algebra over the integers, hence exact.

METHOD (Alexander matrix from the Wirtinger presentation, evaluated at t = -1)
  Every crossing of a knot diagram gives a relation between three arcs:
      over-strand o, incoming under-strand u_ein, outgoing under-strand u_aus.
  The Alexander matrix has one row per crossing:
      Column u_ein:  -t      Column u_aus:  +1      Column o:  t - 1
  (Sign convention drops out of the absolute value of the determinant.)
  det(K) = |Delta_K(-1)| = |determinant of the matrix after deleting ONE column|, at t = -1.

PD CONVENTION (spherogram): X(a, b, c, d) lists the four edge labels counter-clockwise,
starting at the INCOMING UNDER-STRAND. So: a = u_ein, c = u_aus, and {b, d} = the over-strand.
"""
from fractions import Fraction

def determinante_aus_pd(pd, t=-1):
    kanten = sorted({e for x in pd for e in x})
    idx = {e: i for i, e in enumerate(kanten)}
    n = len(pd)
    zeilen = []
    for (a, b, c, d) in pd:
        z = [Fraction(0)] * len(kanten)
        z[idx[a]] += Fraction(-t)        # incoming under-strand
        z[idx[c]] += Fraction(1)         # outgoing under-strand
        z[idx[b]] += Fraction(t - 1)     # over-strand (b and d are the same arc,
        # in a knot diagram in the PD code the over-strand is given by b,d;
        # it is counted ONCE: if b != d, both belong to the same arc only
        # if no further crossing lies between them — therefore corrected below.
        zeilen.append(z)
    return zeilen, kanten

def det_ganzzahlig(m):
    """Exact determinant of a Fraction matrix via Gauss."""
    m = [row[:] for row in m]
    n = len(m)
    if n == 0:
        return Fraction(1)
    det = Fraction(1)
    for i in range(n):
        p = next((r for r in range(i, n) if m[r][i] != 0), None)
        if p is None:
            return Fraction(0)
        if p != i:
            m[i], m[p] = m[p], m[i]; det = -det
        det *= m[i][i]
        inv = Fraction(1) / m[i][i]
        for r in range(i + 1, n):
            f = m[r][i] * inv
            if f:
                for c in range(i, n):
                    m[r][c] -= f * m[i][c]
    return det

def determinante(pd, t=-1):
    """|Delta(t)| via the Alexander matrix; over-strand correctly as ONE arc."""
    # Determine arcs: in the PD code an arc is the chain of edges between two
    # underpasses. We work edge-wise and delete one column at the end —
    # that is permissible because the edges of an over-strand are linked by the relations
    # are identified.
    kanten = sorted({e for x in pd for e in x})
    idx = {e: i for i, e in enumerate(kanten)}
    zeilen = []
    for (a, b, c, d) in pd:
        z = [Fraction(0)] * len(kanten)
        z[idx[a]] += Fraction(-t)
        z[idx[c]] += Fraction(1)
        z[idx[b]] += Fraction(t) - 1
        z[idx[d]] += Fraction(0)      # d belongs to the same over-strand as b
        zeilen.append(z)
    # Over-strand equality b == d as additional rows: edge b and edge d are
    # the same arc, so x_b - x_d = 0
    for (a, b, c, d) in pd:
        if b != d:
            z = [Fraction(0)] * len(kanten)
            z[idx[b]] += Fraction(1); z[idx[d]] -= Fraction(1)
            zeilen.append(z)
    # make square: delete one column, then choose rows so that it fits
    import itertools
    sp = 0
    gestrichen = [[z[i] for i in range(len(kanten)) if i != sp] for z in zeilen]
    m, k = len(gestrichen), len(kanten) - 1
    if m < k:
        return None
    # Rank-preserving square selection: Gauss on the full matrix, product of the pivots
    # -> |det| of the minor results as greatest common value; we take the
    # Smith-like evaluation over the first k independent rows.
    aus, rang = [], 0
    probe = []
    for z in gestrichen:
        probe.append(z)
        if rang_von(probe) > rang:
            rang += 1; aus.append(z)
        else:
            probe.pop()
        if rang == k:
            break
    if rang < k:
        return Fraction(0)
    return abs(det_ganzzahlig(aus))

def rang_von(m):
    m = [r[:] for r in m]
    zeilen, spalten = len(m), len(m[0]) if m else 0
    r = 0
    for c in range(spalten):
        p = next((i for i in range(r, zeilen) if m[i][c] != 0), None)
        if p is None: continue
        m[r], m[p] = m[p], m[r]
        inv = Fraction(1) / m[r][c]
        for i in range(zeilen):
            if i != r and m[i][c] != 0:
                f = m[i][c] * inv
                for j in range(c, spalten):
                    m[i][j] -= f * m[r][j]
        r += 1
        if r == zeilen: break
    return r

# ============================================================================
# SELFTEST — frozen, so that the calibration is not lost.
# The PD codes come from spherogram BRAID WORDS (Link(braid_closure=...)),
# not from a knot table: pure construction, blind-test clean.
# ============================================================================
EICHFAELLE = [
    # (name, braid, PD code, target determinant)
    ("Kleeblatt  = Zopf (1,1,1)",      [(1,5,2,4),(5,3,0,2),(3,1,4,0)], 3),
    ("Achterknoten = Zopf (1,-2,1,-2)", None, 5),   # PD inserted below
]

def selbsttest():
    import json, os
    print("=" * 76)
    print("SELF-TEST — determinant from PD (own route, without Sage/Singular)")
    print("=" * 76)
    faelle = [
        ("Kleeblatt   Zopf (1,1,1)",
         [(1,5,2,4), (5,3,0,2), (3,1,4,0)], 3),
        ("5_1         Zopf (1,1,1,1,1)",
         [(1,7,2,6), (7,3,8,2), (3,9,4,8), (9,5,0,4), (5,1,6,0)], 5),
        ("4_1         Zopf (1,-2,1,-2)",
         [(3,7,4,6), (1,4,2,5), (7,3,0,2), (5,0,6,1)], 5),
        ("det-7-Fall  Zopf (1,1,1,2,-1,2)",
         [(8,2,9,1), (2,10,3,9), (10,4,11,3), (11,7,0,6), (7,4,8,5), (5,1,6,0)], 7),
        ("7_1         Zopf (1,1,1,1,1,1,1)",
         [(1,9,2,8), (9,3,10,2), (3,11,4,10), (11,5,12,4), (5,13,6,12), (13,7,0,6), (7,1,8,0)], 7),
    ]
    ok = True
    for name, pd, soll in faelle:
        d = determinante(pd)
        gut = (d == soll); ok &= gut
        print(f"  [{'OK ' if gut else 'RED'}] {name:<30} det = {d}  (target {soll})")
    # COUNTER-CASE: a crossing change on the trefoil must give the unknot.
    # (The first attempt at a counter-case — setting b := d — DID NOT BITE and was
    #  discarded: it did not change the matrix. A checker that does not notice a falsification
    #  does not notice, measures nothing.)
    pd = [list(x) for x in [(1,5,2,4),(5,3,0,2),(3,1,4,0)]]
    a, b, c, e = pd[0]
    pd[0] = [b, c, e, a]            # rotation by one = under/over-strand swapped
    d = determinante([tuple(x) for x in pd])
    gut = (d == 1); ok &= gut
    print(f"  [{'OK ' if gut else 'RED'}] counter-case: crossing change on the trefoil "
          f"-> det = {d}  (target 1, unknot)")
    print()
    print("All PD codes above are MEASURED (spherogram braid closure), not typed —")
    print("a 5_1 code inserted by hand was wrong on the first attempt and was caught.")
    print()
    print("CALIBRATION:", "ALL GREEN" if ok else "🔴 RED")
    return ok

if __name__ == "__main__":
    import sys
    sys.exit(0 if selbsttest() else 1)
