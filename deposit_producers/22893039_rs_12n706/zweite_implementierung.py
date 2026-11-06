#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""SECOND, INDEPENDENT implementation of the enumeration of reduced singular arc diagrams.

Assignment (Werner, 22.09.2026): GEGENWIND §3 demands a second implementation for R^s_5.
The first is in `enum_bogendiagramme.py` (Euler formula on a self-built rotation system).

This one is built differently in THREE points:

  (A) PLANARITY via a FOREIGN algorithm instead of our own face counting.
      Every crossing is replaced by a quadrilateral whose four corners alternately belong to
      strand 1 and strand 2. In every planar embedding of a quadrilateral
      the four corners in exactly this cyclic order ⟹ the two strands MUST
      cross transversally. Hence:
          Gauss word realizable in S^2  <==>  gadget graph planar
      and the decision is made by `networkx.check_planarity` (LR algorithm, foreign, tested).
      The first version instead tries 2^n arrangements and counts faces.

  (B) THE REDUCEDNESS CONDITIONS ARE NOT FIXED BUT MEASURED.
      The full text of the paper (arXiv:2606.02390v1, section 3) says verbatim:
        (1) "lambda does not contain an embedded arc meeting the same crossing in two points"
        (2) "if c is a simple closed curve meeting lambda in a single point, then c bounds a
             disk containing an aglet and no crossings"
        (3) "no crossing is incident to both aglets"
      (1) and (2) each allow two readings (see below). Which one the authors meant,
      is decided by the calibration 1/2/8 from Lemma 3.1, Lemma 3.2, Proposition 3.3 — not by me.
      Therefore this script computes ALL FOUR combinations and lets the numbers speak.

  (C) Different generation of the words (recursive placement instead of permutations + canonicalization).

Usage:  python3 zweite_implementierung.py [max_c]        (default max_c = 5)
Read and compute only; writes nothing except the result file, if --json is given.
"""
import site_config
import itertools
import json
import sys

import networkx as nx


# ----------------------------------------------------------------- generate words
def woerter(n):
    """All double-occurrence words over 1..n in first-appearance normalization.

    Placed recursively: each position gets either a new symbol (the next
    unused) or the second half of an already open symbol. That generates every
    class under renaming exactly once — without permutations and without post-canonicalization.
    """
    aus = []

    def weiter(wort, naechstes, offen):
        if len(wort) == 2 * n:
            aus.append(tuple(wort))
            return
        if naechstes <= n:                      # open a new symbol
            weiter(wort + [naechstes], naechstes + 1, offen | {naechstes})
        for s in sorted(offen):                 # close an open symbol
            weiter(wort + [s], naechstes, offen - {s})

    weiter([], 1, set())
    return aus


def kanonisch(wort):
    """Normal form under reversal of the arc (swapping the two aglets)."""
    formen = []
    for w in (tuple(wort), tuple(reversed(wort))):
        abb, naechstes, neu = {}, 1, []
        for z in w:
            if z not in abb:
                abb[z] = naechstes
                naechstes += 1
            neu.append(abb[z])
        formen.append(tuple(neu))
    return min(formen)


# ------------------------------------------------------- condition (1): two readings
def schlinge_eng(w):
    """Narrow reading: the loop is EMPTY (adjacent equal symbols) = classical R1."""
    return any(w[i] == w[i + 1] for i in range(len(w) - 1))


def schlinge_weit(w):
    """Wide reading, verbatim: there is an EMBEDDED sub-arc that the same
    crossing twice. Sub-arc w[i..j] with w[i] == w[j]; embedded means:
    no symbol in the interior w[i+1..j-1] occurs twice there (otherwise the
    sub-arc intersects itself). Crossings that the sub-arc traverses only once are allowed —
    it stays embedded, other strands merely cross it.
    """
    m = len(w)
    for i in range(m):
        for j in range(i + 1, m):
            if w[i] != w[j]:
                continue
            inneres = w[i + 1:j]
            if all(inneres.count(z) == 1 for z in set(inneres)):
                return True
    return False


# ------------------------------------------------------- condition (2): two readings
def zerlegbar_eng(w):
    """Narrow reading, verbatim: a circle meets lambda in ONE point. Such a circle
    cuts the arc into two pieces, each with one aglet. Forbidden is that
    BOTH pieces carry crossings ⟹ a prefix/suffix cut with disjoint
    symbol sets and both sides non-empty.
    """
    for k in range(1, len(w)):
        vorn, hinten = set(w[:k]), set(w[k:])
        if vorn and hinten and not (vorn & hinten):
            return True
    return False


def zerlegbar_weit(w):
    """Wide reading (this is how the first implementation computes): EVERY proper contiguous
    subword that contains its symbols completely and disjointly from the rest — hence also
    a bubble IN THE MIDDLE of the arc (circle meets lambda in TWO points).
    """
    n = len(w)
    for i in range(n):
        for j in range(i + 1, n + 1):
            teil = w[i:j]
            if len(teil) in (0, n):
                continue
            sym = set(teil)
            if len(teil) != 2 * len(sym):
                continue
            if sym & (set(w[:i]) | set(w[j:])):
                continue
            return True
    return False


# ------------------------------------------------------------------- condition (3)
def kreuzung_an_beiden_aglets(w):
    """The crossing that one reaches from the aglet without an intervening crossing,
    is the first resp. last of the word. If they are the same, (3) is violated.
    (This reading is fixed by the statement of the paper that a reduced diagram has
    EXACTLY TWO aglet crossings — under any other reading there would be more.)
    """
    return w[0] == w[-1]


# --------------------------------------------------------- planarity (foreign test)
def planar(w):
    """Realizability of the Gauss word in S^2, decided by networkx.

    Edges e_0..e_2n: e_0 from aglet A0 to the crossing at position 0, e_k from position k-1
    to position k, e_2n from the last position to aglet A1.
    Every crossing c becomes the quadrilateral (c,0)-(c,1)-(c,2)-(c,3)-(c,0).
    Strand 1 (positions p1) attaches to the OPPOSITE corners 0 and 2,
    Strand 2 (positions p2) to the corners 1 and 3 ⟹ in every planar embedding
    the strands alternate around the crossing, which is exactly transversality.
    """
    n = len(w) // 2
    vork = {}
    for pos, c in enumerate(w):
        vork.setdefault(c, []).append(pos)

    G = nx.Graph()
    ecke = {}          # (position, 'rein'/'raus') -> gadget corner
    for c, (p1, p2) in vork.items():
        for i in range(4):
            G.add_edge((c, i), (c, (i + 1) % 4))
            G.add_edge((c, 'm'), (c, i))     # centre: fills the quadrilateral so that no
            #                                  foreign edge can run through the interior
        ecke[(p1, 'rein')] = (c, 0)
        ecke[(p1, 'raus')] = (c, 2)
        ecke[(p2, 'rein')] = (c, 1)
        ecke[(p2, 'raus')] = (c, 3)

    m = len(w)
    G.add_edge('A0', ecke[(0, 'rein')])
    for k in range(1, m):
        G.add_edge(ecke[(k - 1, 'raus')], ecke[(k, 'rein')])
    G.add_edge(ecke[(m - 1, 'raus')], 'A1')

    ok, _ = nx.check_planarity(G, counterexample=False)
    return ok


# ------------------------------------------------------------------------ Enumeration
LESARTEN = {
    'eng_eng':   (schlinge_eng,  zerlegbar_eng),
    'eng_weit':  (schlinge_eng,  zerlegbar_weit),
    'weit_eng':  (schlinge_weit, zerlegbar_eng),
    'weit_weit': (schlinge_weit, zerlegbar_weit),
}


def enumerieren(n, lesart):
    schlinge, zerlegbar = LESARTEN[lesart]
    treffer = set()
    for w in woerter(n):
        if schlinge(w) or kreuzung_an_beiden_aglets(w) or zerlegbar(w):
            continue
        if not planar(w):
            continue
        treffer.add(kanonisch(w))
    return treffer


def main():
    max_c = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 5
    soll = {2: 1, 3: 2, 4: 8}          # Lemma 3.1 / Lemma 3.2 / Proposition 3.3
    ergebnis = {}

    print('=' * 88)
    print('SECOND IMPLEMENTATION — four readings of reducedness, planarity via networkx')
    print('Calibration from arXiv:2606.02390v1:  c=2 -> 1 (Lem. 3.1) · c=3 -> 2 (Lem. 3.2) · c=4 -> 8 (Prop. 3.3)')
    print('=' * 88)
    kopf = 'Lesart (Schlinge/Zerlegung)'.ljust(28) + ''.join(f'c={c}'.rjust(8) for c in range(2, max_c + 1))
    print(kopf)
    print('-' * len(kopf))
    for lesart in LESARTEN:
        zeile, ergebnis[lesart] = lesart.ljust(28), {}
        for c in range(2, max_c + 1):
            menge = enumerieren(c, lesart)
            ergebnis[lesart][c] = sorted(''.join(map(str, t)) for t in menge)
            marke = ''
            if c in soll:
                marke = '' if len(menge) == soll[c] else '*'
            zeile += (f'{len(menge)}{marke}').rjust(8)
        print(zeile)
    print('-' * len(kopf))
    print('* = deviates from the calibration of the paper')
    print()

    geeicht = [l for l in LESARTEN
               if all(len(ergebnis[l][c]) == soll[c] for c in soll if c <= max_c)]
    print(f'Readings that hit the calibration 1/2/8: {geeicht if geeicht else "NONE"}')
    for l in geeicht:
        if max_c >= 5:
            print(f'  {l}: c=5 -> {len(ergebnis[l][5])} diagrams')

    if '--json' in sys.argv:
        ziel = (site_config.ROOT + '/maschinen_ideen/RS4_12N706_20260919/zweite_implementierung_ergebnis.json')
        with open(ziel, 'w') as f:
            json.dump({'max_c': max_c, 'soll': soll, 'ergebnis': ergebnis}, f, indent=1)
        print(f'\nwritten: {ziel}')


if __name__ == '__main__':
    main()
