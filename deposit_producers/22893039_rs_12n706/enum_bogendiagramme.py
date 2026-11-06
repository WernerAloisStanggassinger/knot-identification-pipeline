#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Enumerate reduced singular arc diagrams — calibrated on Akash et al. 2606.02390.

GOAL: decide r_s(12n706) = 4. For that a LABELED KNOTOID DIAGRAM with four is needed
crossings whose symmetric union is 12n706 (Lemma 2.3/2.4 there, det(K)=det(cl(k))^2).

STEP 1: enumerate the underlying SINGULAR ARC DIAGRAMS ourselves — an
immersed arc in S^2 with n double points, without crossing information, up to planar
isotopy and symmetry.

CALIBRATION (target from the paper, recomputed independently):
  Lemma 3.1       c=2 -> 1   · Lemma 3.2      c=3 -> 2   · Proposition 3.3 c=4 -> 8

MODEL (combinatorial map, Euler formula)
  Vertices: n crossings (degree 4) + 2 aglets (degree 1).  Edges: e_0 ... e_2n.
  Edge e_k runs from passage k-1 to passage k; e_0 starts at aglet 0, e_2n ends at aglet 1.
  Edge end uniquely as (k, 'A') = start of e_k, (k, 'B') = end of e_k.
    -> At the crossing with word positions p1 < p2 lie exactly the four ends
       (p1,'B'), (p1+1,'A'), (p2,'B'), (p2+1,'A').
       This stays unambiguous even if p2 = p1+1 (then (p1+1,'A') and (p1+1,'B') are
       two different ends of the SAME edge) — exactly this is where the first version
       of this script failed (it returned 0/0/0 instead of 1/2/8).
  TRANSVERSALITY: at a crossing the two passages must cross, not
  touch. In the cyclic order this means: the ends of the two passages
  alternate. Exactly two arrangements satisfy this -> 2^n cases.
  Faces: cycles of sigma(alpha(.)), alpha = other end of the same edge,
  sigma = next end in the rotation of the vertex. Realizable in S^2 <=> V - E + F = 2.

REDUCED (Section 3 of the paper, implemented verbatim)
  (1) no R1 loop            -> no adjacent equal symbols
  (2) prime (no connected sum) -> no proper subword that closes a symbol set
                                         completely and disjointly
  (3) no crossing at both aglets -> first symbol != last symbol
"""
import itertools, sys
from collections import defaultdict

# ------------------------------------------------------------------ Gauss words
def gauss_woerter(n):
    basis = sum([[i, i] for i in range(1, n + 1)], [])
    gesehen = set()
    for perm in set(itertools.permutations(basis)):
        abb, naechste, wort = {}, 1, []
        for z in perm:
            if z not in abb:
                abb[z] = naechste; naechste += 1
            wort.append(abb[z])
        gesehen.add(tuple(wort))
    return sorted(gesehen)

def normalform(wort):
    kandidaten = []
    for w in (tuple(wort), tuple(reversed(wort))):
        abb, naechste, neu = {}, 1, []
        for z in w:
            if z not in abb:
                abb[z] = naechste; naechste += 1
            neu.append(abb[z])
        kandidaten.append(tuple(neu))
    return min(kandidaten)

# ------------------------------------------------------------------ reducedness
def hat_schlinge(w):
    return any(w[i] == w[i + 1] for i in range(len(w) - 1))

def ist_zerlegbar(w):
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

def beide_aglets(w):
    return w[0] == w[-1]

# ------------------------------------------------------- planar realizability
def karten_pruefen(w):
    """True if ONE transversal arrangement yields a sphere map (Euler 2)."""
    n = len(w) // 2
    m = len(w)                      # positions 0..m-1, edges 0..m
    vork = defaultdict(list)
    for pos, c in enumerate(w):
        vork[c].append(pos)
    symbole = sorted(vork)

    # alpha: other end of the same edge
    def alpha(d):
        k, s = d
        return (k, 'B' if s == 'A' else 'A')

    for wahl in itertools.product((0, 1), repeat=n):
        rot = {}                    # vertex -> cyclic list of edge ends
        for idx, c in enumerate(symbole):
            p1, p2 = vork[c]
            e1a, e1b = (p1, 'B'), (p1 + 1, 'A')     # strand 1: in, out
            e2a, e2b = (p2, 'B'), (p2 + 1, 'A')     # strand 2: in, out
            # both arrangements alternate the strands (= transversal)
            rot[('x', c)] = [e1a, e2a, e1b, e2b] if wahl[idx] == 0 else [e1a, e2b, e1b, e2a]
        rot[('aglet', 0)] = [(0, 'A')]
        rot[('aglet', 1)] = [(m, 'B')]

        # map edge end -> vertex and index in its rotation
        wo = {}
        for knoten, liste in rot.items():
            for i, d in enumerate(liste):
                if d in wo:                     # must not occur
                    return False
                wo[d] = (knoten, i)
        # completeness: 2*(m+1) ends must all be placed
        alle = {(k, s) for k in range(m + 1) for s in ('A', 'B')}
        if set(wo) != alle:
            continue

        def sigma(d):
            knoten, i = wo[d]
            liste = rot[knoten]
            return liste[(i + 1) % len(liste)]

        besucht, F = set(), 0
        for d in alle:
            if d in besucht:
                continue
            F += 1
            cur = d
            while cur not in besucht:
                besucht.add(cur)
                cur = sigma(alpha(cur))
        V, E = n + 2, m + 1
        if V - E + F == 2:
            return True
    return False

def enumerieren(n, laut=True):
    treffer = set()
    for w in gauss_woerter(n):
        if hat_schlinge(w) or beide_aglets(w) or ist_zerlegbar(w):
            continue
        if not karten_pruefen(w):
            continue
        treffer.add(normalform(w))
    if laut:
        print(f"  c={n}: {len(treffer)} reduced singular arc diagrams "
              f"-> {[''.join(map(str,t)) for t in sorted(treffer)]}")
    return treffer

if __name__ == "__main__":
    print("=" * 84)
    print("CALIBRATION against Akash-Corrado-Placke-Sanketh-Starns-Timothy-Zupan, arXiv:2606.02390")
    print("  Target: c=2 -> 1 (Lemma 3.1) · c=3 -> 2 (Lemma 3.2) · c=4 -> 8 (Proposition 3.3)")
    print("=" * 84)
    soll, ist = {2: 1, 3: 2, 4: 8}, {}
    for n in (2, 3, 4):
        ist[n] = len(enumerieren(n))
    print()
    alles = True
    for n in (2, 3, 4):
        gut = ist[n] == soll[n]; alles &= gut
        print(f"[{'OK ' if gut else 'RED'}] c={n}: ours {ist[n]} vs paper {soll[n]}")
    print()
    print("✅ ENUMERATOR CALIBRATED — step 2 may be built" if alles
          else "🔴 NOT CALIBRATED — do not build further")
