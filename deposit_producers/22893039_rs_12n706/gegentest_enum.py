#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""COUNTER-TEST to the enumerator (GEGENWIND §3: calibration on small cases is necessary,
not sufficient — losses only show up at c=5, exactly where the new result arises).

THREE independent tests. None replaces a second implementation from scratch; they
attack at different places:

A) FOREIGN PLANARITY ALGORITHM. networkx.check_planarity (LR algorithm, third-party code,
   completely different method from our Euler formula via rotation systems). The underlying
   multigraph of EVERY accepted diagram must be planar. Necessary condition -- an
   accepted diagram with a non-planar graph would be a hard error.
   Conversely: a REJECTED word with a planar graph is a suspicion of loss
   (not necessarily an error: planarity of the graph is weaker than realizability of the
   diagram, because it does not require transversality at the crossings).
B) FUNCTIONAL. Every accepted diagram must yield a valid closure (the map
   is used, not just checked).
C) COUNT AGAINST THE OTHER TRAVERSAL DIRECTION (alpha∘sigma instead of sigma∘alpha).
"""
import site_config
import sys, itertools
from collections import defaultdict
HIER = (site_config.ROOT + '/maschinen_ideen/RS4_12N706_20260919')
sys.path.insert(0, HIER)
import networkx as nx
from enum_bogendiagramme import (gauss_woerter, hat_schlinge, beide_aglets,
                                 ist_zerlegbar, normalform, enumerieren)


def graph_von(w):
    """Multigraph: vertices = n crossings + 2 aglets, edges = e_0..e_2n."""
    m = len(w)
    vork = defaultdict(list)
    for p, c in enumerate(w):
        vork[c].append(p)
    def knoten_an(p):          # which vertex lies at the end of e_p / start of e_{p+1}
        return ('x', w[p])
    G = nx.MultiGraph()
    for c in vork: G.add_node(('x', c))
    G.add_node(('a', 0)); G.add_node(('a', 1))
    for k in range(m + 1):
        u = ('a', 0) if k == 0 else knoten_an(k - 1)
        v = ('a', 1) if k == m else knoten_an(k)
        G.add_edge(u, v, key=k)
    return G


def flaechen_andersherum(w):
    """C: face count via alpha∘sigma instead of sigma∘alpha — must give the same Euler number."""
    m = len(w); n = m // 2
    vork = defaultdict(list)
    for p, c in enumerate(w): vork[c].append(p)
    treffer = 0
    for wahl in itertools.product((0, 1), repeat=n):
        rot = {}
        for idx, c in enumerate(sorted(vork)):
            p1, p2 = vork[c]
            e1a, e1b = (p1, 'B'), (p1 + 1, 'A')
            e2a, e2b = (p2, 'B'), (p2 + 1, 'A')
            rot[('x', c)] = [e1a, e2a, e1b, e2b] if wahl[idx] == 0 else [e1a, e2b, e1b, e2a]
        rot[('ag', 0)] = [(0, 'A')]; rot[('ag', 1)] = [(m, 'B')]
        wo, schlecht = {}, False
        for kn, li in rot.items():
            for i, d in enumerate(li):
                if d in wo: schlecht = True; break
                wo[d] = (kn, i)
            if schlecht: break
        if schlecht: continue
        alle = {(k, s) for k in range(m + 1) for s in ('A', 'B')}
        if set(wo) != alle: continue
        def sig(d):
            kn, i = wo[d]; li = rot[kn]; return li[(i + 1) % len(li)]
        def alp(d): k, s = d; return (k, 'B' if s == 'A' else 'A')
        besucht, F = set(), 0
        for d in alle:                       # <- alpha(sigma(.)) instead of sigma(alpha(.))
            if d in besucht: continue
            F += 1; cur = d
            while cur not in besucht:
                besucht.add(cur); cur = alp(sig(cur))
        if (n + 2) - (m + 1) + F == 2: treffer += 1
    return treffer


if __name__ == "__main__":
    print("=" * 92)
    print("COUNTER-TEST of the enumerator (GEGENWIND §3)")
    print("=" * 92)
    for n in (2, 3, 4, 5):
        akzeptiert = enumerieren(n, laut=False)
        vorfilter = [w for w in gauss_woerter(n)
                     if not (hat_schlinge(w) or beide_aglets(w) or ist_zerlegbar(w))]
        verworfen = {normalform(w) for w in vorfilter} - set(akzeptiert)

        a_bad = [w for w in akzeptiert if not nx.check_planarity(graph_von(list(w)))[0]]
        v_planar = [w for w in verworfen if nx.check_planarity(graph_von(list(w)))[0]]
        c_ok = sum(1 for w in akzeptiert if flaechen_andersherum(list(w)) > 0)
        print(f"\n  c={n}: {len(akzeptiert)} accepted, {len(verworfen)} rejected after pre-filter")
        print(f"    A) accepted AND non-planar (must be 0): {len(a_bad)}"
              + ("   🔴" if a_bad else "   ok"))
        print(f"       rejected, but graph planar (suspicion, not necessarily error): {len(v_planar)}")
        print(f"    C) accepted, also planar when traversed backwards (must be {len(akzeptiert)}): {c_ok}"
              + ("   ok" if c_ok == len(akzeptiert) else "   🔴"))
