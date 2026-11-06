#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Step 2a: For every reduced arc diagram determine the planar maps and
ask whether the two aglets border a COMMON face.

WHY THIS IS THE DECISIVE QUESTION
The closure cl(kappa) connects the aglets by an arc that only has over-crossings
(paper p. 5). If both aglets lie on the same face, this arc can be drawn WITHOUT CROSSINGS
draw — then cl(kappa) is simply the cyclically closed Gauss word, and the determinant
of the partial knot is directly computable. That is the filter from Lemma 2.4:
      det(K) = det(cl(kappa))^2,  so for 12n706 det(cl(kappa)) = 7.
If they lie on different faces, the arc must cross and the closure becomes more laborious.
"""
import site_config
import itertools, sys
from collections import defaultdict
sys.path.insert(0, (site_config.ROOT + '/maschinen_ideen/RS4_12N706_20260919'))
from enum_bogendiagramme import enumerieren

def karten(w):
    """All planar maps (Euler 2) for a Gauss word; per map the faces as
    sets of edge ends, plus the faces at the two aglets."""
    n = len(w) // 2
    m = len(w)
    vork = defaultdict(list)
    for pos, c in enumerate(w):
        vork[c].append(pos)
    symbole = sorted(vork)
    def alpha(d):
        k, s = d
        return (k, 'B' if s == 'A' else 'A')
    gefunden = []
    for wahl in itertools.product((0, 1), repeat=n):
        rot = {}
        for idx, c in enumerate(symbole):
            p1, p2 = vork[c]
            e1a, e1b = (p1, 'B'), (p1 + 1, 'A')
            e2a, e2b = (p2, 'B'), (p2 + 1, 'A')
            rot[('x', c)] = [e1a, e2a, e1b, e2b] if wahl[idx] == 0 else [e1a, e2b, e1b, e2a]
        rot[('aglet', 0)] = [(0, 'A')]
        rot[('aglet', 1)] = [(m, 'B')]
        wo = {}
        schlecht = False
        for knoten, liste in rot.items():
            for i, d in enumerate(liste):
                if d in wo: schlecht = True; break
                wo[d] = (knoten, i)
            if schlecht: break
        if schlecht: continue
        alle = {(k, s) for k in range(m + 1) for s in ('A', 'B')}
        if set(wo) != alle: continue
        def sigma(d):
            kn, i = wo[d]; li = rot[kn]
            return li[(i + 1) % len(li)]
        besucht, flaechen = set(), []
        for d in sorted(alle):
            if d in besucht: continue
            zyklus, cur = [], d
            while cur not in besucht:
                besucht.add(cur); zyklus.append(cur); cur = sigma(alpha(cur))
            flaechen.append(set(zyklus))
        if (n + 2) - (m + 1) + len(flaechen) != 2: continue
        # which face does each aglet lie on? (edge end at the aglet)
        f0 = [i for i, f in enumerate(flaechen) if (0, 'A') in f]
        f1 = [i for i, f in enumerate(flaechen) if (m, 'B') in f]
        gefunden.append({"wahl": wahl, "flaechen": len(flaechen),
                         "aglet0": set(f0), "aglet1": set(f1),
                         "gemeinsam": bool(set(f0) & set(f1))})
    return gefunden

if __name__ == "__main__":
    print("=" * 86)
    print("MAP ANALYSIS: do both aglets border a common face?")
    print("=" * 86)
    for n in (2, 3, 4):
        print(f"\n--- c = {n} ---")
        for w in sorted(enumerieren(n, laut=False)):
            ks = karten(w)
            mit = sum(1 for k in ks if k["gemeinsam"])
            wort = "".join(map(str, w))
            print(f"  {wort:<10} planar maps: {len(ks):>2}   of these with a common "
                  f"aglet face: {mit:>2}   -> closure without crossings possible: "
                  f"{'YES' if mit else 'NO'}")
