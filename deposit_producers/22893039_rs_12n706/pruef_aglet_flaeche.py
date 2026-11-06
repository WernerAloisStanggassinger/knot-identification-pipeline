#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""CALIBRATION of the aglet-face question via TWO independent routes (step 2a).

QUESTION: do both aglets of an arc diagram border a common face? Equivalent:
can the closure cl(kappa) be drawn WITHOUT CROSSINGS?

🔴 OCCASION (20.09.2026): I wrongly took "both aglets on one face" to be NECESSARY
(argument: kappa = half of D # -D, hence lies in a half-plane with both aglets on the
axis). The map code reported ZERO for '1212' = 2_1, and I suspected the code.
Wrong — the EXPECTATION was wrong. The paper itself says it (Lemma 5.2): the closure may
have more crossings than kappa. For 2_1 it even MUST: rs = 2 demands det(K) = 9
(Corollary 5.1), hence a partial knot with det 3 — and a diagram with only 2 crossings
is the unknot.
  ==> Lesson: a calibration probe needs the right expectation, otherwise it measures the checker.

ROUTE A (combinatorial map): faces as orbits of sigma∘alpha, Euler check V-E+F=2,
        then check whether both aglet ends lie in the same orbit.
ROUTE B (Gauss' parity condition, completely independent): an arc can be closed
        without crossings exactly if the CYCLICALLY closed Gauss word is realizable as a closed
        curve in S^2. Necessary for that is Gauss' parity: between the
        two occurrences of each symbol there is an EVEN number of symbols.
The two routes must never contradict. That is the calibration.
"""
import itertools, sys
from collections import defaultdict


# --------------------------------------------------------------- ROUTE A: maps
def karten(w, laut=False):
    n, m = len(w) // 2, len(w)
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
        wo, schlecht = {}, False
        for knoten, liste in rot.items():
            for i, d in enumerate(liste):
                if d in wo:
                    schlecht = True; break
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
        f0 = {i for i, f in enumerate(flaechen) if (0, 'A') in f}
        f1 = {i for i, f in enumerate(flaechen) if (m, 'B') in f}
        gefunden.append({"wahl": wahl, "flaechen": flaechen, "rot": rot, "wo": wo,
                         "aglet0": f0, "aglet1": f1, "gemeinsam": bool(f0 & f1)})
        if laut:
            print(f"    map {wahl}: F={len(flaechen)} aglet0={f0} aglet1={f1}"
                  f" -> {'SHARED' if f0 & f1 else 'separate'}")
    return gefunden


def weg_a(w):
    return any(k["gemeinsam"] for k in karten(w))


# ------------------------------------------------- ROUTE B: Gauss' parity
def weg_b(w):
    """Necessary condition for the realizability of the cyclically closed word."""
    vork = defaultdict(list)
    for pos, c in enumerate(w):
        vork[c].append(pos)
    for c, (p1, p2) in vork.items():
        if (p2 - p1 - 1) % 2 != 0:
            return False
    return True


# ------------------------------------------------------------------- calibration
WOERTER = {
    2: ["1212"],
    3: ["123123", "123132"],
    4: ["12314234", "12314324", "12324143", "12341234",
        "12341243", "12341423", "12341432", "12342143"],
}


def selbsttest(laut=False):
    print("=" * 92)
    print("CALIBRATION: aglet faces via TWO independent routes")
    print("  ROUTE A = combinatorial map (Euler)   ROUTE B = Gauss' parity condition")
    print("=" * 92)
    ok, n_probe = True, 0

    print("\n-- Positive controls: cases whose answer is settled independently --")
    # A cut-open CLOSED diagram always has the aglets next to each other.
    for name, wort, soll in [
        ("leeres Wort (eingebetteter Bogen)", "", True),
        ("'11'   R1-Schlinge, beide Enden aussen", "11", True),
        ("'1122' zwei R1-Schlingen", "1122", True),
        ("'123123' = Kleeblatt aufgeschnitten", "123123", True),
    ]:
        w = [int(c) for c in wort]
        a, b = weg_a(w), weg_b(w)
        gut = (a == soll) and (b == soll); ok &= gut; n_probe += 1
        print(f"  [{'OK ' if gut else 'RED'}] {name:<42} A={str(a):<5} B={str(b):<5} target={soll}")

    print("\n-- COUNTER-CASE: '1212' MUST be separate (Gauss parity violated) --")
    w = [1, 2, 1, 2]
    a, b = weg_a(w), weg_b(w)
    gut = (a is False) and (b is False); ok &= gut; n_probe += 1
    print(f"  [{'OK ' if gut else 'RED'}] '1212' = 2_1                               "
          f"A={str(a):<5} B={str(b):<5} target=False")
    print("       Receipt: between the two 1s there is exactly ONE symbol (odd).")
    print("       Confirmed on the substance by Corollary 5.1: rs=2 ==> det(K)=9 ==> partial knot")
    print("       with det 3 ==> at least 3 crossings ==> the closure MUST cross.")

    print("\n-- FULL COMPARISON of both routes over all reduced diagrams c=2,3,4 --")
    for n in (2, 3, 4):
        for wort in WOERTER[n]:
            w = [int(c) for c in wort]
            a, b = weg_a(w), weg_b(w)
            gut = (a == b); ok &= gut; n_probe += 1
            print(f"  [{'OK ' if gut else 'RED'}] c={n}  {wort:<10} "
                  f"ROUTE A = {str(a):<5}  ROUTE B = {str(b):<5}"
                  f"{'' if gut else '   <-- CONTRADICTION'}")
            if laut:
                karten(w, laut=True)

    print()
    print(f"Base set: {n_probe} probes (not 0 — footgun 32).")
    print("CALIBRATION:", "ALL GREEN — both routes agree" if ok
          else "🔴 RED — the routes contradict each other")
    return ok


if __name__ == "__main__":
    sys.exit(0 if selbsttest("--laut" in sys.argv) else 1)
