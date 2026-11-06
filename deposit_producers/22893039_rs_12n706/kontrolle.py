#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""COUNTER-CASE CHECK for the 12n706 finding. A procedure that finds every target knot,
measures nothing. So measure both:

NEGATIVE: all knots from Table 3 whose r_s >= 5 is PROVEN (column "lower" = Thm. 1.5),
         must NOT be found -- their Delta must not be among our 27.
POSITIVE: knots from Table 3 with r_s = 4 PROVEN must be found -- otherwise
         our procedure finds only the one case we want to see.
"""
import site_config
import itertools, sys, signal
HIER = (site_config.ROOT + '/maschinen_ideen/RS4_12N706_20260919')
sys.path.insert(0, HIER)
from snappy import *
from spherogram import Link
from eiche_union import koeffs, koeffs_von_liste, faelle
from band_rand import alexander_von, band_rand_sequenz
from abschluss import aglet_ueber
from pruef_aglet_flaeche import karten, WOERTER
from dt_bau import link_aus_sequenz

# Table 3, copied: (name, det, Delta coefficients, r_s status)
NEGATIV = [   # r_s >= 5 proven via Thm. 1.5  -> must NOT show up
    ("12n4",   [1,-7,19,-27,19,-7,1]), ("12n43",  [1,-5,10,-15,19,-15,10,-5,1]),
    ("12n49",  [4,-20,33,-20,4]),      ("12n106", [1,-4,10,-16,19,-16,10,-4,1]),
    ("12n170", [6,-20,29,-20,6]),      ("12n380", [2,-9,18,-23,18,-9,2]),
    ("12n504", [1,-6,15,-24,29,-24,15,-6,1]), ("12n676", [2,-2,-4,9,-4,-2,2]),
    ("12n702", [2,-12,28,-37,28,-12,2]),("12n782", [2,-8,18,-25,18,-8,2]),
    ("12n802", [1,-5,14,-25,31,-25,14,-5,1]), ("12n657", [1,-4,9,-16,21,-16,9,-4,1]),
]
POSITIV = [   # r_s = 4 proven -> must be found
    ("12n48",  [2,-12,21,-12,2]), ("12n87",  [4,-12,17,-12,4]),
    ("12n708", [1,-3,6,-9,11,-9,6,-3,1]), ("12n817", [2,-6,10,-13,10,-6,2]),
]
GRENZE = 3


def mit_frist(f, sek=120):
    def h(s, fr): raise TimeoutError()
    alt = signal.signal(signal.SIGALRM, h); signal.alarm(sek)
    try:    return f()
    except Exception as e: return f"({type(e).__name__})"
    finally: signal.alarm(0); signal.signal(signal.SIGALRM, alt)


if __name__ == "__main__":
    print("=" * 98)
    print("Step 1 — compute all 256 mod-2 diagrams once, Delta -> remember sites")
    print("=" * 98)
    orte = {}
    for wort in WOERTER[4]:
        for w, ueber, lab, innere, n_str in faelle(wort):
            for s in alexander_von(w, ueber, lab, True, True):
                k = koeffs(s)
                if k:
                    orte.setdefault(k, []).append((wort, dict(ueber), dict(lab)))
    print(f"  {len(orte)} distinct Delta (Table 1 has 27)")

    print("\n" + "=" * 98)
    print("NEGATIVE CONTROL — r_s >= 5 proven, must NOT show up")
    print("=" * 98)
    nok = 0
    for name, d in NEGATIV:
        drin = koeffs_von_liste(d) in orte
        nok += (not drin)
        print(f"  [{'OK ' if not drin else 'RED'}] {name:<8} Delta in our 27? "
              f"{'YES — CONTRADICTION' if drin else 'no'}")
    print(f"  ==> {nok}/{len(NEGATIV)} correctly not found")

    print("\n" + "=" * 98)
    print(f"POSITIVE CONTROL — r_s = 4 proven, MUST be found (|a| <= {GRENZE})")
    print("=" * 98)
    ung = [a for a in range(-GRENZE, GRENZE + 1) if a % 2]
    ger = [a for a in range(-GRENZE, GRENZE + 1) if a % 2 == 0]
    pok = 0
    for name, d in POSITIV + [("12n706", [1,-4,6,-8,11,-8,6,-4,1])]:
        k = koeffs_von_liste(d)
        stellen = orte.get(k, [])
        ziel = Link('K' + name)
        ze = ziel.exterior()
        gefunden = None
        for wort, ueber, lab0 in stellen:
            if gefunden: break
            w = [int(c) for c in wort]
            frei = sorted(lab0)
            for karte in karten(w):
                if gefunden: break
                for werte in itertools.product(*[(ung if lab0[j] % 2 else ger) for j in frei]):
                    seq = band_rand_sequenz(w, ueber, karte, dict(zip(frei, werte)), True, True)
                    L, dt = link_aus_sequenz(seq)
                    if L is None or len(L.link_components) != 1: continue
                    L.simplify('global')
                    if len(L.crossings) > 16: continue
                    E = L.exterior()
                    if mit_frist(lambda: bool(E.is_isometric_to(ze)), 90) is True:
                        gefunden = (wort, dict(zip(frei, werte)), len(L.crossings), dt)
                        break
        pok += bool(gefunden)
        if gefunden:
            print(f"  [OK ] {name:<8} FOUND: {gefunden[0]} labels={gefunden[1]} "
                  f"c={gefunden[2]}  (is_isometric_to confirmed)")
            print(f"           DT {gefunden[3]}")
        else:
            print(f"  [RED] {name:<8} NOT found at {len(stellen)} Delta sites")
    print(f"  ==> {pok}/{len(POSITIV)+1} found")
