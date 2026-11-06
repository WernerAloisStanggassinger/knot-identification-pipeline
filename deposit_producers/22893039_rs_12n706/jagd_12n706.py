#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""THE CLOSURE: is the knot at the found site really 12n706?

Delta proves NOTHING (not a complete invariant) -- therefore three stages:
  1. Jones bit-identical (sees chirality, unlike Delta)
  2. volume of the complement
  3. is_isometric_to against 12n706 == PROOF (Gordon-Luecke), as in the corpus on 17./18.09.

Site (from rs4_voll.py): arc diagram 12314234, inner crossings {2:1, 3:6},
Labels (1,2,3) = (odd, even, odd). Sought are the REAL labels.
flag_seite is varied along: it produces the mirror image, and Jones sees the difference.
"""
import site_config
import itertools, sys, signal
HIER = (site_config.ROOT + '/maschinen_ideen/RS4_12N706_20260919')
sys.path.insert(0, HIER)
from snappy import *
from spherogram import Link
from abschluss import aglet_ueber
from band_rand import band_rand_sequenz
from pruef_aglet_flaeche import karten
from dt_bau import link_aus_sequenz

GRENZE = 5


def mit_frist(f, sek=120, *a, **k):
    def h(s, fr): raise TimeoutError()
    alt = signal.signal(signal.SIGALRM, h); signal.alarm(sek)
    try:
        return f(*a, **k)
    except Exception as e:
        return f"({type(e).__name__})"
    finally:
        signal.alarm(0); signal.signal(signal.SIGALRM, alt)


if __name__ == "__main__":
    ziel = Link('K12n706')
    zj = ziel.jones_polynomial()
    zjm = ziel.mirror().jones_polynomial()
    ze = ziel.exterior()
    zv = float(ze.volume())
    print("=" * 98)
    print(f"TARGET 12n706: crossings {len(ziel.crossings)}  vol {zv:.10f}")
    print(f"  Jones        {zj}")
    print(f"  Jones mirror {zjm}")
    print("=" * 98)

    w = [int(c) for c in '12314234']
    fest, vork = aglet_ueber(w)
    ueber = dict(fest); ueber.update({2: 1, 3: 6})
    ung = [a for a in range(-GRENZE, GRENZE + 1) if a % 2]
    ger = [a for a in range(-GRENZE, GRENZE + 1) if a % 2 == 0]

    n = 0; jones_treffer = []; bestand = set()
    for flag_seite in (True, False):
        for karte in karten(w):
            for a1, a2, a3 in itertools.product(ung, ger, ung):
                n += 1
                seq = band_rand_sequenz(w, ueber, karte, {1: a1, 2: a2, 3: a3}, True, flag_seite)
                L, dt = link_aus_sequenz(seq)
                if L is None or len(L.link_components) != 1:
                    continue
                L.simplify('global')
                bestand.add(len(L.crossings))
                j = mit_frist(L.jones_polynomial, 90)
                if j == zj or j == zjm:
                    jones_treffer.append((flag_seite, a1, a2, a3, len(L.crossings),
                                          'direkt' if j == zj else 'Spiegel', L))
    print(f"\ncomputed: {n} label combinations (|a| <= {GRENZE}) x 2 maps x 2 side conventions")
    print(f"crossing numbers after simplify: {sorted(bestand)}")
    print(f"JONES HITS: {len(jones_treffer)}")

    for flag_seite, a1, a2, a3, nc, art, L in jones_treffer[:12]:
        E = L.exterior()
        v = mit_frist(lambda: float(E.volume()), 120)
        iso = mit_frist(lambda: bool(E.is_isometric_to(ze)), 240)
        ids = mit_frist(lambda: [str(x) for x in E.identify()], 120)
        print(f"  side={flag_seite} labels=({a1},{a2},{a3}) c={nc} Jones {art}")
        print(f"      vol {v}   is_isometric_to(12n706) = {iso}   identify {ids}")
    if not jones_treffer:
        print("  NONE — then this site carries Delta(12n706) but not the knot.")
        print("  That would be NO proof of r_s >= 5 as long as the label bound |a| <= "
              f"{GRENZE} is not justified.")
