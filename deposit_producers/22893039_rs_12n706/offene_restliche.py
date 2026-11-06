#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""The remaining three open cases from Table 3: 12n23 {4,5} · 12n49 {5,6,7,8} · 12n504 {5,6}.
12n49 (det 81) and 12n504 (det 121) cannot have r_s=4 at all: Corollary 5.1 caps
det <= 49. That leaves 12n23 (det 9, Delta = 2-5t+2t^2, lies in R^s_4) -- decidable.
"""
import site_config
import itertools, sys, signal
HIER = (site_config.ROOT + '/maschinen_ideen/RS4_12N706_20260919')
sys.path.insert(0, HIER)
from snappy import *
from spherogram import Link
from eiche_union import koeffs, koeffs_von_liste, faelle
from band_rand import alexander_von, band_rand_sequenz
from pruef_aglet_flaeche import karten, WOERTER
from dt_bau import link_aus_sequenz

GRENZE = 5
def mit_frist(f, sek=90):
    def h(s, fr): raise TimeoutError()
    alt = signal.signal(signal.SIGALRM, h); signal.alarm(sek)
    try: return f()
    except Exception as e: return f"({type(e).__name__})"
    finally: signal.alarm(0); signal.signal(signal.SIGALRM, alt)

if __name__ == "__main__":
    orte = {}
    for wort in WOERTER[4]:
        for w, ueber, lab, innere, n_str in faelle(wort):
            for s in alexander_von(w, ueber, lab, True, True):
                k = koeffs(s)
                if k: orte.setdefault(k, []).append((wort, dict(ueber), dict(lab)))
    print("=" * 96)
    for name, d, det in [("12n23", [2,-5,2], 9), ("12n49", [4,-20,33,-20,4], 81),
                         ("12n504", [1,-6,15,-24,29,-24,15,-6,1], 121)]:
        k = koeffs_von_liste(d)
        stellen = orte.get(k, [])
        print(f"\n--- {name}  det {det}  Delta in R^s_4? {'YES' if stellen else 'NO'}"
              f"  ({len(stellen)} sites) ---")
        if not stellen:
            print(f"    ==> r_s = 4 EXCLUDED (Corollary 5.1: r_s=4 ==> det <= 49)")
            continue
        ziel = Link('K' + name); ze = ziel.exterior()
        ung = [a for a in range(-GRENZE, GRENZE+1) if a % 2]
        ger = [a for a in range(-GRENZE, GRENZE+1) if a % 2 == 0]
        gefunden, geprueft = None, 0
        for wort, ueber, lab0 in stellen:
            if gefunden: break
            w = [int(c) for c in wort]; frei = sorted(lab0)
            for karte in karten(w):
                if gefunden: break
                for werte in itertools.product(*[(ung if lab0[j]%2 else ger) for j in frei]):
                    seq = band_rand_sequenz(w, ueber, karte, dict(zip(frei, werte)), True, True)
                    L, dt = link_aus_sequenz(seq)
                    if L is None or len(L.link_components) != 1: continue
                    L.simplify('global'); geprueft += 1
                    if len(L.crossings) > 16: continue
                    E = L.exterior()
                    if mit_frist(lambda: bool(E.is_isometric_to(ze))) is True:
                        gefunden = (wort, dict(zip(frei, werte)), len(L.crossings), dt); break
        if gefunden:
            print(f"    ✅ FOUND: {gefunden[0]} labels={gefunden[1]} c={gefunden[2]}")
            print(f"       DT {gefunden[3]}")
            print(f"    ==> r_s({name}) <= 4")
        else:
            print(f"    not found ({geprueft} diagrams checked, |a| <= {GRENZE})")
            print(f"    ==> NO proof of r_s >= 5 — only: not present in this search space")
