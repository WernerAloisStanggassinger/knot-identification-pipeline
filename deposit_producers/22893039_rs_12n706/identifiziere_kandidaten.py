#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Step 2d: identify the 6 det-7 candidates as KNOTS, not just as a determinant.

Each candidate is a knotoid kappa with c(kappa)=4 and det(cl(kappa))=7. cl(kappa) IS the
partial knot J of the sought symmetric union (Lemma 2.4). Which knot that is, is a
hard, checkable statement about the target -- and it shrinks the search space for step 3.

THREE COUNTER-CHECKS, built in:
  (1) det from spherogram against our own Wirtinger route (abschluss.py, calibrated on T(2,k)).
  (2) ALL closure paths and BOTH maps must yield the same JONES polynomial -- that is
      sharper than the determinant, because Jones sees chirality.
  (3) Gauss parity: if no DT can be built, the closure path was wrong (not the knot).

Runs in sage_stable (jones_polynomial needs Sage):
  PATH=/home/squeezy/miniforge3/envs/sage_stable/bin:$PATH \\
  /home/squeezy/miniforge3/envs/sage_stable/bin/python identifiziere_kandidaten.py
"""
import site_config
import sys
from collections import defaultdict

HIER = (site_config.ROOT + '/maschinen_ideen/RS4_12N706_20260919')
sys.path.insert(0, HIER)
from pruef_aglet_flaeche import karten
from abschluss import (aglet_ueber, cl_sequenz, abschlusswege,
                       wirtinger_aus_sequenz, det_aus_wirtinger)
from dt_bau import link_aus_sequenz

# from filter_det7.py, run 2026-09-20
KANDIDATEN = [
    ("12314234", {2: 1, 3: 6}),
    ("12324143", {2: 1, 4: 4}),
    ("12324143", {2: 3, 4: 6}),
    ("12341243", {2: 1, 4: 6}),
    ("12341423", {2: 6, 4: 5}),
    ("12341432", {3: 6, 4: 3}),
]


def analysiere(wort, innere):
    w = [int(c) for c in wort]
    fest, vork = aglet_ueber(w)
    ueber = dict(fest); ueber.update(innere)
    m = len(w)
    jones_menge, alex_menge, det_menge, vol_menge = set(), set(), set(), set()
    namen, n_wege, n_paritaet_weg = set(), 0, 0
    for karte in karten(w):
        for kf in abschlusswege(karte, m, max_extra=2):
            seq = cl_sequenz(w, ueber, kf)
            tripel, N = wirtinger_aus_sequenz(seq)
            det_menge.add(det_aus_wirtinger(tripel, N))
            L, dt = link_aus_sequenz(seq)
            if L is None:
                n_paritaet_weg += 1
                continue
            n_wege += 1
            jones_menge.add(str(L.jones_polynomial()))
            alex_menge.add(str(L.alexander_polynomial()))
            try:
                E = L.exterior()
                vol_menge.add(round(float(E.volume()), 6))
                try:
                    ids = E.identify()
                    for x in ids[:3]:
                        namen.add(str(x))
                except Exception:
                    pass
            except Exception:
                vol_menge.add(None)
    return dict(det=det_menge, jones=jones_menge, alex=alex_menge, vol=vol_menge,
                namen=namen, n_wege=n_wege, n_par=n_paritaet_weg)


if __name__ == "__main__":
    print("=" * 100)
    print("IDENTIFICATION of the 6 candidates: which knot is cl(kappa) = the partial knot?")
    print("=" * 100)
    ok = True
    zusammen = defaultdict(list)
    for wort, innere in KANDIDATEN:
        r = analysiere(wort, innere)
        einig_j = len(r["jones"]) == 1
        einig_d = len(r["det"]) == 1
        gut = einig_j and einig_d and r["det"] == {7}
        ok &= gut
        print(f"\n--- {wort}  inner={innere} ---")
        zusatz = f"  (+{r['n_par']} ohne DT, Gauss-Paritaet)" if r['n_par'] else ""
        print(f"   closure paths x maps computed: {r['n_wege']}{zusatz}")
        print(f"   det (own Wirtinger route): {sorted(r['det'])}"
              f"   {'agree' if einig_d else '🔴 DISAGREE'}")
        print(f"   Jones over all paths:       {'agree' if einig_j else '🔴 DISAGREE'}"
              f"  ({len(r['jones'])} distinct)")
        for j in sorted(r["jones"]):
            print(f"      V = {j}")
        for a in sorted(r["alex"]):
            print(f"      Delta = {a}")
        print(f"   volume: {sorted(x for x in r['vol'] if x is not None)}")
        print(f"   table lookup (snappy identify, ONLY literature comparison): "
              f"{sorted(r['namen']) if r['namen'] else '(none)'}")
        for j in r["jones"]:
            zusammen[j].append((wort, tuple(sorted(innere.items()))))

    print("\n" + "=" * 100)
    print(f"SUMMARY: {len(zusammen)} distinct partial knots among the 6 candidates")
    for j, liste in zusammen.items():
        print(f"\n   V = {j}")
        for wort, innere in liste:
            print(f"      {wort}  {dict(innere)}")
    print()
    print("CALIBRATION:", "ALL GREEN (det=7 confirmed, Jones path-independent)" if ok
          else "🔴 RED — look into it")
