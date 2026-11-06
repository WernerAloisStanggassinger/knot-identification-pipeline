#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""THE SHARPEST CALIBRATION: run through all 256 mod-2-labeled knotoid diagrams with c=4 and
hold it against Table 1 (Theorem 1.5: R^s_4 = 27 polynomials). If the constructor hits the set
EXACT, it is completely proven -- and then the question about 12n706 is cleanly posed.

Afterwards: does one of the 6 det-7 candidates carry the Alexander polynomial of 12n706?
  Delta(12n706) = 1 - 4t + 6t^2 - 8t^3 + 11t^4 - 8t^5 + 6t^6 - 4t^7 + t^8   (Table 1, det 49)
If NO, r_s(12n706) = 4 is EXCLUDED and thereby r_s(12n706) >= 5 is proven.
If YES, these are the sites where the disc must lie.
"""
import site_config
import itertools, sys
HIER = (site_config.ROOT + '/maschinen_ideen/RS4_12N706_20260919')
sys.path.insert(0, HIER)
from eiche_union import koeffs, koeffs_von_liste, faelle
from band_rand import alexander_von
from pruef_aglet_flaeche import WOERTER

# Table 1, copied verbatim (27 polynomials, ordered by determinant).
# The misprint in the line "1 - 5t + 11t^2 - 15t^3 + 11t^4 - 5t^6 + t^6" (t^6 twice) is
# read as -5t^5 + t^6 -- only so does the determinant 49 stated there result.
TAB1 = [
    [1], [1,0,-3,0,1], [1,-1,-1,3,-1,-1,1], [1,-2,1,2,-5,2,1,-2,1],
    [1,0,-2,0,3,0,-2,0,1], [1,-3,-1,7,-1,-3,1], [1,-3,2,3,-7,3,2,-3,1],
    [2,-3,-2,7,-2,-3,2], [2,0,-5,0,2],
    [1,-1,0,-1,3,-1,0,-1,1], [1,-1,1,-3,1,-1,1], [1,-1,-3,7,-3,-1,1],
    [1,-2,3,-2,1], [1,-2,0,4,-7,4,0,-2,1], [2,-5,2],
    [1,-2,3,-4,5,-4,3,-2,1], [1,-3,5,-7,5,-3,1], [1,-6,11,-6,1],
    [2,-6,9,-6,2], [6,-13,6],
    [1,-3,6,-9,11,-9,6,-3,1], [1,-4,6,-8,11,-8,6,-4,1],
    [1,-5,11,-15,11,-5,1], [2,-6,10,-13,10,-6,2],
    [2,-12,21,-12,2], [3,-12,19,-12,3], [4,-12,17,-12,4],
]
ZIEL = [1,-4,6,-8,11,-8,6,-4,1]          # Delta(12n706)

KANDIDATEN = [("12314234", {2:1,3:6}), ("12324143", {2:1,4:4}), ("12324143", {2:3,4:6}),
              ("12341243", {2:1,4:6}), ("12341423", {2:6,4:5}), ("12341432", {3:6,4:3})]

if __name__ == "__main__":
    soll = {koeffs_von_liste(v) for v in TAB1}
    ziel = koeffs_von_liste(ZIEL)
    print("=" * 98)
    print(f"R^s_4 FULL CHECK: Table 1 has {len(TAB1)} entries, {len(soll)} distinct")
    print(f"Target polynomial Delta(12n706) in Table 1? {'YES' if ziel in soll else 'NO'}")
    print("=" * 98)
    det = lambda v: abs(sum(c * (-1)**i for i, c in enumerate(v)))
    print("Determinant counter-check over Table 1:",
          {d: sum(1 for v in TAB1 if det(v) == d) for d in sorted({det(v) for v in TAB1})})

    ist, n, ohne_dt, links = set(), 0, 0, 0
    treffer = []
    for wort in WOERTER[4]:
        for w, ueber, lab, innere, n_str in faelle(wort):
            n += 1
            for s in alexander_von(w, ueber, lab, True, True):
                if s is None:
                    ohne_dt += 1
                elif isinstance(s, str) and s.startswith("LINK"):
                    links += 1
                else:
                    k = koeffs(s)
                    ist.add(k)
                    if k == ziel:
                        treffer.append((wort, dict(innere), dict(lab)))
    print(f"\ncomputed: {n} mod-2 diagrams  (expected 8 diagrams x 4 crossing infos x 8 labels = 256)")
    print(f"of these without DT: {ohne_dt} · multi-component: {links}")
    print(f"generated polynomials: {len(ist)}   Table 1: {len(soll)}")
    print(f"  generated \\ Table1 : {sorted(ist - soll)}")
    print(f"  Table1 \\ generated : {sorted(soll - ist)}")
    print(f"  ==> {'HIT EXACTLY' if ist == soll else 'DEVIATION'}")

    print("\n" + "=" * 98)
    print(f"DOES A CANDIDATE CARRY Delta(12n706)?  -> {len(treffer)} hits")
    for t in treffer:
        print(f"   {t[0]}  inner={t[1]}  labels={t[2]}")
    if not treffer:
        print("   NONE. Then r_s(12n706) >= 5 would hold — but only if the full check above")
        print("   hit exactly; otherwise the constructor is to blame, not the mathematics.")
