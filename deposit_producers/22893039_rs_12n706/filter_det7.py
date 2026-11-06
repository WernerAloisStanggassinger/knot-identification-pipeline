#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Step 2c: the sharp filter. Which c=4 knotoids have det(cl(kappa)) = 7?

Lemma 2.4: det(K) = det(cl(kappa))^2. For 12n706 det = 49 ==> det(cl) = 7.
A knotoid whose closure does NOT have determinant 7 cannot carry 12n706 as symmetric
union -- no matter which twist labels one chooses, because cl(kappa) depends on the labels
at all.

POSITIVE CONTROL: Table 1 (Theorem 1.5) lists for Rs4 a polynomial with det 49, namely
Delta(12n706) = 1-4t+6t^2-8t^3+11t^4-8t^5+6t^6-4t^7+t^8. So there MUST be at least one
hit -- if this run finds none, the run is broken, not the mathematics.
"""
import site_config
import itertools, sys
from collections import defaultdict
sys.path.insert(0, (site_config.ROOT + '/maschinen_ideen/RS4_12N706_20260919'))
from abschluss import aglet_ueber, det_cl
from pruef_aglet_flaeche import WOERTER

if __name__ == "__main__":
    print("=" * 96)
    print("FILTER det(cl(kappa)) = 7  —  candidates for rs(12n706) = 4")
    print("=" * 96)
    treffer, bilanz, gesamt = [], defaultdict(int), 0
    for wort in WOERTER[4]:
        w = [int(c) for c in wort]
        fest, vork = aglet_ueber(w)
        innere = sorted(c for c in vork if c not in fest)
        aglet = sorted(fest)
        print(f"\n--- {wort}   aglet crossings {aglet} (fixed, under-crossing) "
              f"· inner {innere} (free) ---")
        for wahl in itertools.product(*[[vork[c][0], vork[c][1]] for c in innere]):
            ueber = dict(fest); ueber.update(dict(zip(innere, wahl)))
            d, nw = det_cl(w, ueber)
            gesamt += 1; bilanz[d] += 1
            marke = "  <== TREFFER det 7" if d == 7 else ""
            print(f"   over(inner) = {dict(zip(innere, wahl))}   det(cl) = {d:<4}"
                  f" [{nw} paths agree]{marke}")
            if d == 7:
                treffer.append((wort, dict(zip(innere, wahl)), ueber))

    print("\n" + "=" * 96)
    print(f"Base set: {gesamt} crossing informations over {len(WOERTER[4])} diagrams")
    print("Distribution of det(cl):", dict(sorted(bilanz.items())))
    print(f"==> det(cl) = 7 in {len(treffer)} cases")
    for wort, wahl, _ in treffer:
        print(f"      {wort}  inner {wahl}")
    print()
    if treffer:
        print("✅ POSITIVE CONTROL passed: Table 1 demands at least one det-49 case.")
    else:
        print("🔴 NO HIT — contradicts Table 1 (Delta(12n706) in Rs4). Check the run!")
    print("The remaining cases are EXCLUDED for 12n706, label-independently")
    print("(cl(kappa) does not depend on the twist labels).")
