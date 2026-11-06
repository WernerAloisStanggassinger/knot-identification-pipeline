#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""STAGE 2 — Theta on the freshly built mutant pairs. Runs in `sage_stable`.

First DE-DUPLICATE: the 30 pairs from stage 1 contain repetitions. Pretzels
are invariant under cyclic permutation of their parameters, so P(2,3,5,7),
P(2,3,7,5) and P(2,5,7,3) are THE SAME knot — with equal determinant, equal
volume, equal isometry_signature. Whoever counts them one by one reports 30 cases,
where there are ten. De-duplication is done on the pair of isometry_signatures, not on the name.

Then theta at three points (`Theta_numeric`), not symbolically: a different
VALUE proves different functions and thereby the separation. Symbolic is
needlessly expensive at c=29 and not necessary for the statement.
[[feedback_point_eval_beats_symbolic_form]]

Counter-check built in: the pair P(2,3,5,7) is already in the corpus and is separated there by
theta. If it does NOT come out as separated here, the chain is wrong.
"""
import site_config
import json
import os
import sys

sys.path.insert(0, (site_config.ROOT + '/freeze_manifest/pipeline_scripts'))
from pd_to_theta_input import pd_to_cs_phi

# The way the night-run worker does it: sage_cmdline + load() with ONE argument.
# (`from sage.all import ...` gives a load() that insists on globals().)
from sage.all_cmdline import *   # noqa: F401,F403

BN = (site_config.ROOT + '/external_tools/bar_natan_theta')
load(os.path.join(BN, "Theta_numeric.sage"))

PUNKTE = [(QQ(22) / 7, QQ(21) / 13), (QQ(5) / 3, QQ(-7) / 4), (QQ(13) / 8, QQ(17) / 9)]


def theta_werte(pd):
    cs, phi = pd_to_cs_phi([tuple(c) for c in pd])
    Cs = matrix(ZZ, [list(r) for r in cs])
    ph = vector(ZZ, phi)
    out = []
    for a, b in PUNKTE:
        out.append(str(Theta_numeric([Cs, ph], a, b)[2]))
    return tuple(out)


def main():
    roh = json.load(open("mutanten_kandidaten.json"))
    gut = [x for x in roh if x.get("brauchbar")]
    print("Stage-1 pairs with verdict ECHTES PAAR: %d" % len(gut))

    # de-duplicate on the pair of isometry_signatures
    gesehen, eindeutig = set(), []
    for x in gut:
        schluessel = frozenset((x["K"]["isosig"], x["Kmut"]["isosig"]))
        if schluessel in gesehen:
            continue
        gesehen.add(schluessel)
        eindeutig.append(x)
    print("after de-duplication (isometry_signature pair): %d" % len(eindeutig))
    if not eindeutig:
        print("🔴 n=0 is no green."); sys.exit(1)

    ergebnisse, getrennt = [], 0
    for x in sorted(eindeutig, key=lambda y: y["K"]["crossings"]):
        nK, nM = x["K"]["name"], x["Kmut"]["name"]
        try:
            wK, wM = theta_werte(x["K"]["pd"]), theta_werte(x["Kmut"]["pd"])
        except Exception as e:
            print("   %-15s vs %-15s  🔴 THETA ERROR %s: %s"
                  % (nK, nM, type(e).__name__, str(e)[:60]))
            ergebnisse.append({"K": nK, "Kmut": nM, "fehler": str(e)[:200]})
            continue
        trennt = wK != wM
        getrennt += 1 if trennt else 0
        print("   %-15s vs %-15s c=%2d det=%-5d  Theta %s"
              % (nK, nM, x["K"]["crossings"], x["K"]["determinant"],
                 "SEPARATES" if trennt else "equal"))
        ergebnisse.append({"K": nK, "Kmut": nM, "c": x["K"]["crossings"],
                           "det": x["K"]["determinant"], "vol": x["K"]["volume"],
                           "trennt": bool(trennt),
                           "theta_K": list(wK), "theta_Kmut": list(wM)})

    gerechnet = [e for e in ergebnisse if "trennt" in e]
    print("\n== RESULT: theta separates %d of %d computed pairs ==" % (getrennt, len(gerechnet)))

    eich = [e for e in gerechnet if e["K"] == "pz4_2_3_5_7"]
    if eich:
        print("   counter-check pz4_2_3_5_7 (already separated in the corpus): %s"
              % ("SEPARATES ✓" if eich[0]["trennt"] else "🔴 EQUAL — the chain is wrong!"))
    json.dump(ergebnisse, open("ergebnis_mutanten_theta.json", "w"),
              ensure_ascii=False, indent=1)
    print("   written: ergebnis_mutanten_theta.json")


if __name__ == "__main__":
    main()
