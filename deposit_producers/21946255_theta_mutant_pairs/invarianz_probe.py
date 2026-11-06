#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""INVARIANCE PROBE: is the compared theta value a knot invariant
or does it hang on the diagram?

Without this probe every "separation" is worthless: different values could
simply come from different diagrams. Nine of the built knots occur
under TWO parameter sets of (cyclic rotation of the pretzel) — same
isometry_signature, i.e. the same knot, but a different diagram.
Theta MUST give the same value there.
"""
import site_config
import json, os, sys
from collections import defaultdict
sys.path.insert(0, (site_config.ROOT + '/freeze_manifest/pipeline_scripts'))
from pd_to_theta_input import pd_to_cs_phi
from sage.all_cmdline import *   # noqa
BN = (site_config.ROOT + '/external_tools/bar_natan_theta')
load(os.path.join(BN, "Theta_numeric.sage"))
PUNKTE = [(QQ(22)/7, QQ(21)/13), (QQ(5)/3, QQ(-7)/4)]

def werte(pd):
    cs, phi = pd_to_cs_phi([tuple(c) for c in pd])
    return tuple(str(Theta_numeric([matrix(ZZ, [list(r) for r in cs]),
                                    vector(ZZ, phi)], a, b)[2]) for a, b in PUNKTE)

roh = json.load(open("mutanten_kandidaten.json"))
seiten = []
for x in roh:
    if x.get("brauchbar"):
        seiten += [(x["K"]["isosig"], x["K"]["name"], x["K"]["pd"], x["K"]["crossings"]),
                   (x["Kmut"]["isosig"], x["Kmut"]["name"], x["Kmut"]["pd"], x["Kmut"]["crossings"])]
g = defaultdict(list)
for iso, name, pd, c in seiten:
    if not any(name == n for _, n, _, _ in g[iso]):
        g[iso].append((iso, name, pd, c))
mehrfach = {k: v for k, v in g.items() if len(v) > 1}
print("Knots that occur under several diagrams: %d" % len(mehrfach))
if not mehrfach:
    print("🔴 no probe possible — n=0 is no green."); sys.exit(1)
gut = schlecht = 0
for iso, eintraege in sorted(mehrfach.items()):
    w = {}
    for _, name, pd, c in eintraege:
        w.setdefault(werte(pd), []).append(name)
    ok = len(w) == 1
    gut += ok; schlecht += (not ok)
    print("   %-9s c=%2d  %-38s %s" % (iso[:9], eintraege[0][3],
          " = ".join(n for _, n, _, _ in eintraege), "EQUAL ✓" if ok else "🔴 DIFFERENT"))
    if not ok:
        for wert, namen in w.items():
            print("        %s -> %s" % (namen, wert[0][:40]))
print("\n== %d/%d knots give the SAME theta value under different diagrams ==" % (gut, gut + schlecht))
if schlecht:
    print("🔴 Theta depends on the diagram here — every separation statement is thereby void.")
    sys.exit(2)
print("✓ The compared value behaves like an invariant. The separations hold.")
