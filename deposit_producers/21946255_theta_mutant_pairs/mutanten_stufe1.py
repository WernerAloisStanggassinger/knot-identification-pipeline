#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""STAGE 1 — build and sieve mutant pairs. Runs in `sage_stable`.

⚠️ ENV: NOT gnottheorie_env. There snappy throws `SageNotAvailable` for
determinant/alexander/jones — the polynomials need SnapPy INSIDE Sage.
In sage_stable both are present (snappy 3.3.2 + spherogram + Sage 10.7).

BLUEPRINT — and the calibration that saved it:
The first draft took `RationalTangle(p)` (horizontal twist). The sum of those
is ONE long twist band, the closure therefore the torus knot T(2, sum p) —
P(2,3,5,7) came out with **det 17** instead of 247. No error was visible, the
knot even had the right crossing number 17.
Correct is the VERTICAL twist, i.e. the fraction 1/p: `RationalTangle(1, p)`.
Calibrated on the two pairs that already lie in the corpus:
    P(2,3,5,7) -> det 247 ✓     P(2,3,5,3) -> det 123 ✓
Both bit-exact against `corpus_db`. Without this calibration the whole series would have been
silently wrong. [[project_calibrate_convention_against_proven_range]]

SIEVE, in this order — every stage may reject:
 1. both hyperbolic with a geometric solution?     otherwise: no test case
 2. isometry_signature DIFFERENT?                     otherwise: the same knot
 3. Alexander + Jones + det EQUAL?                    otherwise: already separated classically
Only what passes all three is a theta test case.
"""
import itertools
import json
import sys

import snappy          # 🔴 MANDATORY: without this import spherogram throws a
                      # RuntimeError in exterior() ("SnapPy doesn't seem to be
                      # available"). On 14.08. my own except block turned that into
                      # 31x "no geometric solution" — an error disguised as a
                      # finding. The tell was the clean 0 of 31.
from spherogram import RationalTangle

ZIEL = "mutanten_kandidaten.json"
EICHUNG = {(2, 3, 5, 7): 247, (2, 3, 5, 3): 123}


def pretzel(ps):
    t = RationalTangle(1, ps[0])
    for p in ps[1:]:
        t = t + RationalTangle(1, p)
    L = t.numerator_closure()
    L.simplify("global")
    return L


def profil(L, name):
    d = {"name": name, "crossings": len(L.crossings),
         "komponenten": len(L.link_components)}
    try:
        d["pd"] = [list(map(int, c)) for c in L.PD_code(min_strand_index=1)]
    except Exception as e:
        d["pd"] = None
        d["pd_fehler"] = "%s: %s" % (type(e).__name__, e)
    for feld, ruf in (("determinant", lambda: int(abs(L.determinant()))),
                      ("alexander", lambda: str(L.alexander_polynomial())),
                      ("jones", lambda: str(L.jones_polynomial()))):
        try:
            d[feld] = ruf()
        except Exception as e:
            d[feld] = None
            d[feld + "_fehler"] = "%s: %s" % (type(e).__name__, e)
    try:
        M = L.exterior()
        M.simplify()
        # A non-geometric solution does NOT mean "not hyperbolic" — it means
        # "this triangulation is no good". Re-roll until it is geometric.
        versuche = 0
        while not str(M.solution_type()).startswith("all tetrahedra") and versuche < 12:
            M.randomize()
            versuche += 1
        d["solution_type"] = str(M.solution_type())
        d["randomize_versuche"] = versuche
        d["volume"] = float(M.volume())
        d["isosig"] = (M.isometry_signature()
                       if d["solution_type"].startswith("all tetrahedra") else None)
    except Exception as e:
        d.update({"volume": None, "isosig": None,
                  "geo_fehler": "%s: %s" % (type(e).__name__, e)})
    return d


def main():
    print("== CALIBRATION (without it the whole series is worthless) ==")
    for ps, soll in EICHUNG.items():
        L = pretzel(ps)
        got = int(abs(L.determinant()))
        print("   P%s  c=%d  det=%d  (corpus_db: %d)  %s"
              % (ps, len(L.crossings), got, soll, "✓" if got == soll else "🔴"))
        if got != soll:
            print("   🔴 Calibration failed — abort, nothing is measured here.")
            sys.exit(2)

    kand = []
    for b, c, e in itertools.permutations([3, 5, 7, 9, 11], 3):
        if b > c:                      # (b,c) and (c,b) would be the same pair
            continue
        kand.append((2, b, c, e))
    for k in EICHUNG:
        if k not in kand:
            kand.append(k)

    print("\n== %d quadruples P(2,b,c,d), each against P(2,c,b,d) ==" % len(kand))
    erg = []
    for p in kand:
        a, b, c, e = p
        q = (a, c, b, e)
        nK, nM = "pz4_%d_%d_%d_%d" % p, "pz4_%d_%d_%d_%d" % q
        A, B = profil(pretzel(p), nK), profil(pretzel(q), nM)

        if A["komponenten"] != 1 or B["komponenten"] != 1:
            urteil = "kein Knoten (Verschlingung) — verworfen"
        elif not (A["isosig"] and B["isosig"]):
            urteil = "keine geometrische Loesung — verworfen"
        elif A["isosig"] == B["isosig"]:
            urteil = "SELBER Knoten (isosig gleich) — kein Prueffall"
        elif not (A["alexander"] == B["alexander"] and A["jones"] == B["jones"]
                  and A["determinant"] == B["determinant"]):
            urteil = "klassisch schon getrennt — als Theta-Prueffall wertlos"
        else:
            urteil = "ECHTES PAAR"
        print("   %-15s vs %-15s c=%2d det=%-6s vol=%-10s %s"
              % (nK, nM, A["crossings"], A["determinant"],
                 ("%.5f" % A["volume"]) if A["volume"] else "-", urteil))
        erg.append({"params": list(p), "params_mut": list(q), "urteil": urteil,
                    "brauchbar": urteil == "ECHTES PAAR", "K": A, "Kmut": B})

    fehler = [x for x in erg if x["K"].get("geo_fehler")]
    if len(fehler) == len(erg) and erg:
        print("\n   🔴 ALL %d pairs fail at the same place: %s"
              % (len(erg), erg[0]["K"]["geo_fehler"]))
        print("      This is a TOOL error, not a finding about the knots. Abort.")
        sys.exit(3)

    gut = [x for x in erg if x["brauchbar"]]
    print("\\n== %d of %d are real test cases ==" % (len(gut), len(erg)))
    for x in gut:
        print("   P%s vs P%s  c=%d det=%d vol=%.5f"
              % (tuple(x["params"]), tuple(x["params_mut"]),
                 x["K"]["crossings"], x["K"]["determinant"], x["K"]["volume"]))
    json.dump(erg, open(ZIEL, "w"), ensure_ascii=False, indent=1)
    print("\n   written: %s" % ZIEL)
    if not gut:
        print("   🔴 n=0 is no green — then this family carries no test cases.")
        sys.exit(1)


if __name__ == "__main__":
    main()
