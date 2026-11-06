#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""
NEW UPPER BOUNDS, second run — fat budget AND an honest witness test.

Why a second run (11.08. 20:5x):
  Run 1 (17:46) reported 6 "u EXACT". The counter-check REFUTED 5 of them: ist_unknoten()
  accepted the unknot already at abs(det)==1, and det=1 does not mean unknot
  (11n34). Calibrated on the proven range: trefoil with one flipped crossing delivers
  Jones == 1 exactly ⟹ Jones != 1 REFUTES the unknot.

Three stages instead of one — sieve cheaply, judge expensively:
  1. SIEVE    simplify('global'); 0 crossings -> witness at once. Otherwise det != 1 -> next.
              (det is the cheapest filter and lets no real witness through.)
  2. CHECK Jones != 1  -> REFUTED, keep searching. (Costs only for the few
              candidates that survive stage 1.)
  3. PROOF    group of the complement = Z (1 generator, 0 relations), otherwise backtrack rounds
              down to 0 crossings. Only THAT is called "witness".
  If Jones==1 stays, group not reduced and backtrack unsuccessful: 'kandidat_unbewiesen'
  — is NOT counted as a witness, but written along, so that it is not lost.

What is carried over: the 37 'erschoepfend' verdicts from run 1 remain valid (the old
test was TOO GENEROUS, not too strict — it would have let every real witness through,
so the absence is still receipted). The 5 refuted knots come back into the
pot; their search had stopped at the false witness.

Usage (one worker of N):
  PATH=/home/squeezy/miniforge3/envs/sage_stable/bin:$PATH \
  /home/squeezy/miniforge3/envs/sage_stable/bin/python neue_obere_schranken2.py <w> <N> <budget_s>
Every worker writes atomically to lauf2/teil_<w>.json (footgun 28).
"""
import ast
import itertools
import json
import math
import os
import sys
import time

import spherogram
import snappy  # must be imported, otherwise spherogram does not know .exterior()

HIER = os.path.dirname(os.path.abspath(__file__))
LAUF = os.path.join(HIER, "lauf2")
ALT = os.path.join(HIER, "neue_obere_schranken.json")

BEWEIS_BUDGET = 90.0          # per candidate, for backtrack rounds
# sys.argv is read ONLY in main() — otherwise an import (pruefe_die_fuenf.py) swallows
# the caller's arguments and dies of the IndexError.


def kippe(pd, idx):
    return [((k[1], k[2], k[3], k[0]) if i in idx else tuple(k)) for i, k in enumerate(pd)]


def urteile_ueber(pd_gekippt):
    """
    Returns (art, info):
      'zeuge'                -> proven the unknot
      'widerlegt'            -> Jones != 1, certainly NOT the unknot
      'kandidat_unbewiesen'  -> Jones == 1, but no proof succeeded
    """
    L = spherogram.Link(pd_gekippt)
    L.simplify(mode="global")
    if len(L.crossings) == 0:
        return "zeuge", "0 Kreuzungen nach simplify"
    if abs(L.determinant()) != 1:
        return "widerlegt", "det != 1"

    J = L.jones_polynomial()
    if J != 1:
        return "widerlegt", f"Jones = {J}"

    G = L.exterior().fundamental_group()
    if G.num_generators() == 1 and G.num_relators() == 0:
        return "zeuge", "Gruppe des Aussenraums = Z"

    t0 = time.time()
    for s in (5, 10, 20, 40, 80, 160, 320):
        if time.time() - t0 > BEWEIS_BUDGET:
            break
        K = spherogram.Link(pd_gekippt)
        K.simplify(mode="global")
        K.backtrack(steps=s)
        K.simplify(mode="global")
        K.simplify(mode="pickup")
        K.simplify(mode="global")
        if len(K.crossings) == 0:
            return "zeuge", f"backtrack(steps={s}) auf 0 Kreuzungen"
    return "kandidat_unbewiesen", f"Jones=1, Gruppe ({G.num_generators()},{G.num_relators()}), backtrack erfolglos"


def main():
    W = int(sys.argv[1])
    N = int(sys.argv[2])
    BUDGET = float(sys.argv[3]) if len(sys.argv) > 3 else 7200.0

    os.makedirs(LAUF, exist_ok=True)
    ziel = os.path.join(LAUF, f"teil_{W}.json")

    alt = {x["stem"]: x for x in json.load(open(ALT))}
    # validly decided: only the exhaustive ones. The 6 find reports are to be re-evaluated,
    # 5 of them are refuted -> back into the pot. The one receipted stands in zeugen_urteil.json.
    fertig = {s for s, x in alt.items() if x["urteil"].startswith("erschoepfend")}
    if os.path.exists(os.path.join(HIER, "zeugen_urteil.json")):
        for u in json.load(open(os.path.join(HIER, "zeugen_urteil.json"))):
            if u["urteil"].startswith("✅"):
                fertig.add(u["stem"])
    # Restart over ALL partial files, not only the own: otherwise a worker, after
    # a restart with a different worker count, recomputes foreign results once more.
    for f in sorted(os.listdir(LAUF)) if os.path.isdir(LAUF) else []:
        if f.startswith("teil_") and f.endswith(".json"):
            fertig |= {x["stem"] for x in json.load(open(os.path.join(LAUF, f)))}

    b = json.load(open(os.path.join(HIER, "bilanz_daten.json")))["zeilen"]
    kand = [z for z in b if z["ist_vertreter"] and z["u_obere"] is None
            and z["u_lower"] is not None and 1 <= z["u_lower"] <= 4
            and z["stem"] not in fertig]

    # estimate cost: measured rate from run 1, otherwise conservatively 5/s
    def kosten(z):
        n = len(ast.literal_eval(json.load(open(z["pflicht_path"]))["pd_list_str"]))
        anzahl = math.comb(n, z["u_lower"])
        a = alt.get(z["stem"])
        rate = (a["geprueft"] / a["sekunden"]) if a and a.get("sekunden") else 5.0
        return anzahl / rate, anzahl, n

    mit = sorted(((kosten(z), z) for z in kand), key=lambda t: t[0][0])
    meine = [(k, z) for i, (k, z) in enumerate(mit) if i % N == W]
    print(f"[w{W}] {len(meine)} of {len(mit)} knots, budget {BUDGET:.0f}s per knot", flush=True)

    stand = json.load(open(ziel)) if os.path.exists(ziel) else []
    for (schaetz, anzahl, n), z in meine:
        pd = [tuple(x) for x in ast.literal_eval(json.load(open(z["pflicht_path"]))["pd_list_str"])]
        k = z["u_lower"]
        t1, geprueft, zeuge, kandidaten, widerlegt = time.time(), 0, None, [], 0
        for idx in itertools.combinations(range(n), k):
            geprueft += 1
            art, info = urteile_ueber(kippe(pd, set(idx)))
            if art == "zeuge":
                zeuge = (list(idx), info)
                break
            if art == "kandidat_unbewiesen":
                kandidaten.append([list(idx), info])
            else:
                widerlegt += 1
            if time.time() - t1 > BUDGET:
                break

        if zeuge:
            urteil = f"u = {k} EXAKT (Zeuge BEWIESEN: {zeuge[1]})"
        elif geprueft >= anzahl:
            urteil = f"erschoepfend: KEIN Zeuge bei k={k} im gespeicherten Diagramm"
        else:
            urteil = f"kein Zeuge in {BUDGET:.0f}s ({geprueft} von {anzahl}) — offen"
        e = dict(stem=z["stem"], c=z["c"], u_lower=k, n_pd=n, kombinationen=anzahl,
                 geprueft=geprueft, gefunden=zeuge[0] if zeuge else None,
                 beweis=zeuge[1] if zeuge else None, kandidaten_unbewiesen=kandidaten,
                 urteil=urteil, schaetzung_s=round(schaetz, 1),
                 sekunden=round(time.time() - t1, 1))
        stand.append(e)
        json.dump(stand, open(ziel + ".tmp", "w"), indent=1)
        os.replace(ziel + ".tmp", ziel)
        print(f"[w{W}] {urteil:56} c={z['c']:3} n={n:3} {z['stem'][:34]}", flush=True)

    print(f"[w{W}] done: {len(stand)} knots", flush=True)


if __name__ == "__main__":
    main()
