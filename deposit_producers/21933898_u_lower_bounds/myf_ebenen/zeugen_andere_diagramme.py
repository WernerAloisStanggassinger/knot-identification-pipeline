#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""
WITNESS SEARCH IN OTHER DIAGRAMS OF THE SAME KNOT (lauf3).

Why (11.08. 23:5x, after the probe `probe_vereinfacht.py`):
  lauf1+lauf2 reported "exhaustive: NO witness at k=u_lower" for ~70 knots.
  This verdict is honest, but it is a statement about THE STORED DIAGRAM,
  not about the knot: u(K) is the minimum over ALL diagrams. Another diagram
  of the same knot may have a witness at the same k.
  The probe showed that mere simplifying changes hardly anything (median x1.2 fewer
  combinations, same crossing number) — so here the diagram SPACE is searched:
  backtrack(steps=s) generates, via Reidemeister moves, another diagram of the SAME knot.

What a find is worth, and what not:
  FIND      witness at k=u_lower  ==>  u <= k, and with the proven lower bound k
            stands  u = k  EXACT. A find is a PROOF, no matter in which diagram it falls.
  NO FIND still says NOTHING. That is why completeness is not computed here
            (it would be paid per diagram and worthless per diagram), but scattered:
            many diagrams, one sample of k-subsets each. Chance instead of exhaustiveness
            is the right currency here, because only the positive case proves anything.

Verdict on a candidate: taken over from neue_obere_schranken2.urteile_ueber
  (sieve det -> check Jones != 1 refuted -> proof group=Z or backtrack to 0).
  Exactly this test overturned 5 of 6 false witnesses on 11.08.; it stays unchanged.

Counter-check against Frankenstein: every generated diagram must have the same determinant as
  the stored one. If it deviates, the diagram is discarded and counted (must never
  occur — backtrack preserves the knot type).

Receipt: on a find, the GENERATING diagram (PD) and the flipped indices are
  written along, so the witness can be recomputed independently.

Usage:
  <w> <N> <budget_s>   one worker of N, seconds per knot
  --probe <stem>       positive control: a single knot, logging aloud
"""
import ast
import itertools
import json
import math
import os
import random
import sys
import time

import spherogram
import snappy  # noqa: F401

from neue_obere_schranken2 import kippe, urteile_ueber  # reads sys.argv only in main()

HIER = os.path.dirname(os.path.abspath(__file__))
LAUF = os.path.join(HIER, "lauf3")
PROBEN_JE_DIAGRAMM = 400      # sample of k-subsets per generated diagram
BACKTRACK_SCHRITTE = (20, 40, 80, 160, 320)


def lade_pd(z):
    return [tuple(x) for x in ast.literal_eval(json.load(open(z["pflicht_path"]))["pd_list_str"])]


def anderes_diagramm(pd, schritte, rng):
    """Another diagram of the SAME knot: Reidemeister moves, then tidy up."""
    L = spherogram.Link([tuple(x) for x in pd])
    L.backtrack(steps=schritte)
    L.simplify(mode="global")
    L.simplify(mode="pickup")
    L.simplify(mode="global")
    return L


def teilmengen(n, k, wieviel, rng):
    """how many distinct k-subsets out of n — full count if cheap, otherwise scattered."""
    gesamt = math.comb(n, k)
    if gesamt <= wieviel:
        alle = list(itertools.combinations(range(n), k))
        rng.shuffle(alle)
        return alle
    gesehen, raus = set(), []
    while len(raus) < wieviel:
        s = tuple(sorted(rng.sample(range(n), k)))
        if s not in gesehen:
            gesehen.add(s)
            raus.append(s)
    return raus


def offene_knoten():
    """All representatives whose witness search in the stored diagram came up empty."""
    b = {z["stem"]: z for z in json.load(open(os.path.join(HIER, "bilanz_daten.json")))["zeilen"]}
    offen = {}
    quellen = [os.path.join(HIER, "neue_obere_schranken.json")]
    for lauf in ("lauf2", "lauf3"):
        d = os.path.join(HIER, lauf)
        if os.path.isdir(d):
            quellen += [os.path.join(d, f) for f in sorted(os.listdir(d))
                        if f.startswith("teil_") and f.endswith(".json")]
    erledigt = set()
    for q in quellen:
        for x in json.load(open(q)):
            if "EXAKT" in x["urteil"]:
                erledigt.add(x["stem"])          # already proven — not once more
            elif x["urteil"].startswith("erschoepfend") or "offen" in x["urteil"]:
                offen[x["stem"]] = x
    for s in erledigt:
        offen.pop(s, None)
    return [(s, b[s], x) for s, x in offen.items() if s in b and b[s]["u_lower"]]


def jage(stem, z, budget, laut=False):
    pd = lade_pd(z)
    k = z["u_lower"]
    det_soll = int(abs(spherogram.Link([tuple(x) for x in pd]).determinant()))
    rng = random.Random(hash(stem) & 0xFFFFFFFF)   # reproducible per knot
    t0 = time.time()
    diagramme, geprueft, verworfen, kandidaten = 0, 0, 0, []

    while time.time() - t0 < budget:
        s = rng.choice(BACKTRACK_SCHRITTE)
        L = anderes_diagramm(pd, s, rng)
        n = len(L.crossings)
        if n < k:
            continue
        if int(abs(L.determinant())) != det_soll:
            verworfen += 1
            continue
        pd_alt = [tuple(x) for x in L.PD_code()]
        diagramme += 1
        for idx in teilmengen(n, k, PROBEN_JE_DIAGRAMM, rng):
            geprueft += 1
            art, info = urteile_ueber(kippe(pd_alt, set(idx)))
            if art == "zeuge":
                return dict(stem=stem, c=z["c"], k=k, urteil=f"u = {k} EXAKT (Zeuge BEWIESEN: {info})",
                            gefunden=list(idx), beweis=info, backtrack_schritte=s,
                            n_diagramm=n, diagramm_pd=[list(x) for x in pd_alt],
                            diagramme=diagramme, geprueft=geprueft, verworfen=verworfen,
                            det=det_soll, sekunden=round(time.time() - t0, 1))
            if art == "kandidat_unbewiesen":
                kandidaten.append([list(idx), info, [list(x) for x in pd_alt]])
            if time.time() - t0 > budget:
                break
        if laut:
            print(f"  [{stem[:28]}] {diagramme} diagrams, {geprueft} samples, "
                  f"{time.time()-t0:.0f}s", flush=True)

    return dict(stem=stem, c=z["c"], k=k,
                urteil=f"kein Zeuge in {diagramme} anderen Diagrammen ({geprueft} Proben) — offen",
                gefunden=None, beweis=None, n_diagramm=None, diagramme=diagramme,
                geprueft=geprueft, verworfen=verworfen, kandidaten_unbewiesen=kandidaten[:20],
                det=det_soll, sekunden=round(time.time() - t0, 1))


def main():
    os.makedirs(LAUF, exist_ok=True)

    if sys.argv[1] == "--probe":
        stem = sys.argv[2]
        budget = float(sys.argv[3]) if len(sys.argv) > 3 else 600.0
        b = {z["stem"]: z for z in json.load(open(os.path.join(HIER, "bilanz_daten.json")))["zeilen"]}
        z = b[stem]
        print(f"PROBE {stem}  c={z['c']}  k={z['u_lower']}  Budget {budget:.0f}s", flush=True)
        e = jage(stem, z, budget, laut=True)
        print(json.dumps({k: v for k, v in e.items() if k != "diagramm_pd"}, indent=1), flush=True)
        ziel = os.path.join(LAUF, f"probe_{stem[:40]}.json")
        json.dump(e, open(ziel + ".tmp", "w"), indent=1)
        os.replace(ziel + ".tmp", ziel)
        return

    W, N = int(sys.argv[1]), int(sys.argv[2])
    budget = float(sys.argv[3]) if len(sys.argv) > 3 else 1800.0
    # Second layer (when run2 frees the cores): unroll the same heap from the back.
    # No re-sharding — that would tear apart the allocation of the workers already running.
    rueckwaerts = len(sys.argv) > 4 and sys.argv[4].startswith("rueck")
    ziel = os.path.join(LAUF, f"teil_{'r' if rueckwaerts else ''}{W}.json")
    stand = json.load(open(ziel)) if os.path.exists(ziel) else []

    # cheap first: small diagrams yield more diagrams per second
    pool = sorted(offene_knoten(), key=lambda t: t[1]["c"], reverse=rueckwaerts)
    # skip what ANY worker has already filed (not only its own file) —
    # otherwise the second layer recomputes the first.
    schon = {x["stem"] for x in stand}
    for f in sorted(os.listdir(LAUF)):
        if f.startswith("teil_") and f.endswith(".json"):
            schon |= {x["stem"] for x in json.load(open(os.path.join(LAUF, f)))}
    meine = [(s, z, x) for i, (s, z, x) in enumerate(pool) if i % N == W and s not in schon]
    print(f"[w{W}] {len(meine)} of {len(pool)} open knots, budget {budget:.0f}s per knot",
          flush=True)

    for stem, z, _alt in meine:
        e = jage(stem, z, budget)
        stand.append(e)
        json.dump(stand, open(ziel + ".tmp", "w"), indent=1)
        os.replace(ziel + ".tmp", ziel)
        print(f"[w{W}] {e['urteil'][:64]:64} c={e['c']:3} {stem[:32]}", flush=True)

    print(f"[w{W}] done: {len(stand)} knots", flush=True)


if __name__ == "__main__":
    main()
