#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""
PROOF of the upper bounds ON THE STORED DIAGRAM.

Occasion (11.08.2026): the balance reports 63 distinct knots with u(K) EXACTLY determined
(lower == upper bound). The audit shows: **29 of them** rest on a witness
that was computed on a DIFFERENT diagram than the one that stands in the record.
For all 29 the search diagram is LARGER — the benign direction, because a witness on
any diagram of the same knot is a valid upper bound. **But it is not verifiable from
the record.** Exactly that is healed here:

  For each of these knots a search is run exhaustively over all subsets of size k = u_lower
  in the STORED PD. If one is found, the exact determination of u(K) is traceable
  **from the record alone**.

If none is found, that is NO counter-proof — the other representation can very well
realize the bound. It only means: not reproducible from the record alone.
This separation is the whole point.

Unknot test: simplify -> 0 crossings OR |det| == 1 (det is an invariant,
|det| != 1 safely excludes the unknot). NO except around the computation
(Disc #613 / the blunder of 14:3x) — must run in sage_stable.

Usage:
  PATH=/home/squeezy/miniforge3/envs/sage_stable/bin:$PATH \
  /home/squeezy/miniforge3/envs/sage_stable/bin/python nachweis_obere_schranken.py
"""
import site_config
import ast
import itertools
import json
import math
import os
import time

import spherogram

HIER = os.path.dirname(os.path.abspath(__file__))
DB = (site_config.ROOT + '/freeze_manifest/pipeline_scripts/corpus_db.json')
ZIEL = os.path.join(HIER, "nachweis_obere_schranken.json")
DECKEL = 400_000          # combinations per knot; above that the run is honestly aborted


def kippe(pd, idx):
    return [((k[1], k[2], k[3], k[0]) if i in idx else tuple(k)) for i, k in enumerate(pd)]


def ist_unknoten(pd):
    L = spherogram.Link(pd)
    L.simplify(mode="global")
    if len(L.crossings) == 0:
        return True
    return abs(L.determinant()) == 1          # no except: wrong env should crash


def main():
    bil = {z["stem"]: z for z in json.load(open(os.path.join(HIER, "bilanz_daten.json")))["zeilen"]}
    audit = {x["stem"]: x for x in json.load(open(os.path.join(HIER, "audit_unknotting_diagramm.json")))}
    pfade = {r["stem"]: r["pflicht_path"] for r in json.load(open(DB))["knots"]}

    ziel = [z for z in bil.values()
            if z["ist_vertreter"] and z["u_lower"] is not None and z["u_obere"] is not None
            and z["u_lower"] == z["u_obere"]
            and audit.get(z["stem"], {}).get("urteil") != "OK"]
    ziel.sort(key=lambda z: (z["u_lower"], z["c"] or 0))
    print(f"{len(ziel)} exact determinations whose witness arose on a different diagram", flush=True)

    aus, t0 = [], time.time()
    for z in ziel:
        pd = [tuple(x) for x in ast.literal_eval(json.load(open(pfade[z["stem"]]))["pd_list_str"])]
        n, k = len(pd), z["u_lower"]
        anzahl = math.comb(n, k)
        e = dict(stem=z["stem"], c=z["c"], u=k, n_pd=n, kombinationen=anzahl,
                 gefunden=None, urteil=None, sekunden=None)
        t1 = time.time()
        if anzahl > DECKEL:
            e["urteil"] = "NICHT-GEPRUEFT (zu viele Kombinationen)"
        else:
            for idx in itertools.combinations(range(n), k):
                if ist_unknoten(kippe(pd, set(idx))):
                    e["gefunden"] = list(idx)
                    break
            e["urteil"] = "AM RECORD NACHWEISBAR" if e["gefunden"] else "aus dem Record nicht reproduzierbar"
        e["sekunden"] = round(time.time() - t1, 1)
        aus.append(e)
        print(f"  {e['urteil']:38} u={k} c={z['c']:3} n={n:3} ({anzahl} comb., {e['sekunden']}s)  "
              f"{z['stem'][:34]}  {e['gefunden']}", flush=True)
        json.dump(aus, open(ZIEL + ".tmp", "w"), indent=1)
        os.replace(ZIEL + ".tmp", ZIEL)

    n_ok = sum(1 for x in aus if x["urteil"] == "AM RECORD NACHWEISBAR")
    n_nein = sum(1 for x in aus if x["urteil"].startswith("aus dem Record"))
    n_gross = sum(1 for x in aus if x["urteil"].startswith("NICHT-GEPRUEFT"))
    print(f"\nRESULT after {time.time()-t0:.0f}s: {n_ok} verifiable on the record, "
          f"{n_nein} not reproducible, {n_gross} too large")
    print("⚠️ 'not reproducible' is NO counter-proof — the other representation can")
    print("   realize the bound. It only means: not traceable from the record alone.")


if __name__ == "__main__":
    main()
