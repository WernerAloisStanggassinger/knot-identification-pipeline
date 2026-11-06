#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Braid-index repair (AUFTRAG_ZOPFINDEX_SUNDOWN_20260831, pattern: double-6 retraction 25.08.).

Disease: stage 5 copies the STRAND COUNT of the braid found by Sage into the scalar field
`braid_index` — an upper bound that appears as a value. 130 records thus carry the upper
end of an open MFW interval, 192 a scalar with no bounds at all.

Repair = make the provenance honest, at KNOT level:
  per distinct_knot_id:  lower bound = max(MFW-lower of all members)
                         upper bound = min(strand count reached by all members)
  (valid, because all members are the same knot; every braid word is a witness of an
   upper bound, every HOMFLY-MFW a witness of a lower bound.)
  lower == upper  ->  braid_index exact (status 'exact', witnesses in the block)
  otherwise       ->  braid_index = None, interval remains (status 'interval' resp.
                      'upper_bound_only' where no MFW exists)
  lower > upper   ->  CONTRADICTION: write nothing, report loudly (never smooth over).

Every changed record: backup `.bak_20260905_zopfindex` + block
`braid_index_reconciled_20260905` with the old values, reason and witnesses.

COUNTER-CASE (mandatory, otherwise exit 3): sundown must become exactly 7 (witness trainstation),
trainstation must stay 7, and at least 300 records must change
(dry-run measurement: 130+192+1). Call without flag = dry run; --schreiben writes.
"""
import site_config
import json, os, shutil, sys
from collections import defaultdict

ROOT = site_config.ROOT
DB = site_config.PIPELINE_DIR + "/corpus_db.json"
DATUM = "20260905"


def main(schreiben=False):
    db = json.load(open(DB))
    recs = db["knots"]
    info = {}   # stem -> (path, dia, did)
    for r in recs:
        pf = r.get("pflicht_path")
        if not pf or not os.path.exists(pf):
            continue
        d = json.load(open(pf))
        dia = d.get("pflichtliste", {}).get("DIAGRAMM", {})
        info[r["stem"]] = (pf, d, dia, r.get("distinct_knot_id"))

    # Knot level: collect bounds + witnesses
    gruppen = defaultdict(list)
    for stem, (pf, d, dia, did) in info.items():
        lo = dia.get("braid_index_lower_bound")
        up = dia.get("braid_index_upper_bound")
        if up is None:
            up = dia.get("braid_n_strands")
        if lo is None and up is None and dia.get("braid_index") is None:
            continue
        gruppen[did].append((stem, lo, up))

    knoten = {}
    widersprueche = []
    for did, mem in gruppen.items():
        los = [(l, s) for s, l, u in mem if l is not None]
        ups = [(u, s) for s, l, u in mem if u is not None]
        klo = max(los) if los else (None, None)
        kup = min(ups) if ups else (None, None)
        if klo[0] is not None and kup[0] is not None and klo[0] > kup[0]:
            widersprueche.append((did, klo, kup, mem))
            continue
        knoten[did] = {"lower": klo[0], "lower_zeuge": klo[1],
                       "upper": kup[0], "upper_zeuge": kup[1]}

    geaendert, exakt_neu, intervall, nur_ober = [], 0, 0, 0
    for stem, (pf, d, dia, did) in info.items():
        if did not in knoten:
            continue
        k = knoten[did]
        neu_scalar = k["lower"] if (k["lower"] is not None and k["lower"] == k["upper"]) else None
        if neu_scalar is not None:
            status = "exact"
        elif k["lower"] is not None:
            status = "interval"
        else:
            status = "upper_bound_only"
        alt = {f: dia.get(f) for f in ("braid_index", "braid_index_interval",
                                       "braid_index_lower_bound", "braid_index_upper_bound",
                                       "braid_index_status")}
        neu = {"braid_index": neu_scalar,
               "braid_index_interval": [k["lower"], k["upper"]],
               "braid_index_lower_bound": k["lower"],
               "braid_index_upper_bound": k["upper"],
               "braid_index_status": status}
        wertfelder = ("braid_index", "braid_index_interval",
                      "braid_index_lower_bound", "braid_index_upper_bound")
        if all(alt.get(f) == neu[f] for f in wertfelder):
            continue  # values unchanged — no intervention, no status field either
        geaendert.append(stem)
        if status == "exact": exakt_neu += 1
        elif status == "interval": intervall += 1
        else: nur_ober += 1
        if schreiben:
            shutil.copy2(pf, pf + f".bak_{DATUM}_zopfindex")
            dia.update(neu)
            dia[f"braid_index_reconciled_{DATUM}"] = {
                "grund": ("Scalar was the strand count of the braid found = only an "
                          "upper bound (AUFTRAG_ZOPFINDEX_SUNDOWN_20260831). Reconciled at "
                          "knot level over all realisations of the same distinct_knot_id."),
                "alt": alt,
                "zeuge_untere_schranke": k["lower_zeuge"],
                "zeuge_obere_schranke": k["upper_zeuge"],
            }
            with open(pf, "w") as f:
                json.dump(d, f, indent=2, default=str)

    print(f"{'WRITING' if schreiben else 'DRY-RUN'}: {len(geaendert)} records to change "
          f"(of which newly exact: {exakt_neu} · interval: {intervall} · upper bound only: {nur_ober})")
    print(f"Knot groups: {len(gruppen)} · contradictions (NOT touched): {len(widersprueche)}")
    for did, klo, kup, mem in widersprueche:
        print(f"   🔴 {did}: lower {klo} > upper {kup} — {mem}")

    # COUNTER-CASE
    sd = knoten.get("sundown")
    ok = (sd and sd["lower"] == sd["upper"] == 7
          and "sundown" in geaendert and "trainstation" not in geaendert
          and len(geaendert) >= 300)
    print(f"Counter-case: sundown->7 exact: {bool(sd and sd['lower']==sd['upper']==7)} · "
          f"sundown changed: {'sundown' in geaendert} · trainstation untouched: "
          f"{'trainstation' not in geaendert} · >=300 changes: {len(geaendert)>=300}")
    if not ok:
        print("🔴 COUNTER-CASE NOT MET — nothing counts.")
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main(schreiben="--schreiben" in sys.argv))
