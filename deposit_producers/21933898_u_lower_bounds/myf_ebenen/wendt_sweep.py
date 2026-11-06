#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""
SWEEP: Wendt bound u(K) >= mu(H_1(Sigma_2(K))) over all 788 corpus records.

Method + counter-checks see wendt_probe.py (run small, 13/13 det check green).
Here only the full run, incremental and atomic (footgun 28: save after every block,
tmp + os.replace, never directly into the target). Restart recognizes finished work BY THE ARTIFACT,
not at the log (footgun 25) — and counts only DECIDED ones as done.

Usage:  /home/squeezy/miniforge3/envs/gnottheorie_env/bin/python wendt_sweep.py [--neu]
"""
import site_config
import ast
import json
import os
import sys
import time

from sympy import factorint

import spherogram

WURZEL = site_config.ROOT
DB = WURZEL + "/freeze_manifest/pipeline_scripts/corpus_db.json"
HIER = os.path.dirname(os.path.abspath(__file__))
ZIEL = os.path.join(HIER, "wendt_sweep.json")

sys.path.insert(0, HIER)
from wendt_probe import det_ganz, mu_von  # noqa: E402  (the same computation, no second code)


def sichern(stand):
    json.dump(stand, open(ZIEL + ".tmp", "w"), indent=1)
    os.replace(ZIEL + ".tmp", ZIEL)


def main():
    neu = "--neu" in sys.argv
    db = json.load(open(DB))
    stand = {}
    if os.path.exists(ZIEL) and not neu:
        alt = json.load(open(ZIEL))
        # only DECIDED ones count as done (footgun 25)
        stand = {k: v for k, v in alt.get("zeilen", {}).items() if v.get("mu") is not None}
        print(f"↻ restart: {len(stand)} count as decided, "
              f"{len(alt.get('zeilen', {})) - len(stand)} as UNDECIDED")

    # Cheap ones first (ascending by c): then the bulk is in after minutes and only
    # the expensive tail runs on afterwards. The box hits the heat ceiling (theta runs 24 workers)
    # — more processes would bring no throughput, a better order would.
    offen = sorted([r for r in db["knots"] if r["stem"] not in stand],
                   key=lambda r: (r.get("c") or 0))
    print(f"open: {len(offen)} of {len(db['knots'])} "
          f"(c from {offen[0].get('c') if offen else '-'} to {offen[-1].get('c') if offen else '-'})",
          flush=True)
    t0 = time.time()
    for i, r in enumerate(offen, 1):
        stem = r["stem"]
        eintrag = dict(stem=stem, c=r.get("c"), det_record=r.get("det"),
                       category=r.get("category"), folder=r.get("folder"),
                       distinct_knot_id=r.get("distinct_knot_id"),
                       ist_vertreter=bool(r.get("is_distinct_representative")),
                       knotinfo_id=r.get("knotinfo_id"), sig=r.get("sig"),
                       mu=None, fehler=None)
        try:
            d = json.load(open(r["pflicht_path"]))
            pd = d.get("pd_list_str")
            if not pd:
                raise ValueError("no pd_list_str")
            L = spherogram.Link([tuple(x) for x in ast.literal_eval(pd)])
            V = [[int(x) for x in zeile] for zeile in L.seifert_matrix()]
            n = len(V)
            M = [[V[a][b] + V[b][a] for b in range(n)] for a in range(n)]
            dM = det_ganz(M)
            mu, pro_p = mu_von(M, dM)
            eintrag.update(n=n, det_VVT=abs(dM), mu=mu,
                           pro_primzahl={str(k): v for k, v in pro_p.items()},
                           det_stimmt=(abs(dM) == abs(int(r["det"]))) if r.get("det") is not None else None)
        except Exception as ex:
            eintrag["fehler"] = f"{type(ex).__name__}: {str(ex)[:120]}"
        stand[stem] = eintrag
        if i % 10 == 0 or i == len(offen):
            sichern({"gebaut_am": "2026-08-11", "n": len(stand), "zeilen": stand})
            rate = i / max(1e-9, time.time() - t0)
            print(f"[{i}/{len(offen)}] {rate:.1f}/s  last: {stem[:34]} "
                  f"mu={eintrag['mu']} error={eintrag['fehler']}", flush=True)

    sichern({"gebaut_am": "2026-08-11", "n": len(stand), "zeilen": stand})
    ok = [v for v in stand.values() if v.get("mu") is not None]
    fehler = [v for v in stand.values() if v.get("fehler")]
    dstimmt = [v for v in ok if v.get("det_stimmt") is True]
    dfalsch = [v for v in ok if v.get("det_stimmt") is False]
    print(f"\ndone: {len(ok)} computed, {len(fehler)} errors, {time.time()-t0:.0f}s")
    print(f"G1  |det(V+V^T)| == det(Record): {len(dstimmt)} ok, {len(dfalsch)} RED")
    for v in dfalsch[:10]:
        print(f"    RED {v['stem']}: {v.get('det_VVT')} != {v.get('det_record')}")
    from collections import Counter
    print("error kinds:", Counter(v["fehler"].split(":")[0] for v in fehler))


if __name__ == "__main__":
    main()
