#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""
COUNTER-CHECK with a SECOND, INDEPENDENT construction.

wendt_sweep.py computes mu(H_1(Sigma_2)) via the SEIFERT surface:  M = V + V^T.
Here the same mu is computed via the GOERITZ matrix (checkerboard coloring,
Gordon-Litherland way) — a different path from the same diagram to the same group.
Additionally the invariant divisors come exactly from the Smith normal form (Sage),
instead of from my rank-mod-p computation.

If both agree, the checker is not measuring itself (footgun 15).

Run (Sage env, PATH preamble is mandatory — TOOL_ENV_MAP):
  PATH=/home/squeezy/miniforge3/envs/sage_stable/bin:$PATH \
  /home/squeezy/miniforge3/envs/sage_stable/bin/python gegenprobe_goeritz.py [c_max]
"""
import site_config
import ast
import json
import os
import sys
import time

import spherogram
from sage.all import ZZ, matrix

HIER = os.path.dirname(os.path.abspath(__file__))
DB = (site_config.ROOT + '/freeze_manifest/pipeline_scripts/corpus_db.json')
ZIEL = os.path.join(HIER, "gegenprobe_goeritz.json")


def mu_smith(G):
    """Invariant divisors from the Smith normal form; mu = number of divisors != 1."""
    S = G.smith_form()[0]
    d = [ZZ(S[i, i]) for i in range(min(S.nrows(), S.ncols()))]
    nicht_eins = [abs(x) for x in d if abs(x) not in (0, 1)]
    nullen = sum(1 for x in d if x == 0)
    return len(nicht_eins) + nullen, [int(x) for x in nicht_eins], nullen


def main():
    c_max = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    db = json.load(open(DB))
    wendt = json.load(open(os.path.join(HIER, "wendt_sweep.json")))["zeilen"]

    aus, t0 = {}, time.time()
    kandidaten = [r for r in db["knots"]
                  if (r.get("c") or 999) <= c_max and r["stem"] in wendt
                  and wendt[r["stem"]].get("mu") is not None]
    print(f"Counter-check over {len(kandidaten)} records with c <= {c_max}", flush=True)

    for i, r in enumerate(kandidaten, 1):
        stem = r["stem"]
        e = dict(stem=stem, c=r.get("c"), det_record=r.get("det"),
                 mu_seifert=wendt[stem]["mu"], mu_goeritz=None, fehler=None)
        try:
            d = json.load(open(r["pflicht_path"]))
            L = spherogram.Link([tuple(x) for x in ast.literal_eval(d["pd_list_str"])])
            G = L.goeritz_matrix()
            G = matrix(ZZ, G)
            mu, teiler, nullen = mu_smith(G)
            e.update(mu_goeritz=mu, invariantenteiler=teiler, freier_rang=nullen,
                     det_goeritz=int(abs(G.det())),
                     det_stimmt=(int(abs(G.det())) == abs(int(r["det"]))))
        except Exception as ex:
            e["fehler"] = f"{type(ex).__name__}: {str(ex)[:100]}"
        aus[stem] = e
        if i % 25 == 0:
            json.dump(aus, open(ZIEL + ".tmp", "w"), indent=1)
            os.replace(ZIEL + ".tmp", ZIEL)
            print(f"[{i}/{len(kandidaten)}] {i/(time.time()-t0):.1f}/s", flush=True)

    json.dump(aus, open(ZIEL + ".tmp", "w"), indent=1)
    os.replace(ZIEL + ".tmp", ZIEL)

    ok = [v for v in aus.values() if v["mu_goeritz"] is not None]
    gleich = [v for v in ok if v["mu_goeritz"] == v["mu_seifert"]]
    ungleich = [v for v in ok if v["mu_goeritz"] != v["mu_seifert"]]
    dfalsch = [v for v in ok if v.get("det_stimmt") is False]
    fehler = [v for v in aus.values() if v["fehler"]]
    print()
    print(f"COUNTER-CHECK  mu(Goeritz+Smith) == mu(V+V^T, rank mod p) : "
          f"{len(gleich)}/{len(ok)} equal, {len(ungleich)} DEVIATION")
    for v in ungleich[:15]:
        print(f"    RED {v['stem']} c={v['c']}: Goeritz {v['mu_goeritz']} != Seifert {v['mu_seifert']} "
              f"(divisors {v.get('invariantenteiler')})")
    print(f"det(Goeritz) == det(Record): {len(ok)-len(dfalsch)}/{len(ok)} ok, {len(dfalsch)} RED")
    for v in dfalsch[:10]:
        print(f"    RED {v['stem']}: {v['det_goeritz']} != {v['det_record']}")
    print(f"not computed: {len(fehler)}")
    if fehler:
        from collections import Counter
        print("   ", dict(Counter(v["fehler"].split(":")[0] for v in fehler)))


if __name__ == "__main__":
    main()
