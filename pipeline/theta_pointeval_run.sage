# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# theta_pointeval_run.sage — evaluates all literals from _pointeval_literals.json at
# (t1,t2)=(22/7,21/13), ONE session. Writes theta_pointeval_RESULTS.json.
import json, time

import site_config
ROOT = site_config.ROOT
load(ROOT + "/external_tools/bar_natan_theta/Theta_numeric.sage")
LITS = ROOT + "/claude_setup/tmp_durable/_pointeval_literals.json"
OUT  = ROOT + "/claude_setup/tmp_durable/theta_pointeval_RESULTS.json"

lits = json.load(open(LITS))
t1 = QQ(22)/7; t2 = QQ(21)/13
out = {}
for stem in sorted(lits, key=lambda s: lits[s].get("c") or 0):
    e = lits[stem]
    if "lit" not in e:
        out[stem] = {"ok": False, "error": e.get("error")}
        continue
    try:
        Knot = sage_eval(e["lit"], locals={'matrix': matrix})
        t0 = time.time()
        Delta1, theta, theta_total, expo = Theta_numeric(Knot, t1, t2)
        dt = time.time() - t0
        is_zero = (theta_total == 0)
        out[stem] = {"ok": True, "c": e.get("c"), "src": e.get("src"),
                     "theta_at_22_7": str(theta_total),
                     "theta_core": str(theta),
                     "theta_zero": bool(is_zero),
                     "secs": float(dt)}
        tag = " θ=0 amphi" if is_zero else ""
        print("  %-46s c=%-4s %.3fs%s" % (stem[:46], str(e.get("c")), dt, tag), flush=True)
    except Exception as ex:
        out[stem] = {"ok": False, "c": e.get("c"), "error": "%s: %s" % (type(ex).__name__, ex)}
        print("  %-46s c=%-4s FEHLER %s" % (stem[:46], str(e.get("c")), type(ex).__name__), flush=True)
    json.dump(out, open(OUT, "w"), indent=2)

okv = sum(1 for v in out.values() if v.get("ok"))
z = sum(1 for v in out.values() if v.get("ok") and v.get("theta_zero"))
print("\n== %d/%d point-eval ok | %d θ=0 (amphi) -> %s ==" % (okv, len(out), z, OUT))
