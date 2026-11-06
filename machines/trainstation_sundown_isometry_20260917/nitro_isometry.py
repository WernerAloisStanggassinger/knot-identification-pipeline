# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Scene 7: nitro ≡ machine curve kf191_sh+1.375_tw+1.75 ? — isometry of the complements, each side from its own .fseries.
# In addition the two further cluster members kf201/kf202; counter-cases: ping-pong (S_24_d81_D6_amphi) and trainstation.
# Werner Stanggassinger + Claude Code, 17.09.2026.
import site_config
import sys, types, resource, json, glob
resource.setrlimit(resource.RLIMIT_AS, (6 << 30, 6 << 30))
import numpy as np
for _n, _b in (("complex", complex), ("float", float), ("int", int), ("bool", bool)):
    if not hasattr(np, _n): setattr(np, _n, _b)
sys.modules.setdefault("regina", types.ModuleType("regina"))
ROOT = site_config.ROOT
sys.path.insert(0, ROOT + "/claude_setup/tools"); sys.path.insert(0, site_config.PIPELINE_DIR)
from build_curve_json import parse_fseries, reconstruct_curve
from r3_polygon_regina_det import rot, pd_from_polygon
import spherogram, snappy

def aus_kurve(path, N=2048, seeds=(0, 9001, 9002, 9003)):
    off, a, b = parse_fseries(path); pts = reconstruct_curve(off, a, b, N); out = []
    for sd in seeds:
        try:
            pd, C = pd_from_polygon(pts @ rot(sd).T)
            L = spherogram.Link(pd); L.simplify(mode="global"); E = L.exterior()
            out.append(dict(seed=sd, roh=C, simpl=len(L.crossings), vol=float(E.volume()), sol=E.solution_type(), E=E, L=L))
        except Exception as e:
            out.append(dict(seed=sd, fehler=str(e)[:120]))
    return out

def iso(A, B, versuche=8):
    last = ''
    for i in range(versuche):
        try: return bool(A.is_isometric_to(B))
        except Exception as e:
            A = A.copy(); A.randomize(); last = str(e)[:80]
    return 'unentschieden: ' + last

C = site_config.GENERATED_KNOTS_DIR + "/new_found_knots/S_24_d961_D6_amphi"
pfade = {"nitro": ROOT + "/own_knots/nitro/nitro.fseries",
         "kf191": glob.glob(C + "/*kf191*/*.fseries")[0],
         "kf201": glob.glob(C + "/*kf201*/*.fseries")[0],
         "kf202": glob.glob(C + "/*kf202*/*.fseries")[0],
         "trainstation(GEGENFALL)": ROOT + "/own_knots/trainstation/trainstation.fseries"}
pp = glob.glob(ROOT + "/own_knots/ping-pong/*.fseries") or glob.glob(ROOT + "/own_knots/ping_pong/*.fseries")
if pp: pfade["ping-pong(GEGENFALL)"] = [p for p in pp if not p.split('/')[-1].startswith('_')][0] if any(not p.split('/')[-1].startswith('_') for p in pp) else pp[0]
K = {}
for nm, p in pfade.items():
    X = aus_kurve(p); K[nm] = [x for x in X if 'E' in x]
    for x in X: print(nm, {k: v for k, v in x.items() if k not in ("E", "L")}, flush=True)
n = K["nitro"]
res = {}
for nm in pfade:
    if nm == "nitro" or not K[nm]: continue
    res[f"nitro~{nm}"] = iso(n[0]["E"], K[nm][0]["E"])
if len(n) > 1 and len(K["kf191"]) > 1:
    res["nitro(seed2)~kf191(seed2)"] = iso(n[1]["E"], K["kf191"][1]["E"])
print("ERGEBNIS", json.dumps(res), flush=True)
import sympy
for nm in ("nitro", "kf191"):
    E = K[nm][0]["E"]; G = E.symmetry_group()
    print(nm, "H1", E.homology(), "Sym", G, "amphicheiral", G.is_amphicheiral(),
          "det", abs(sympy.Matrix((lambda V: [[V[i][j]+V[j][i] for j in range(len(V))] for i in range(len(V))])([list(r) for r in K[nm][0]["L"].seifert_matrix()])).det()), "min simpl crossings", min(x["simpl"] for x in K[nm]), flush=True)
