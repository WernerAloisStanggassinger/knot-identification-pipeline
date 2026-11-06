# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# trainstation ≡ sundown? — isometry of the complements, each side INDEPENDENTLY from its own .fseries (not from record PD:
# trainstation's diagram fields were copied from sundown on 28.06., Disc #354 C — those would be circular).
# Plus: each side individually against Burton 18ah_5090151 (DT qfapcjbmgeiokhrldn). Counter-case: a foreign knot must fail.
# Werner Stanggassinger + Claude Code, 17.09.2026.
import site_config
import sys, types, resource, json
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

def aus_kurve(stem, N=2048, seeds=(0, 9001, 9002, 9003)):
    off, a, b = parse_fseries(f"{ROOT}/own_knots/{stem}/{stem}.fseries")
    pts = reconstruct_curve(off, a, b, N)
    out = []
    for sd in seeds:
        try:
            pd, C = pd_from_polygon(pts @ rot(sd).T)
            L = spherogram.Link(pd); L.simplify(mode="global")
            E = L.exterior(); v = E.volume()
            out.append(dict(seed=sd, roh=C, simpl=len(L.crossings), vol=float(v), sol=E.solution_type(), E=E, L=L))
        except Exception as e:
            out.append(dict(seed=sd, fehler=str(e)[:120]))
    return out

def iso(A, B, versuche=8):
    for i in range(versuche):
        try:
            r = A.is_isometric_to(B)
            return bool(r)
        except Exception as e:
            A = A.copy(); A.randomize(); last = str(e)[:80]
    return 'unentschieden: ' + last

T = aus_kurve("trainstation"); S = aus_kurve("sundown")
for name, X in (("trainstation", T), ("sundown", S)):
    for x in X: print(name, {k: v for k, v in x.items() if k not in ("E", "L")})
t = next(x for x in T if 'E' in x); s = next(x for x in S if 'E' in x)
sys.path.insert(0, site_config.MACHINE_YARD + "/corpus_vs_maguire/work"); import knotlib as K
B = snappy.Manifold("DT:" + repr([tuple(K.decode_dt("qfapcjbmgeiokhrldn"))]))  # as in corpus_vs_burton/stage_b_v4.py
print("Burton 18ah_5090151 vol", float(B.volume()), B.solution_type())
res = {}
res["trainstation~sundown"] = iso(t["E"], s["E"])
res["trainstation~burton"] = iso(t["E"], B)
res["sundown~burton"] = iso(s["E"], B)
# second, independent projection per side (other seed) against the first one of the other side
t2 = [x for x in T if 'E' in x][1] if len([x for x in T if 'E' in x]) > 1 else t
s2 = [x for x in S if 'E' in x][1] if len([x for x in S if 'E' in x]) > 1 else s
res["trainstation(seed2)~sundown(seed2)"] = iso(t2["E"], s2["E"])
# Counter-cases: nitro (another D6-amphi knot) and a Burton neighbor must fail
N = aus_kurve("nitro", seeds=(0,)); n = next(x for x in N if 'E' in x)
print("nitro", {k: v for k, v in n.items() if k not in ("E", "L")})
res["GEGENFALL nitro~sundown"] = iso(n["E"], s["E"])
res["GEGENFALL figure8~sundown"] = iso(snappy.Manifold("4_1"), s["E"])
# homological check that they are knot complements + symmetry groups
for nm, E in (("trainstation", t["E"]), ("sundown", s["E"]), ("burton", B)):
    try:
        G = E.symmetry_group(); print(nm, "H1", E.homology(), "Sym", G, "amphicheiral", G.is_amphicheiral())
    except Exception as e: print(nm, "sym error", str(e)[:80])
print("ERGEBNIS", json.dumps(res))
