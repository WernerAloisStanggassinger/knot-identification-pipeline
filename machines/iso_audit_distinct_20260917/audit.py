# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Isometry audit of all multi-groups of the distinctness key (17.09.2026, Werner: "website iso check").
# Occasion Disc #795: label_distinct_knots.py merges by (det, sig, alex, θ@1 point, vol@6) — that is NO identity proof.
# Per group: each member from ITS served .fseries (folder/stem.fseries) -> complement; is_isometric_to against the representative.
# Counter-case per group: representative against the volume-nearest representative of ANOTHER group (must be False).
# Non-hyperbolic groups: isometry not applicable -> only flagged (proof route proven_identities / braid conjugation).
# Latch: every task in its own subprocess, 3 GB RLIMIT_AS, 900 s.
import site_config
import json, os, sys, subprocess, time, math
from concurrent.futures import ThreadPoolExecutor, as_completed
HIER = os.path.dirname(os.path.abspath(__file__))
ROOT = site_config.ROOT
PY = site_config.PYTHON_MAIN

WORKER = r'''
import sys, types, resource, json
resource.setrlimit(resource.RLIMIT_AS, (3 << 30, 3 << 30))
import numpy as np
for _n, _b in (("complex", complex), ("float", float), ("int", int), ("bool", bool)):
    if not hasattr(np, _n): setattr(np, _n, _b)
sys.modules.setdefault("regina", types.ModuleType("regina"))
import site_config
ROOT = site_config.ROOT
sys.path.insert(0, ROOT + "/claude_setup/tools"); sys.path.insert(0, site_config.PIPELINE_DIR)
from build_curve_json import parse_fseries, reconstruct_curve
from r3_polygon_regina_det import rot, pd_from_polygon
import spherogram, snappy
def ext(fs, vol_soll):
    off, a, b = parse_fseries(fs); pts = reconstruct_curve(off, a, b, 2048); best = None
    for sd in (0, 9001, 9002, 9003):
        try:
            pd, C = pd_from_polygon(pts @ rot(sd).T)
            L = spherogram.Link(pd); L.simplify(mode="global"); E = L.exterior()
            v = float(E.volume()); ok = str(E.solution_type()).startswith("all tetrahedra positively")
            info = dict(seed=sd, vol=v, geom=ok, simpl=len(L.crossings))
            if ok and (vol_soll is None or abs(v - vol_soll) < 1e-4): return E, info
            best = best or (E, info)
        except Exception as e:
            continue
    return best if best else (None, dict(fehler="keine Projektion"))
def iso(A, B):
    last = ""
    for i in range(8):
        try: return bool(A.is_isometric_to(B))
        except Exception as e:
            A = A.copy(); A.randomize(); last = str(e)[:80]
    return "unentschieden: " + last
job = json.loads(sys.argv[1])
EA, ia = ext(job["fa"], job["va"]); EB, ib = ext(job["fb"], job["vb"])
out = dict(job, info_a=ia, info_b=ib)
out["iso"] = iso(EA, EB) if (EA is not None and EB is not None) else "fehler: " + ("a " if EA is None else "") + ("b" if EB is None else "")
print("REC" + json.dumps(out))
'''

def fseries(r):
    # served .fseries lies next to the mandatory JSON (folder is only the folder name, no path)
    p = os.path.join(os.path.dirname(r.get("pflicht_path") or ""), r["stem"] + ".fseries")
    return p if os.path.exists(p) else None

def lauf(job):
    try:
        # OpenBLAS reserves large virtual buffers per thread and otherwise fails on RLIMIT_AS (measured 17.09.: 313/313 aborted)
        env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
        pr = subprocess.run([PY, "-c", WORKER, json.dumps(job)], capture_output=True, text=True, timeout=900, env=env)
        for z in pr.stdout.splitlines():
            if z.startswith("REC"): return json.loads(z[3:])
        return dict(job, iso="abbruch rc=%s %s" % (pr.returncode, (pr.stderr or "")[-160:]))
    except subprocess.TimeoutExpired:
        return dict(job, iso="timeout 900s")

def main():
    db = json.load(open(site_config.PIPELINE_DIR + "/corpus_db.json"))
    recs = list(db["knots"].values()) if isinstance(db["knots"], dict) else db["knots"]
    G = {}
    for r in recs:
        if r.get("distinct_knot_id"): G.setdefault(r["distinct_knot_id"], []).append(r)
    multi = {k: v for k, v in G.items() if len(v) > 1}
    jobs, nichthyp, ohne = [], [], []
    reps = {}
    for k, v in multi.items():
        rep = next((r for r in v if r.get("is_distinct_representative")), v[0])
        hyp = bool(rep.get("is_hyperbolic")) and rep.get("vol_valid") is not False and (rep.get("vol") or 0) > 0.5
        if not hyp: nichthyp.append(dict(gruppe=k, n=len(v), klasse=rep.get("realizes_class"))); continue
        reps[k] = rep
        for r in v:
            if r is rep: continue
            fa, fb = fseries(rep), fseries(r)
            if not fa or not fb: ohne.append(dict(gruppe=k, mitglied=r["stem"], rep=rep["stem"])); continue
            jobs.append(dict(art="mitglied", gruppe=k, rep=rep["stem"], mitglied=r["stem"], fa=fa, fb=fb, va=rep.get("vol"), vb=r.get("vol"),
                             own=("own_knots" in r["folder"]) or ("own_knots" in rep["folder"])))
    ks = sorted(reps, key=lambda k: reps[k]["vol"])

    def spiegelverdacht(a, b):
        """Is b the MIRROR of a? Then b is no good as counter-case.

        A chiral knot and its mirror image have ISOMETRIC complements —
        the isometry may be orientation-reversing, and is_isometric_to()
        does not distinguish that. Gordon-Luecke needs orientation-preserving
        homeomorphy. So a mirror pair reports 'wrongly merged',
        where the labeler is right: there the SIGNATURE separates.
        Exactly this way two spurious findings arose on 17.09.2026
        (11n_133/inv_11n_133 sig -4/+4 and 12n_839/12n_839m sig +4/-4).
        Tell: same det, same volume, sig(b) == -sig(a) != 0.
        """
        try:
            sa, sb = int(a.get("sig")), int(b.get("sig"))
        except (TypeError, ValueError):
            return False
        if sa == 0 or sa != -sb:
            return False
        if str(a.get("det")) != str(b.get("det")):
            return False
        va, vb = a.get("vol"), b.get("vol")
        if not (isinstance(va, (int, float)) and isinstance(vb, (int, float))):
            return False
        return abs(va - vb) < 1e-6

    for i, k in enumerate(ks):
        # volume-nearest group that is NOT a mirror partner
        nb = None
        for j in list(range(i + 1, len(ks))) + list(range(i - 1, -1, -1)):
            if not spiegelverdacht(reps[k], reps[ks[j]]):
                nb = ks[j]
                break
        if nb is None:
            continue
        a, b = reps[k], reps[nb]
        if fseries(a) and fseries(b):
            jobs.append(dict(art="GEGENFALL", gruppe=k, rep=a["stem"], mitglied=b["stem"], gegen_gruppe=nb, fa=fseries(a), fb=fseries(b),
                             va=a.get("vol"), vb=b.get("vol"), own=False))
    jobs.sort(key=lambda j: (not j["own"], j["art"] != "mitglied"))
    json.dump(dict(gruppen=len(multi), hyperbolisch=len(reps), nichthyperbolisch=nichthyp, ohne_fseries=ohne, jobs=len(jobs)),
              open(os.path.join(HIER, "plan.json"), "w"), indent=1)
    print(f"groups {len(multi)} · hyperbolic {len(reps)} · non-hyperbolic {len(nichthyp)} · without fseries {len(ohne)} · tasks {len(jobs)}", flush=True)
    fertig = set()
    aus = os.path.join(HIER, "ergebnisse.jsonl")
    if os.path.exists(aus):
        for z in open(aus): d = json.loads(z); fertig.add((d["art"], d["gruppe"], d["mitglied"]))
    offen = [j for j in jobs if (j["art"], j["gruppe"], j["mitglied"]) not in fertig]
    with open(aus, "a") as f, ThreadPoolExecutor(int(sys.argv[1]) if len(sys.argv) > 1 else 6) as ex:
        for i, fu in enumerate(as_completed([ex.submit(lauf, j) for j in offen]), 1):
            d = fu.result(); f.write(json.dumps(d) + "\n"); f.flush()
            if (d["art"] == "mitglied" and d["iso"] is not True) or (d["art"] == "GEGENFALL" and d["iso"] is not False):
                print("!!", d["art"], d["gruppe"], d["rep"], d["mitglied"], d["iso"], flush=True)
            if i % 25 == 0: print(time.strftime("%H:%M:%S"), i, "/", len(offen), flush=True)
    print("DONE", time.strftime("%F %T"), flush=True)

if __name__ == "__main__":
    main()
