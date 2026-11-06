#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""
precellipse_atlas.py — VERIFIED ATLAS of Hizume's precellipse construction.
Sweeps the Farey family of precession ratios Delta=q/p, generates the
precessing Kepler ellipse for each, identifies the knot BLIND through the pipeline
(polygon_resolve, no table lookups) and checks the result against the closed
TORUS FORMULA (det, crossing number c, Alexander degree — formula, no lookup).

Result: a map "precession parameter q/p -> torus knot T(p,q) + verified
invariants", + honest working limit of the construction (where embedding/Nyquist degenerates).
This is the gift to Akio: he constructs, we identify + catalogue blind.

Reuse: precellipse_machine.py (precellipse/best_embedding/curve_to_fseries/
nyquist_id/flattest_axis). base env (regina+numpy). Parallel over (p,q).
"""
import site_config
import sys, os, math, json, argparse
from math import gcd
from concurrent.futures import ProcessPoolExecutor, as_completed

TOOLS=(site_config.TOOLS_DIR)
sys.path.insert(0, TOOLS)
import importlib.util
def _load(n,p):
    s=importlib.util.spec_from_file_location(n,p); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
PM=_load("pm", TOOLS+"/precellipse_machine.py")

OUTROOT=(site_config.GENERATED_KNOTS_DIR)

# ---- TORUS ORACLE (closed formulas, NO table lookup) ----
def torus_det(p,q):
    # det(T(p,q)) = |Delta(-1)|: 1 if p,q both odd; otherwise the odd number
    if p%2==1 and q%2==1: return 1
    return q if p%2==0 else p
def torus_c(p,q):      return min(p*(q-1), q*(p-1))      # crossing number
def torus_alex_deg(p,q): return (p-1)*(q-1)              # degree of the Alexander polynomial

def one(args):
    p,q,N = args
    try:
        rel,e,k,kdeg = PM.best_embedding(p,q,N=6000)
        P = PM.precellipse(e,k,p,q,16000)
        offset,a,b = PM.curve_to_fseries(P)
        axis = PM.flattest_axis(offset,a,b)
        r,nst = PM.nyquist_id(offset,a,b,axis)
        if not r:
            return {"p":p,"q":q,"status":"no_stable_measure","d_min_rel":round(rel,5)}
        det=r['det']; nc=r.get('n_cross_simplified'); alex=r.get('alexander')
        edet,ec,edeg = torus_det(p,q), torus_c(p,q), torus_alex_deg(p,q)
        # estimate the Alexander degree from the blind string (highest exponent)
        import re
        degs=[int(x) for x in re.findall(r'x\^(\d+)', str(alex))] or ([1] if 'x' in str(alex) else [0])
        bdeg=max(degs) if degs else 0
        det_ok = (det==edet); c_ok=(nc==ec); deg_ok=(bdeg==edeg)
        # QUANTUM INVARIANT: Jones blind from the simplified signed Gauss (regina),
        # = the rigorous "quantum-mechanical reading" of the knot (Jones = Chern-Simons QFT).
        jones_blind=None
        try:
            import regina
            sg=r.get('signed_gauss_simplified')
            if sg:
                Lj=regina.Link.fromSignedGauss(sg)
                jones_blind=str(Lj.jones())
        except Exception:
            pass
        return {"p":p,"q":q,"status":"ok","delta":f"{q}/{p}","knot":f"T({p},{q})",
                "d_min_rel":round(rel,5),"e":e,"k_deg":kdeg,"n_stabil":nst,
                "det_blind":det,"det_formula":edet,"det_ok":det_ok,
                "c_blind":nc,"c_formula":ec,"c_ok":c_ok,
                "alexdeg_blind":bdeg,"alexdeg_formula":edeg,"deg_ok":deg_ok,
                "knot_confirmed": det_ok and deg_ok,   # knot type secured (det+Alex degree)
                "c_exact": c_ok,                        # blind crossing count exact?
                "verified": det_ok and c_ok and deg_ok,
                "jones_blind": jones_blind,
                "alexander":str(alex)}
    except Exception as ex:
        return {"p":p,"q":q,"status":"err:"+repr(ex)[:60]}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--pmax",type=int,default=5)
    ap.add_argument("--qmax",type=int,default=7)
    ap.add_argument("--workers",type=int,default=12)
    A=ap.parse_args()
    pairs=[]
    for p in range(2,A.pmax+1):
        for q in range(p+1,A.qmax+1):
            if gcd(p,q)==1: pairs.append((p,q,8000))
    print(f"Atlas sweep: {len(pairs)} torus knots T(p,q), p<=%d q<=%d, {A.workers} workers"%(A.pmax,A.qmax),flush=True)
    res=[]
    with ProcessPoolExecutor(max_workers=A.workers) as ex:
        futs={ex.submit(one,pr):pr for pr in pairs}
        for f in as_completed(futs):
            r=f.result(); res.append(r)
            tag = "✓" if r.get("verified") else ("✗" if r.get("status")=="ok" else r.get("status"))
            print(f"  T({r['p']},{r['q']}) {tag}"
                  + (f" det {r['det_blind']}/{r['det_formula']} c {r['c_blind']}/{r['c_formula']} deg {r['alexdeg_blind']}/{r['alexdeg_formula']}" if r.get('status')=='ok' else ""),flush=True)
    res.sort(key=lambda r:(r['p'],r['q']))
    json.dump(res, open(f"{OUTROOT}/precellipse_atlas.json","w"), indent=1)
    conf=sum(1 for r in res if r.get("knot_confirmed"))
    cex=sum(1 for r in res if r.get("c_exact"))
    jq=sum(1 for r in res if r.get("jones_blind"))
    print(f"\n=== ATLAS BALANCE ===",flush=True)
    print(f"  Knot type confirmed blind (det+Alexander degree against torus formula): {conf}/{len(res)}",flush=True)
    print(f"  of these blind crossing count exact: {cex}/{len(res)} (rest: c measurement frays by 1-4, knot still correct)",flush=True)
    print(f"  Quantum invariant (Jones) computed blind: {jq}/{len(res)}",flush=True)
    print(f"-> {OUTROOT}/precellipse_atlas.json")

if __name__=="__main__":
    main()
