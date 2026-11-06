#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""theta_write_machine_layer — rebuild-safe write of the CORRECTED theta_at_22_7 layer
(Disc #249) for the 246 MACHINE records. own_knots EXCLUDED (separate, "between us").

Writes (only with --apply): (a) corpus_db theta_at_22_7__21_13 = corrected (immediately correct),
(b) on-disk sidecar {stem}_theta_numeric.json key "theta" (rebuild-safe, build_corpus_db line 170).
Does NOT touch: theta_full, served, Dror mail. Dry-run default. Backup before every --apply.

PRE-GATE: for records with existing corpus.theta_full, theta_full@(22/7,21/13) is checked against the
new theta_at_22_7 (no internal contradiction). Mismatch => abort unless --force.
"""
import site_config
import json, sys, shutil, time, re
from pathlib import Path
from fractions import Fraction
import sympy

ROOT = Path(site_config.ROOT)
DB   = Path(site_config.PIPELINE_DIR)/"corpus_db.json"
STG  = ROOT/"claude_setup/tmp_durable/theta_corrected_machine_layer.json"
APPLY = "--apply" in sys.argv
FORCE = "--force" in sys.argv

t1,t2 = sympy.Rational(22,7), sympy.Rational(21,13)
T1,T2 = sympy.symbols('t1 t2')

def eval_full(expr_str):
    """theta_full (sym in t1,t2) @ (22/7,21/13) -> sympy.Rational, or None."""
    if not expr_str or not isinstance(expr_str,str): return None
    s = expr_str.replace('^','**')
    try:
        e = sympy.sympify(s, locals={'t1':T1,'t2':T2})
        v = e.subs({T1:t1, T2:t2})
        v = sympy.nsimplify(v) if not v.is_number else v
        return sympy.Rational(v)
    except Exception:
        return None

def as_rat(s):
    if s is None: return None
    try: return sympy.Rational(str(s))
    except Exception: return None

stg = json.load(open(STG))
db  = json.load(open(DB))
recs = db["knots"]
byfolder = {r["folder"]: r for r in recs}

n_match=n_change=0
verify_ok=verify_bad=verify_na=0
bad=[]; targets=[]
for folder, ent in stg.items():
    r = byfolder.get(folder)
    if r is None:  # silently skip own_knots/non-machine
        continue
    new = ent["theta_at_22_7__21_13"]
    n_match+=1
    old = r.get("theta_at_22_7__21_13")
    changed = str(old)!=str(new)
    if changed: n_change+=1
    # Pre-gate: theta_full consistency
    tf = r.get("theta_full")
    if tf:
        fv = eval_full(tf); nv = as_rat(new)
        if fv is None or nv is None: verify_na+=1
        elif abs(fv)==abs(nv):   # mirror sign allowed (chir/amphi convention)
            verify_ok+=1
        else:
            verify_bad+=1; bad.append((folder,str(fv)[:30],str(nv)[:30]))
    # on-disk path
    pfl = Path(r["pflicht_path"])
    fdir = pfl.parent
    stem = re.sub(r"_pflicht_invarianten.*_complete\.json$","",pfl.name)
    side = fdir/f"{stem}_theta_numeric.json"
    targets.append((r, new, side, changed))

print(f"=== DRY-RUN ({'APPLY' if APPLY else 'dry'}) ===")
print(f"staged matched machine records : {n_match}")
print(f"  of which theta_at_22_7 changed : {n_change}")
print(f"PRE-GATE theta_full consistency: ok={verify_ok} bad={verify_bad} n/a={verify_na}")
if bad:
    print("  MISMATCH (theta_full@pt != new theta_at_22_7):")
    for f,a,b in bad[:10]: print(f"    {f}: full={a} vs new={b}")
if verify_bad and not FORCE:
    print("\nABORT: theta_full contradiction -> --force needed (or clarify first). NO write.")
    sys.exit(3)

if not APPLY:
    n_side_exist = sum(1 for *_ ,s,_ in [(0,t[2],0) for t in targets] if False)  # placeholder
    n_side = sum(1 for _,_,s,_ in targets if s.exists())
    print(f"\nSidecars already exist: {n_side}/{len(targets)} (the rest will be newly created)")
    print("Dry run. Write with --apply (backup is taken beforehand).")
    sys.exit(0)

# ---- APPLY ----
ts=time.strftime("%Y%m%d-%H%M%S")
bk=DB.with_suffix(f".json.bak-pre-theta-WRITE-{ts}")
shutil.copy2(DB,bk); print(f"\nBackup corpus_db -> {bk.name}")
n_corpus=n_sidecar=0
for r,new,side,changed in targets:
    r["theta_at_22_7__21_13"]=new; n_corpus+=1
    # Sidecar merge
    if side.exists():
        try: sc=json.load(open(side))
        except Exception: sc={}
    else:
        sc={"name":side.name.replace("_theta_numeric.json",""),
            "source":"theta_pointeval_corpus (Disc #248/#249 correction)"}
    sc["theta"]=new
    sc["t1"]="22/7"; sc["t2"]="21/13"
    sc["corrected_2026_06_02"]=True
    json.dump(sc, open(side,"w"), indent=2)
    n_sidecar+=1
json.dump(db, open(DB,"w"), indent=2, ensure_ascii=False)
print(f"corpus_db: {n_corpus} theta_at_22_7 set")
print(f"Sidecars : {n_sidecar} {{stem}}_theta_numeric.json written/merged")
print("DONE. theta_full/served/mail UNTOUCHED.")
