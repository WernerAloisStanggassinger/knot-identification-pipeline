#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""COUNTER-CASE BEFORE BUILDING for the alex normal form (Disc #809).

The fix must merge EXACTLY ONE pair — the proven one. If it merges more,
the corpus loses real knots; if it merges none, it has no effect.
"""
import site_config
import json, sys
sys.path.insert(0, (site_config.TOOLS_DIR))
from alexander_normal_form import canonical, terms
from collections import defaultdict

DB = (site_config.PIPELINE_DIR + '/corpus_db.json')
ks = json.load(open(DB))['knots']
f = 0

print("="*76); print("COUNTER-CASE alex normal form (Disc #809)"); print("="*76)

# (1) Full coverage: no record may become unparsable
mit = [m for m in ks if m.get('alex')]
unparsbar = [m for m in mit if terms(m['alex']) is None]
print(f"\n(1) records with alex: {len(mit)}   unparsable: {len(unparsbar)} "
      f"{'✅' if not unparsbar else '🔴'}")
for m in unparsbar[:5]: print(f"      {m['stem']}: {str(m['alex'])[:60]}")
f += len(unparsbar)

# (2) POSITIVE: the proven pair must have the same normal form
A = [m for m in ks if m['stem']=='poly_s08'][0]
B = [m for m in ks if m['stem']=='kf332_polys08_sh+0.5_tw-1'][0]
gleich = canonical(A['alex']) == canonical(B['alex'])
print(f"\n(2) POSITIVE poly_s08 / kf332_polys08 same normal form: {gleich} "
      f"{'✅' if gleich else '🔴 fix has no effect'}")
f += (not gleich)

# (3) NEGATIVE: the normal form must not make DIFFERENT knots equal.
#     Count: how many distinct groups arise NEWLY through the normal form?
def schluessel(m, kanon):
    v = m.get('vol')
    volk = round(v,6) if isinstance(v,(int,float)) and v>0.5 else None
    a = canonical(m.get('alex')) if kanon else str(m.get('alex'))
    return (str(m.get('det')), str(m.get('sig')), a,
            str(m.get('theta_at_22_7__21_13')), volk)
for kanon,lab in ((False,'ALT (String)'),(True,'NEU (Normalform)')):
    g = defaultdict(list)
    for m in ks: g[schluessel(m,kanon)].append(m)
    print(f"    {lab:22s} -> {len(g)} groups")
    if kanon: neu = g
    else: alt = g
delta = len(alt) - len(neu)
print(f"\n(3) groups OLD {len(alt)} -> NEW {len(neu)}   merged: {delta} "
      f"{'✅ exactly one' if delta==1 else '🔴 expected exactly 1'}")
f += (delta != 1)

# (4) WHICH one was merged? Must be the proven pair.
alt_von_neu = defaultdict(set)
for m in ks:
    alt_von_neu[schluessel(m,True)].add(schluessel(m,False))
verschmolzen = {k:v for k,v in alt_von_neu.items() if len(v)>1}
print(f"\n(4) merged groups: {len(verschmolzen)}")
for k,v in verschmolzen.items():
    stems = sorted({m['stem'] for m in ks if schluessel(m,True)==k})
    print(f"      det={k[0]} sig={k[1]} vol={k[4]}  <- {len(v)} old groups")
    print(f"      stems: {stems}")
    soll = {'poly_s08','kf332_polys08_sh+0.5_tw-1','kf336_polys08_sh+1_tw-1','kf333_polys08_sh+0.5_tw0'}
    ok = set(stems)==soll
    print(f"      == the receipted pair (4 records): {ok} {'✅' if ok else '🔴'}")
    f += (not ok)

# (5) NEGATIVE: proven_distinctions must NOT be lifted
import pathlib
pd = pathlib.Path((site_config.PIPELINE_DIR + '/proven_distinctions.json'))
if pd.exists():
    d = json.loads(pd.read_text()).get('distinctions', {})
    print(f"\n(5) proven_distinctions: {len(d)} case(s)")
    for fall,e in d.items():
        a,b = e.get('abgespalten'), e.get('verbleibt')
        ma = [m for m in ks if m['stem']==a]; mb = [m for m in ks if m['stem']==b]
        if not ma or not mb: print(f"      {fall}: stems missing"); continue
        kollidiert = schluessel(ma[0],True)==schluessel(mb[0],True)
        print(f"      {fall}: {a} / {b} -> normal-form key equal: {kollidiert} "
              f"(separated by proven_distinctions anyway)")

print("\n"+"="*76)
print("ALL PROBES PASSED" if f==0 else f"{f} ERRORS")
print("="*76)
sys.exit(1 if f else 0)
