#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Per-entry receipts for the 40 Maguire retractions (inventory v1 item C1, 2026-10-14).
For every retraction key: our corpus record(s) (isosig of the complement, stored signature, stored Jones) against the
Maguire knot rebuilt from its DT code: SnapPy is_isometric_to (identity up to mirror, Gordon-Luecke) + chirality via
signature (if non-zero) or Jones V(q) vs V(1/q). Writes receipts.json. Read-only on the corpus."""
import site_config
import json, sys, signal, sympy, re
import snappy, spherogram
sys.path.insert(0,(site_config.MACHINE_YARD + '/corpus_vs_maguire/work')); import knotlib as K  # decode_dt = the July decoder (Disc #512/#516)
ROOT=site_config.ROOT
db=json.load(open(site_config.PIPELINE_DIR + '/corpus_db.json'))
mg=json.load(open(site_config.PIPELINE_DIR + '/maguire_retractions.json'))['retractions']
q=sympy.symbols('q')
def jones_sym(s):
    if not s: return None
    t=str(s).replace('^','**'); t=re.sub(r'(\d)\s+([A-Za-z])', r'\1*\2', t)
    return sympy.sympify(t)
def alarm(*a): raise TimeoutError
signal.signal(signal.SIGALRM, alarm)
out={}
for key,m in mg.items():
    recs=[r for r in db['knots'] if r.get('folder')==key or r.get('folder','').split('__')[0]==key or r.get('cluster_folder')==key or (r.get('folder','').startswith(key+'__'))]
    if not recs: recs=[r for r in db['knots'] if 'maguire' in json.dumps(r).lower() and (key in json.dumps(r))]
    try:
        signal.alarm(120)
        num=K.decode_dt(m['maguire_dt']); L=spherogram.Link('DT:'+repr([tuple(num)]))
        MM=L.exterior(); vol_m=float(MM.volume()); sig_m=L.signature(); jm=L.jones_polynomial()
        jm_s=sympy.sympify(str(jm).replace('^','**'))
        signal.alarm(0)
    except Exception as ex:
        out[key]=dict(error='maguire side: '+repr(ex)[:120]); continue
    rows=[]
    for r in recs:
        row=dict(folder=r.get('folder'), stored_sig=r.get('sig'), stored_vol=r.get('vol') or r.get('volume'))
        try:
            signal.alarm(300)
            iso=r.get('isosig')
            if not iso: row['isometric']=None; row['note']='no isosig stored'
            else:
                N=snappy.Manifold(iso)
                row['isometric']=bool(N.is_isometric_to(MM))
                row['our_vol']=float(N.volume())
            signal.alarm(0)
        except Exception as ex:
            signal.alarm(0); row['isometric']=None; row['note']='isometry: '+repr(ex)[:100]
        # chirality
        sr=r.get('sig')
        if sr is not None and sig_m!=0:
            row['chirality']='same' if sr==sig_m else ('mirror' if sr==-sig_m else 'sig mismatch')
        else:
            js=jones_sym(r.get('jones'))
            if js is None: row['chirality']='undecided (sig 0, no stored Jones)'
            else:
                x=list(js.free_symbols)
                if len(x)!=1: row['chirality']='undecided (jones symbols %s)'%x
                else:
                    x=x[0]; jms=jm_s.subs(list(jm_s.free_symbols)[0], x) if jm_s.free_symbols else jm_s
                    if sympy.simplify(js-jms)==0: row['chirality']='same (Jones)'
                    elif sympy.simplify(js-jms.subs(x,1/x))==0: row['chirality']='mirror (Jones)'
                    else: row['chirality']='Jones mismatch'
        rows.append(row)
    out[key]=dict(maguire_dt=m['maguire_dt'], maguire_sig=sig_m, maguire_vol=vol_m, true_min_crossing=m['true_min_crossing'], n_records=len(rows), records=rows)
    print(key, 'sig_m=%s'%sig_m, [(x['folder'], x.get('isometric'), x.get('chirality')) for x in rows], flush=True)
json.dump(out, open(site_config.MACHINE_YARD+'/corpus_vs_burton/maguire_receipts_20261014/receipts.json','w'), indent=1)
iso_all=sum(1 for v in out.values() for x in v.get('records',[]) if x.get('isometric') is True)
n=sum(v.get('n_records',0) for v in out.values())
print('DONE records=%d isometric=%d errors=%d'%(n, iso_all, sum(1 for v in out.values() if 'error' in v)))
