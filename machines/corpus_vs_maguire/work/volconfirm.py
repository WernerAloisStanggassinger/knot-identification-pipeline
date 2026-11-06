# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
import site_config
import json, importlib, sys, time
import knotlib as K
sys.path.insert(0,(site_config.MACHINE_YARD + '/corpus_vs_maguire/work'))
import snappy
WORK=(site_config.MACHINE_YARD + '/corpus_vs_maguire/work')

raw=json.load(open(WORK+'/match_raw.json'))
cands=raw['candidates']
targets={r['distinct_knot_id']:r for r in json.load(open(WORK+'/targets_novelhyp.json'))}

def dt_volume(dtstr):
    num=K.decode_dt(dtstr)
    M=snappy.Manifold('DT:'+repr([tuple(num)]))
    try:
        v=M.volume()
    except Exception:
        return None
    # try to firm up
    return float(v)

volcache={}
results=[]
t0=time.time()
for i,c in enumerate(cands):
    dt=c['dt']
    if dt not in volcache:
        try:
            volcache[dt]=dt_volume(dt)
        except Exception as e:
            volcache[dt]=None
    mvol=volcache[dt]
    for tid in c['target_ids']:
        tr=targets[tid]
        ourvol=tr.get('vol')
        vol_ok = (mvol is not None and ourvol is not None and abs(mvol-ourvol)<1e-3)
        results.append({
            'target_id':tid,'our_vol':ourvol,'our_det':tr.get('det'),'our_c_upperbound':tr.get('c'),
            'maguire_cross':c['cross'],'maguire_dt':dt,'maguire_vol':mvol,
            'maguire_det':c['maguire_det'],
            'vol_confirmed':bool(vol_ok),
            'vol_diff':(None if mvol is None or ourvol is None else abs(mvol-ourvol)),
        })
    if (i+1)%20==0:
        print("  vol %d/%d (%.0fs)"%(i+1,len(cands),time.time()-t0),flush=True)

json.dump(results, open(WORK+'/volconfirm_results.json','w'), indent=1)
conf=[r for r in results if r['vol_confirmed']]
coll=[r for r in results if not r['vol_confirmed']]
retracted_ids=sorted(set(r['target_id'] for r in conf))
collision_ids=sorted(set(r['target_id'] for r in coll) - set(retracted_ids))
print("=== VOLUME CONFIRM ===",flush=True)
print("candidate pairs:",len(results),flush=True)
print("volume-CONFIRMED retractions (distinct targets):",len(retracted_ids),flush=True)
print("collision-only (jones match, vol mismatch) distinct targets:",len(collision_ids),flush=True)
for r in conf:
    print("  RETRACT %s  our_c<=%s -> TRUE c=%s  vol %.6f==%.6f  det %s"%(
        r['target_id'],r['our_c_upperbound'],r['maguire_cross'],r['maguire_vol'] or -1,r['our_vol'] or -1,r['our_det']),flush=True)
json.dump({'retracted':retracted_ids,'collision_only':collision_ids},open(WORK+'/adjudication_ids.json','w'))
