# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
import site_config
import json, os, glob, time, knotlib as K
import snappy
MD=(site_config.MACHINE_YARD + '/corpus_vs_maguire/maguire_data')
WORK=(site_config.MACHINE_YARD + '/corpus_vs_maguire/work')
res=json.load(open((site_config.MACHINE_YARD + '/corpus_vs_maguire/RESULT.json')))
byid={r['distinct_knot_id']:r for r in json.load(open(WORK+'/targets_novelhyp.json'))}
unres=[u for u in res['unresolved_no_poly'] if u.get('vol') and u.get('det')]
undets={u['det'] for u in unres}
# targets by det for quick lookup
from collections import defaultdict
by_det=defaultdict(list)
for u in unres: by_det[u['det']].append(u)
print("refine: %d unresolved with vol+det, %d distinct dets"%(len(unres),len(undets)),flush=True)

# stream once, collect det-matching rows
rows=[]
files=sorted(glob.glob(MD+'/*_hyp_jones.tar.xz'))
t0=time.time(); n=0
for p in files:
    cross=os.path.basename(p).split('_')[0]
    for dt,poly in K.stream_maguire(p):
        n+=1
        d=K.det_from_jones(K.parse_jones(poly))
        if d in undets: rows.append((d,cross,dt))
    print("scanned %s total=%d rows_hit=%d %.0fs"%(cross,n,len(rows),time.time()-t0),flush=True)
print("det-matching maguire rows: %d"%len(rows),flush=True)

# compute volumes, match to unresolved targets by (det, vol +-1e-4)
resolved=[]; t1=time.time()
for i,(d,cross,dt) in enumerate(rows):
    try: v=float(snappy.Manifold('DT:'+repr([tuple(K.decode_dt(dt))])).volume())
    except: continue
    for u in by_det[d]:
        if abs(v-float(u['vol']))<1e-4:
            resolved.append({'target_id':u['target_id'],'det':d,'maguire_cross':cross,
                             'maguire_dt':dt,'maguire_vol':v,'our_vol':u['vol'],'c_ub':u['c_ub']})
    if (i+1)%2000==0: print("  vol %d/%d resolved=%d %.0fs"%(i+1,len(rows),len(resolved),time.time()-t1),flush=True)
# dedup resolved by target
seen={}
for r in resolved:
    if r['target_id'] not in seen or r['maguire_cross']<seen[r['target_id']]['maguire_cross']:
        seen[r['target_id']]=r
json.dump(list(seen.values()), open(WORK+'/unresolved_resolved.json','w'), indent=1)
print("=== UNRESOLVED REFINEMENT: %d of %d now volume-CONFIRMED retractions ==="%(len(seen),len(unres)),flush=True)
for r in seen.values():
    print("  RETRACT(no-poly) %-28s det=%s -> TRUE c=%s vol=%.5f"%(r['target_id'],r['det'],r['maguire_cross'],r['maguire_vol']),flush=True)
print("REFINE DONE",flush=True)
