# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
import site_config
import json, importlib
import knotlib as K
WORK=(site_config.MACHINE_YARD + '/corpus_vs_maguire/work')
OUT=(site_config.MACHINE_YARD + '/corpus_vs_maguire')

targets=json.load(open(WORK+'/targets_novelhyp.json'))
byid={r['distinct_knot_id']:r for r in targets}
raw=json.load(open(WORK+'/match_raw.json'))
volres=json.load(open(WORK+'/volconfirm_results.json'))
det_set=set(json.load(open(WORK+'/maguire_det_set.json')))
max_det=max(det_set)

# per-target best volume-confirmed record
retracted={}      # tid -> record (min crossing)
collision={}      # tid -> best collision record
for r in volres:
    tid=r['target_id']
    if r['vol_confirmed']:
        cur=retracted.get(tid)
        # prefer smallest crossing / smallest vol_diff
        if cur is None or r['maguire_cross']<cur['maguire_cross']:
            retracted[tid]=r
    else:
        cur=collision.get(tid)
        if cur is None or (r['vol_diff'] or 9e9)<(cur['vol_diff'] or 9e9):
            collision[tid]=r

retracted_ids=set(retracted)
collision_ids=set(collision)-retracted_ids

# classify all 395 targets
cats={'retracted':[], 'collision_only':[], 'confirmed_beyond_jones':[],
      'confirmed_beyond_det':[], 'unresolved_no_poly':[]}
for tid,r in byid.items():
    has_j=bool(r.get('jones'))
    det=r.get('det')
    if tid in retracted_ids:
        cats['retracted'].append(tid); continue
    if tid in collision_ids:
        cats['collision_only'].append(tid); continue
    if has_j:
        # jones present, not volume-confirmed -> its exact canon-jones is absent from Maguire<=17
        # (Maguire hyp census complete) => rigorously beyond 17
        cats['confirmed_beyond_jones'].append(tid)
    else:
        if det is None:
            cats['unresolved_no_poly'].append(tid)
        elif det not in det_set:
            cats['confirmed_beyond_det'].append(tid)   # no <=17 hyp knot has this det
        else:
            cats['unresolved_no_poly'].append(tid)

summary={k:len(v) for k,v in cats.items()}
n_conf_beyond = summary['confirmed_beyond_jones']+summary['confirmed_beyond_det']

# build retracted detail
retr_detail=[]
for tid in sorted(retracted_ids):
    r=retracted[tid]; t=byid[tid]
    retr_detail.append({
        'target_id':tid,'our_c_upperbound':t.get('c'),'true_min_crossing':r['maguire_cross'],
        'det':t.get('det'),'vol':t.get('vol'),'maguire_vol':r['maguire_vol'],
        'vol_diff':r['vol_diff'],'maguire_dt':r['maguire_dt'],'n_realisierungen':t.get('n_realisierungen'),
        'stem':t.get('stem')})

coll_detail=[]
for tid in sorted(collision_ids):
    r=collision[tid]; t=byid[tid]
    coll_detail.append({'target_id':tid,'our_c_upperbound':t.get('c'),'det':t.get('det'),
        'our_vol':t.get('vol'),'nearest_maguire_vol':r['maguire_vol'],'vol_diff':r['vol_diff'],
        'maguire_cross':r['maguire_cross'],'maguire_dt':r['maguire_dt']})

# unresolved detail (small set, list them)
unres_detail=[{'target_id':tid,'det':byid[tid].get('det'),'c_ub':byid[tid].get('c'),
               'vol':byid[tid].get('vol')} for tid in sorted(cats['unresolved_no_poly'],
               key=lambda x:(byid[x].get('det') or 0))]

result={
 'n_targets_novelhyp_this_run':len(targets),
 'meta_distinct_novel_hyperbolic':373,
 'meta_mirror_collapsed_headline':368,
 'maguire_max_det_le17_hyp':max_det,
 'maguire_det_set_size':len(det_set),
 'maguire_rows_scanned':raw['nrows'],
 'summary':summary,
 'n_confirmed_beyond_17':n_conf_beyond,
 'retracted':retr_detail,
 'collision_only':coll_detail,
 'unresolved_no_poly':unres_detail,
}
json.dump(result, open(OUT+'/RESULT.json','w'), indent=1)
print(json.dumps({'summary':summary,'n_conf_beyond':n_conf_beyond,
      'n_retracted':len(retracted_ids),'n_collision':len(collision_ids),
      'max_det':max_det}, indent=1))
print("\nRETRACTED:")
for d in retr_detail:
    print("  %-34s c<=%-4s TRUE=%s det=%s vdiff=%s"%(d['target_id'],d['our_c_upperbound'],
          d['true_min_crossing'],d['det'],('%.2e'%d['vol_diff']) if d['vol_diff'] is not None else 'NA'))
print("\nCOLLISION-ONLY (jones match, vol mismatch):")
for d in coll_detail:
    print("  %-34s c<=%s det=%s our_vol=%.4f nearMag=%.4f vdiff=%.3f cross=%s"%(d['target_id'],
          d['our_c_upperbound'],d['det'],d['our_vol'] or -1, d['nearest_maguire_vol'] or -1,
          d['vol_diff'] or -1, d['maguire_cross']))
print("\nUNRESOLVED no-poly (det achievable <=17, no corpus polynomial):",len(unres_detail))
for d in unres_detail[:60]:
    print("  %-30s det=%s c_ub=%s vol=%s"%(d['target_id'],d['det'],d['c_ub'],round(d['vol'],3) if d['vol'] else None))
