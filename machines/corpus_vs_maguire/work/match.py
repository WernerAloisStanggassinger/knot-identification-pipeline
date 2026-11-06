# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
import site_config
import json, os, glob, time, importlib
import knotlib as K
MD=(site_config.MACHINE_YARD + '/corpus_vs_maguire/maguire_data')
WORK=(site_config.MACHINE_YARD + '/corpus_vs_maguire/work')

targets=json.load(open(WORK+'/targets_novelhyp.json'))
withj=[r for r in targets if r.get('jones')]
# fingerprint: canon_jones -> list of target records
fp={}
for r in withj:
    cj=K.canon_jones(K.parse_jones(r['jones']))
    fp.setdefault(cj,[]).append(r['distinct_knot_id'])
print("fingerprints (canon jones) from %d with-jones targets: %d distinct"%(len(withj),len(fp)),flush=True)

target_dets={ (r.get('det') or -1) for r in targets }

# files: c=4..17 hyp jones, both a and n
files=sorted(glob.glob(MD+'/*_hyp_jones.tar.xz'))
candidates=[]         # (target_ids, crossfile, dt, maguire_poly)
maguire_det_set=set()
max_det_a=0
t0=time.time(); nrows=0
for p in files:
    base=os.path.basename(p)
    cross=base.split('_')[0]      # e.g. '17n'
    is_a = 'a' in cross
    cnt=0
    for dt,poly in K.stream_maguire(p):
        nrows+=1; cnt+=1
        pj=K.parse_jones(poly)
        d=K.det_from_jones(pj)
        maguire_det_set.add(d)
        if is_a and d>max_det_a: max_det_a=d
        cj=K.canon_jones(pj)
        if cj in fp:
            candidates.append({'target_ids':fp[cj],'cross':cross,'dt':dt,'maguire_jones':poly,'maguire_det':d})
        if cnt % 500000==0:
            print("  %s: %d rows (total %d, %.0fs) cands=%d"%(cross,cnt,nrows,time.time()-t0,len(candidates)),flush=True)
    print("DONE %s: %d rows (total %d, %.0fs) cands=%d maxdet_a=%d"%(cross,cnt,nrows,time.time()-t0,len(candidates),max_det_a),flush=True)
    # checkpoint
    json.dump({'candidates':candidates,'max_det_a':max_det_a,'ndet':len(maguire_det_set),
               'files_done':p,'nrows':nrows},
              open(WORK+'/match_checkpoint.json','w'))

out={'candidates':candidates,'max_det_a':max_det_a,
     'maguire_det_set_size':len(maguire_det_set),
     'nrows':nrows,
     'maguire_det_sorted_top':sorted(maguire_det_set)[-5:],
     # which target dets are achievable in maguire <=17
     'target_dets_in_maguire':sorted(d for d in target_dets if d in maguire_det_set)}
json.dump(out, open(WORK+'/match_raw.json','w'), indent=1)
# also persist full det set for no-jones adjudication
json.dump(sorted(maguire_det_set), open(WORK+'/maguire_det_set.json','w'))
print("MATCH DONE. rows=%d candidates=%d max_det_a=%d det_set=%d"%(nrows,len(candidates),max_det_a,len(maguire_det_set)),flush=True)
