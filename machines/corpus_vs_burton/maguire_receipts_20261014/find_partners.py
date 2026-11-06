# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""For the five retractions whose stored Maguire DT is NOT isometric to our record (receipts.json, retry_five.json):
find the true census partner = a row of the same Maguire file with the same (mirror-canonical) Jones polynomial as the
stored DT (mutants share Jones) that IS isometric to our record. Writes partners.json."""
import site_config
import json, sys, signal
import snappy, spherogram
sys.path.insert(0,(site_config.MACHINE_YARD + '/corpus_vs_maguire/work')); import knotlib as K
ROOT=site_config.ROOT; MD=site_config.MACHINE_YARD+'/corpus_vs_maguire/maguire_data/'
db=json.load(open(site_config.PIPELINE_DIR+'/corpus_db.json'))
mg=json.load(open(site_config.PIPELINE_DIR+'/maguire_retractions.json'))['retractions']
five={'S_17_d247_Z2_chiral_mut':'pz4_2_3_5_7__vs__2_5_3_7_Kmut','conway_kt':'11_58__m127_conway_kt','fleet_n4_2_3_3_7':'fleet_n4_2_3_3_7','fleet_n4_2_3_5_5':'fleet_n4_2_3_5_5','fleet_n5_2_3_3_3_5':'fleet_n5_2_3_3_3_5'}
def iso(N, MM, tries=8):
    for i in range(tries):
        try:
            if N.is_isometric_to(MM): return True
        except Exception: pass
        N.randomize(); MM.randomize()
    return False
out={}
for key,folder in five.items():
    r=[x for x in db['knots'] if x.get('folder')==folder][0]
    N=snappy.Manifold(r['isosig']); our_vol=float(N.volume())
    cfile=mg[key]['true_min_crossing']; stored=mg[key]['maguire_dt']
    rows=list(K.stream_maguire(MD+cfile+'_hyp_jones.tar.xz'))
    jref=[K.canon_jones(K.parse_jones(p)) for d,p in rows if d==stored]
    if not jref: out[key]=dict(error='stored DT not in '+cfile); print(key,'stored DT not found'); continue
    cands=[d for d,p in rows if K.canon_jones(K.parse_jones(p))==jref[0]]
    res=[]
    for d in cands:
        L=spherogram.Link('DT:'+repr([tuple(K.decode_dt(d))])); MM=L.exterior()
        v=float(MM.volume()); same_vol=abs(v-our_vol)<1e-9
        ok=iso(N.copy(),MM) if same_vol else False
        res.append(dict(dt=d, vol=v, same_vol=same_vol, isometric=ok, sig=L.signature(), is_stored=(d==stored)))
        print(key, d, 'vol',round(v,9), 'iso',ok, 'sig',L.signature(), '(stored)' if d==stored else '', flush=True)
    hits=[x['dt'] for x in res if x['isometric']]
    out[key]=dict(folder=folder, file=cfile, stored_dt=stored, our_vol=our_vol, our_sig=r.get('sig'), same_jones_rows=len(cands), candidates=res, partner_dts=hits)
json.dump(out, open('partners.json','w'), indent=1)
print('DONE', {k:(v.get('partner_dts'), v.get('same_jones_rows')) for k,v in out.items()})
