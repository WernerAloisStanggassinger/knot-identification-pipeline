# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Retry the five non-isometric verdicts with retriangulation (is_isometric_to is a search, False is not a proof).
Also report whether the pair is a mutant situation (same HOMFLY) so a wrong-mutant match can be told apart."""
import site_config
import json, sys, signal
import snappy, spherogram
sys.path.insert(0,(site_config.MACHINE_YARD + '/corpus_vs_maguire/work')); import knotlib as K
ROOT=site_config.ROOT
db=json.load(open(site_config.PIPELINE_DIR + '/corpus_db.json'))
mg=json.load(open(site_config.PIPELINE_DIR + '/maguire_retractions.json'))['retractions']
five=['S_17_d247_Z2_chiral_mut','conway_kt','fleet_n4_2_3_3_7','fleet_n4_2_3_5_5','fleet_n5_2_3_3_3_5']
folders={'S_17_d247_Z2_chiral_mut':'pz4_2_3_5_7__vs__2_5_3_7_Kmut','conway_kt':'11_58__m127_conway_kt','fleet_n4_2_3_3_7':'fleet_n4_2_3_3_7','fleet_n4_2_3_5_5':'fleet_n4_2_3_5_5','fleet_n5_2_3_3_3_5':'fleet_n5_2_3_3_3_5'}
out={}
for key in five:
    r=[x for x in db['knots'] if x.get('folder')==folders[key]][0]
    num=K.decode_dt(mg[key]['maguire_dt']); L=spherogram.Link('DT:'+repr([tuple(num)])); MM=L.exterior()
    N=snappy.Manifold(r['isosig'])
    res=None; tries=0
    for i in range(20):
        tries=i+1
        try:
            if N.is_isometric_to(MM): res=True; break
        except Exception as ex: pass
        N.randomize(); MM.randomize()
    # second route: canonical triangulation + isometry signature
    try:
        s1=N.isometry_signature(); s2=MM.isometry_signature(); sig_eq=(s1==s2)
    except Exception as ex: s1=s2=None; sig_eq='err '+repr(ex)[:60]
    out[key]=dict(folder=folders[key], isometric_after_retriangulation=res, tries=tries, isometry_signature_equal=sig_eq, our_isosig_sig=s1, maguire_isosig_sig=s2,
                  stored_sig=r.get('sig'), maguire_sig=L.signature(), vol_ours=float(N.volume()), vol_maguire=float(MM.volume()),
                  record_knotinfo_id=r.get('knotinfo_id'), record_crosscheck=r.get('crosscheck_knotinfo'), record_distinct=r.get('distinct_knot_id'), retraction=mg[key])
    print(key, res, tries, 'isosig_equal=',sig_eq, 'stored_sig',r.get('sig'),'maguire_sig',L.signature(), flush=True)
json.dump(out, open('retry_five.json','w'), indent=1); print('DONE')
