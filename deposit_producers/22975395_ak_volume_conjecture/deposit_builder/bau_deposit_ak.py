# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Builds the Zenodo deposit 'Andersen-Kashaev volume conjecture for 213 knots outside the HT census' (26.09.2026).
Source: maschinen_ideen/QG_PRIORART_20260925/famed/deposit_frisch.jsonl (computed fresh on 26.09.).
Fresh HERE: HT exclusion (identify on ordered + 5 geometric retriangulations) for all without c_lb>=15;
Isometry duplicates (volume intervals overlap -> is_isometric_to) + orientation on PD complements.
Checking code for reproduce.py is cut VERBATIM from famed_unabhaengig.py / meridian_q.py (ast), never retyped.
Runs with the sage_stable Python. Writes HOCHLADEN/{README.md, ak_data.json, reproduce.py, SHA256SUMS}, .zenodo.json,
METADATEN_ZUM_EINFUEGEN.txt. Aborts (assert) if any knot fails a check."""
import site_config
import os, sys, json, ast, re, hashlib, signal, resource, datetime
resource.setrlimit(resource.RLIMIT_AS, (8000 << 20, 8000 << 20))
HIER = os.path.dirname(os.path.abspath(__file__))
FAMED = (site_config.ROOT + '/maschinen_ideen/QG_PRIORART_20260925/famed')
UP = os.path.join(HIER, 'HOCHLADEN'); os.makedirs(UP, exist_ok=True)
import snappy, regina, sage.version
class Frist(Exception): pass
def _alarm(*a): raise Frist()
signal.signal(signal.SIGALRM, _alarm)

R = [json.loads(l) for l in open(os.path.join(FAMED, 'deposit_frisch.jsonl'))]
R = {r['did']: r for r in R}
N_ERWARTET = int(sys.argv[1]) if len(sys.argv) > 1 else 213   # v2 (27.09.2026): 220, counted against deposit_frisch.jsonl
assert len(R) == N_ERWARTET, len(R)
for d, r in R.items():
    assert 'fehler' not in r and 'q_fehler' not in r, (d, r.get('fehler'), r.get('q_fehler'))
    assert r['famed'] and r['geometric_verified'] and r['ordered'], d
    assert r['n1'] == '0' and r['eindeutig'] is True and r['mu'] in ('1', '-1'), (d, r['n1'], r['mu'])
    assert r['H1'] == 'Z', d

# --- families, sculptures ---
db = json.load(open((site_config.ROOT + '/freeze_manifest/pipeline_scripts/corpus_db.json')))
skulptur = {}
for k in db['knots']:
    if k.get('distinct_knot_id') in R and 'own_knots' in str(k.get('pflicht_path', '')):
        skulptur.setdefault(k['distinct_knot_id'], []).append(k['stem'])
del db
STEEL = {'trainstation', 'number-5', 'xyz-free-view', 'double-3', 'celtic-trefoil'}
def familie(d):
    m = re.fullmatch(r'T_(\d+)_(\d+)_2_(\d+)', d)
    if m: return f'twisted torus knot T({m[1]},{m[2]};2,{m[3]})'
    m = re.fullmatch(r'pz_(\d+)_(\d+)_(\d+)', d)
    if m: return f'pretzel knot P({m[1]},{m[2]},{m[3]})'
    return 'other'

# --- HT exclusion, fresh ---
H = snappy.HTLinkExteriors
for probe in ['K14n1207', 'K12n706']:
    M = H[probe]; M.randomize(); assert H.identify(M), probe   # counter-check: hyperbolic HT knots must be found
for d, r in R.items():
    if (r['c_lb'] or 0) >= 15:
        r['outside_HT'] = dict(method='proven crossing number lower bound', c_lower_bound=r['c_lb']); continue
    M = snappy.Manifold(regina.Triangulation3.tightDecoding(r['ordered_triangulation_tight'])); M.set_name('k')
    geo = 0; hits = set()
    for i in range(30):                      # until 5 geometric ones seen, at most 30 attempts
        if geo >= 5: break
        if i: M.randomize()
        if M.solution_type() != 'all tetrahedra positively oriented': continue
        geo += 1; signal.alarm(300)
        h = H.identify(M) or H.identify(M, extends_to_link=False); signal.alarm(0)
        if h: hits.add(str(h))
    assert not hits, (d, hits)
    assert geo >= 3, (d, geo)
    r['outside_HT'] = dict(method='SnapPy HTLinkExteriors.identify on geometric triangulations: no match', geometric_triangulations_tried=geo)
    print('HT', d, geo, flush=True)

# --- duplicates, fresh ---
from sage.all import RealIntervalField
RIF = RealIntervalField(212)
vol = {d: RIF(r['volume_interval']) for d, r in R.items()}   # Sage reads the ?-notation as an interval
ds = sorted(R)
paare = []
for i, a in enumerate(ds):
    for b in ds[i+1:]:
        if not vol[a].overlaps(vol[b]): continue
        A = snappy.Manifold(R[a]['isosig']); B = snappy.Manifold(R[b]['isosig'])
        signal.alarm(600); iso = A.is_isometric_to(B); signal.alarm(0)
        if not iso: continue
        LA, LB = snappy.Link(R[a]['pd']), snappy.Link(R[b]['pd'])
        EA, EB = LA.exterior(), LB.exterior()
        dets = sorted({int(x.cusp_maps()[0].det()) for x in EA.is_isometric_to(EB, return_isometries=True)})
        paare.append(dict(knot_a=a, knot_b=b, isometric_complements=True, cusp_map_determinants=dets,
                          chern_simons_a=float(EA.chern_simons()), chern_simons_b=float(EB.chern_simons()),
                          verdict='mirror pair (orientation-reversing isometry)' if dets == [-1] else 'SAME KNOT?'))
        print('PAAR', paare[-1], flush=True)
assert all(p['verdict'].startswith('mirror') for p in paare), paare
L = snappy.Link(R[paare[0]['knot_a']]['pd'])   # counter-check of orientation
gp = (sorted({int(x.cusp_maps()[0].det()) for x in L.exterior().is_isometric_to(L.mirror().exterior(), return_isometries=True)}),
      sorted({int(x.cusp_maps()[0].det()) for x in L.exterior().is_isometric_to(L.copy().exterior(), return_isometries=True)}))
assert gp == ([-1], [1]), gp
n_kompl = len(R) - len(paare)

# --- data ---
knots = []
for d in ds:
    r = R[d]
    knots.append(dict(
        id=d, family=familie(d), sculpture_records=sorted(skulptur.get(d, [])) or None,
        crossing_number_lower_bound=r['c_lb'], diagram_crossings=len(r['pd']), outside_HT_census=r['outside_HT'],
        ordered_triangulation=dict(regina_tight_encoding=r['ordered_triangulation_tight'], tetrahedra=r['tet'],
                                   isosig_unordered_for_reference=r['ordered_triangulation_isosig'],
                                   origin='authors\' search code (A. Guilloux, FAMEDexploration) + ordering; verdict below is ours'),
        famed=dict(verdict=r['famed'], checker='independent implementation of Ben Aribi-Wong Def. 1.1 (1)-(4), exact over Q', global_sign_of_E=r['famed_note'].split()[-1]),
        geometric=dict(verified=r['geometric_verified'], method='SnapPy verify_hyperbolicity(bits_prec=212) on the same ordered triangulation',
                       min_imaginary_part_lower_bound=r['min_im_shape_lower']),
        volume_verified_interval=r['volume_interval'], H1=r['H1'],
        meridian=dict(n1=r['n1'], mu=r['mu'], unique=r['eindeutig'], meridian_in_triangulation_basis=r['meridian'],
                      homological_longitude=r['long'], isometries_used=r['n_isos'], retriangulation_attempt=r['versuch']),
        diagram=dict(pd=r['pd'], source='project corpus (knot-structures.stainlesssteel4u.de)')))
fam_n = {}
for k in knots:
    f = k['family'].split(' T(')[0].split(' P(')[0]; fam_n[f] = fam_n.get(f, 0) + 1
tet = [k['ordered_triangulation']['tetrahedra'] for k in knots]
n_clb = sum(1 for k in knots if k['outside_HT_census']['method'].startswith('proven'))
tt_clb = [k['crossing_number_lower_bound'] for k in knots if k['family'].startswith('twisted')]
skul = sorted(s for k in knots for s in (k['sculpture_records'] or []))
meta = dict(
    title=f'The Andersen-Kashaev volume conjecture for {len(knots)} knots outside the Hoste-Thistlethwaite census',
    date=datetime.date.today().isoformat(), knots=len(knots), distinct_complements=n_kompl, mirror_pairs=paare,
    orientation_counter_check=dict(link_vs_mirror=gp[0], link_vs_itself=gp[1]),
    counts=dict(outside_HT_by_proven_c_ge_15=n_clb, outside_HT_by_identify=len(knots) - n_clb,
                tetrahedra_min=min(tet), tetrahedra_max=max(tet), tetrahedra_gt_23=sum(t > 23 for t in tet),
                families=fam_n, twisted_torus_c_lower_bound_min=min(tt_clb), twisted_torus_c_lower_bound_max=max(tt_clb),
                n1_zero=sum(k['meridian']['n1'] == '0' for k in knots)),
    sculpture_records=skul, sculpture_records_in_steel=sorted(set(skul) & STEEL),
    tools=dict(snappy=snappy.__version__, regina=regina.versionString(), sage=sage.version.version),
    theorem_used='F. Ben Aribi, K. H. Wong, arXiv:2410.10776v2, Corollary 1.9: "Conjecture 1.5 is true for all FAMED geometric triangulations"',
    census_covered_before='F. Ben Aribi, A. Guilloux, K. H. Wong, arXiv:2512.17437v1, Theorem 3.1 (HT census, <=14 crossings, <=23 tetrahedra, all but 6 of 42 123)')
json.dump(dict(meta=meta, knots=knots), open(os.path.join(UP, 'ak_data.json'), 'w'), indent=1, ensure_ascii=False)
json.dump(meta, open(os.path.join(HIER, 'meta_fuer_readme.json'), 'w'), indent=1, ensure_ascii=False)

# --- reproduce.py: cut the checking code verbatim ---
def schneide(pfad, namen):
    src = open(pfad).read(); t = ast.parse(src); zeilen = src.splitlines()
    out = []
    for n in t.body:
        if isinstance(n, ast.FunctionDef) and n.name in namen:
            out.append('\n'.join(zeilen[n.lineno - 1:n.end_lineno]))
    assert len(out) == len(namen), (pfad, namen, len(out))
    return '\n\n'.join(out)
code_famed = schneide(os.path.join(FAMED, 'famed_unabhaengig.py'), ['eigen_famed'])
code_q = schneide(os.path.join(FAMED, 'meridian_q.py'), ['nz_zz2', 'q_fuer'])
code_st = schneide(os.path.join(FAMED, 'gegenfall_4.py'), ['mit_stoerung'])
vorlage = open(os.path.join(HIER, 'reproduce_vorlage.py')).read()
assert vorlage.count('#<<EIGEN_FAMED>>') == 1 and vorlage.count('#<<Q_FUER>>') == 1 and vorlage.count('#<<MIT_STOERUNG>>') == 1
open(os.path.join(UP, 'reproduce.py'), 'w').write(vorlage.replace('#<<EIGEN_FAMED>>', code_famed).replace('#<<Q_FUER>>', code_q).replace('#<<MIT_STOERUNG>>', code_st))
print('BUILD OK', json.dumps(meta['counts']), 'complements', n_kompl, 'pairs', len(paare), flush=True)
