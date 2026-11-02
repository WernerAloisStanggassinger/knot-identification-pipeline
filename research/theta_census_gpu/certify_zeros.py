# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Rigorous chirality / amphichirality certificate for HYPERBOLIC census knots with theta = 0 at (22/7, 21/13) (first written for the
# 352 such knots with c <= 14 and 15nh_008268 of our theta deposit, DOI 10.5281/zenodo.23856098). Different signatures => chiral: by
# Mostow-Prasad rigidity every self-homeomorphism of a hyperbolic exterior is homotopic to an isometry. SnapPy isometry_signature(verified=True,
# ignore_orientation=False) is a complete invariant of the oriented cusped hyperbolic manifold (canonical cell decomposition certified by
# interval arithmetic, or by exact arithmetic where it is not tetrahedral — SnapPy switches automatically). If sig(M) == sig(-M), the exterior has an orientation-reversing isometry, hence an orientation-reversing
# self-homeomorphism; by Gordon–Luecke (a homeomorphism of knot exteriors takes meridian to meridian, so it extends over the solid torus)
# the knot is equivalent to its mirror image, i.e. amphichiral. Counter-cases first: 4_1 must give equal, 5_2 different signatures.
# Usage: python certify_zeros.py <zeros.tsv> <census_dir> <out.jsonl> <workers> [<file pattern>]
#   zeros.tsv: rows name<TAB>…<TAB>0 for θ-zeros; codes are read from <census_dir>/<file pattern> (Burton census csv.bz2 with
#   columns name, dt_code; default '*.csv.bz2'; the authors' runs used one small candidates file per crossing number).
import sys, csv, bz2, json, time, glob, os
from multiprocessing import Pool
import snappy
def dt(s): return [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]
def cert(item):
    name, M = item
    t0 = time.time()
    if isinstance(M, str): M = snappy.Manifold('DT:' + str([tuple(dt(M))]))
    r = dict(name=name); errs = []
    # attempt 0 on the given triangulation; on failure (e.g. FindExactShapesError) retry on randomized triangulations of the same
    # manifold — the signature is an invariant of the manifold, so any triangulation that certifies gives the same answer
    for k in range(0, 21):
        try:
            N = M.copy()
            for _ in range(k): N.randomize()
            Nm = N.copy(); Nm.reverse_orientation()
            s = N.isometry_signature(verified=True, ignore_orientation=False)
            sm = Nm.isometry_signature(verified=True, ignore_orientation=False)
            r.update(sig=s, sig_mirror=sm, oriented_equal=(s == sm), randomizations=k, earlier_errors=errs,
                     verdict='AMPHICHIRAL (verified: oriented signatures equal)' if s == sm else 'CHIRAL (verified: oriented signatures differ)')
            break
        except Exception as e:
            errs.append(f'{k}: {type(e).__name__}')
    else:
        r.update(error='; '.join(errs), verdict='NOT CERTIFIED')
    r['seconds'] = round(time.time() - t0, 1)
    return r
if __name__ == '__main__':
    tsv, raw, out, w = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    pattern = sys.argv[5] if len(sys.argv) > 5 else '*.csv.bz2'
    for nm, exp in (('4_1', True), ('5_2', False)):
        r = cert((nm, snappy.Manifold(nm))); assert r.get('oriented_equal') is exp, r
        print('counter-case', nm, r['verdict'], flush=True)
    zeros = [row[0] for row in csv.reader(open(tsv), delimiter='\t') if row and row[2] == '0']
    need = set(zeros); codes = {}
    for f in sorted(glob.glob(os.path.join(raw, pattern))):
        for row in csv.DictReader(bz2.open(f, 'rt')):
            if row['name'] in need: codes[row['name']] = row['dt_code']
        if len(codes) == len(need): break
    assert len(codes) == len(need), sorted(need - set(codes))[:5]
    print('theta-zeros to certify:', len(zeros), flush=True)
    res = []
    with Pool(w) as p, open(out, 'w') as fo:
        for r in p.imap_unordered(cert, [(n, codes[n]) for n in zeros]):
            fo.write(json.dumps(r) + '\n'); fo.flush(); res.append(r)
    from collections import Counter
    print('verdicts:', dict(Counter(r['verdict'] for r in res)), flush=True)
