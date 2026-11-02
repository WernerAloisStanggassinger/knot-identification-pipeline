# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Chirality certificate v2 for hyperbolic knots (02.11.2026). One knot per call (run under timeout + memory latch), result = one JSON line.
#
# AMPHICHIRAL (exact): P = SnapPea proto-canonical triangulation of the exterior M (canonize(), optionally after randomize()), read into Regina
#   and oriented; P has a combinatorial automorphism with odd facet permutations (orientation-reversing) ⟹ the manifold of P has an
#   orientation-reversing self-homeomorphism. P ≅ M is shown EXACTLY, independent of SnapPea: Regina simplify() of P and of M, then a common
#   isomorphism signature in their height-1 Pachner neighbourhoods (height 2 from M as fallback). Gordon–Luecke: meridian ↦ meridian, so
#   K is isotopic to its mirror image.
#   Why this route: an earlier version used the canonical RE-triangulation, which has finite vertices; its homeomorphism check to M
#   exploded at height 2. The proto-canonical one is ideal. homeo() searches M's simplified triangulation B at retriangulation
#   height hB = 1, then 2 (height = extra tetrahedra over B, Regina 2-3/3-2 moves), and P's simplified triangulation A at height 1;
#   the log field homeo_pachner_distance_bound is 0 (equal after simplify), hB, or hB + 1 (A's height added) — not a move count.
# CHIRAL (exact, verified): canonical_retriangulation(verified=True) by interval arithmetic only (exact shapes are out of reach for these
#   knots the authors tried; no log of that probe is part of this release) succeeds and has NO orientation-reversing automorphism ⟹ no
#   orientation-reversing isometry (the canonical decomposition is isometry-invariant) ⟹ by Mostow no orientation-reversing homeomorphism
#   ⟹ chiral.
# Otherwise OPEN. Counter-cases: 4_1 amphichiral; 5_2 and 15nh_048139 (census; chiral) must not come out AMPHICHIRAL.
import sys, json, bz2, csv, time
import snappy, regina
dt = lambda s: [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]
def autos(t):
    t = regina.Triangulation3(t); assert t.isOrientable(); t.orient(); n = [0, 0]
    def act(iso):
        s = {iso.facetPerm(i).sign() for i in range(t.size())}; assert len(s) == 1; n[0 if s == {1} else 1] += 1; return False
    t.findAllIsomorphisms(t, act); return n
def homeo(P, M):
    A = regina.Triangulation3.fromSnapPea(P._to_string()); A.simplify(); B = regina.Triangulation3.fromSnapPea(M._to_string()); B.simplify()
    if A.isoSig() == B.isoSig(): return 0
    for hB in (1, 2):
        S = {B.isoSig()}; B.retriangulate(hB, 4, lambda s, tr: (S.add(s), False)[1])
        if A.isoSig() in S: return hB
        h = []; A.retriangulate(1, 4, lambda s, tr: (h.append(s), True)[1] if s in S else False)
        if h: return hB + 1
    return None
def certify(M, tries=8):
    r = dict(tets=M.num_tetrahedra())
    for k in range(tries):
        P = M.copy()
        if k: P.randomize()
        P.canonize(); a = autos(regina.Triangulation3.fromSnapPea(P._to_string()))
        if a[1]:
            h = homeo(P, M)
            if h is not None:
                r.update(proto_tets=P.num_tetrahedra(), automorphisms=a, try_=k, homeo_pachner_distance_bound=h,
                         verdict='AMPHICHIRAL (exact: orientation-reversing automorphism of a triangulation shown homeomorphic to M)'); return r
            r.setdefault('notes', []).append(f'try {k}: reversing automorphism, homeomorphism to M not shown')
    try:
        C = M.canonical_retriangulation(verified=True, exact_bits_prec_and_degrees=None)
        a = autos(regina.Triangulation3.fromSnapPea(C._to_string())); r.update(verified_canonical_tets=C.num_tetrahedra(), verified_canonical_automorphisms=a)
        r['verdict'] = 'CHIRAL (exact: verified canonical retriangulation has no orientation-reversing automorphism)' if a[1] == 0 else \
                       'AMPHICHIRAL (verified canonical retriangulation has an orientation-reversing automorphism)'
    except Exception as e:
        r['verified_canonical_error'] = repr(e)[:120]; r['verdict'] = 'OPEN'
    return r

if __name__ == '__main__':
    import os as _os
    sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
    from census_codes import main
    from census_codes import codes
    test = certify; from_dt = lambda code: snappy.Link('DT:' + str([tuple(dt(code))])).exterior()
    ok_ch = lambda r: (not r['verdict'].startswith('AMPHICHIRAL')) and (r['verdict'].startswith('CHIRAL') or r['verdict'] == 'OPEN')
    counter = [('4_1', (lambda: snappy.Manifold('4_1')), (lambda r: r['verdict'].startswith('AMPHICHIRAL'))),
               ('5_2', (lambda: snappy.Manifold('5_2')), ok_ch),
               ('15nh_048139', (lambda: from_dt(codes(['15nh_048139'])['15nh_048139'])), ok_ch)]
    main(test, counter, from_dt)
