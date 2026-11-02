# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Σn linking-form chirality witness, n = 3..6 (01.11.2026). Σn(K) = n-fold cyclic branched cover: SnapPy n-fold cyclic cover of the
# orbifold obtained by (n,0)-filling the exterior (see cover_inv); converted to Regina; Regina 7.4 HomologicalData → Kawauchi–Kojima invariants of the torsion
# linking form (torsionRankVector, torsionSigmaVector [2-primary], torsionLegendreSymbolVector [odd p]) = invariants of the ORIENTED Σn.
# Σn(mK) = -Σn(K). K amphichiral ⟹ Σn(K) ≅ Σn(mK) orientation-preservingly ⟹ equal invariants. Differ (for some n) ⟹ chiral.
# Σn(K) and Σn(mK) are built separately from K's and mK's own exteriors. (Viro: invariant under mutation for n = 2; for n ≥ 3 we do not
# rely on mutation invariance.)
# FILLING CHECK: H1(Σn) from the cover must equal coker(Γ^n − (Γ−I)^n), Γ = (V^T − V)^{-1} V^T (Seifert matrix, Sage) — else result void.
# Counter-cases: amphichiral 4_1, 8_17, 10_99 must be equal for all n; chiral 3_1, 5_2, 9_42, 10_125 are reported as measured.
import sys, json, bz2, csv, time
import sage.all
from sage.all import ZZ, QQ, matrix, identity_matrix
import snappy, regina
dt = lambda s: [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]
def seifert_h1(L, n):
    V = L.sage_link().seifert_matrix().change_ring(QQ); A = V.transpose() - V
    G = A.inverse() * V.transpose(); assert all(x in ZZ for x in G.list())
    P = (G ** n - (G - identity_matrix(QQ, G.nrows())) ** n).change_ring(ZZ)
    return sorted(int(abs(d)) for d in P.elementary_divisors() if abs(d) != 1)
def cover_inv(E, n):
    # orbifold route: fill the exterior
    # with (n,0), take the n-fold cyclic cover of the orbifold — SnapPy lifts the filling itself.
    F = E.copy(); F.dehn_fill((n, 0)); C = [c for c in F.covers(n, cover_type='cyclic')]; assert len(C) == 1; C = C[0]
    assert C.num_cusps() == 1 and tuple(round(x) for x in C.cusp_info('filling')[0]) == (1, 0); T = C.filled_triangulation()
    R = regina.Triangulation3(T._to_string()); assert R.isOrientable() and R.isClosed()
    H = regina.HomologicalData(R)
    h1 = R.homology()
    tors = sorted(int(str(h1.invariantFactor(i))) for i in range(h1.countInvariantFactors()))
    return dict(rank=int(h1.rank()), torsion=tors, rk=H.torsionRankVectorString(), sigma=H.torsionSigmaVectorString(), legendre=H.torsionLegendreSymbolVectorString())
def test(L, ns=(3, 4, 5, 6), check=True):
    E = L.exterior(); Em = L.mirror().exterior(); out = {}; differ = []
    for n in ns:
        a = cover_inv(E, n); b = cover_inv(Em, n)
        if check:
            s = seifert_h1(L, n)
            if a['rank'] == 0: assert a['torsion'] == s, ('FILLING CHECK FAILED', n, a['torsion'], s)
        out[n] = dict(K=a, mK=b, equal=(a == b))
        if a != b: differ.append(n)
    return dict(differ_at=differ, detail=out, verdict=('CHIRAL (Σn linking form ≠ mirror at n=%s)' % differ) if differ else 'Σ3..Σ6 blind')
def from_dt(code): return snappy.Link('DT:' + str([tuple(dt(code))]))

if __name__ == '__main__':
    import os as _os
    sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
    from census_codes import main
    counter = [(nm, (lambda nm=nm: snappy.Link(nm)), (lambda r, a=a: (not r['differ_at']) if a else True)) for nm, a in
               (('4_1', True), ('8_17', True), ('10_99', True), ('3_1', False), ('5_2', False), ('9_42', False), ('10_125', False))]
    main(test, counter, from_dt)
