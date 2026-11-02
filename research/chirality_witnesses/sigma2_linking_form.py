# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Linking-form chirality witness (01.11.2026). Double branched cover Σ2(K): H1 = coker(B), B = V + V^T (V Seifert matrix, Sage), linking
# form λ(x,y) = x^T B^{-1} y mod 1 (one fixed sign convention, used for K and mK alike). Mirror: V -> -V^T ⟹ λ(mK) = -λ(K).
# K amphichiral ⟹ Σ2(K) has an orientation-reversing self-homeomorphism ⟹ λ ≅ -λ ⟹ the multiset {q(x) = λ(x,x) : x ∈ H1} equals
# {-q(x)}. Multisets differ ⟹ chiral. Exact rational arithmetic; independent of SnapPy's geometry and of Jones/HOMFLY/Kh/HFK.
# Counter-cases: 3_1 (det 3) differ; 4_1, 8_17, 10_99, 8_18 (amphichiral) equal.
import sys, json, bz2, csv, itertools, collections
import sage.all
from sage.all import QQ, ZZ, matrix
import snappy
dt = lambda s: [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]
def qdist(L):
    V = L.sage_link().seifert_matrix(); B = (V + V.transpose()).change_ring(ZZ)
    D, U, W = B.smith_form()                       # U*B*W = D
    d = [D[i, i] for i in range(D.nrows())]
    Ui = U.inverse(); Binv = B.change_ring(QQ).inverse()
    idx = [i for i in range(len(d)) if abs(d[i]) > 1]
    gens = [Ui.column(i) for i in idx]; ords = [int(abs(d[i])) for i in idx]
    Gm = matrix(QQ, [[(gens[a] * Binv * gens[b]) for b in range(len(gens))] for a in range(len(gens))])
    cnt = collections.Counter()
    for x in itertools.product(*[range(o) for o in ords]):
        v = matrix(QQ, [x]); q = (v * Gm * v.transpose())[0, 0]; q = q - q.floor(); cnt[q] += 1
    neg = collections.Counter({((-q) - (-q).floor()): c for q, c in cnt.items()})
    return ords, cnt == neg, int(abs(B.det()))
def test(L):
    ords, sym, det = qdist(L)
    return dict(H1=ords, det=det, symmetric=sym, verdict='CHIRAL (linking form of Σ2 not ≅ its negative)' if not sym else 'linking form blind')
def from_dt(code): return snappy.Link('DT:' + str([tuple(dt(code))]))

if __name__ == '__main__':
    import os as _os
    sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
    from census_codes import main
    counter = [(nm, (lambda nm=nm: snappy.Link(nm)), (lambda r, s=s: r['symmetric'] is s)) for nm, s in
               (('3_1', False), ('4_1', True), ('8_17', True), ('10_99', True), ('8_18', True))]
    main(test, counter, from_dt)
