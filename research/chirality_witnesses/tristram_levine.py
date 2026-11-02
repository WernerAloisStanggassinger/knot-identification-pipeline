# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Tristram–Levine signature witness (01.11.2026). sigma_w(K) = sign((1-w)V + (1-conj w)V^T), V = Seifert matrix (Sage Link.seifert_matrix).
# Mirror: V -> -V^T ⟹ sigma_w(mK) = -sigma_w(K); amphichiral ⟹ sigma_w ≡ 0 for all w on the unit circle.
# sigma_w is locally constant away from the unit-circle roots of Delta; we evaluate at the midpoint of every arc between consecutive
# roots (numeric, mpmath 60 digits; each evaluation records the smallest |eigenvalue| as margin — margin < 1e-20 ⟹ result NOT used).
# Any nonzero value with margin ⟹ chiral (numerical witness, independent of SnapPy's isometry signature and of Jones/HOMFLY/Kh).
# Counter-cases: 3_1 nonzero, 9_42 nonzero, 4_1 and 8_3 all zero.
import sys, json, bz2, csv
import sage.all
from sage.all import PolynomialRing, QQ, CC
import mpmath, snappy
mpmath.mp.dps = 60
dt = lambda s: [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]
def tl(S):
    V = S.seifert_matrix(); n = V.nrows()
    if n == 0: return [], 1.0
    Vm = mpmath.matrix([[int(V[i, j]) for j in range(n)] for i in range(n)])
    D = S.alexander_polynomial(); R = PolynomialRing(QQ, 't'); t = R.gen()
    num = R(D * t ** (-min(D.exponents()))) if hasattr(D, 'exponents') else R(D)
    roots = sorted({float(mpmath.arg(complex(r))) % (2 * mpmath.pi) for r in num.roots(CC, multiplicities=False) if abs(abs(complex(r)) - 1) < 1e-9})
    roots = [x for x in roots if 1e-12 < x < float(2 * mpmath.pi) - 1e-12]
    pts = [float(mpmath.pi)] if not roots else [ (a + b) / 2 for a, b in zip([0.0] + roots, roots + [float(2 * mpmath.pi)]) ]
    out, margin = [], 1e99
    for th in pts:
        w = mpmath.exp(1j * th); A = (1 - w) * Vm + (1 - mpmath.conj(w)) * Vm.T
        ev = mpmath.eighe(A, eigvals_only=True) if hasattr(mpmath, 'eighe') else mpmath.eigh(A, eigvals_only=True)
        margin = min(margin, min(abs(float(e)) for e in ev)); out.append((round(float(th), 6), int(sum(1 for e in ev if e > 0) - sum(1 for e in ev if e < 0))))
    return out, margin
def run(name, code=None):
    L = snappy.Link('DT:' + str([tuple(dt(code))])) if code else snappy.Link(name)
    vals, margin = tl(L.sage_link()); nz = [v for v in vals if v[1] != 0]
    ok = margin > 1e-20
    return dict(name=name, tl_values=vals, margin=float(margin), verdict=('CHIRAL (Tristram-Levine signature nonzero, numeric)' if nz and ok else
                ('TL all zero' if ok else 'MARGIN TOO SMALL')))

if __name__ == '__main__':
    import os as _os
    sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
    from census_codes import main
    def test(L):                       # same verdict rule as run() above, on a given link
        vals, margin = tl(L.sage_link()); nz = [v for v in vals if v[1] != 0]; ok = margin > 1e-20
        return dict(tl_values=vals, margin=float(margin), verdict=('CHIRAL (Tristram-Levine signature nonzero, numeric)' if nz and ok else
                    ('TL all zero' if ok else 'MARGIN TOO SMALL')))
    from_dt = lambda code: snappy.Link('DT:' + str([tuple(dt(code))]))
    counter = [(nm, (lambda nm=nm: snappy.Link(nm)), (lambda r, e=e: r['verdict'].startswith('CHIRAL') is e)) for nm, e in
               (('3_1', True), ('9_42', True), ('4_1', False), ('8_3', False))]
    main(test, counter, from_dt)
