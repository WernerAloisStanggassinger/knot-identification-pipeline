# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Chern–Simons chirality witness (02.11.2026). SnapPy complex_volume(verified_modulo_2_torsion=True): vol + i*CS verified by interval
# arithmetic, CS modulo pi^2/2. Orientation reversal negates CS, so amphichiral => CS = -CS mod pi^2/2 => CS in {0, pi^2/4} mod pi^2/2.
# Verified interval of CS mod pi^2/2 excluding both 0 and pi^2/4 => CHIRAL (exact). Otherwise no information.
# Counter-cases (must hold): 4_1, 8_17 amphichiral -> no info; 5_2, 9_42 chiral (9_42 HOMFLY-blind) -> CHIRAL.
import sys, json, bz2, csv, os, snappy
from sage.all import RealIntervalField, pi
dt = lambda s: [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]
def test(M, prec=212):
    z = M.complex_volume(verified_modulo_2_torsion=True, bits_prec=prec); R = RealIntervalField(prec)
    m = R(pi) ** 2 / 2; cs = z.imag(); k = (cs / m).floor() if (cs / m).floor().is_exact() else None
    # distance of CS to the lattice (pi^2/4) Z: chiral iff the interval CS/(pi^2/4) contains no integer
    q = cs / (m / 2); contains_int = q.floor().upper() != q.ceil().lower() or q.contains_integer() if hasattr(q, 'contains_integer') else True
    lo, hi = q.lower().floor(), q.upper().ceil()
    hit = any(R(n) in q for n in range(int(lo), int(hi) + 1))
    return dict(cs=str(cs), cs_over_pi2_4=str(q), chiral=not hit, verdict='CHIRAL (CS not in {0, pi^2/4} mod pi^2/2)' if not hit else 'CS consistent with amphichiral')

if __name__ == '__main__':
    import os as _os
    sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
    from census_codes import main
    from_dt = lambda code: snappy.Link('DT:' + str([tuple(dt(code))])).exterior()
    counter = [(nm, (lambda nm=nm: snappy.Manifold(nm)), (lambda r, e=e: r['chiral'] is e)) for nm, e in
               (('4_1', False), ('8_17', False), ('5_2', True), ('9_42', True))]
    main(test, counter, from_dt)
