# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Writhe witness for alternating knots, from the authors' writhe_alt.py (same rule, now a module with the shared driver).
# Murasugi / Thistlethwaite: the writhe of a REDUCED alternating diagram is a knot invariant, and
# the mirror has writhe -w  ⟹  w != 0 ⟹ chiral. Checked per knot: the diagram is alternating and reduced, i.e. has no nugatory crossing
# (census_codes.is_reduced: every chord of the Gauss word is interlaced with another one). Counter-cases: 4_1 (w = 0), 3_1 (|w| = 3),
# and a 7-crossing alternating but NOT reduced diagram of the amphichiral square knot 3_1 # m3_1 (writhe 1), which must give no verdict.
# (The authors' original writhe_alt.py tested reducedness with SnapPy simplify('basic'), which misses that case.)
import sys, os, snappy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from census_codes import is_reduced, SQUARE_KNOT_7
dt = lambda s: [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]
def test(L):
    n = len(L.crossings); alt = L.is_alternating(); red = is_reduced(L.PD_code()); w = L.writhe()
    usable = alt and red
    return dict(crossings=n, alternating_diagram=alt, reduced=red, writhe=w,
                verdict=('CHIRAL (writhe of a reduced alternating diagram != 0)' if usable and w != 0 else
                         ('writhe 0 (no information)' if usable else 'not a reduced alternating diagram (no information)')))
if __name__ == '__main__':
    import os as _os
    sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
    from census_codes import main
    from_dt = lambda code: snappy.Link('DT:' + str([tuple(dt(code))]))
    def square7():
        import regina
        return snappy.Link([tuple(x - 1 for x in c) for c in regina.Link.fromSignedGauss(SQUARE_KNOT_7).pdData()])
    counter = [('4_1', (lambda: snappy.Link('4_1')), (lambda r: r['writhe'] == 0 and r['alternating_diagram'])),
               ('3_1', (lambda: snappy.Link('3_1')), (lambda r: abs(r['writhe']) == 3 and r['verdict'].startswith('CHIRAL'))),
               ('square knot, 7-crossing non-reduced diagram', (lambda: square7()), (lambda r: not r['reduced'] and not r['verdict'].startswith('CHIRAL')))]
    main(test, counter, from_dt)
