# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Knot Floer homology chirality witness (01.11.2026): SnapPy knot_floer_homology (Szabó's program) for K and for its mirror, computed
# separately from each diagram. Amphichiral ⟹ HFK(K) ≅ HFK(mK) ⟹ same ranks dict, same tau, same epsilon. Differ ⟹ chiral.
# Counter-cases: 3_1 differ, 4_1 / 8_17 / 10_17 equal.
import sys, json, bz2, csv, glob
import snappy
dt = lambda s: [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]
def hfk(L):
    h = L.knot_floer_homology(); return dict(ranks=sorted((tuple(k), v) for k, v in h['ranks'].items()), tau=h['tau'], epsilon=h['epsilon'],
                                              total_rank=h['total_rank'], seifert_genus=h['seifert_genus'], fibered=h['fibered'])
def test(L):
    a, b = hfk(L), hfk(L.mirror()); eq = a == b
    return dict(tau=a['tau'], epsilon=a['epsilon'], tau_mirror=b['tau'], total_rank=a['total_rank'], genus=a['seifert_genus'], fibered=a['fibered'],
                hfk_equal_mirror=eq, verdict='CHIRAL (HFK differs from mirror)' if not eq else 'HFK blind (equal to mirror)')

if __name__ == '__main__':
    import os as _os
    sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
    from census_codes import main
    from_dt = lambda code: snappy.Link('DT:' + str([tuple(dt(code))]))
    counter = [(nm, (lambda nm=nm: snappy.Link(nm)), (lambda r, e=e: r['hfk_equal_mirror'] is e)) for nm, e in
               (('3_1', False), ('4_1', True), ('8_17', True), ('10_17', True))]
    main(test, counter, from_dt)
