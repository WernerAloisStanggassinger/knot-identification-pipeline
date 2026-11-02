# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# HOMFLY chirality witness. Sage 10.9 Link.homfly_polynomial
# (libhomfly). Amphichiral ⟹ H(K) == H(mK); so H(K) != H(mK) ⟹ chiral, independent of SnapPy's isometry signature.
# H(K) == H(mK) proves nothing (HOMFLY is chirality-blind for e.g. 9_42). Counter-cases: 3_1 differ, 5_2 differ, 4_1 equal, 9_42 equal (chiral!).
import sys, json, bz2, csv, signal
from multiprocessing import Pool
dt = lambda s: [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]
def init():
    global snappy; import sage.all, snappy
def work(item):
    name, code = item
    signal.signal(signal.SIGALRM, lambda *a: (_ for _ in ()).throw(TimeoutError())); signal.alarm(900)
    try:
        L = snappy.Link('DT:' + str([tuple(dt(code))])) if code else snappy.Link(name)
        S = L.sage_link(); H = S.homfly_polynomial(); Hm = S.mirror_image().homfly_polynomial()
        return dict(name=name, homfly=str(H), homfly_mirror=str(Hm), mirror_equal=(H == Hm),
                    verdict='CHIRAL (HOMFLY != mirror)' if H != Hm else 'HOMFLY blind (equal to mirror)')
    except Exception as e:
        return dict(name=name, verdict='ERROR', error=type(e).__name__ + ' ' + str(e)[:100])
    finally: signal.alarm(0)

if __name__ == '__main__':
    import os as _os
    sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
    from census_codes import main
    init()
    def test(L):
        S = L.sage_link(); H = S.homfly_polynomial(); Hm = S.mirror_image().homfly_polynomial()
        return dict(homfly=str(H), homfly_mirror=str(Hm), mirror_equal=(H == Hm),
                    verdict='CHIRAL (HOMFLY != mirror)' if H != Hm else 'HOMFLY blind (equal to mirror)')
    from_dt = lambda code: snappy.Link('DT:' + str([tuple(dt(code))]))
    counter = [(nm, (lambda nm=nm: snappy.Link(nm)), (lambda r, e=e: r['mirror_equal'] is e)) for nm, e in
               (('3_1', False), ('5_2', False), ('4_1', True), ('9_42', True))]
    main(test, counter, from_dt)
