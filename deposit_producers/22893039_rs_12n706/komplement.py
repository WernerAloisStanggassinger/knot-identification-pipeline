#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Build the UNFILLED complement M: knot + twist circles instead of twist crossings.

WHY: a twist with a crossings is a DEHN FILLING at the twist circle. If one knows M, then for
EVERY filling Thurston holds: vol(M(a)) < vol(M) -- and Futer-Kalfagianni-Purcell gives the
quantitative version   vol(M(a)) >= (1 - (2pi/l)^2)^(3/2) * vol(M),  l = length of the filling slope.
This turns "not found in this search space" into a PROOF: as soon as l is large enough,
vol(M(a)) lies above vol(target knot), and all larger labels are excluded.

THE TWIST CIRCLE: an oval that encloses the two band sides. It crosses each side twice
(4 crossings); when running through the oval it goes over both strands and then under both:
     Oval:    O(A) O(B) U(B) U(A)
     Strand A: U (where the oval is on top) ... O (where it is at the bottom)
Without twist crossings the band side no longer changes -- only the cap turns it.

CALIBRATION (built in, not claimed): M filled with the matching slope must have the SAME
volume as the directly built K_a. If not, the circle is inserted wrongly.
"""
import site_config
import sys
from collections import defaultdict

HIER = (site_config.ROOT + '/maschinen_ideen/RS4_12N706_20260919')
sys.path.insert(0, HIER)
from pruef_aglet_flaeche import karten
from band_rand import strands


# ------------------------------------------------- DT code for LINKS
def sequenzen_zu_dt(komponenten):
    """komponenten = list of sequences [(kreuzungs_id, 'O'|'U'), ...], one per component.
    Returns the DT tuple list, or None if the parity does not work out."""
    pos, lauf = {}, 0
    komp_von = {}
    for ci, seq in enumerate(komponenten):
        for c, t in seq:
            lauf += 1
            pos.setdefault(c, []).append((lauf, t, ci))
    for c, e in pos.items():
        if len(e) != 2:
            return None
        (p1, _, _), (p2, _, _) = e
        if p1 % 2 == p2 % 2:            # one even, one odd -- otherwise no DT
            return None
    tupel = defaultdict(list)
    for c, e in pos.items():
        u = [x for x in e if x[0] % 2 == 1][0]
        g = [x for x in e if x[0] % 2 == 0][0]
        vz = -1 if g[1] == 'O' else 1   # convention as in dt_bau.py
        tupel[u[2]].append((u[0], vz * g[0]))
    aus = []
    for ci in range(len(komponenten)):
        aus.append(tuple(v for _, v in sorted(tupel[ci])))
    return aus


def dreh(seq, k):
    return seq[k:] + seq[:k]


def dt_mit_drehung(komponenten, versuche=400):
    """The parity depends on the starting points. Instead of guessing: rotate until it works out."""
    import itertools
    laengen = [len(s) for s in komponenten]
    for offsets in itertools.islice(itertools.product(*[range(l) for l in laengen]), versuche):
        dt = sequenzen_zu_dt([dreh(s, o) for s, o in zip(komponenten, offsets)])
        if dt is not None:
            return dt, offsets
    return None, None


# ------------------------------------------------- band boundary WITH twist circles
def band_mit_kreisen(w, ueber, karte, kreis_strands, flag_sing, flag_seite):
    """Like band_rand_sequenz, but instead of twist crossings at every strand from
    kreis_strands there is a twist CIRCLE (own component). All labels are 0."""
    m = len(w)
    vork = defaultdict(list)
    for p, c in enumerate(w):
        vork[c].append(p)
    rot = karte["rot"]
    sv, n_str = strands(w, ueber)

    kreis_auf_kante = {}
    for j in kreis_strands:
        k = next(k for k in range(m + 1) if sv[k] == j)
        kreis_auf_kante[k] = j

    def andere_pos(c, p):
        a, b = vork[c]
        return b if p == a else a

    def reihenfolge(c, p):
        q = andere_pos(c, p)
        liste = rot[('x', c)]
        ix = liste.index((p, 'B'))
        zuerst_plus = (liste[(ix + 1) % 4] == (q, 'B'))
        if not flag_seite:
            zuerst_plus = not zuerst_plus
        return (0, 1) if zuerst_plus else (1, 0)

    def durch(c):
        a, b = vork[c]
        u = ueber[c]
        return u if flag_sing else (b if u == a else a)

    def kr(c, p, eigene, fremde):
        a, b = vork[c]
        return (c, eigene, fremde) if p == a else (c, fremde, eigene)

    def rolle(c, p, eigene, fremde):
        if durch(c) == p:
            return 'U' if fremde == reihenfolge(c, p)[0] else 'O'
        return 'O' if eigene == reihenfolge(c, andere_pos(c, p))[0] else 'U'

    haupt, seite = [], 0
    # pass 1 (band side A): at the circle first U (oval on top), then O (oval at the bottom)
    for k in range(m + 1):
        if k in kreis_auf_kante:
            j = kreis_auf_kante[k]
            haupt += [((f"ov{j}", 1), 'U'), ((f"ov{j}", 4), 'O')]
        if k < m:
            c = w[k]
            for jj in reihenfolge(c, k):
                haupt.append((kr(c, k, seite, jj), rolle(c, k, seite, jj)))
    seite ^= 1                                   # cap
    # pass 2 (band side B), backwards
    for k in range(m, -1, -1):
        if k in kreis_auf_kante:
            j = kreis_auf_kante[k]
            haupt += [((f"ov{j}", 3), 'O'), ((f"ov{j}", 2), 'U')]
        if k > 0:
            p = k - 1
            c = w[p]
            for jj in reversed(reihenfolge(c, p)):
                haupt.append((kr(c, p, seite, jj), rolle(c, p, seite, jj)))
    # per circle one own component: over both strands, then under both
    komponenten = [haupt]
    for j in kreis_strands:
        komponenten.append([((f"ov{j}", 1), 'O'), ((f"ov{j}", 2), 'O'),
                            ((f"ov{j}", 3), 'U'), ((f"ov{j}", 4), 'U')])
    return komponenten
