#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Step 3/4: build the SYMMETRIC UNION from a labeled knotoid diagram kappa.

THE PICTURE (Figure 4 of the paper, looked at — not reconstructed from the prose):
the symmetric ribbon disc is a BAND along kappa; its BOUNDARY is the knot K.
  · every aglet becomes a cap  -> there the boundary switches from one side of the band to the other
  · every crossing of kappa becomes a ribbon singularity -> one band goes through a
    slit in the other; in the diagram the TWO boundary lines each cross there ==> 4 crossings

🔴 WHAT THE FIRST DRAFT GOT WRONG (20.09., flipped by the calibration): it treated the
singularity like a CROSSING -- band A completely under band B, all four crossings with A below.
Result: Delta = 1 for ALL conventions and both labels, i.e. the unknot (Footgun 29:
a unanimous result means that a quantity does not arrive at all).
Figure 2 of the paper lists "singularity" and "crossing" as DIFFERENT building blocks, and the
difference is exactly this: at a singularity A PIERCES the surface B -- it enters
at one edge of B from below and exits at the other above. This also matches Lemma 2.3 verbatim:
whoever runs on the SLIT band meets "two over-crossings"; whoever runs on the PUSHED-THROUGH
band meets "an over- and under-crossing". From this it also follows which band goes through:
  "two over-crossings (corresponding to an under-crossing in kappa)"  -> slit = under-passage
  "an over- and under-crossing (corresponding to an over-crossing)"   -> pushed through = over-passage
  · every label a on a strand becomes a twists of the band ==> a crossings
    between the two boundary lines
Crossing number of K: 4*c(kappa) + sum|labels| (not minimal, but that does not matter).

TWO CONVENTIONS THAT THE FIGURES DO NOT GIVE BIT-EXACTLY — therefore CALIBRATED, not guessed
(feedback_calibrate_each_convention_alone_not_the_product):
  flag_sing    : does the band at the OVER position go through the slit, or the one at the under position?
                 (Lemma 2.3 says: the one at the over position -- the calibration must confirm that)
  flag_seite   : which boundary line of the crossed band is hit first?
(flag_eintritt dropped: the entry line is NOT free, but the one that the pushed-through
 band hits FIRST -- there it enters from below. That falls out of `reihenfolge` and needs no flag.)

🔴 SECOND ERROR, likewise flipped by the calibration: an axis crossing (twist) SWAPS the
two band sides. The first draft kept computing on "side 0" afterwards; the sequence became
unrealizable as a result (Gauss parity violated), link_aus_sequenz returned None -- and an `if x`
in the calibration script silently filtered that away, so that both labels seemed to yield the same polynomial.
Footgun 32, built in by ourselves. Since then the side is CARRIED ALONG and turned at every twist.
The referee is the calibration ladder from the source:
  R^s_2 = {1-2t+3t^2-2t^3+t^4, 2-5t+2t^2}                       (Lemma 4.3, 2 polynomials)
  R^s_3 = 6 polynomials                                          (Lemma 4.4)
A convention choice only counts if it produces EXACTLY these sets -- and the three others
choices must fail on it, otherwise the calibration has distinguished nothing (Footgun 16).

STRAND = piece between underpasses (n crossings -> n+1 strands; paper p. 5 + Remark 2.5).
The aglet strand carries 0 (Figure 6); labels suffice mod 2 (Lemma 4.2).
"""
import site_config
import itertools, sys
from collections import defaultdict

HIER = (site_config.ROOT + '/maschinen_ideen/RS4_12N706_20260919')
sys.path.insert(0, HIER)
from pruef_aglet_flaeche import karten
from abschluss import aglet_ueber
from dt_bau import link_aus_sequenz


def strands(w, ueber):
    """strand_von_kante[k] = index of the strand in which edge e_k lies.
    strand boundaries are the UNDERPASSES."""
    m = len(w)
    sv, j = [], 0
    for k in range(m + 1):
        sv.append(j)
        if k < m and ueber[w[k]] != k:      # position k is an underpass
            j += 1
    return sv, j + 1                         # j+1 = number of strands = n+1


def band_rand_sequenz(w, ueber, karte, lab, flag_sing, flag_seite, flag_eintritt=None):
    """Gauss sequence of the boundary K. lab[j] in {0,1} = twist label of strand j (mod 2).
    The boundary line ("side") is CARRIED ALONG: every twist turns it, every cap as well."""
    m = len(w)
    vork = defaultdict(list)
    for p, c in enumerate(w):
        vork[c].append(p)
    rot = karte["rot"]
    sv, n_str = strands(w, ueber)

    twist_auf_kante = {}          # edge -> label a (integer, |a| crossings)
    for j in range(n_str):
        a = lab.get(j, 0)
        if a:
            k = next(k for k in range(m + 1) if sv[k] == j)
            twist_auf_kante[k] = a

    def andere_pos(c, p):
        a, b = vork[c]
        return b if p == a else a

    def reihenfolge(c, p):
        """(first, then): the boundary lines of the CROSSED band in the order they are hit.
        The one hit first is the ENTRY LINE -- there the band dives under."""
        q = andere_pos(c, p)
        liste = rot[('x', c)]
        ix = liste.index((p, 'B'))
        zuerst_plus = (liste[(ix + 1) % 4] == (q, 'B'))
        if not flag_seite:
            zuerst_plus = not zuerst_plus
        return (0, 1) if zuerst_plus else (1, 0)

    def durch(c):
        """Position of the band that is PUSHED THROUGH the slit."""
        a, b = vork[c]
        u = ueber[c]
        return u if flag_sing else (b if u == a else a)

    def kr(c, p, eigene, fremde):
        a, b = vork[c]
        return (c, eigene, fremde) if p == a else (c, fremde, eigene)

    def rolle(c, p, eigene, fremde):
        if durch(c) == p:                                  # pushed-through band
            return 'U' if fremde == reihenfolge(c, p)[0] else 'O'
        return 'O' if eigene == reihenfolge(c, andere_pos(c, p))[0] else 'U'

    seq, gesehen, seite = [], set(), 0

    def twist(k, vorwaerts):
        """A label a is a 2-strand twist with |a| crossings between the two
        boundary lines. The roles ALTERNATE (as in the braid sigma^a), and every crossing
        turns the side. Backwards they are passed in reverse order."""
        nonlocal seite
        a = twist_auf_kante.get(k)
        if not a:
            return
        n = abs(a)
        basis = 'O' if a > 0 else 'U'
        for i in (range(n) if vorwaerts else range(n - 1, -1, -1)):
            nm = (f"tw{k}", i)
            r = basis if i % 2 == 0 else ('U' if basis == 'O' else 'O')
            if nm in gesehen:
                r = 'U' if r == 'O' else 'O'               # second pass: complementary
            gesehen.add(nm)
            seq.append((nm, r))
            seite ^= 1                                     # <- the side swap

    # ---- pass 1: forwards, aglet0 -> aglet1
    for k in range(m + 1):
        twist(k, True)
        if k < m:
            c = w[k]
            for j in reihenfolge(c, k):
                seq.append((kr(c, k, seite, j), rolle(c, k, seite, j)))
    seite ^= 1                                             # cap at aglet1
    # ---- pass 2: backwards, aglet1 -> aglet0
    for k in range(m, -1, -1):
        twist(k, False)
        if k > 0:
            p = k - 1
            c = w[p]
            for j in reversed(reihenfolge(c, p)):
                seq.append((kr(c, p, seite, j), rolle(c, p, seite, j)))
    return seq


def alexander_von(w, ueber, lab, flag_sing, flag_seite, flag_eintritt=None, nur_erste_karte=True):
    """Alexander polynomials over all (or the first) planar map.
    Failures are NOT concealed: None = no DT, LINKn = multi-component."""
    erg = set()
    ks = karten(w)
    for karte in (ks[:1] if nur_erste_karte else ks):
        seq = band_rand_sequenz(w, ueber, karte, lab, flag_sing, flag_seite)
        L, dt = link_aus_sequenz(seq)
        if L is None:
            erg.add(None); continue
        if len(L.link_components) != 1:
            erg.add(f"LINK{len(L.link_components)}"); continue
        erg.add(str(L.alexander_polynomial()))
    return erg
