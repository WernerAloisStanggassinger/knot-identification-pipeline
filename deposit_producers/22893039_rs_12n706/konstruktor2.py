#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""SECOND, INDEPENDENT implementation of the union constructor (labeled knotoid -> Delta_K).

Assignment from Werner, 22.09.2026: "also second-implement the constructor" (orig. „konstruktor auch zweit-implementieren").

FIRST ATTEMPT, FAILED AND RECORDED: the construction from Lamm's side (knotoid on the left,
mirror image on the right, twists on the axis) is WRONG for the general case. It presupposes
that both aglets lie in the same face. Already the only diagram with two crossings,
`1212`, violates that (both crossings are at positions of equal parity), and 8 of the 11
diagrams up to c=4 likewise. The calibration against R^s_2 immediately turned red. Exactly for that reason the first
implementation goes via the BOUNDARY OF THE RIBBON SURFACE — that is necessary, not roundabout.

THIS VERSION goes the same way, but is built differently in four points:
  (1) The planar embedding comes from `networkx` (gadget graph from `zweite_implementierung.py`),
      not from our own map code `pruef_aglet_flaeche.karten`.
  (2) The boundary traversal is newly derived (two rails, caps, slit passage),
      without a line from `band_rand.py`.
  (3) Delta is computed by **Sage**, not spherogram.
  (4) The conventions are fixed ONLY on R^s_2 (2 polynomials) and then VALIDATED on R^s_3 (6) and
      R^s_4 (27) — 33 withheld foreign data points. Whoever uses all 35 for tuning
      has checked nothing, but fitted.
REMAINS SHARED exactly one step: spherogram builds the planar embedding of the
FINISHED knot. An error there would hit both versions.

Model (paper arXiv:2606.02390v1, section 2, Figure 4, read verbatim 22.09.2026):
aglet -> cap, arc -> strip, crossing -> ribbon singularity (one band through a slit
in the other). The boundary of the strip consists of two rails; at every cap one changes the
rail. At the slit every rail of the pushed-through band crosses both rails of the
slit band — four crossings per singularity, each traversed exactly twice.

Usage (sage_stable env):
  PATH=/home/squeezy/miniforge3/envs/sage_stable/bin:$PATH \
  /home/squeezy/miniforge3/envs/sage_stable/bin/python konstruktor2.py eichen
  ... konstruktor2.py 5 --json
"""
import itertools
import json
import os
import sys

ALLE_KARTEN = bool(os.environ.get('ALLE_KARTEN'))
HIER = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HIER)


# ----------------------------------------------------------- map from networkx
def karte_nx(w):
    """The planar map of the arc diagram, read from the embedding of `networkx`.

    Every crossing becomes a quadrilateral with a centre point (wheel W4). The wheel is 3-fold
    connected, so its embedding is unique up to mirroring ⟹ the four
    edge ends necessarily lie alternately, that is transversality. Which of the
    two mirror images is present at a crossing is told by the found embedding itself;
    it is globally compatible, and exactly that is the map.

    ⚠️ Measured and recorded: the PATTERN CHOICE when building the gadget decides nothing —
    every quadrilateral can be mirrored individually, therefore all 2^n pattern choices are planar.
    The map stands only in the EMBEDDING, not in the construction.

    Returns: {crossing: True | False}, True = cyclic order (in@p1, in@p2, out@p1, out@p2).
    """
    import networkx as nx

    vork = {}
    for pos, c in enumerate(w):
        vork.setdefault(c, []).append(pos)
    G = nx.Graph()
    ecke = {}
    for c, (p1, p2) in vork.items():
        for i in range(4):
            G.add_edge((c, i), (c, (i + 1) % 4))
            G.add_edge((c, 'm'), (c, i))
        ecke[(p1, 'rein')], ecke[(p2, 'rein')] = (c, 0), (c, 1)
        ecke[(p1, 'raus')], ecke[(p2, 'raus')] = (c, 2), (c, 3)
    m = len(w)
    G.add_edge('A0', ecke[(0, 'rein')])
    for k in range(1, m):
        G.add_edge(ecke[(k - 1, 'raus')], ecke[(k, 'rein')])
    G.add_edge(ecke[(m - 1, 'raus')], 'A1')

    ok, emb = nx.check_planarity(G)
    if not ok:
        return None
    aus = {}
    for c in vork:
        L = list(emb.neighbors_cw_order((c, 0)))
        i1 = L.index((c, 1))
        aus[c] = (L[(i1 + 1) % len(L)] == (c, 'm'))    # clockwise: 1, then centre
    return aus


# ------------------------------------------------------------------ strands (section 2)
def straenge(w, ueber):
    """sv[k] = strand number of edge e_k; boundaries are the UNDERPASSES (n+1 strands)."""
    m, sv, j = len(w), [], 0
    for k in range(m + 1):
        sv.append(j)
        if k < m and ueber[w[k]] != k:
            j += 1
    return sv, j + 1


# --------------------------------------------------------- boundary traversal of the ribbon surface
def rand_sequenz(w, ueber, lab, rot, konv):
    """Gauss sequence of the boundary K = (name, 'O'|'U') in traversal order.

    konv = (durch_ist_ueber, durchgestecktes_band_unten, ordnung_gedreht) — three switches,
    which are fixed on R^s_2 alone.
    """
    durch_ist_ueber, durch_unten, ordnung_gedreht = konv
    m = len(w)
    vork = {}
    for p, c in enumerate(w):
        vork.setdefault(c, []).append(p)
    sv, n_str = straenge(w, ueber)

    twist_kante = {}
    for j in range(n_str):
        if lab.get(j, 0) % 2:
            twist_kante[next(k for k in range(m + 1) if sv[k] == j)] = j

    def durchgesteckt(c, p):
        """Is the band at position p the one pushed through the slit?"""
        return (ueber[c] == p) if durch_ist_ueber else (ueber[c] != p)

    def reihenfolge(c, p):
        """In which order does our rail hit the two rails of the other
        band? From the embedding: if we run 'with' the rotation, we hit 0 first."""
        p1, p2 = vork[c]
        mit = rot[c] if p == p1 else not rot[c]
        if ordnung_gedreht:
            mit = not mit
        return (0, 1) if mit else (1, 0)

    def rolle(c, p, eigene, fremde):
        """Role of our rail at the crossing with the foreign rail `fremde`.

        The pushed-through band dives under at the FIRST-hit edge of the slit band
        and comes out above it at the second — it goes through the surface, after all.
        The slit band sees the opposite role at the same crossing.
        """
        p1, p2 = vork[c]
        t = p1 if durchgesteckt(c, p1) else p2          # position of the pushed-through band
        zuerst = reihenfolge(c, t)[0]                   # rail of the other band hit first
        if p == t:                                      # we are the pushed-through band
            unten = (fremde == zuerst)
        else:                                           # we are the slit band
            unten = not (eigene == zuerst)
        if not durch_unten:
            unten = not unten
        return 'U' if unten else 'O'

    def name(c, p, eigene, fremde):
        """Canonical name: (crossing, rail of the band at p1, rail of the band at p2)."""
        p1, _ = vork[c]
        return (c, eigene, fremde) if p == p1 else (c, fremde, eigene)

    seq, schiene = [], 0

    def twists(k, vorwaerts):
        nonlocal schiene
        j = twist_kante.get(k)
        if j is None:
            return
        seq.append((('tw', j), 'O' if (schiene == 0) == vorwaerts else 'U'))
        schiene ^= 1

    def passage(c, p):
        """Two crossings: our rail against both rails of the other band."""
        for fremde in reihenfolge(c, p):
            seq.append((name(c, p, schiene, fremde), rolle(c, p, schiene, fremde)))

    for k in range(m + 1):                       # pass 1: forwards on rail 0
        twists(k, True)
        if k < m:
            passage(w[k], k)
    schiene ^= 1                                 # cap at aglet 1
    for k in range(m, -1, -1):                   # pass 2: backwards on the other
        twists(k, False)
        if k > 0:
            c, p = w[k - 1], k - 1
            for fremde in reversed(reihenfolge(c, p)):
                seq.append((name(c, p, schiene, fremde), rolle(c, p, schiene, fremde)))
    return seq


# --------------------------------------------------------------- Gauss sequence -> DT code
def dt_code(seq):
    pos = {}
    for i, (nm, rolle) in enumerate(seq, start=1):
        pos.setdefault(nm, []).append((i, rolle))
    if any(len(v) != 2 for v in pos.values()):
        return None
    paare = []
    for nm, ((i1, r1), (i2, r2)) in pos.items():
        if i1 % 2 == i2 % 2:
            return None
        (u, _), (g, rg) = ((i1, r1), (i2, r2)) if i1 % 2 else ((i2, r2), (i1, r1))
        paare.append((u, -g if rg == 'O' else g))
    paare.sort()
    return tuple(g for _, g in paare)


# ------------------------------------------------------------------------ DT -> Delta (Sage)
def alexander(dt):
    from spherogram import Link as SLink
    from sage.all import Link as SageLink
    try:
        L = SLink(f'DT: [{dt}]')
    except Exception:
        return 'KEIN_LINK'
    if len(L.link_components) != 1:
        return f'LINK{len(L.link_components)}'
    try:
        pd = [[int(e) + 1 for e in c] for c in L.PD_code(KnotTheory=False)]
        p = SageLink(pd).alexander_polynomial()
    except Exception as e:
        return f'FEHLER:{type(e).__name__}'
    return normiere(p)


def normiere(p):
    c, e = list(p.coefficients()), list(p.exponents())
    if not c:
        return (0,)
    lo, hi = min(e), max(e)
    v = [0] * (hi - lo + 1)
    for ci, ei in zip(c, e):
        v[ei - lo] = int(ci)
    return kanon(v)


def kanon(v):
    while v and v[0] == 0:
        v = v[1:]
    while v and v[-1] == 0:
        v = v[:-1]
    if not v:
        return (0,)
    return min(tuple(v), tuple(v[::-1]), tuple(-x for x in v), tuple(-x for x in v[::-1]))


# ------------------------------------------------------------------------- task building
def aufgaben(n):
    from zweite_implementierung import enumerieren
    aus = []
    for t in sorted(enumerieren(n, 'eng_weit')):
        w = list(t)
        vork = {}
        for p, c in enumerate(w):
            vork.setdefault(c, []).append(p)
        # aglet arcs dive under (paper section 4) => the over-passage lies elsewhere
        fest = {}
        p1, p2 = vork[w[0]]
        fest[w[0]] = p2 if p1 == 0 else p1
        p1, p2 = vork[w[-1]]
        fest[w[-1]] = p1 if p2 == len(w) - 1 else p2
        innere = sorted(c for c in vork if c not in (w[0], w[-1]))
        rot = karte_nx(w)
        for wahl in itertools.product(*[vork[c] for c in innere]):
            ueber = dict(fest)
            ueber.update(dict(zip(innere, wahl)))
            _, n_str = straenge(w, ueber)
            frei = [j for j in range(n_str) if j not in (0, n_str - 1)]
            for bits in itertools.product((0, 1), repeat=len(frei)):
                aus.append((w, ueber, dict(zip(frei, bits)), rot))
    return aus


def menge_fuer(n, konv, laut=False):
    menge, kaputt = set(), {}
    for w, ueber, lab, rot in aufgaben(n):
        seq = rand_sequenz(w, ueber, lab, rot, konv)
        dt = dt_code(seq)
        if dt is None:
            kaputt['kein DT'] = kaputt.get('kein DT', 0) + 1
            continue
        r = alexander(dt)
        if isinstance(r, tuple):
            menge.add(r)
        else:
            kaputt[r] = kaputt.get(r, 0) + 1
    if laut:
        print(f'   c={n}: {len(menge)} polynomials, not usable {kaputt}', flush=True)
    return menge, kaputt


# ------------------------------------------------------------- calibration sets from the paper
def paper_mengen():
    """R^s_2 (Lemma 4.3) and R^s_3 (Lemma 4.4), verbatim from the arXiv full text."""
    rs2 = [[1, -2, 3, -2, 1], [2, -5, 2]]
    rs3 = [[1],
           [1, -1, -1, 3, -1, -1, 1],
           [1, 0, -3, 0, 1],
           [1, -3, 5, -7, 5, -3, 1],
           [1, -6, 11, -6, 1],
           [2, -6, 9, -6, 2]]
    return {2: {kanon(v) for v in rs2}, 3: {kanon(v) for v in rs3}}


def tabelle1():
    """The 27 polynomials from Table 1, parsed from the downloaded full text."""
    import re
    pfad = os.path.join(HIER, 'paper_tabelle1.txt')
    if not os.path.exists(pfad):
        return None
    aus = set()
    for zeile in open(pfad):
        zeile = zeile.strip()
        if not zeile:
            continue
        v = {}
        for vz, koef, exp in re.findall(r'([+-]?)\s*(\d*)\s*t?\^?\{?(\d*)\}?', zeile):
            pass
        aus.add(zeile)
    return aus


def main():
    arg = sys.argv[1] if len(sys.argv) > 1 else 'eichen'
    soll = paper_mengen()

    if arg == 'eichen':
        print('FIXING CONVENTIONS — only on R^s_2 (2 polynomials), then validate on R^s_3')
        treffer = []
        for konv in itertools.product((True, False), repeat=3):
            m2, k2 = menge_fuer(2, konv)
            gut = (m2 == soll[2])
            print(f'  {konv}: {len(m2)} polynomials  {"HITS R^s_2" if gut else ""} {k2 if k2 else ""}')
            if gut:
                treffer.append(konv)
        print(f'\nConventions that hit R^s_2: {treffer}')
        if not treffer:
            print('🔴 NONE — the model is wrong.')
            return
        if len(treffer) == len(list(itertools.product((True, False), repeat=3))):
            print('⚠️ ALL hit — the calibration distinguishes nothing (Footgun 16).')
        print('\nVALIDATION on R^s_3 (6 polynomials, withheld):')
        for konv in treffer:
            m3, k3 = menge_fuer(3, konv)
            fehlt, zuviel = soll[3] - m3, m3 - soll[3]
            print(f'  {konv}: {len(m3)} polynomials · '
                  f'{"GREEN" if not fehlt and not zuviel else "RED"}'
                  f'{"" if not fehlt else f" missing {sorted(fehlt)}"}'
                  f'{"" if not zuviel else f" too many {sorted(zuviel)}"} {k3 if k3 else ""}')
        return

    n = int(arg)
    konv = eval(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].startswith('(') else (True, True, False)
    m, k = menge_fuer(n, konv, laut=True)
    print(f'c={n} with convention {konv}: {len(m)} polynomials, not usable {k}')
    if '--json' in sys.argv:
        ziel = os.path.join(HIER, f'konstruktor2_rs{n}.json')
        json.dump({'c': n, 'konvention': list(konv), 'polynome': sorted(map(list, m)), 'kaputt': k},
                  open(ziel, 'w'), indent=1)
        print('written:', ziel)


if __name__ == '__main__':
    main()
