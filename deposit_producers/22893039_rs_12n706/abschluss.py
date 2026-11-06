#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Step 2b: build the closure cl(kappa) and compute its determinant.

GOAL (Lemma 2.4, Akash et al. arXiv:2606.02390): det(K) = det(cl(kappa))^2.
For 12n706 det(K) = 49, so the sharp filter is  det(cl(kappa)) = 7.

WHY A DUAL PATH IS NEEDED
cl(kappa) connects the aglets "with an arc that contains only over-crossings" (p. 5).
If the aglets lie on different faces -- and that is the case for 8 of the 11 reduced
diagrams with c<=4 (pruef_aglet_flaeche.py, both routes agree) -- this must
arc must cross edges of the diagram. Which ones is told by a path in the DUAL graph from the
aglet-1 face to the aglet-0 face. Every crossed edge gives a crossing in which the
closure arc lies OVER.
  ⚠️ cl(kappa) does NOT depend on the chosen path (an arc lying over everything
  can be moved freely in the half-space above the plane). That is not an assumption here,
  but a built-in COUNTER-CHECK: we compute several paths and both maps through
  and demand the same value.

WHY THE DETERMINANT NEEDS NO SIGNS
The Wirtinger relation at a crossing is x_out = x_over^(+-1) x_in x_over^(-+1),
depending on the sign of the crossing. Fox-differentiated this gives the row
      -t*x_ein + x_aus + (t-1)*x_ueber   or   -t^-1*x_ein + x_aus + (t^-1-1)*x_ueber.
At t = -1 we have t = t^-1, so both cases coincide:
      x_ein + x_aus - 2*x_ueber.
==> for det = |Delta(-1)| NO crossing signs are needed, only the over/under
    structure -- and that is fully contained in the Gauss word. (For the full Alexander polynomial
    later this does NOT hold; there the signs come in.)

CROSSING INFORMATION (p. 16, verbatim)
  "we can assume that all aglet arcs meet an aglet crossing as an under-crossing (or else we
   could perform a simplifying aglet move), and so need only consider crossing information
   for interior crossings."
==> Passage at position 0 and at position 2n-1 is UNDER. Only the inner ones are free.
"""
import site_config
import itertools, sys
from collections import defaultdict, deque
from fractions import Fraction

sys.path.insert(0, (site_config.ROOT + '/maschinen_ideen/RS4_12N706_20260919'))
from pruef_aglet_flaeche import karten


# =================================================== determinant from Wirtinger
def det_aus_wirtinger(tripel, n_boegen):
    """tripel = list of (over, in, out) of arc indices, one per crossing.
    Returns |Delta(-1)| = det(K), exact over Fraction."""
    N = len(tripel)
    if N != n_boegen:
        raise ValueError(f"#crossings {N} != #arcs {n_boegen}")
    if N == 0:
        return 1
    zeilen = []
    for (o, ein, aus) in tripel:
        z = [Fraction(0)] * N
        z[ein] += 1; z[aus] += 1; z[o] -= 2
        zeilen.append(z)

    def minor(zeile_weg, spalte_weg):
        m = [[z[c] for c in range(N) if c != spalte_weg]
             for r, z in enumerate(zeilen) if r != zeile_weg]
        return abs(_det(m))

    d = minor(0, 0)
    # COUNTER-CHECK: all deletion minors are equal up to a unit; at t=-1 thus in absolute value.
    d2 = minor(min(1, N - 1), min(1, N - 1))
    if N > 1 and d != d2:
        raise AssertionError(f"minors disagree: {d} vs {d2} -- presentation broken")
    return int(d)


def _det(m):
    m = [r[:] for r in m]
    n = len(m)
    if n == 0:
        return Fraction(1)
    det = Fraction(1)
    for i in range(n):
        p = next((r for r in range(i, n) if m[r][i] != 0), None)
        if p is None:
            return Fraction(0)
        if p != i:
            m[i], m[p] = m[p], m[i]; det = -det
        det *= m[i][i]
        inv = Fraction(1) / m[i][i]
        for r in range(i + 1, n):
            f = m[r][i] * inv
            if f:
                for c in range(i, n):
                    m[r][c] -= f * m[i][c]
    return det


def wirtinger_aus_sequenz(seq):
    """seq = cyclic list (crossing_id, 'O'|'U') along the closed knot.
    Arcs = maximal segments between consecutive under-passages."""
    u_pos = [i for i, (c, t) in enumerate(seq) if t == 'U']
    N = len(u_pos)
    if N == 0:
        return [], 0
    rang = {p: i for i, p in enumerate(u_pos)}   # the i-th U opens arc i

    def bogen_von(i):
        """Index of the arc in which position i of the sequence lies."""
        vorher = [p for p in u_pos if p < i]
        return rang[max(vorher)] if vorher else N - 1

    o_bogen, u_index = {}, {}
    for i, (c, t) in enumerate(seq):
        if t == 'O':
            o_bogen[c] = bogen_von(i)
        else:
            u_index[c] = rang[i]
    tripel = []
    for c in sorted(o_bogen, key=str):   # IDs are int (from kappa) and str (closure)
        j = u_index[c]
        tripel.append((o_bogen[c], (j - 1) % N, j))
    return tripel, N


# =================================================== dual graph / closure paths
def dual_kanten(karte, m):
    """Edge e -> the (up to two) faces bounding it."""
    fl = karte["flaechen"]
    res = {}
    for e in range(m + 1):
        a = next(i for i, f in enumerate(fl) if (e, 'A') in f)
        b = next(i for i, f in enumerate(fl) if (e, 'B') in f)
        res[e] = (a, b)
    return res


def abschlusswege(karte, m, max_extra=2):
    """All simple face paths from the aglet-1 face to the aglet-0 face,
    up to max_extra edges longer than the shortest. Returns the edge sequence per path."""
    dk = dual_kanten(karte, m)
    nachbarn = defaultdict(list)
    for e, (a, b) in dk.items():
        if a != b:
            nachbarn[a].append((b, e)); nachbarn[b].append((a, e))
    start = next(iter(karte["aglet1"]))
    ziel = next(iter(karte["aglet0"]))
    if start == ziel:
        return [[]]
    # shortest length via BFS
    dist, q = {start: 0}, deque([start])
    while q:
        u = q.popleft()
        for v, e in nachbarn[u]:
            if v not in dist:
                dist[v] = dist[u] + 1; q.append(v)
    if ziel not in dist:
        return []
    grenze = dist[ziel] + max_extra
    wege, stapel = [], [(start, [], {start})]
    while stapel:
        u, kanten, gesehen = stapel.pop()
        if u == ziel:
            wege.append(kanten); continue
        if len(kanten) >= grenze:
            continue
        for v, e in nachbarn[u]:
            if v not in gesehen:
                stapel.append((v, kanten + [e], gesehen | {v}))
    return wege


def cl_sequenz(w, ueber, kantenfolge):
    """Gauss sequence of the closed knot cl(kappa).
    ueber[c] = position at which crossing c is crossed over.
    kantenfolge = edges crossed by the closure arc, in its direction of travel
                  (from aglet 1 back to aglet 0)."""
    m = len(w)
    ab = {e: f"ab{i}" for i, e in enumerate(kantenfolge)}   # one crossing per edge
    seq = []
    for p in range(m):
        if p in ab:                      # edge e_p is crossed by the closure arc
            seq.append((ab[p], 'U'))
        c = w[p]
        seq.append((c, 'O' if ueber[c] == p else 'U'))
    if m in ab:
        seq.append((ab[m], 'U'))
    for e in kantenfolge:                # the closure arc itself: over everywhere
        seq.append((ab[e], 'O'))
    return seq


def aglet_ueber(w):
    """Crossing information of the two aglet crossings: position 0 and 2n-1 are UNDER."""
    m = len(w)
    vork = defaultdict(list)
    for pos, c in enumerate(w):
        vork[c].append(pos)
    fest = {}
    for pos in (0, m - 1):
        c = w[pos]
        p1, p2 = vork[c]
        fest[c] = p2 if pos == p1 else p1
    return fest, vork


def det_cl(w, ueber, nur_erste_karte=False, max_extra=2):
    """det(cl(kappa)) -- over ALL maps and ALL closure paths; demands agreement."""
    m = len(w)
    werte, n_wege = set(), 0
    ks = karten(w)
    for karte in (ks[:1] if nur_erste_karte else ks):
        for kf in abschlusswege(karte, m, max_extra):
            seq = cl_sequenz(w, ueber, kf)
            tripel, N = wirtinger_aus_sequenz(seq)
            werte.add(det_aus_wirtinger(tripel, N)); n_wege += 1
    if not werte:
        return None, 0
    if len(werte) > 1:
        raise AssertionError(f"closure not path-independent: {sorted(werte)}")
    return werte.pop(), n_wege


# =================================================== SELF-TEST
def _torus_sequenz(k):
    """T(2,k) as closure of the braid sigma_1^k: the over/under role alternates."""
    seq = []
    for runde in range(2):
        for i in range(k):
            seq.append((i + 1, 'O' if (i + runde) % 2 == 0 else 'U'))
    return seq


def selbsttest():
    print("=" * 92)
    print("SELF-TEST abschluss.py")
    print("=" * 92)
    ok, n = True, 0

    print("\n-- (1) determinant from Wirtinger, calibrated on T(2,k): det = k --")
    for k in (3, 5, 7, 9, 11):
        tripel, N = wirtinger_aus_sequenz(_torus_sequenz(k))
        d = det_aus_wirtinger(tripel, N)
        gut = (d == k); ok &= gut; n += 1
        print(f"  [{'OK ' if gut else 'RED'}] T(2,{k:<2}) det = {d:<3} (target {k})")

    print("\n-- (2) COUNTER-CASE: one crossing flipped turns T(2,3) into the unknot --")
    seq = _torus_sequenz(3)
    seq = [(c, ('U' if t == 'O' else 'O') if c == 1 else t) for c, t in seq]
    tripel, N = wirtinger_aus_sequenz(seq)
    d = det_aus_wirtinger(tripel, N)
    gut = (d == 1); ok &= gut; n += 1
    print(f"  [{'OK ' if gut else 'RED'}] trefoil with crossing change det = {d} (target 1)")
    print("       (if the counter-case did not bite, the checker would measure nothing)")

    print("\n-- (3) closure calibrated on the paper: c=2 ==> det(cl) MUST be 3 --")
    print("       Corollary 5.1: rs(K)=2 ==> det(K)=9; Lemma 2.4: det(K)=det(cl)^2 ==> 3.")
    w = [1, 2, 1, 2]
    fest, vork = aglet_ueber(w)
    innere = [c for c in vork if c not in fest]
    gut_alle = True
    for wahl in itertools.product(*[[vork[c][0], vork[c][1]] for c in innere]):
        ueber = dict(fest); ueber.update(dict(zip(innere, wahl)))
        d, nw = det_cl(w, ueber)
        gut = (d == 3); gut_alle &= gut; n += 1
        print(f"  [{'OK ' if gut else 'RED'}] '1212' inner={innere or '(none)'}  "
              f"det(cl) = {d}  (target 3)   [{nw} closure paths x maps, all agree]")
    ok &= gut_alle

    print("\n-- (4) closure calibrated on the paper: c=3 ==> det(cl) in {1,5} --")
    print("       Corollary 5.1: rs(K)=3 ==> det(K) = 1 or 25 ==> det(cl) = 1 or 5.")
    for wort in ("123123", "123132"):
        w = [int(c) for c in wort]
        fest, vork = aglet_ueber(w)
        innere = sorted(c for c in vork if c not in fest)
        for wahl in itertools.product(*[[vork[c][0], vork[c][1]] for c in innere]):
            ueber = dict(fest); ueber.update(dict(zip(innere, wahl)))
            d, nw = det_cl(w, ueber)
            gut = (d in (1, 5)); ok &= gut; n += 1
            print(f"  [{'OK ' if gut else 'RED'}] '{wort}' over={dict(zip(innere,wahl))}  "
                  f"det(cl) = {d}  (target 1 or 5)   [{nw} paths agree]")

    print()
    print(f"Base set: {n} probes.")
    print("CALIBRATION:", "ALL GREEN" if ok else "🔴 RED — do not build further")
    return ok


if __name__ == "__main__":
    sys.exit(0 if selbsttest() else 1)
