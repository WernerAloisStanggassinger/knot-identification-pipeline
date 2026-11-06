#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""skein_pass.py — HOMFLY/Jones by TRANSFER MATRIX over a cut,
where BOTH existing tools hit the wall (Disc #714 (d)).

── WHY THIS MUST EXIST ──────────────────────────────────────────────────────
Two tools, two different walls, measured:
  · regina dies of MEMORY — on the flip_ family already at c=70 with 9 GB.
    The code difference raw/`Algorithm.Treewidth` is measured DEAD as a cause
    (26.08.: c=70 falls at 6 GB in ALL three variants).
  · our Hecke trace (`homfly_interpolation.py`) dies of the STRAND COUNT,
    effort ~ n!. The open records have n = 23..58 strands.
The number neither of them sees is the CUT WIDTH: 10..14 arcs in the field
(`bandwidth.py`, 26.08.) — while c goes up to 133 and n up to 58.

── THE IDEA, IN ONE SENTENCE ────────────────────────────────────────────────
Work through the diagram crossing by crossing and carry along between two steps only the
state of the CUT: which open arcs in the already processed
part are connected to each other. The price is the number of states per
cut — not c, not n.

⚠️ **HOW MANY STATES THAT IS, IS MEASURED HERE, NOT CLAIMED.**
`bandwidth.py` prints k! as an UPPER BOUND. This tool counts the
states actually reached (`max_zustaende` in the result) — that is the
measurement that was explicitly left open there.

── STAGE 1 (this file): KAUFFMAN BRACKET ⟹ JONES ────────────────────────────
Deliberately first, although HOMFLY is the goal. The bracket needs NO
orientation and NO Hecke algebra — it checks the whole scaffold
(ordering, cut management, merging, loop factor) against something that
regina computes independently: `Link.bracket()`. A HOMFLY run on an
uncalibrated scaffold would be two unchecked claims on top of each other.

── WHAT IS NOT GUESSED HERE ─────────────────────────────────────────────────
The assignment "which of the two planar smoothings is the A-smoothing"
depends on regina's PD convention AND on the sign of the crossing. It is
NOT derived but CALIBRATED: all four combinations are run against
`Link.bracket()` on the right- and left-handed trefoil, exactly one
may hit both (`--eiche-konvention`). Exactly the error that cost the winding
direction in the Seifert builder (Hard Rule 9).

── COUNTER-CASES (--selbsttest; each one must BITE) ─────────────────────────
  A trefoil right  : bracket == regina.bracket()
  B trefoil left   : bracket == regina.bracket()   (separates the conventions)
  C figure-eight   : bracket == regina.bracket()
  D unknot         : bracket == 1
  E WRONG convention must FAIL  (otherwise A/B/C check nothing)
  F WRONG delta must FAIL       (otherwise the loop counting is blind)
  G ORDER INVARIANCE: 12 random crossing orders, one result
  H Jones == regina.jones() on 4 knots
  I 3_1#3_1#3_1: max_zustaende stays small although c=9 (does it secretly measure c?)
"""
import site_config
import sys, os, json, time, random, argparse, math
from collections import defaultdict

sys.path.insert(0, (site_config.FREEZE_DIR))
from pipeline_seal import seal_all
seal_all(verbose=False)
# 🔴 REGINA STAYS A DIRECT IMPORT HERE — with reason (Disc #750, 2026-08-29).
# Regina is GPL v2+, and the pipeline is switching calls over to subprocesses. THIS file
# is not: it passes a DIAGRAM on after `simplify()` (signedGauss/PD), and Regina's
# simplify() on larger diagrams picks, depending on the process, a different one of the many
# equivalent results. Measured on polygon_resolve.py: same invariants, but a
# different Gauss code and in one resolution even a different crossing number.
# Whoever wants to switch it anyway must put ALL calls into ONE subprocess — not one
# process per call. Bridge: pipeline_scripts/regina_bridge.py
import regina

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bandwidth import bogen_inzidenz, schnittbreite, breite_fuer

ROOT = site_config.ROOT

# ---------------------------------------------------------------------------
# Laurent polynomials in A as dict {exponent: coefficient}. Exact, integer.
# ---------------------------------------------------------------------------
def lp_mono(e, k=1):
    return {e: k} if k else {}

def lp_add(p, q):
    r = dict(p)
    for e, k in q.items():
        n = r.get(e, 0) + k
        if n:
            r[e] = n
        elif e in r:
            del r[e]
    return r

def lp_mul(p, q):
    r = {}
    for e1, k1 in p.items():
        for e2, k2 in q.items():
            e = e1 + e2
            n = r.get(e, 0) + k1 * k2
            if n:
                r[e] = n
            elif e in r:
                del r[e]
    return r

def lp_str(p):
    if not p:
        return '0'
    return ' + '.join(f'{k}*A^{e}' for e, k in sorted(p.items()))

DELTA = {2: -1, -2: -1}          # delta = -A^2 - A^-2
DELTA_INV_CACHE = {}

def lp_teile_durch_delta(p):
    """Polynomial division by delta = -(A^2 + A^-2). Raises if remainder != 0.

    ⚠️ The LATCH here is not cosmetics: without it the loop runs forever on
    a NON-divisible polynomial — exactly so counter-case F hung on the
    first run (wrong delta) instead of failing. A test that hangs is
    no test."""
    if not p:
        return {}
    rest = dict(p)
    erg = {}
    grenze = (max(p) - min(p)) + 8
    for _ in range(grenze):
        if not rest:
            return {e: k for e, k in erg.items() if k}
        e = max(rest)
        k = rest[e]
        q_e, q_k = e - 2, -k            # leading term of delta is -A^2
        erg[q_e] = erg.get(q_e, 0) + q_k
        rest = lp_add(rest, lp_mul({q_e: q_k}, {2: 1, -2: 1}))
    raise ValueError('not divisible by delta — loop bookkeeping broken')


# ---------------------------------------------------------------------------
# The run-through
# ---------------------------------------------------------------------------
def _komponenten(kanten, knoten):
    """Graph with degree<=2 -> list of (type, endpoints). type: 'zyklus'|'pfad'."""
    nach = defaultdict(list)
    for u, v in kanten:
        nach[u].append(v)
        nach[v].append(u)
    gesehen = set()
    out = []
    # paths first (starting from nodes of degree 1)
    for s in knoten:
        if s in gesehen or len(nach[s]) != 1:
            continue
        pfad = [s]
        gesehen.add(s)
        vor, akt = None, s
        while True:
            weiter = [w for w in nach[akt] if w != vor]
            if not weiter:
                break
            vor, akt = akt, weiter[0]
            if akt in gesehen:
                break
            gesehen.add(akt)
            pfad.append(akt)
        out.append(('pfad', (pfad[0], pfad[-1])))
    # the rest are cycles
    for s in knoten:
        if s in gesehen:
            continue
        gesehen.add(s)
        vor, akt = None, s
        while True:
            weiter = [w for w in nach[akt] if w != vor]
            if not weiter:
                break
            vor, akt = akt, weiter[0]
            if akt in gesehen:
                break
            gesehen.add(akt)
        out.append(('zyklus', None))
    return out


def _repr_ende(x, boegen4):
    """End node of a path -> the arc number it represents after the step."""
    typ, w = x
    return w if typ == 'b' else boegen4[w]


def klammer(pd, ordnung, glaettung, delta=DELTA, laut=False):
    """Kauffman bracket by transfer matrix.
    pd: list of 4-tuples (arc numbers, cyclic ccw).
    glaettung(sign_unknown, t) -> ((pairingA, weightA), (pairingB, weightB)),
      pairing = ((i,j),(k,l)) with port indices 0..3.
    Returns: (bracket_polynomial, max_states, sum_states)
    """
    m = len(pd)
    if m == 0:
        return {0: 1}, 1, 1
    zaehler = defaultdict(int)
    for t in pd:
        for b in t:
            zaehler[b] += 1
    # State: frozenset of frozenset pairs of open arcs -> coefficient
    zustaende = {frozenset(): {0: 1}}
    offen_cnt = defaultdict(int)
    max_z, summe_z = 1, 0
    for schritt_i, ci in enumerate(ordnung):
        boegen4 = list(pd[ci])
        vorher_offen = {b for b, v in offen_cnt.items() if 0 < v < zaehler[b]}
        for b in boegen4:
            offen_cnt[b] += 1
        nachher_offen = {b for b, v in offen_cnt.items() if 0 < v < zaehler[b]}
        neu = {}
        for zus, koeff in zustaende.items():
            paare = {}
            for pr in zus:
                u, v = tuple(pr)
                paare[u] = v
                paare[v] = u
            for paarung, gewicht in glaettung(ci):
                kanten, knoten = [], set()
                # 1) existing pairings
                for pr in zus:
                    u, v = tuple(pr)
                    kanten.append((('b', u), ('b', v)))
                    knoten.add(('b', u)); knoten.add(('b', v))
                # 2)/3) attach ports to the arcs
                erst = {}
                for i, b in enumerate(boegen4):
                    knoten.add(('p', i))
                    if b in vorher_offen:
                        kanten.append((('p', i), ('b', b)))
                        knoten.add(('b', b))
                    else:
                        if b in erst:
                            kanten.append((('p', erst[b]), ('p', i)))
                        else:
                            erst[b] = i
                # 4) smoothing
                for i, j in paarung:
                    kanten.append((('p', i), ('p', j)))
                komp = _komponenten(kanten, knoten)
                kreise = sum(1 for t, _ in komp if t == 'zyklus')
                neu_paare = []
                ok = True
                for t, ends in komp:
                    if t != 'pfad':
                        continue
                    a, b_ = ends
                    ra, rb = _repr_ende(a, boegen4), _repr_ende(b_, boegen4)
                    if ra not in nachher_offen or rb not in nachher_offen:
                        ok = False
                        break
                    neu_paare.append(frozenset((ra, rb)))
                if not ok:
                    raise RuntimeError('Path end is not an open arc — '
                                       'bookkeeping broken')
                nz = frozenset(neu_paare)
                f = dict(gewicht)
                for _ in range(kreise):
                    f = lp_mul(f, delta)
                bei = lp_mul(koeff, f)
                neu[nz] = lp_add(neu.get(nz, {}), bei)
        zustaende = {k: v for k, v in neu.items() if v}
        max_z = max(max_z, len(zustaende))
        summe_z += len(zustaende)
        if laut and (schritt_i % 10 == 0 or schritt_i == m - 1):
            groesse = sum(len(v) for v in zustaende.values())
            print(f'    [{schritt_i+1}/{m}] states {len(zustaende):>7} '
                  f'terms {groesse:>9}', flush=True)
    if list(zustaende) != [frozenset()]:
        raise RuntimeError(f'open arcs at the end: {list(zustaende)[:2]}')
    roh = zustaende[frozenset()]
    # every closed circle has received a delta; convention <unknot> = 1
    return lp_teile_durch_delta(roh), max_z, summe_z


def zaehle_zustaende(pd, ordnung, riegel=None):
    """Propagate only the STATE SET, without coefficients.

    The states do not depend on the coefficients — so this is the
    honest measurement of the state count, and it costs nothing. Exactly this number
    was left open by `bandwidth.py` (there k! stands as UPPER BOUND, with the
    explicit remark that the real number is only delivered by the built run-through).
    riegel: max allowed state count, otherwise abort with value.
    """
    m = len(pd)
    if m == 0:
        return {'max_zustaende': 1, 'abgebrochen': False, 'bei': None}
    zaehler = defaultdict(int)
    for t in pd:
        for b in t:
            zaehler[b] += 1
    zustaende = {frozenset()}
    offen_cnt = defaultdict(int)
    max_z = 1
    for schritt_i, ci in enumerate(ordnung):
        boegen4 = list(pd[ci])
        vorher_offen = {b for b, v in offen_cnt.items() if 0 < v < zaehler[b]}
        for b in boegen4:
            offen_cnt[b] += 1
        nachher_offen = {b for b, v in offen_cnt.items() if 0 < v < zaehler[b]}
        neu = set()
        for zus in zustaende:
            for paarung in (PAARUNG_1, PAARUNG_2):
                kanten, knoten = [], set()
                for pr in zus:
                    u, v = tuple(pr)
                    kanten.append((('b', u), ('b', v)))
                    knoten.add(('b', u)); knoten.add(('b', v))
                erst = {}
                for i, b in enumerate(boegen4):
                    knoten.add(('p', i))
                    if b in vorher_offen:
                        kanten.append((('p', i), ('b', b)))
                        knoten.add(('b', b))
                    elif b in erst:
                        kanten.append((('p', erst[b]), ('p', i)))
                    else:
                        erst[b] = i
                for i, j in paarung:
                    kanten.append((('p', i), ('p', j)))
                paare = []
                for t, ends in _komponenten(kanten, knoten):
                    if t != 'pfad':
                        continue
                    ra = _repr_ende(ends[0], boegen4)
                    rb = _repr_ende(ends[1], boegen4)
                    paare.append(frozenset((ra, rb)))
                neu.add(frozenset(paare))
        zustaende = neu
        max_z = max(max_z, len(zustaende))
        if riegel and len(zustaende) > riegel:
            return {'max_zustaende': max_z, 'abgebrochen': True,
                    'bei': schritt_i + 1}
    return {'max_zustaende': max_z, 'abgebrochen': False, 'bei': None}


# ---------------------------------------------------------------------------
# Smoothing conventions: NOT derived, but calibrated.
# ---------------------------------------------------------------------------
PAARUNG_1 = ((0, 1), (2, 3))
PAARUNG_2 = ((1, 2), (3, 0))

def mach_glaettung(vorzeichen, wahl):
    """wahl: (for_positive, for_negative) each 1 or 2 = which pairing is A."""
    def g(ci):
        s = vorzeichen[ci]
        w = wahl[0] if s > 0 else wahl[1]
        pa = PAARUNG_1 if w == 1 else PAARUNG_2
        pb = PAARUNG_2 if w == 1 else PAARUNG_1
        return ((pa, {1: 1}), (pb, {-1: 1}))
    return g


def vorzeichen_liste(L):
    return [int(L.crossing(i).sign()) for i in range(L.size())]


def rechne(L, wahl=(1, 1), neustarts=40, ordnung=None, laut=False, delta=DELTA):
    pd = L.pdData()
    if ordnung is None:
        k, ordnung = schnittbreite(pd, neustarts)
    else:
        m, inz = bogen_inzidenz(pd)
        k = breite_fuer(inz, ordnung)
    g = mach_glaettung(vorzeichen_liste(L), wahl)
    t0 = time.time()
    br, mz, sz = klammer(pd, ordnung, g, delta=delta, laut=laut)
    return {'bracket': br, 'schnittbreite': k, 'max_zustaende': mz,
            'summe_zustaende': sz, 'sek': round(time.time() - t0, 3),
            'c': len(pd)}


def jones_aus_klammer(br, writhe):
    """V(A) = (-A^3)^-w * <D>, then t = A^-4.
    Returns: dict {exponent_in_A: coeff} — regina's jones() is in x = t^1/2,
    so x^k <-> A^-2k. We return both, comparison is in A."""
    f = {-3 * writhe: (-1) ** (writhe % 2)}
    return lp_mul(br, f)


def regina_bracket_als_lp(L):
    """regina.bracket() -> dict. Laurent in A."""
    p = L.bracket()
    d = {}
    for e in range(p.minExp(), p.maxExp() + 1):
        k = int(str(p[e]))
        if k:
            d[e] = k
    return d


def regina_jones_als_lp_in_A(L):
    """regina.jones() is in x with t = x^2; our Jones is in A with
    t = A^-4  ⟹  x^k corresponds to A^-2k."""
    p = L.jones()
    d = {}
    for e in range(p.minExp(), p.maxExp() + 1):
        k = int(str(p[e]))
        if k:
            d[-2 * e] = d.get(-2 * e, 0) + k
    return {e: k for e, k in d.items() if k}


# ---------------------------------------------------------------------------
# CALIBRATED CONVENTION — measured 2026-08-26, not derived.
#
#   Against regina.bracket() on three knots, all four combinations run:
#       knot            (1,1)   (1,2)   (2,1)   (2,2)
#       trefoil right   HITS    HITS      --      --
#       trefoilLeft     HITS      --    HITS      --
#       figureEight     HITS      --      --    HITS
#   ⟹ only (1,1) carries all three; every wrong choice fails on at least one.
#   Chirality separates column 3, the figure-eight separates column 2 — the
#   test is thus not blind (a test that lets everything through cannot be told apart from one
#   that rejects everything).
#   Finding: the A-smoothing in regina's pdData is, INDEPENDENT OF SIGN, the
#   pairing (0,1)(2,3). That was NOT predictable — the PD order is
#   defined via the orientation.
# ---------------------------------------------------------------------------
WAHL = (1, 1)


def selbsttest():
    ok = True

    def probe(name, bed, extra=''):
        nonlocal ok
        ok &= bool(bed)
        print(f"  {name:<52} {'OK' if bed else '🔴 FAIL'} {extra}")

    # A/B/C — bracket against regina, three knots
    for nm, L in [('A trefoil right', regina.ExampleLink.trefoil()),
                  ('B trefoil left  ', regina.ExampleLink.trefoilLeft()),
                  ('C figure-eight   ', regina.ExampleLink.figureEight())]:
        r = rechne(L, WAHL, neustarts=10)
        probe(f'{nm} bracket == regina.bracket()',
              r['bracket'] == regina_bracket_als_lp(L),
              f"(width {r['schnittbreite']}, states {r['max_zustaende']})")

    # D — unknot
    U = regina.ExampleLink.unknot()
    r = rechne(U, WAHL, neustarts=4)
    probe('D unknot: bracket == 1', r['bracket'] == {0: 1})

    # E COUNTER-CASE — wrong convention MUST fail
    L = regina.ExampleLink.trefoil()
    soll = regina_bracket_als_lp(L)
    schlecht = [w for w in ((1, 2), (2, 1), (2, 2))
                if rechne(L, w, neustarts=10)['bracket'] == soll]
    probe('E wrong convention fails (on at least one knot)',
          len(schlecht) < 3, f'(still equal at {schlecht})')
    L2 = regina.ExampleLink.figureEight()
    soll2 = regina_bracket_als_lp(L2)
    probe('E2 (1,2)/(2,1) fail on the figure-eight',
          rechne(L2, (1, 2), neustarts=10)['bracket'] != soll2 and
          rechne(L2, (2, 1), neustarts=10)['bracket'] != soll2)

    # F COUNTER-CASE — wrong delta MUST fail
    try:
        falsch = rechne(regina.ExampleLink.figureEight(), WAHL, neustarts=10,
                        delta={2: -1, -2: 1})
        f_ok = falsch['bracket'] != soll2
        f_wie = 'different polynomial'
    except ValueError:
        f_ok, f_wie = True, 'not divisible by delta'
    probe('F wrong delta fails', f_ok, f'({f_wie})')

    # G — order invariance: the result must not depend on the order
    L3 = regina.ExampleLink.torus(2, 7)
    L3.simplify()
    pd = L3.pdData()
    basis = rechne(L3, WAHL, neustarts=20)['bracket']
    rng = random.Random(4711)
    gleich = 0
    for _ in range(12):
        o = list(range(len(pd)))
        rng.shuffle(o)
        gleich += (rechne(L3, WAHL, ordnung=o)['bracket'] == basis)
    probe('G 12 random orders, one result', gleich == 12,
          f'({gleich}/12)')

    # H — Jones against regina
    treffer = 0
    proben = [regina.ExampleLink.trefoil(), regina.ExampleLink.trefoilLeft(),
              regina.ExampleLink.figureEight(), regina.ExampleLink.torus(3, 5)]
    for L4 in proben:
        L4.simplify()
        r = rechne(L4, WAHL, neustarts=20)
        v = jones_aus_klammer(r['bracket'], int(L4.writhe()))
        treffer += (v == regina_jones_als_lp_in_A(L4))
    probe('H Jones == regina.jones()', treffer == len(proben),
          f'({treffer}/{len(proben)})')

    # I COUNTER-CASE — does it secretly measure the crossing number?
    K = regina.ExampleLink.trefoil()
    K.composeWith(regina.ExampleLink.trefoil())
    K.composeWith(regina.ExampleLink.trefoil())
    rk = rechne(K, WAHL, neustarts=20)
    r1 = rechne(regina.ExampleLink.trefoil(), WAHL, neustarts=10)
    probe('I 3_1#3_1#3_1 (c=9): states stay small',
          rk['max_zustaende'] <= 4 * r1['max_zustaende'],
          f"(c={rk['c']}, width {rk['schnittbreite']}, "
          f"states {rk['max_zustaende']} against {r1['max_zustaende']})")

    # J — the witnesses themselves must bite
    L5 = regina.ExampleLink.trefoil()
    j5 = jones_aus_klammer(rechne(L5, WAHL, neustarts=10)['bracket'],
                           int(L5.writhe()))
    gut, wie = zeugen_pruefe(j5, det=3)
    probe('J witnesses hold on the trefoil (det=3)', gut, f'({wie})')
    schlecht, _ = zeugen_pruefe(j5, det=5)
    probe('J2 witness fails on WRONG det', not schlecht)
    verbogen = dict(j5)
    e0 = max(verbogen)
    verbogen[e0] = verbogen[e0] + 1
    s3, _ = zeugen_pruefe(verbogen, det=3)
    probe('J3 witness fails on a bent polynomial', not s3)
    # J4: a negative exponent must not slip into floating point
    tst = {-8: 1, -4: -1, 0: 1}
    probe('J4 negative exponents stay integer',
          all(isinstance(v, int) for v in
              (zeugen(tst)['v1'], zeugen(tst)['v_minus1'])))

    print('  ⇒ ' + ('✅ the scaffold holds' if ok else '🔴 do NOT use'))
    return 0 if ok else 1


def zeugen(jones_A):
    """Two independent, CONVENTION-FREE witnesses for a Jones polynomial.

      V(1)   == 1        (every knot)
      |V(-1)| == det     (standard; exactly integer here)

    t = A^-4, so t=1 ⟹ sum of the coefficients; t=-1 ⟹ A^4=-1 ⟹
    A^e = (-1)^(e/4) for e divisible by 4.

    ⚠️ COMPUTED IN INTEGERS, on purpose. `(-1)**(e//4)` gives in Python for
    a NEGATIVE exponent a float (-1.0) — exactly the error that already sat in the
    own witness once on 24.08.2026 (`s1**i`, error log).
    At det = 1.4e13 this is not noticed, at the 23-digit det of
    k6_s3 it would be wrong. Hence a sign via modulo here.
    """
    ungerade = [e for e in jones_A if e % 4]
    if ungerade:
        return {'v1': None, 'v_minus1': None,
                'grund': f'{len(ungerade)} exponents not divisible by 4 '
                         '(link, not a knot?)'}
    v1 = sum(jones_A.values())
    vm1 = 0
    for e, k in jones_A.items():
        vm1 += k if (e // 4) % 2 == 0 else -k
    return {'v1': v1, 'v_minus1': vm1, 'grund': None}


def zeugen_pruefe(jones_A, det=None):
    z = zeugen(jones_A)
    if z['grund']:
        return False, z['grund']
    if z['v1'] != 1:
        return False, f"V(1) = {z['v1']}, not 1"
    if det is not None and abs(z['v_minus1']) != int(det):
        return False, f"|V(-1)| = {abs(z['v_minus1'])} != det {int(det)}"
    return True, f"V(1)=1 · |V(-1)|={abs(z['v_minus1'])}" + (
        ' == det' if det is not None else ' (no det to compare)')


def _det_aus_record(stem, pfad):
    f = f'{ROOT}/{pfad}/{stem}_pflicht_invarianten_complete.json'
    if not os.path.exists(f):
        f = f'{ROOT}/{pfad}/{stem}_pflicht_complete.json'
    if not os.path.exists(f):
        return None
    d = json.load(open(f))
    return d.get('pflichtliste', {}).get('ALGEBRA', {}).get('determinant')


def link_von_stem(stem, pfad):
    f = f'{ROOT}/{pfad}/{stem}_pflicht_invarianten_complete.json'
    if not os.path.exists(f):
        f = f'{ROOT}/{pfad}/{stem}_pflicht_complete.json'
    d = json.load(open(f))
    sg = d.get('signed_gauss') or (d.get('raw_subprocess_outputs', {})
                                   .get('regina', {}).get('signed_gauss_reduced'))
    return regina.Link.fromSignedGauss(sg) if sg else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--selbsttest', action='store_true')
    ap.add_argument('--liste', help='JSON with name/pfad (like homfly_wand.py)')
    ap.add_argument('--stem'); ap.add_argument('--pfad')
    ap.add_argument('--neustarts', type=int, default=40)
    ap.add_argument('--nur-messen', action='store_true',
                    help='only cut width + state count, without bracket')
    ap.add_argument('--riegel', type=int, default=2_000_000,
                    help='max states in count mode')
    ap.add_argument('--laut', action='store_true')
    ap.add_argument('--out')
    a = ap.parse_args()
    if a.selbsttest:
        sys.exit(selbsttest())
    ziele = []
    if a.stem:
        ziele = [{'name': a.stem, 'pfad': a.pfad}]
    elif a.liste:
        ziele = json.load(open(a.liste))
    else:
        sys.exit('🔴 --selbsttest, --stem or --liste required')
    erg = []
    for e in ziele:
        L = link_von_stem(e['name'], e['pfad'])
        if L is None:
            print(f"{e['name']:22} 🔴 no signed_gauss"); continue
        L.simplify()
        try:
            if a.nur_messen:
                pd = L.pdData()
                k, ordn = schnittbreite(pd, a.neustarts)
                t0 = time.time()
                z = zaehle_zustaende(pd, ordn, riegel=a.riegel)
                r = {'name': e['name'], 'c': len(pd), 'schnittbreite': k,
                     'max_zustaende': z['max_zustaende'],
                     'abgebrochen': z['abgebrochen'], 'bei': z['bei'],
                     'hecke_schranke': math.factorial(k),
                     'sek': round(time.time() - t0, 2)}
                print(f"{e['name']:22} c={r['c']:<4} width {k:<3} "
                      f"states {r['max_zustaende']:<9,} "
                      f"(k! would be {r['hecke_schranke']:>18,}) "
                      f"{'ABORTED at '+str(r['bei']) if z['abgebrochen'] else ''}"
                      f" {r['sek']:>7}s")
                erg.append(r)
                continue
            r = rechne(L, WAHL, neustarts=a.neustarts, laut=a.laut)
            r['name'] = e['name']
            r['jones_A'] = {str(k): v for k, v in
                            jones_aus_klammer(r['bracket'], int(L.writhe())).items()}
            r['bracket'] = {str(k): v for k, v in r['bracket'].items()}
            det = e.get('det')
            if det is None:
                det = _det_aus_record(e['name'], e['pfad'])
            jA = {int(k): v for k, v in r['jones_A'].items()}
            r['zeugen_ok'], r['zeugen'] = zeugen_pruefe(jA, det)
            print(f"{e['name']:22} c={r['c']:<4} width {r['schnittbreite']:<3} "
                  f"states {r['max_zustaende']:<8} {r['sek']:>8}s  "
                  f"terms {len(r['jones_A']):<4} "
                  f"{'✅' if r['zeugen_ok'] else '🔴'} {r['zeugen']}")
            erg.append(r)
        except Exception as ex:
            print(f"{e['name']:22} 🔴 {type(ex).__name__}: {ex}")
    if a.out:
        json.dump(erg, open(a.out, 'w'), indent=1)
        print(f'\n  written: {a.out}')


if __name__ == '__main__':
    main()
