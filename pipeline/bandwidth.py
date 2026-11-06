#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""bandwidth.py — measures the CUT WIDTH of a diagram in STRANDS.

WHAT FOR (Disc #714 (d), 26.08.2026). The way out proposed there — "compute
at the point, interpolate afterwards" — stands or falls with a number that nobody has
measured yet. Not the crossing number (#714 (b): irrelevant), not regina's treewidth
(which counts crossings), but: **how many strands does the narrowest cut through the
diagram cross?**

WHY EXACTLY THIS NUMBER. Whoever works through the diagram as a transfer matrix carries
between two steps the state of a cut with k arcs. At a numerical point every basis state
holds ONE number instead of a polynomial — exactly the explosion that kills regina
disappears (Disc #714 (d)). The price is then solely the NUMBER OF STATES per cut.

⚠️ **HOW MANY STATES THAT IS, IS NOT PROVEN HERE.** The obvious bound
is the rank of the Hecke algebra, k!, because the HOMFLY skein of a rectangle
with n inputs and n outputs is H_n (rank n!). Whether a cut with k
arcs and mixed orientation really carries k! states or considerably
fewer, I have **not recomputed**. Therefore k! is printed below as an
**UPPER BOUND** and never as "that is how many there are".
The true number is a MEASUREMENT that only the built run delivers
(it counts the basis elements reached) — not a formula from the head
(`feedback_guessing_is_strictly_forbidden`).
For the order of magnitude, if the bound were sharp:
    k=8  -> 40,320       (0.3 MB)   k=12 ->    479,001,600   (3.8 GB)
    k=10 -> 3,628,800    (29 MB)    k=13 ->  6,227,020,800   (50 GB, dead)
    k=11 -> 39,916,800   (319 MB)   k=14 -> 87,178,291,200   (dead)
⟹ The route carries up to **k=11**, with effort k=12. Hence it is measured first
   and only then built. A number that fits nothing is reported, not
   presented (same doctrine as `breite.py`, #714).

WHAT IS MEASURED. Run over an order of the crossings; after each
prefix the "open" arcs are those that have exactly ONE end in the prefix.
Cut width = maximum over all prefixes. Minimised via:
  · greedy from each start point (up to 40 start points)
  · greedy with random noise (restarts), seed FIXED -> reproducible
This is an UPPER bound: found means found, not optimal.
It is therefore reported as "<= k", never as "the cut width is k".

COUNTER-CASES (--selbsttest, all must bite):
  A trefoil        -> cut width 4 (every cut separates 2 strands = 4 ends)
  B figure-eight   -> small (<= 6)
  C random order is WORSE than the greedy one (otherwise nothing optimises)
  D a diagram with double the crossing number but the same structure must NOT
    double the width (otherwise it measures c and not the width)
"""
import site_config
import json, os, random, sys, time, argparse

sys.path.insert(0, (site_config.FREEZE_DIR))
from pipeline_seal import seal_all
seal_all(verbose=False)
# 🔴 REGINA STAYS A DIRECT IMPORT HERE — with reason (Disc #750, 2026-08-29).
# Regina is GPL v2+, and the pipeline switches calls over to subprocesses. THIS file
# does not: it passes on a DIAGRAM after `simplify()` (signedGauss/PD), and regina's
# simplify() picks, for larger diagrams, a different one of the many equivalent
# results depending on the process. Measured on polygon_resolve.py: same invariants, but
# different Gauss code and in one resolution even a different crossing number.
# Whoever wants to switch it anyway must put ALL calls into ONE subprocess — not one
# process per call. Bridge: pipeline_scripts/regina_bridge.py
import regina

ROOT = site_config.ROOT


def bogen_inzidenz(pd):
    """pdData() -> list of 4-tuples (arc numbers per crossing).
    Returns: (m, inz) with inz[i] = set of arcs at crossing i."""
    return len(pd), [set(t) for t in pd]


def breite_fuer(inz, ordnung):
    """Maximum number of open arcs over all prefixes of this order."""
    zaehler = {}
    for s in inz:
        for b in s:
            zaehler[b] = zaehler.get(b, 0) + 1
    offen, best = {}, 0
    for i in ordnung:
        for b in inz[i]:
            offen[b] = offen.get(b, 0) + 1
        # an arc is closed when ALL its ends are in the prefix
        k = sum(1 for b, v in offen.items() if v < zaehler[b])
        best = max(best, k)
    return best


def gierig(inz, start, rng=None, rausch=0.0):
    m = len(inz)
    zaehler = {}
    for s in inz:
        for b in s:
            zaehler[b] = zaehler.get(b, 0) + 1
    drin, offen, ordnung = set(), {}, []
    akt = start
    while len(drin) < m:
        drin.add(akt); ordnung.append(akt)
        for b in inz[akt]:
            offen[b] = offen.get(b, 0) + 1
        if len(drin) == m:
            break
        # candidates: touch an open arc; otherwise any
        kand = [j for j in range(m) if j not in drin
                and any(b in offen for b in inz[j])]
        if not kand:
            kand = [j for j in range(m) if j not in drin]
        def kosten(j):
            o = dict(offen)
            for b in inz[j]:
                o[b] = o.get(b, 0) + 1
            k = sum(1 for b, v in o.items() if v < zaehler[b])
            if rng is not None and rausch:
                k += rng.random() * rausch
            return k
        akt = min(kand, key=kosten)
    return ordnung


def schnittbreite(pd, neustarts=40, saat=20260826):
    m, inz = bogen_inzidenz(pd)
    if m == 0:
        return 0, []
    rng = random.Random(saat)
    best, beste_ordnung = None, None
    starts = list(range(m))
    rng.shuffle(starts)
    for st in starts[:min(neustarts, m)]:
        for rausch in (0.0, 0.7):
            o = gierig(inz, st, rng, rausch)
            w = breite_fuer(inz, o)
            if best is None or w < best:
                best, beste_ordnung = w, o
    return best, beste_ordnung


def hecke_rang(k):
    """UPPER BOUND for the number of states of a cut with k arcs: k!.
    Not the proven rank — see head. Only computed, not asserted."""
    import math
    return math.factorial(k)


def link_von_stem(stem, pfad):
    f = f'{ROOT}/{pfad}/{stem}_pflicht_invarianten_complete.json'
    if not os.path.exists(f):
        f = f'{ROOT}/{pfad}/{stem}_pflicht_complete.json'
    d = json.load(open(f))
    sg = d.get('signed_gauss') or (d.get('raw_subprocess_outputs', {})
                                   .get('regina', {}).get('signed_gauss_reduced'))
    if not sg:
        return None
    return regina.Link.fromSignedGauss(sg)


def messe(L, neustarts=40):
    roh = int(L.size())
    L.simplify()
    n_simpl = int(L.size())
    t0 = time.time()
    b_tree = int(regina.TreeDecomposition(L).width())
    L.improveTreewidth()
    b_tree2 = int(regina.TreeDecomposition(L).width())
    n_nach = int(L.size())
    pd = L.pdData()
    k, _ = schnittbreite(pd, neustarts)
    return {'roh': roh, 'simpl': n_simpl, 'nach_improve': n_nach,
            'baumweite': b_tree, 'baumweite_improve': b_tree2,
            'schnittbreite': k, 'hecke_rang': hecke_rang(k),
            'sek': round(time.time() - t0, 2)}


def selbsttest():
    ok = True
    kleeblatt = regina.ExampleLink.trefoil()
    a = messe(kleeblatt, 20)
    a_ok = a['schnittbreite'] <= 4
    print(f"  A trefoil  cut width {a['schnittbreite']} (<=4)      : "
          f"{'OK' if a_ok else '🔴 FAIL'}")
    ok &= a_ok
    acht = regina.ExampleLink.figureEight()
    b = messe(acht, 20)
    b_ok = b['schnittbreite'] <= 6
    print(f"  B figure-eight cut width {b['schnittbreite']} (<=6)    : "
          f"{'OK' if b_ok else '🔴 FAIL'}")
    ok &= b_ok
    # C COUNTER-CASE: a RANDOM order must be worse than the greedy one.
    #   Without this probe `schnittbreite` could constantly deliver anything.
    L = regina.ExampleLink.torus(3, 7); L.simplify()
    pd = L.pdData(); m, inz = bogen_inzidenz(pd)
    g, _ = schnittbreite(pd, 20)
    rng = random.Random(7)
    schlecht = []
    for _ in range(30):
        o = list(range(m)); rng.shuffle(o)
        schlecht.append(breite_fuer(inz, o))
    c_ok = g <= min(schlecht)
    print(f"  C greedy {g} <= best of 30 random orders {min(schlecht)} : "
          f"{'OK' if c_ok else '🔴 FAIL — it optimises nothing'}")
    ok &= c_ok
    # D COUNTER-CASE: connected sum trefoil#trefoil has DOUBLE
    #   crossing number, but the cut stays narrow. If it measures the number c, it shows.
    K = regina.ExampleLink.trefoil()
    K.composeWith(regina.ExampleLink.trefoil())
    K.composeWith(regina.ExampleLink.trefoil())
    d = messe(K, 20)
    d_ok = d['schnittbreite'] <= a['schnittbreite'] + 2
    print(f"  D 3_1#3_1#3_1: c={d['simpl']}, width {d['schnittbreite']} "
          f"(<= {a['schnittbreite']}+2)  : "
          f"{'OK' if d_ok else '🔴 FAIL — measures the crossing number'}")
    ok &= d_ok
    print('  ⇒ ' + ('✅ the measurement separates' if ok else '🔴 DO NOT use'))
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--liste', help='JSON with name/pfad (like homfly_wand.py)')
    ap.add_argument('--neustarts', type=int, default=40)
    ap.add_argument('--out')
    ap.add_argument('--selbsttest', action='store_true')
    a = ap.parse_args()
    if a.selbsttest:
        sys.exit(selbsttest())
    if not a.liste:
        sys.exit('🔴 --liste needed')
    erg = []
    print(f"{'Record':22} {'c':>4} {'Treew':>6} {'Cut':>8} {'k! (upper bound)':>20}  {'s':>6}")
    for e in json.load(open(a.liste)):
        L = link_von_stem(e['name'], e['pfad'])
        if L is None:
            print(f"{e['name']:22} 🔴 no signed_gauss"); continue
        r = messe(L, a.neustarts); r['name'] = e['name']; r['c_db'] = e.get('c')
        erg.append(r)
        print(f"{e['name']:22} {r['c_db']:>4} {r['baumweite_improve']:>6} "
              f"{r['schnittbreite']:>8} {r['hecke_rang']:>20,}  {r['sek']:>6}")
    if a.out:
        json.dump(erg, open(a.out, 'w'), indent=1)
        print(f"\n  written: {a.out}")


if __name__ == '__main__':
    main()
