#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""CALIBRATION of the union constructor against the source: R^s_2 (Lemma 4.3) and R^s_3 (Lemma 4.4).

Delta is only determined up to a unit +-t^k (Footgun 12) ==> what is compared is the
COEFFICIENT SEQUENCE up to sign and reversal, not the string.
A convention choice only counts if it hits R^s_2 EXACTLY -- and there must be choices
that FAIL on it, otherwise the calibration distinguishes nothing (Footgun 16).
"""
import site_config
import itertools, re, sys
HIER = (site_config.ROOT + '/maschinen_ideen/RS4_12N706_20260919')
sys.path.insert(0, HIER)
from abschluss import aglet_ueber
from band_rand import alexander_von, strands

RS2 = [[1, -2, 3, -2, 1], [2, -5, 2]]
RS3 = [[1], [1, -1, -1, 3, -1, -1, 1], [1, 0, -3, 0, 1],
       [1, -3, 5, -7, 5, -3, 1], [1, -6, 11, -6, 1], [2, -6, 9, -6, 2]]


def koeffs(s):
    """Polynomial string -> normalized coefficient sequence (up to +-t^k and reversal)."""
    if s is None or not isinstance(s, str) or s.startswith("LINK"):
        return None
    s = s.replace(" ", "").replace("-", "+-")
    d = {}
    for term in [t for t in s.split("+") if t]:
        m = re.fullmatch(r"(-?\d*)\*?t\^?(-?\d*)", term)
        if m:
            k = m.group(1); e = m.group(2)
            k = -1 if k == "-" else (1 if k == "" else int(k))
            e = 1 if e == "" else int(e)
        else:
            m2 = re.fullmatch(r"(-?\d+)", term)
            if not m2:
                return None
            k, e = int(m2.group(1)), 0
        d[e] = d.get(e, 0) + k
    if not d:
        return None
    lo, hi = min(d), max(d)
    v = [d.get(e, 0) for e in range(lo, hi + 1)]
    while v and v[0] == 0: v.pop(0)
    while v and v[-1] == 0: v.pop()
    if not v: return None
    kand = [v, v[::-1], [-x for x in v], [-x for x in v][::-1]]
    return min(tuple(x) for x in kand)


def normmenge(liste): return {koeffs_von_liste(x) for x in liste}
def koeffs_von_liste(v):
    kand = [v, v[::-1], [-x for x in v], [-x for x in v][::-1]]
    return min(tuple(x) for x in kand)


def faelle(wort):
    """All mod-2-labeled knotoid diagrams for one arc diagram."""
    w = [int(c) for c in wort]
    fest, vork = aglet_ueber(w)
    innere = sorted(c for c in vork if c not in fest)
    for wahl in itertools.product(*[[vork[c][0], vork[c][1]] for c in innere]):
        ueber = dict(fest); ueber.update(dict(zip(innere, wahl)))
        _, n_str = strands(w, ueber)
        frei = [j for j in range(n_str) if j not in (0, n_str - 1)]
        for bits in itertools.product((0, 1), repeat=len(frei)):
            lab = dict(zip(frei, bits))
            yield w, ueber, lab, dict(zip(innere, wahl)), n_str


if __name__ == "__main__":
    soll2 = {koeffs_von_liste(v) for v in RS2}
    soll3 = {koeffs_von_liste(v) for v in RS3}
    print("=" * 96)
    print("CALIBRATION of the union constructor against R^s_2 (Lemma 4.3) and R^s_3 (Lemma 4.4)")
    print("=" * 96)
    gewinner = []
    for flag_sing, flag_seite, flag_eintritt in itertools.product((True, False), (True, False), (0, 1)):
        if True:
            ist2, n2 = set(), 0
            for w, ueber, lab, innere, n_str in faelle("1212"):
                n2 += 1
                for s in alexander_von(w, ueber, lab, flag_sing, flag_seite, flag_eintritt):
                    ist2.add(koeffs(s))
            treffer2 = (ist2 == soll2)
            print(f"\n--- flag_sing={flag_sing}  flag_seite={flag_seite}  flag_eintritt={flag_eintritt} ---")
            print(f"  c=2: {n2} mod-2 diagrams (paper says 2)")
            fehl = sum(1 for x in ist2 if x is None)
            print(f"       generated: {sorted(x for x in ist2 if x)}"
                  + (f"   [+{fehl} without DT — NOT concealed]" if fehl else ""))
            print(f"       target : {sorted(soll2)}")
            print(f"       ==> R^s_2 {'HIT' if treffer2 else 'missed'}")
            if treffer2:
                ist3, n3 = set(), 0
                for wort in ("123123", "123132"):
                    for w, ueber, lab, innere, n_str in faelle(wort):
                        n3 += 1
                        for s in alexander_von(w, ueber, lab, flag_sing, flag_seite, flag_eintritt):
                            ist3.add(koeffs(s))
                drin = {x for x in ist3 if x} <= soll3
                print(f"  c=3: {n3} mod-2 diagrams, {len({x for x in ist3 if x})} polynomials")
                print(f"       all in R^s_3? {'YES' if drin else 'NO'}")
                fehlt = soll3 - {x for x in ist3 if x}
                print(f"       not generated: {sorted(fehlt) if fehlt else '(nothing — complete)'}")
                gewinner.append((flag_sing, flag_seite, flag_eintritt, drin))
    print("\n" + "=" * 96)
    if len(gewinner) == 1:
        print(f"✅ EXACTLY ONE convention hits R^s_2: flag_sing={gewinner[0][0]}, "
              f"flag_seite={gewinner[0][1]}, flag_eintritt={gewinner[0][2]} — the calibration DISTINGUISHES.")
    elif not gewinner:
        print("🔴 NO convention hits R^s_2 — the constructor is still wrong.")
    else:
        print(f"⚠️ {len(gewinner)} conventions hit R^s_2 — the calibration does not distinguish "
              f"sharply enough; c=3 must decide: {gewinner}")
