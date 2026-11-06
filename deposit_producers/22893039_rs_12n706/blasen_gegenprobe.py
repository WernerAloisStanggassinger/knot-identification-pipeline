#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""COUNTER-CHECK on the two diagrams that only the LITERAL reading of condition (2) retains.

The second implementation (`zweite_implementierung.py`) shows: at c=5 the
two readings of condition (2) differ by exactly two arc diagrams, `1213453452` and `1234534512`.
Both carry the bubble `345345` in the middle — a circle meets lambda there in TWO points,
not in one. The prose of the paper says that in such a case the knotoid is "not prime", and
R^s_n is explicitly declared for PRIME knots. Then the two may be missing.

Measured instead of derived, exactly this is done here:
  · Which Alexander polynomials do the two diagrams yield?
  · Do they split (product of two smaller Alexander polynomials = composite knot)?
  · Does one of them NOT lie in the 145 of R^s_5? Only then would our set be too small.

Usage: python3 blasen_gegenprobe.py
"""
import itertools
import json
import os
import sys
import time

HIER = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HIER)

BLASEN = ["1213453452", "1234534512"]
FRIST = 90
ARBEITER = 6


def aufgaben_fuer(wort):
    from abschluss import aglet_ueber
    from band_rand import strands
    w = [int(c) for c in wort]
    fest, vork = aglet_ueber(w)
    innere = sorted(c for c in vork if c not in fest)
    aus = []
    for wahl in itertools.product(*[[vork[c][0], vork[c][1]] for c in innere]):
        ueber = dict(fest)
        ueber.update(dict(zip(innere, wahl)))
        _, n_str = strands(w, ueber)
        frei = [j for j in range(n_str) if j not in (0, n_str - 1)]
        for bits in itertools.product((0, 1), repeat=len(frei)):
            aus.append((wort, {str(a): b for a, b in ueber.items()},
                        {str(a): b for a, b in zip(frei, bits)}))
    return aus


def arbeite(auf):
    from eiche_union import koeffs
    from band_rand import alexander_von
    wort, ueber, lab = auf
    w = [int(c) for c in wort]
    ueber = {int(k): v for k, v in ueber.items()}
    lab = {int(k): v for k, v in lab.items()}
    aus = []
    for s in alexander_von(w, ueber, lab, True, True):
        k = koeffs(s)
        aus.append(list(k) if k else None)
    return aus


def in_kind(auf):
    r, w = os.pipe()
    pid = os.fork()
    if pid == 0:
        try:
            os.close(r)
            os.write(w, json.dumps(arbeite(auf)).encode())
            os.close(w)
        except BaseException:
            pass
        finally:
            os._exit(0)
    os.close(w)
    return pid, r


def norm(v):
    k = [list(v), list(v)[::-1], [-x for x in v], [-x for x in v][::-1]]
    return min(tuple(x) for x in k)


def det(v):
    return abs(sum(c * (-1) ** j for j, c in enumerate(v)))


def faktorisiert(v):
    """Does the polynomial over Q split into more than one nontrivial factor?"""
    try:
        from sympy import Poly, factor_list, symbols
        t = symbols('t')
        p = Poly(list(reversed(list(v))), t)
        _, faktoren = factor_list(p.as_expr(), t)
        grade = []
        for f, e in faktoren:
            g = Poly(f, t).degree()
            grade += [g] * e
        return sorted(g for g in grade if g > 0)
    except Exception as e:                      # sympy missing or special case
        return f"nicht bestimmt ({e})"


def main():
    alle = []
    for wort in BLASEN:
        alle += aufgaben_fuer(wort)
    print(f"{len(alle)} tasks from {len(BLASEN)} bubble diagrams", flush=True)

    menge, laufend, i, erledigt, tot = {}, {}, 0, 0, 0
    t0 = time.time()
    while i < len(alle) or laufend:
        while len(laufend) < ARBEITER and i < len(alle):
            pid, r = in_kind(alle[i])
            laufend[pid] = (r, time.time(), alle[i][0])
            i += 1
        for pid in list(laufend):
            r, ts, wort = laufend[pid]
            fertig, _ = os.waitpid(pid, os.WNOHANG)
            if fertig == 0:
                if time.time() - ts > FRIST:
                    os.kill(pid, 9)
                    os.waitpid(pid, 0)
                    os.close(r)
                    del laufend[pid]
                    erledigt += 1
                    tot += 1
                continue
            daten = b""
            while True:
                st = os.read(r, 65536)
                if not st:
                    break
                daten += st
            os.close(r)
            del laufend[pid]
            erledigt += 1
            try:
                for k in json.loads(daten.decode()):
                    if k:
                        menge.setdefault(norm(k), set()).add(wort)
            except Exception:
                tot += 1
        time.sleep(0.01)

    print(f"done: {erledigt} tasks, {tot} without result, {len(menge)} polynomials, "
          f"{time.time()-t0:.0f}s\n", flush=True)

    rs5 = {tuple(p) for p in json.load(open(os.path.join(HIER, 'rs5_ergebnis.json')))['polynome']}
    neu = [p for p in menge if p not in rs5]
    print(f"R^s_5 (first computation) has {len(rs5)} polynomials.")
    print(f"Of these, {len(neu)} of the bubble polynomials are NOT contained.\n")
    for p in sorted(menge):
        print(f"  det {det(list(p)):4d}  degrees of the factors {faktorisiert(p)}  "
              f"{'NEW (not in R^s_5)' if p in neu else 'already in R^s_5'}  {list(p)}")


if __name__ == '__main__':
    main()
