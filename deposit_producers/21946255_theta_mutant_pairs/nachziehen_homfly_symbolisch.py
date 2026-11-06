#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Compute HOMFLY + carry the 30 new mutant pairs over SYMBOLICALLY. `sage_stable`.

TWO GAPS FROM THE FIRST RUN ARE CLOSED HERE:
 (a) HOMFLY was NOT computed for the new pairs — "classically indistinguishable"
     rested there partly on theory instead of measurement.
 (b) The separation stood only at three points. For a letter it must be as
     an inequality of FUNCTIONS.

HOMFLY via Sage's `Link.homfly_polynomial()` (libhomfly), calibrated both ways:
  · the corpus pair P(2,3,5,7) must yield EQUAL HOMFLY  -> ✓
  · two different knots must yield DIFFERENT                  -> ✓
(The pipeline otherwise uses regina.homfly(); it explodes exponentially at high c,
`aggressive_simplify.py` line 110. Here it is only about EQUAL/DIFFERENT, for that is
every correct implementation is good — as long as it is calibrated.)

SYMBOLIC in three stages, as with the corpus pairs:
  1. f1 - f2 == 0                      -> identical
  2. f1/f2 == unit ±t1^a*t2^b       -> equal up to normalisation, NO separation
  3. otherwise                          -> real separation
Stage 2 is the important one: theta_full is Delta-multiplied, Delta only up to ±t^k
determined (footgun 12).

⚠️ Ascending by crossing number, and **saved atomically after every pair**
(tmp + os.replace) — a run that writes only at the end is on this box a
total-loss risk (footgun 28). Aborting is allowed at any time: what is through, stands.
"""
import site_config
import json
import os
import sys
import time

sys.path.insert(0, (site_config.ROOT + '/freeze_manifest/pipeline_scripts'))
from pd_to_theta_input import pd_to_cs_phi

from sage.all_cmdline import *   # noqa: F401,F403
from sage.knots.link import Link as SageLink

BN = (site_config.ROOT + '/external_tools/bar_natan_theta')
load(os.path.join(BN, "Theta.sage"))

HIER = os.path.dirname(os.path.abspath(__file__))
ZIEL = os.path.join(HIER, "ergebnis_nachziehen.json")
T1, T2 = var("t1"), var("t2")


def homfly(pd):
    return str(SageLink([list(map(int, c)) for c in pd]).homfly_polynomial())


def theta_symbolisch(pd):
    cs, phi = pd_to_cs_phi([tuple(c) for c in pd])
    return Theta([matrix(ZZ, [list(r) for r in cs]), vector(ZZ, phi)])


def vergleich(f1, f2):
    """(verdict, note). Stage 1 -> 2 -> 3."""
    d = (f1 - f2).simplify_full() if hasattr(f1, "simplify_full") else (f1 - f2)
    if d == 0:
        return "identisch", "f1-f2 == 0"
    # 🔴 do NOT use .factor(): it returns a `Factorization` object, and that has
    # no .numerator() — on 14.08. the first attempt ran 6x into an AttributeError with it,
    # which was booked as "theta_fehler" and looked like a computation problem.
    # Theta returns an element of Frac(QQ[t1,t2]); there numerator/denominator exist
    # directly, and those are polynomials with is_monomial().
    q = f1 / f2
    try:
        num, den = q.numerator(), q.denominator()
        einheit = (num.is_monomial() and den.is_monomial()
                   and abs(num.lc() / den.lc()) == 1)
    except Exception as ex:
        einheit = False
        return "unentschieden", "Einheiten-Test nicht durchfuehrbar: %s" % type(ex).__name__
    if einheit:
        return "gleich_bis_auf_einheit", "f1/f2 = %s" % str(qs)[:60]
    return "echte_trennung", "Quotient kein Monom"


def atomar(obj):
    tmp = ZIEL + ".tmp"
    with open(tmp, "w") as f:
        # Belt and braces: in case a Sage type slips through after all.
        json.dump(obj, f, ensure_ascii=False, indent=1,
                  default=lambda o: float(o) if hasattr(o, "__float__") else str(o))
    os.replace(tmp, ZIEL)


def main():
    kand = json.load(open(os.path.join(HIER, "mutanten_kandidaten.json")))
    gut = sorted((x for x in kand if x.get("brauchbar")),
                 key=lambda y: y["K"]["crossings"])
    print("== %d pairs, ascending by c ==" % len(gut), flush=True)

    erg = []
    for i, x in enumerate(gut, 1):
        nK, nM, c = x["K"]["name"], x["Kmut"]["name"], x["K"]["crossings"]
        e = {"K": nK, "Kmut": nM, "c": c, "det": x["K"]["determinant"]}
        t0 = time.time()
        try:
            hK, hM = homfly(x["K"]["pd"]), homfly(x["Kmut"]["pd"])
            e["homfly_gleich"] = (hK == hM)
        except Exception as ex:
            e["homfly_gleich"] = None
            e["homfly_fehler"] = "%s: %s" % (type(ex).__name__, str(ex)[:120])
        e["homfly_sekunden"] = float(round(time.time() - t0, 1))   # float(): Sages round()
        # returns a RealDoubleElement, and that breaks json.dump (footgun 2, #358->#373).

        t1 = time.time()
        try:
            f1, f2 = theta_symbolisch(x["K"]["pd"]), theta_symbolisch(x["Kmut"]["pd"])
            v, notiz = vergleich(f1, f2)
            e["symbolisch"] = v
            e["notiz"] = notiz
        except Exception as ex:
            e["symbolisch"] = None
            e["theta_fehler"] = "%s: %s" % (type(ex).__name__, str(ex)[:120])
        e["theta_sekunden"] = float(round(time.time() - t1, 1))

        erg.append(e)
        atomar(erg)
        print("  [%2d/%d] c=%2d %-15s vs %-15s  HOMFLY %-12s  symbolic %-24s (%.0fs+%.0fs)"
              % (i, len(gut), c, nK, nM,
                 {True: "EQUAL", False: "DIFFERENT", None: "ERROR"}[e["homfly_gleich"]],
                 e.get("symbolisch") or e.get("theta_fehler", "?")[:24],
                 e["homfly_sekunden"], e["theta_sekunden"]), flush=True)

    hg = sum(1 for e in erg if e.get("homfly_gleich") is True)
    tr = sum(1 for e in erg if e.get("symbolisch") == "echte_trennung")
    ok = sum(1 for e in erg if e.get("symbolisch") is not None)
    print("\n== HOMFLY equal: %d of %d ==" % (hg, len(erg)), flush=True)
    print("== computed symbolically: %d · of these real separation: %d ==" % (ok, tr), flush=True)
    if hg != len(erg):
        print("   ⚠️ not all HOMFLY equal — those affected are NO LONGER pure "
              "theta test cases (classically already separable).", flush=True)


if __name__ == "__main__":
    main()
