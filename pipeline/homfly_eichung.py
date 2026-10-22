#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""homfly_eichung.py — calibrates our OWN HOMFLY (Hecke trace, at a point) against regina.

WHY: regina's homfly is the one call that builds the wall for high-crossing knots
(measured at c=84: everything else together 20 s, homfly > 120 s). Our own
tool hecke_homfly_fast.py computes the same quantity in fractions of a second —
but in ITS normalisation ("×common norm K", author's comment). Uncalibrated,
a fast number is worthless.

Model: alex_eichung.py (Disc #570). There variants 0-5 failed completely
and only 6 hit — so the test COULD reject. The same here:
several conventions are offered, the data decide, and a
variant that fits only partially counts as failed.
"""
import sys, os, re, json, argparse, ast
import site_config
from fractions import Fraction

sys.path.insert(0, site_config.HECKE_DIR)


# ── evaluate a regina HOMFLY string ───────────────────────────────────────
_TERM = re.compile(r'([+-]?)\s*(\d+)?\s*(x(?:\^(-?\d+))?)?\s*(y(?:\^(-?\d+))?)?')

def parse_homfly(s):
    """'x^-2 y^2 + 2 x^-2 - x^-4' → {(i,j): coef} with x^i y^j."""
    s = s.strip()
    if not s:
        return None
    # split terms at binary +/- (not at an exponent minus)
    s = re.sub(r'\s*([+-])\s+', r' \1 ', s)
    stuecke = re.split(r'\s(?=[+-]\s)', s)
    out = {}
    for st in stuecke:
        st = st.strip()
        if not st:
            continue
        vz = -1 if st.startswith('-') else 1
        st = st.lstrip('+-').strip()
        m = re.match(r'^(\d+)?\s*(?:x(?:\^(-?\d+))?)?\s*(?:y(?:\^(-?\d+))?)?$', st)
        if not m:
            return None
        koef = int(m.group(1)) if m.group(1) else 1
        hat_x = 'x' in st
        hat_y = 'y' in st
        i = int(m.group(2)) if m.group(2) else (1 if hat_x else 0)
        j = int(m.group(3)) if m.group(3) else (1 if hat_y else 0)
        out[(i, j)] = out.get((i, j), 0) + vz * koef
    return out


def eval_homfly(poly, x, y):
    return sum(c * (x ** i) * (y ** j) for (i, j), c in poly.items())


# ── conventions under which (a,q) is mapped to (x,y) ──────────────────────
VARIANTEN = {
    0: ('x=a,    y=q-1/q',   lambda a, q: (a,     q - 1 / q)),
    1: ('x=1/a,  y=q-1/q',   lambda a, q: (1 / a, q - 1 / q)),
    2: ('x=a,    y=1/q-q',   lambda a, q: (a,     1 / q - q)),
    3: ('x=1/a,  y=1/q-q',   lambda a, q: (1 / a, 1 / q - q)),
    4: ('x=a*q^0 (a<->q)',   lambda a, q: (q,     a - 1 / a)),
    5: ('x=1/q,  y=a-1/a',   lambda a, q: (1 / q, a - 1 / a)),
}


def unsere(w, n, a, q, procs):
    from hecke_homfly_fast import homfly_hecke
    return homfly_hecke(list(w), n, complex(a), complex(q), procs)


def lauf(faelle, a, q, procs, nur=None):
    erg = {}
    for vid, (name, sub) in VARIANTEN.items():
        if nur is not None and vid != nur:
            continue
        x, y = sub(complex(a), complex(q))
        zeilen = []
        for f in faelle:
            p = parse_homfly(f['homfly'])
            if p is None:
                zeilen.append((f, None, None, None)); continue
            r_reg = eval_homfly(p, x, y)
            r_own = unsere(f['w'], f['n'], a, q, procs)
            quot = (r_reg / r_own) if abs(r_own) > 1e-14 else None
            zeilen.append((f, r_reg, r_own, quot))
        erg[vid] = (name, zeilen)
    return erg


# ── normalisation K: the author left it out on purpose ("×common norm K") ──
# FOUND from the data on 23.08., not derived: for variant 0 the
# quotient regina/own is CONSTANT over all 14 calibration cases — also across n=2,3,4,
# so not dependent on the strand count at all. The value 0.690821 at
# (a=1.3, q=1.2) is exactly (q-1/q)/(a-1/a).
def norm_K(a, q):
    return (q - 1 / q) / (a - 1 / a)


def pruefe(faelle, punkte, procs, toleranz=1e-9):
    """Calibration at SEVERAL points. A constant fitted at one point
    proves no functional form — only agreement at further
    points does that (feedback-agreement-is-no-proof-below-resolution)."""
    gesamt_ok = True
    for vid, (vname, sub) in sorted(VARIANTEN.items()):
        treffer = fehl = 0
        schlimmster = 0.0
        for (a, q) in punkte:
            x, y = sub(complex(a), complex(q))
            K = norm_K(complex(a), complex(q))
            for f in faelle:
                p = parse_homfly(f['homfly'])
                if p is None:
                    fehl += 1; continue
                r_reg = eval_homfly(p, x, y)
                r_own = unsere(f['w'], f['n'], a, q, procs) * K
                nen = max(abs(r_reg), abs(r_own), 1e-12)
                rel = abs(r_reg - r_own) / nen
                schlimmster = max(schlimmster, rel)
                if rel < toleranz: treffer += 1
                else: fehl += 1
        ok = (fehl == 0)
        gesamt_ok &= (ok if vid == 0 else True)
        marke = '✅' if ok else '🔴'
        print(f"  {marke} variant {vid} ({vname:<20}): {treffer} hits, {fehl} misses, "
              f"largest rel. deviation {schlimmster:.2e}")
        if vid == 0 and not ok:
            print("     🔴 the calibrated variant does NOT hit — calibration invalid.")
    return 0 if gesamt_ok else 1


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--faelle', required=True)
    ap.add_argument('--a', type=float, default=1.3)
    ap.add_argument('--q', type=float, default=1.2)
    ap.add_argument('--procs', type=int, default=12)
    ap.add_argument('--variante', type=int)
    ap.add_argument('--pruefen', action='store_true',
                    help='check all variants at several points against the norm K')
    a = ap.parse_args()
    faelle = json.load(open(a.faelle))
    if a.pruefen:
        punkte = [(1.3, 1.2), (0.7, 1.45), (2.1, 0.83), (1.11, 1.07)]
        print(f"Calibration at {len(punkte)} points x {len(faelle)} cases "
              f"= {len(punkte)*len(faelle)} probes per variant")
        sys.exit(pruefe(faelle, punkte, a.procs))
    for vid, (name, zeilen) in lauf(faelle, a.a, a.q, a.procs, a.variante).items():
        print(f"\n=== variant {vid}: {name} ===")
        print(f"  {'knot':<32}{'n':>3}  {'regina':>14}  {'own':>14}  {'quotient':>16}")
        for f, rr, ro, qu in zeilen:
            fmt = lambda z: ('—' if z is None else f"{z.real:14.6f}"
                             if abs(getattr(z, 'imag', 0)) < 1e-9 else f"{z:.4f}")
            print(f"  {f['name'][:32]:<32}{f['n']:>3}  {fmt(rr)}  {fmt(ro)}  {fmt(qu)}")
