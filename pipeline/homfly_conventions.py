#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""homfly_conventions.py — make HOMFLY strings of different origin comparable.

OCCASION (2026-08-23): the calibration reported 4 of 89 records "off". They were
NOT off — they carry their HOMFLY in other variables (`L,M` or `v,z`),
while regina prints `x,y`. A text comparison across convention boundaries measures
itself (footgun 15).

The substitution is NOT set from memory, but MEASURED:
all candidates of the form x = c1*u^e, y = c2*w are checked at several random
points against regina's AZ; only what hits everywhere counts, and the
others must fail (memory feedback-calibrate-each-convention-alone).

Documented on regina's own numbers (trefoil):
    AZ : -x^4 + x^2 y^2 + 2 x^2
    LM : -x^4 + x^2 y^2 - 2 x^2      ⟹  LM(x,y) = AZ(i*x, i*y)
"""
import re, cmath, itertools

VARPAARE = [('x', 'y'), ('L', 'M'), ('l', 'm'), ('v', 'z'), ('a', 'z')]


def variablen(s):
    """Which variable pair does the string carry? Measured on the text, not guessed."""
    for u, w in VARPAARE:
        if re.search(rf'(?<![A-Za-z]){re.escape(u)}(?![A-Za-z])', s) and \
           re.search(rf'(?<![A-Za-z]){re.escape(w)}(?![A-Za-z])', s):
            return u, w
    return None


def parse(s, u, w):
    """'-L^12*M^18 + 2*L^10' or 'v**12 - 3*v**10*z**2' → {(i,j): coef}."""
    t = s.strip().replace('**', '^').replace('*', ' ')
    # Split ONLY at binary signs. A '-' directly after '^' belongs to the
    # exponent ('v^-2') — whoever splits there tears the term apart and then reports
    # a deviation that is none (happened exactly so on 23.08. with two lj_ records,
    # footgun 15).
    t = re.sub(r'(?<!\^)([+-])', r'|\1', t)
    out = {}
    for st in t.split('|'):
        st = st.strip()
        if not st:
            continue
        vz = -1 if st.startswith('-') else 1
        st = st.lstrip('+-').strip()
        koef, i, j = 1, 0, 0
        gefunden = False
        for stueck in st.split():
            m = re.fullmatch(rf'{re.escape(u)}(?:\^(-?\d+))?', stueck)
            if m:
                i = int(m.group(1)) if m.group(1) else 1; gefunden = True; continue
            m = re.fullmatch(rf'{re.escape(w)}(?:\^(-?\d+))?', stueck)
            if m:
                j = int(m.group(1)) if m.group(1) else 1; gefunden = True; continue
            m = re.fullmatch(r'(\d+)', stueck)
            if m:
                koef *= int(m.group(1)); gefunden = True; continue
            return None
        if not gefunden:
            return None
        out[(i, j)] = out.get((i, j), 0) + vz * koef
    return out or None


def ausw(p, x, y):
    return sum(c * (x ** i) * (y ** j) for (i, j), c in p.items())


def ausw_str(s, u, w, u0, w0):
    """Value of a HOMFLY string at a complex point.

    First the fast monomial parser. If it fails, SymPy takes over — because
    some records carry SymPy output with FRACTIONS (`z**6/v**2`) and
    bracketed exponents (`v**(-6)`). Writing a parser of one's own for that
    would mean rebuilding SymPy and getting it wrong; on
    23.08. exactly that reported two records as "different" that were not
    (footgun 15, third case in one day)."""
    p = parse(s, u, w)
    if p is not None:
        return ausw(p, u0, w0)
    import sympy
    U, W = sympy.symbols(f'{u} {w}')
    e = sympy.sympify(s.replace('^', '**'), locals={u: U, w: W})
    v = complex(sympy.N(e.subs({U: sympy.sympify(u0), W: sympy.sympify(w0)})))
    return v


# Candidates: x = c1 * u^e , y = c2 * w
EINHEITEN = [(1, '1'), (-1, '-1'), (1j, 'i'), (-1j, '-i')]
KANDIDATEN = [((c1, e, c2), f'x={n1}*u^{e}, y={n2}*w')
              for (c1, n1) in EINHEITEN
              for e in (1, -1)
              for (c2, n2) in EINHEITEN]


def finde_substitution(paare, punkte=None, tol=1e-7):
    """paare = [(record_string, regina_az_string), ...] of the same knot.
    Returns the ONE substitution that hits at all points of all pairs
    — or None if none or more than one hits."""
    if punkte is None:
        punkte = [(1.3 + 0.2j, 0.7 - 0.4j), (0.55, 1.9),
                  (2.1 + 0.3j, 0.83 + 0.11j), (0.77 - 0.5j, 1.45)]
    treffer = []
    for (c1, e, c2), name in KANDIDATEN:
        gut = True
        for rec_s, az_s in paare:
            vp = variablen(rec_s)
            if vp is None:
                gut = False; break
            pa = parse(az_s, 'x', 'y')
            if pa is None:
                gut = False; break
            for (u0, w0) in punkte:
                x = c1 * (u0 ** e)
                y = c2 * w0
                try:
                    lhs = ausw_str(rec_s, vp[0], vp[1], u0, w0)
                except Exception:
                    gut = False; break
                rhs = ausw(pa, x, y)
                nen = max(abs(lhs), abs(rhs), 1e-12)
                if abs(lhs - rhs) / nen > tol:
                    gut = False; break
            if not gut:
                break
        if gut:
            treffer.append(((c1, e, c2), name))
    return treffer


def gleich_ueber_konventionen(rec_s, az_s, punkte=None, tol=1e-7):
    """Does the record carry the same polynomial as the new regina AZ — regardless of
    which convention? Returns (True/False, name_of_substitution)."""
    tr = finde_substitution([(rec_s, az_s)], punkte, tol)
    if not tr:
        return False, None
    return True, tr[0][1]


def selbsttest():
    """Four probes, two must REJECT (footgun 16)."""
    ok = True
    az = '-x^4 + x^2 y^2 + 2 x^2'
    lm = '-x^4 + x^2 y^2 - 2 x^2'
    g, nm = gleich_ueber_konventionen(lm.replace('x', 'L').replace('y', 'M'), az)
    print(f'  Probe 1 LM↔AZ trefoil        : {g} ({nm})',
          '✅' if g else '🔴')
    ok &= bool(g)
    g2, _ = gleich_ueber_konventionen(az, az)
    print(f'  Probe 2 AZ↔AZ (identical)    : {g2}', '✅' if g2 else '🔴')
    ok &= bool(g2)
    # COUNTER-CASE 1: figure-eight knot must not pass as trefoil
    vier = '- x^-2 + 3 - x^2 - y^2'
    g3, _ = gleich_ueber_konventionen(vier.replace('x', 'L').replace('y', 'M'), az)
    print(f'  Probe 3 COUNTER-CASE 4_1 vs 3_1 : {g3}', '✅' if not g3 else '🔴')
    ok &= not g3
    # COUNTER-CASE 2: a falsified coefficient must not pass
    g4, _ = gleich_ueber_konventionen('-L^4 + L^2 M^2 - 3 L^2', az)
    print(f'  Probe 4 COUNTER-CASE coef. 2→3  : {g4}', '✅' if not g4 else '🔴')
    ok &= not g4
    print('  ⇒', '✅ the convention calibrator bites' if ok else '🔴 failed')
    return 0 if ok else 1


if __name__ == '__main__':
    import sys
    sys.exit(selbsttest())
