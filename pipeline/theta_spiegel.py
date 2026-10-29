#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Knot-Identification Pipeline — Werner Alois Stanggassinger & Claude Code (Anthropic).
"""Decompose theta_full — without sympy, CALIBRATED against the corpus.

WHY NO SYMPY: Disc #806 measured that sympy hits the wall on 1.2-MB theta strings
(the same run 9 s instead of >9 min).

WHY CALIBRATED AND NOT BELIEVED: the first draft (19.09. 11:3x) treated
every string that does not begin with '(' as a pure numerator. 149 of the 832 records
however have a PREFACTOR ('2*(...)/(...)'). For those it SILENTLY delivered garbage
— no exception, only wrong monomials. It was caught not by the hand test (six
cases, all beginning with '('), but by the wide run over all 832.

THE CALIBRATION: every parsed polynomial is evaluated at (t1,t2) = (22/7, 21/13)
and checked against the field `theta_at_22_7__21_13` stored in the record. That is
an INDEPENDENT reference — it comes from the pipeline, not from this parser.
If it matches exactly for all 832, the parser is proven instead of asserted.
"""
import re
from fractions import Fraction

_MONOM = re.compile(r'^t(\d+)(?:\^(\d+))?$')


def _monome(ausdruck):
    """'a*t1^2*t2 - 3*t1 + 1' -> {monomial_key: coeff}. No parentheses allowed."""
    terme = {}
    for st in re.findall(r'[+-]?[^+-]+', ausdruck.replace(' ', '')):
        st = st.strip()
        if not st:
            continue
        vz = -1 if st.startswith('-') else 1
        st = st.lstrip('+-')
        if not st:
            continue
        koeff, mon = 1, []
        for t in st.split('*'):
            if not t:
                continue
            if t.lstrip('-').isdigit():
                koeff *= int(t)
            elif _MONOM.match(t):
                mon.append(t)
            else:
                raise ValueError(f"unknown piece {t!r} in {st!r}")
        schluessel = '*'.join(sorted(mon))
        terme[schluessel] = terme.get(schluessel, 0) + vz * koeff
    return {k: v for k, v in terme.items() if v != 0}


def zerlege(s):
    """'K*(NUMERATOR)/(DENOMINATOR)' or '(NUMERATOR)/(DENOMINATOR)' or 'NUMERATOR'
    -> (dict{monom: coeff}, denominator_monomials) or None."""
    if not s or not isinstance(s, str):
        return None
    s = s.strip()
    # Denominator: split at the LAST '/(' (the numerator may itself have parentheses)
    i = s.rfind('/(')
    if i != -1 and s.endswith(')'):
        zaehler_teil, nenner_teil = s[:i], s[i + 2:-1]
    else:
        zaehler_teil, nenner_teil = s, ''
    try:
        nenner = _monome(nenner_teil) if nenner_teil.strip() else {'': 1}
    except ValueError:
        return None
    if len(nenner) != 1:
        return None                      # the denominator is always ONE monomial
    nmon, nkoeff = next(iter(nenner.items()))

    # Numerator: optional prefactor before the parenthesis -> 'K*(...)'
    zaehler_teil = zaehler_teil.strip()
    vorfaktor = 1
    m = re.match(r'^([+-]?\d+)\s*\*\s*\((.*)\)$', zaehler_teil, re.S)
    if m:
        vorfaktor, zaehler_teil = int(m.group(1)), m.group(2)
    elif zaehler_teil.startswith('(') and zaehler_teil.endswith(')'):
        zaehler_teil = zaehler_teil[1:-1]
    try:
        terme = _monome(zaehler_teil)
    except ValueError:
        return None
    if vorfaktor != 1:
        terme = {k: v * vorfaktor for k, v in terme.items()}
    return terme, (nmon, nkoeff)


def auswerten(zerlegt, t1, t2):
    """Evaluate polynomial at (t1, t2) — exactly, via Fraction."""
    terme, (nmon, nkoeff) = zerlegt
    def mon_wert(schluessel):
        w = Fraction(1)
        if schluessel:
            for t in schluessel.split('*'):
                m = _MONOM.match(t)
                basis = t1 if m.group(1) == '1' else t2
                w *= basis ** int(m.group(2) or 1)
        return w
    zaehler = sum(Fraction(v) * mon_wert(k) for k, v in terme.items())
    return zaehler / (Fraction(nkoeff) * mon_wert(nmon))


def ist_negiert(a, b):
    """theta_full a == -b?  ⚠️ SUFFICIENT, NOT NECESSARY: different denominators
    are counted as 'not negated'. A missed case does NOT collapse —
    that is the safe direction (better count two than wrongly one)."""
    za, zb = zerlege(a), zerlege(b)
    if za is None or zb is None:
        return False
    (ta, na), (tb, nb) = za, zb
    if na != nb or not ta or not tb:
        return False
    return ta == {k: -v for k, v in tb.items()}


# ────────────────────────── Mirror detection ──────────────────────────
def theta_punkt(m):
    """Bare theta at the point (22/7, 21/13) as Fraction — or None.

    ⚠️ DIFFERENT NORMALISATION THAN theta_full (Memory project_theta_two_normalizations_in_corpus):
    theta_at_22_7__21_13 is the BARE theta, theta_full is theta·Δ(t1)Δ(t2)Δ(t1t2).
    Harmless here, because it is ALWAYS like against like (point against point,
    full against full) — never crosswise. Whoever changes that produces spurious findings.
    """
    s = m.get('theta_at_22_7__21_13')
    if s in (None, ''):
        return None
    try:
        return Fraction(str(s))
    except (ValueError, ZeroDivisionError):
        return None


def sind_spiegel(d, e):
    """Are d and e mirror images? -> 'sig' | 'theta' | None.  (Disc #808)

    TWO WAYS, because one alone has a hole:
      'sig'   : sigma = -sigma' and sigma != 0. The classical way.
      'theta' : sigma = 0 on BOTH sides — then the signature does not separate and the
                geometry fundamentally never separates mirrors (Gordon-Luecke / Disc #798).
                theta is the only thing left (Disc #807).

    THREE COUNTER-CASES that this rule MUST pass (gegenfall_spiegel.py):
      1. conway_cw vs conway_kt — Conway and Kinoshita-Terasaka. det=1, sigma=0,
         volume equal to 13 digits, but MUTANTS (genus 3 vs 2), not a mirror pair.
         A rule without theta check would merge them.
      2. theta = 0 (amphichiral). 0 = -0 is trivially true; two DIFFERENT amphichiral
         knots with the same (det, vol) would collapse. Hence theta != 0 is demanded.
      3. a knot with itself — theta = -theta holds only for theta = 0, see 2.

    The point evaluation is the cheap FILTER, the full polynomial the expensive CONFIRMATION.
    Point evaluation alone is NOT enough (Memory feedback_agreement_is_no_proof_below_resolution):
    a single point can match by chance.
    """
    try:
        si, ei = int(d.get('sig')), int(e.get('sig'))
    except (TypeError, ValueError):
        return None
    if si == -ei and si != 0:
        return 'sig'
    if si == ei == 0:
        ta, tb = theta_punkt(d), theta_punkt(e)
        if ta is None or tb is None or ta == 0 or ta != -tb:
            return None
        if ist_negiert(d.get('theta_full'), e.get('theta_full')):
            return 'theta'
    return None
