#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""four_genus_verdict.py — the ONE verdict on the smooth 4-genus, pipeline + catch-up sweep.

WHY THIS EXISTS (2026-07-23, right after the slice finding Disc #555):
pflicht_invarianten.py set `smooth_4_genus = abs(tau)`. That is the same overclaim class
as with the slice field: only g_4 >= |tau| holds (Ozsvath-Szabo), NOT equality. A knot with
tau = 0 thus got g_4 = 0 -- i.e. implicitly "slice" -- through the back door.

THE RULE HERE: collect bounds, state a value only when they coincide. This is exactly the
pattern the pipeline already uses for unknotting_number and bridge_number (interval +
lower/upper + note), so nothing foreign.
  LOWER BOUNDS       |tau|            Ozsvath-Szabo
                     |s| / 2          Rasmussen
                     |signature| / 2  Murasugi
                     1                if slice is proven refuted (g_4 = 0 <=> slice)
  UPPER BOUND        g_3              every Seifert surface can be pushed into the 4-ball
  value = lower == upper, otherwise None + interval.

For torus knots and everything positive the interval collapses (Milnor: g_4 = g_3), so a number
still stands there. Where it does not collapse, there is now an interval instead of a claim.

TOPOLOGICAL 4-GENUS: g_4^top <= g_4^smooth (every smooth surface is also topologically flat) and
g_4^top >= |signature|/2 (Murasugi holds topologically). tau and s do NOT bound it. The old
code set topological_4_genus = smooth_4_genus, i.e. equality -- that also becomes an interval.
"""
from fractions import Fraction


def _zahl(x):
    try:
        return Fraction(x)
    except (TypeError, ValueError):
        return None


def four_genus_verdict(a, slice_verdict_value=None):
    """a = ALGEBRA dict. Returns a dict with bounds + value (or None)."""
    hfk = a.get('knot_floer_homology') or {}
    tau = _zahl(a.get('tau_invariant', hfk.get('tau')))
    s = _zahl(a.get('rasmussen_s', a.get('s_invariant')))
    sig = _zahl(a.get('signature'))
    g3 = a.get('genus', hfk.get('seifert_genus'))

    unten = Fraction(0)
    quellen = []
    if tau is not None and abs(tau) > unten:
        unten = abs(tau); quellen = ['|tau|=%s (Ozsvath-Szabo)' % abs(tau)]
    elif tau is not None:
        quellen.append('|tau|=%s' % abs(tau))
    for wert, name in ((s, '|s|/2 (Rasmussen)'), (sig, '|signature|/2 (Murasugi)')):
        if wert is None:
            continue
        b = abs(wert) / 2
        if b > unten:
            unten = b; quellen = ['%s = %s' % (name, b)]
    if slice_verdict_value is False and unten < 1:
        unten = Fraction(1); quellen = ['slice widerlegt => g_4 >= 1']

    # g_4 is an integer -> round up
    unten_ganz = -((-unten.numerator) // unten.denominator)

    # g_3 only works as an UPPER bound if it is a number. Some records hold
    # only a bound there themselves ("≥ 18") -- nothing follows from that for g_4 from above.
    try:
        oben = int(g3) if g3 is not None else None
    except (TypeError, ValueError):
        oben = None

    out = {
        'smooth_4_genus_lower_bound': int(unten_ganz),
        'smooth_4_genus_lower_source': '; '.join(quellen) if quellen else 'keine Schranke berechenbar',
        'smooth_4_genus_upper_bound': oben,
        'smooth_4_genus_upper_source': ('g_3 = %s (Seifert-Flaeche in die 4-Kugel gedrueckt)' % oben)
                                       if oben is not None else 'g_3 unbekannt',
    }
    if oben is not None and int(unten_ganz) == oben:
        out['smooth_4_genus'] = oben
        out['smooth_4_genus_note'] = 'exakt: untere und obere Schranke fallen zusammen (%s)' % oben
        out['smooth_4_genus_interval'] = [oben, oben]
    else:
        out['smooth_4_genus'] = None
        out['smooth_4_genus_interval'] = [int(unten_ganz), oben]
        out['smooth_4_genus_note'] = (
            'unbestimmt: %s <= g_4 <= %s. Die Pipeline kennt nur Schranken; frueher stand hier |tau|, '
            'das ist aber nur die UNTERE Schranke (2026-07-23).' % (int(unten_ganz), oben))
    out['smooth_4_genus_method'] = 'four_genus_verdict.py (2026-07-23): Schranken, Wert nur bei Gleichheit'

    # topological: |signature|/2 <= g_4^top <= g_4^smooth ; tau and s do NOT bound it
    t_unten = int(-((-(abs(sig) / 2).numerator) // (abs(sig) / 2).denominator)) if sig is not None else 0
    t_oben = out['smooth_4_genus'] if out['smooth_4_genus'] is not None else oben
    out['topological_4_genus_lower_bound'] = t_unten
    out['topological_4_genus_upper_bound'] = t_oben
    if t_oben is not None and t_unten == t_oben:
        out['topological_4_genus'] = t_oben
        out['topological_4_genus_note'] = 'exakt: |signature|/2 == obere Schranke'
    else:
        out['topological_4_genus'] = None
        out['topological_4_genus_note'] = (
            'unbestimmt: %s <= g_4^top <= %s (g_4^top <= g_4^smooth; tau und s schranken das '
            'topologische Genus NICHT ein)' % (t_unten, t_oben))
    return out


def apply_to_algebra(a, slice_verdict_value=None, name=None, identitaet=None):
    """With `name`/`identitaet`, a RECEIPTED entry from beweis_register.json wins first
    (key order: see slice_verdict.beweis_eintrag)."""
    from slice_verdict import beweis_eintrag
    e = beweis_eintrag(name, identitaet) if (name or identitaet) else None
    a.update(four_genus_verdict(a, slice_verdict_value))
    if e is not None and 'smooth_4_genus' in e:
        a['smooth_4_genus'] = e['smooth_4_genus']
        a['smooth_4_genus_interval'] = [e['smooth_4_genus'], e['smooth_4_genus']]
        a['smooth_4_genus_note'] = 'BELEGT (beweis_register.json): ' + ' '.join(e['beweis'])
        a['smooth_4_genus_method'] = 'beweis_register.json, belegt am %s' % e.get('belegt_am', '?')
        if 'topological_4_genus' in e:
            a['topological_4_genus'] = e['topological_4_genus']
            a['topological_4_genus_note'] = a['smooth_4_genus_note']
    return a.get('smooth_4_genus')
