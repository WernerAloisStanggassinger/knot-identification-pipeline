#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Knot-Identification Pipeline — Werner Alois Stanggassinger & Claude Code (Anthropic).
"""Bring the Alexander polynomial into a canonical notation — Disc #809.

WHY: Δ(t) is only defined up to ±t^k. `poly_s08` and `kf332_polys08_sh+0.5_tw-1`
carry the same Δ, once as "-x^18 + 13x^17 - ...", once as "x^18 - 13x^17 + ...".
ident() (label_distinct_knots.py) took the STRING into the identity key and
thereby separated ONE knot into TWO — in Jones, HOMFLY, volume, genus, theta and the
isometry of the complements the two cannot be told apart.

WHAT IS NORMALISED HERE — AND WHAT IS NOT:
  ✅ the OVERALL SIGN (leading coefficient positive)
  ✅ the term ORDER (sorted by degree)
  ⛔ NOT the t^k shift. Reason: in the corpus there is not a single such case
     (alex_vorzeichen_sweep.py, 504 representatives, 2 buckets, 1 hit — and that one is
     purely a sign). A more far-reaching normalisation could merge knots
     that are NOT checked. Doctrine: normalise only what is measured.

CALIBRATED against all 843 records (alex_normalform_gegenfall.py): no record becomes unparsable,
and the normal form merges EXACTLY the one documented pair, no further one.
"""
import re

_TERM = re.compile(r'[+-]?\s*[^+-]+')


def terms(s):
    """'-x^18 + 13 x^17 - 1' -> {monomial: coeff}, or None if not parsable."""
    if s is None or not isinstance(s, str) or not s.strip():
        return None
    out = {}
    for t in _TERM.findall(s):
        t = t.strip()
        if not t:
            continue
        vz = -1 if t.startswith('-') else 1
        t = t.lstrip('+- ').strip()
        if not t:
            return None
        m = re.match(r'^(\d+)?\s*(.*)$', t)
        if not m:
            return None
        k = int(m.group(1)) if m.group(1) else 1
        mon = m.group(2).strip()
        # A monomial is '' (constant), 'x' or 'x^<int>'. Everything else is NOT a
        # polynomial — e.g. the string "None" that str(alex) yields for a missing field.
        # The first draft accepted "None" as a monomial and made '0:1' out of it,
        # i.e. a key value out of a MISSING Alexander polynomial. Caught
        # by the Burton-absence assert in the labeler (19.09.2026).
        if mon and not re.match(r'^x(\^-?\d+)?$', mon):
            return None
        out[mon] = out.get(mon, 0) + vz * k
    return {k: v for k, v in out.items() if v != 0}


def _grad(mon):
    """'x^17' -> 17 · 'x' -> 1 · '' -> 0 · 'x^-3' -> -3."""
    if not mon:
        return 0
    m = re.match(r'^x\^(-?\d+)$', mon)
    if m:
        return int(m.group(1))
    return 1 if mon == 'x' else 0


def canonical(s):
    """Canonical notation: leading coefficient positive, sorted by degree.

    Returns the ORIGINAL STRING if not parsable — then the key acts as before
    (no silent degradation, Disc #809).
    """
    t = terms(s)
    if not t:
        # Not parsable (also: alex missing entirely) -> exactly the old behaviour str(alex).
        # MUST be a string: the key is frozen via json.dumps, and
        # None vs "None" are DIFFERENT keys there.
        return str(s)
    highest = max(t, key=_grad)
    if t[highest] < 0:
        t = {k: -v for k, v in t.items()}
    return '|'.join(f"{_grad(k)}:{v}" for k, v in sorted(t.items(), key=lambda kv: -_grad(kv[0])))
