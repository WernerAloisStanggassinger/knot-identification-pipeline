#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""slice_verdict.py — the ONE verdict on slice knots, for the pipeline and the catch-up sweep.

WHY THIS EXISTS (2026-07-23, pipeline_error_log [SLICE-FELD]):
pflicht_invarianten.py set `slice = (tau == 0 and seifert_genus is not None)`. tau = 0 is a
NECESSARY, not a sufficient condition -- so EVERY tau=0 knot was declared a slice knot.
213 records affected, 29 of them provably wrong, 78 live on the website. Noticed on
gtwr_1.9375: tau=0, but det=629 is not a perfect square, so by Fox-Milnor definitely not slice.

THE RULE HERE: the pipeline can REFUTE slice, but not PROVE it (exception: the unknot).
  slice = False  if any obstruction applies              -> proven NOT slice
  slice = True   only for the unknot (genus 0)           -> trivially slice
  slice = None   otherwise                               -> undetermined, NOT "probably"
An empty field is more honest than a plausible-looking wrong one.

THE OBSTRUCTIONS (each one on its own sufficient for "not slice"):
  tau != 0            Ozsvath-Szabo: tau vanishes on slice knots
  s != 0              Rasmussen: s vanishes on slice knots
  signature != 0      Murasugi: the signature vanishes on slice knots
  Arf != 0            Robertello: Arf vanishes on slice knots
  det not a square    Fox-Milnor: Delta(t) = f(t)*f(1/t) => |Delta(-1)| = det is a square
We deliberately do NOT infer anything from smooth_4_genus -- that field is itself derived from |tau|
(same overclaim class, see the error log) and would be circular reasoning.
"""
import json
import os
from math import isqrt

_REGISTER_PFAD = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'beweis_register.json')
_REGISTER = None


def beweis_eintrag(name, identitaet=None):
    """Receipt-backed override from beweis_register.json — ONLY with a non-empty 'beweis' (proof).
    An entry without proof text is ignored; the register must not become a back door.

    Lookup order (2026-08-19):
      1. `name`        — the record stem, e.g. 'double-3'    (specific beats general)
      2. `identitaet`  — the distinct_knot_id, e.g. '12a_1019'
      3. prefix of `name` before '__' — applies under the naming doctrine `<classical>__<original>`,
         i.e. '12a_1019__kf221_sh+1.875_tw+1.5' -> '12a_1019'. A freshly generated record thereby
         inherits the identity receipts automatically, even without a corpus join.
    Reason: sliceness is a property of the KNOT, not of the parametrisation. One entry
    per knot keeps the register small and prevents 39 records from drifting apart."""
    global _REGISTER
    if _REGISTER is None:
        try:
            _REGISTER = json.load(open(_REGISTER_PFAD))
        except Exception:
            _REGISTER = {}
    kandidaten = [name, identitaet]
    if name and '__' in name:
        kandidaten.append(name.split('__')[0])
    for k in kandidaten:
        if not k:
            continue
        e = _REGISTER.get(k)
        if isinstance(e, dict) and e.get('beweis'):
            return e
    return None


def _quadratzahl(n):
    try:
        n = int(n)
    except (TypeError, ValueError):
        return None
    if n < 0:
        return False
    return isqrt(n) * isqrt(n) == n


def slice_verdict(a):
    """a = ALGEBRA dict. Returns (verdict, note, obstructions).
    verdict: True (unknot only) / False (proven not slice) / None (undetermined)."""
    hfk = a.get('knot_floer_homology') or {}
    tau = a.get('tau_invariant', hfk.get('tau'))
    s = a.get('rasmussen_s', a.get('s_invariant'))
    sig = a.get('signature')
    arf = a.get('arf_invariant')
    det = a.get('determinant')
    genus = a.get('genus', hfk.get('seifert_genus'))

    # Do NOT take Arf blindly from the field: Arf is DETERMINED by the determinant
    # (Arf = 0 <=> det ≡ ±1 mod 8) -- that is a theorem, not an estimate, and the determinant
    # comes from three independent paths. Where the stored Arf contradicts it (known for
    # st121_k4, Goeritz parity bug Disc #291), the theorem wins and the contradiction is recorded.
    arf_konflikt = None
    if det is not None:
        try:
            arf_aus_det = 0 if int(det) % 8 in (1, 7) else 1
            if arf is not None and arf != arf_aus_det:
                arf_konflikt = 'stored Arf=%s contradicts det=%s mod 8 (=> Arf=%s); ' \
                               'the theorem holds, the field does not' % (arf, det, arf_aus_det)
            arf = arf_aus_det
        except (TypeError, ValueError):
            pass

    obstr = []
    if tau not in (0, None):
        obstr.append('tau=%s != 0 (Ozsvath-Szabo)' % tau)
    if s not in (0, None):
        obstr.append('s=%s != 0 (Rasmussen)' % s)
    if sig not in (0, None):
        obstr.append('signature=%s != 0 (Murasugi)' % sig)
    if arf not in (0, None):
        obstr.append('Arf=%s != 0 (Robertello%s)' % (arf, ', from det mod 8' if det is not None else ''))
    if det is not None and _quadratzahl(det) is False:
        obstr.append('det=%s not a square (Fox-Milnor, point evaluation at t=-1)' % det)
    # FULL Fox-Milnor: Delta = f(t)*f(1/t) as a factorisation condition. Strictly stronger than
    # the determinant test, which only looks at Delta(-1) (2026-07-23, Disc #558).
    alex = a.get('alexander_polynomial') or a.get('alexander')
    if alex:
        try:
            from fox_milnor import fox_milnor as _fm
            besteht, grund, _d = _fm(alex)
            if besteht is False:
                obstr.append('Fox-Milnor verletzt: %s' % grund)
        except Exception:
            pass

    warn = (' | ⚠ ' + arf_konflikt) if arf_konflikt else ''
    if obstr:
        return False, 'proven NOT slice: ' + '; '.join(obstr) + warn, obstr
    if genus == 0:
        return True, 'unknot (genus 0) — trivially slice' + warn, []
    fehlt = [n for n, v in (('tau', tau), ('s', s), ('signature', sig),
                            ('Arf', arf), ('determinant', det)) if v is None]
    note = ('undetermined: all checked necessary conditions hold '
            '(tau=s=signature=Arf=0, det a square), that is NOT sufficient — '
            'the pipeline cannot prove slice')
    if fehlt:
        note += ' | not checked (field missing): ' + ', '.join(fehlt)
    return None, note + warn, []


def apply_to_algebra(a, name=None, identitaet=None):
    """Sets slice + slice_note + slice_obstructions in place. Returns the verdict.
    With `name` the proof register is asked first: if it holds a RECEIPTED entry, that entry wins,
    and the proof text is copied into the field (so that nobody has to re-derive it later)."""
    e = beweis_eintrag(name, identitaet) if (name or identitaet) else None
    if e is not None and 'slice' in e:
        a['slice'] = e['slice']
        a['slice_note'] = 'RECEIPTED (beweis_register.json): ' + ' '.join(e['beweis'])
        a['slice_obstructions'] = []
        a['slice_method'] = 'beweis_register.json, receipted on %s' % e.get('belegt_am', '?')
        if 'ribbon' in e:
            a['ribbon'] = e['ribbon']
        # Publicly citable source (2026-08-19). The internal slice_method only says
        # "beweis_register.json" -- that is not a source citation for the website. Whoever ships a
        # claim ships the receipt with it; that is why it goes into the record.
        if e.get('quelle_oeffentlich'):
            a['slice_source'] = e['quelle_oeffentlich']
        if e.get('quelle_url'):
            a['slice_source_url'] = e['quelle_url']
        if 'is_amphicheiral' in e:
            a['is_amphicheiral_belegt'] = e['is_amphicheiral']
        return e['slice']
    v, note, obstr = slice_verdict(a)
    a['slice'] = v
    a['slice_note'] = note
    a['slice_obstructions'] = obstr
    a['slice_method'] = 'slice_verdict.py (2026-07-23): refutable, not provable'
    return v
