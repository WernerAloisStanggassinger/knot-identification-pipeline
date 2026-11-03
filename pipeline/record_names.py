#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""record_names.py — lets the corpus layer read records written by the released (English) pipeline.

The corpus layer (build_corpus_db.py, validate_complete_json.py, frankenstein_gate.py, label_distinct_knots.py, …) reads the
authors' records, whose file and field names are German (`<knot>_pflicht_invarianten_complete.json`, `pflichtliste`, …). The
pipeline in this release writes the English names (`<knot>_mandatory_invariants_complete.json`, `mandatory_list`, …), as listed in
tools/translation_table.json. This module maps an English record back to the German names, using that table in reverse:

- load_record(path): the record as a dict. A record with the top-level key `mandatory_list` and without `pflichtliste` is English:
  every dict key that the table lists as a translated key is replaced by its German original (exact string equality), except three
  names that also occur untranslated in the authors' records (_AMBIGUOUS). String values are left as the pipeline wrote them: the
  corpus layer reads values but does not branch on translated texts (no German text literal of the table occurs in its code
  except key names), and mapping values back would turn ordinary words into German ones (a record's name "trefoil" became
  "Kleeblatt" in a first version, because the table also translates that word). Any other record is returned exactly as json.load gives it, so the authors' records are read as before.
- RECORD_GLOBS: the two file-name patterns of a record; record_stem(name): the knot stem of a record file name, for both forms.
"""
import json
import os
import re

_TABLE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'tools', 'translation_table.json')
RECORD_GLOBS = ('*_pflicht_invarianten*_complete.json', '*_mandatory_invariants*_complete.json')
_STEM_RE = re.compile(r'_(?:pflicht_invarianten|mandatory_invariants).*_complete\.json$')
_INV = None
# English names that ALSO occur untranslated in the authors' records (measured 3 November 2026 on the 843 records: key `proof` in
# 777 of them, `route` in 103; `exact` occurs only as a value, which is never mapped); the corpus layer reads `proof`. They are left as they are.
_AMBIGUOUS = {'proof', 'route', 'exact'}


def _inverse():
    global _INV
    if _INV is None:
        t = json.load(open(_TABLE, encoding='utf-8'))
        inv = {}
        for section in ('keys',):
            for de, en in t.get(section, {}).items():
                if isinstance(de, str) and isinstance(en, str) and de != en and en not in _AMBIGUOUS:
                    if en in inv and inv[en] != de:
                        raise ValueError(f'translation table is not invertible at {en!r}')
                    inv[en] = de
        if len(inv) < 50:
            # fail closed: an empty or truncated table would silently leave English records unmapped
            raise RuntimeError(f'{_TABLE}: only {len(inv)} key translations, expected about 118; restore it from git')
        _INV = inv
    return _INV


def to_german(obj):
    """Dict keys back to their German originals; string values are left as written."""
    inv = _inverse()
    if isinstance(obj, dict):
        return {inv.get(k, k) if isinstance(k, str) else k: to_german(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [to_german(v) for v in obj]
    return obj


def is_english(d):
    return isinstance(d, dict) and 'mandatory_list' in d and 'pflichtliste' not in d


def load_record(path):
    with open(path, encoding='utf-8') as f:
        d = json.load(f)
    return to_german(d) if is_english(d) else d


def record_stem(name):
    return _STEM_RE.sub('', name)


def iter_records(root):
    """All record files under root (pathlib.Path), both name forms, each path once: first the German ones in rglob order (so a
    tree of the authors' records is visited exactly as before), then the English ones not already seen."""
    seen = set()
    for pat in RECORD_GLOBS:
        for p in root.rglob(pat):
            if p not in seen:
                seen.add(p)
                yield p
