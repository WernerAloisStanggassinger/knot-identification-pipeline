#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Migration of the frozen identity registry to the alex normal form (Disc #809).

WHY BY HAND AND NOT AUTOMATICALLY: label_distinct_knots.py:415 says explicitly
"migrate registry by hand, do not guess". The frozen key is the ident() tuple as
JSON list; if ONE element of the tuple changes, no entry matches any more, and the labeler
would mint all 504 identities ANEW — including new #-suffixes in places where
nothing is merged. Exactly that happened in the dry run (tb_1237 got '13a_1237#2').

WHAT THIS MIGRATION DOES: in every key it replaces EXACTLY the element with index 2
(the alex string) by its normal form. Nothing else. If two keys collide in doing so,
the ID WITHOUT '#' suffix wins (the one minted first).

CHECKS (aborts if one fails):
  1. exactly ONE collision  — more would mean the normal form merges unproven things
  2. result has 503 entries (504 - 1)
  3. every old ID except 'S_24_d25721_D8_amphi_sig0#2' is kept
  4. no key still has the old alex string after the migration
"""
import site_config
import json
import shutil
import sys
from pathlib import Path
from datetime import date

sys.path.insert(0, (site_config.TOOLS_DIR))
from alexander_normal_form import canonical

REG = Path((site_config.PIPELINE_DIR + '/distinct_id_registry.json'))
VERSCHWINDET = 'S_24_d25721_D8_amphi_sig0#2'
BLEIBT = 'S_24_d25721_D8_amphi_sig0'


def main(apply=False):
    d = json.loads(REG.read_text())
    frozen = d['frozen']
    print(f"old registry: {len(frozen)} entries")

    neu = {}
    kollisionen = []
    for k, e in frozen.items():
        t = json.loads(k)
        assert isinstance(t, list) and len(t) >= 5, f"unerwarteter Key: {k[:80]}"
        t[2] = canonical(t[2])
        assert isinstance(t[2], str), (
            f"Normalform lieferte {t[2]!r} statt eines Strings — json.dumps macht daraus "
            f"null, and null != \"None\" as key. Exactly this killed the first migration "
            f"(Burton-Absenz-Assert, 19.09.2026).")
        nk = json.dumps(t, ensure_ascii=False)
        if nk in neu:
            kollisionen.append((neu[nk]['id'], e['id'], nk))
            # ID without '#' wins
            a, b = neu[nk], e
            gewinner = a if '#' not in a['id'] else (b if '#' not in b['id'] else
                                                     min([a, b], key=lambda x: x['id']))
            verlierer = b if gewinner is a else a
            gewinner = dict(gewinner)
            gewinner['merged_from'] = verlierer['id']
            gewinner['merged_on'] = date.today().isoformat()
            gewinner['merged_why'] = (
                'Disc #809: alex-Gesamtvorzeichen trennte denselben Knoten. Jones und '
                'HOMFLY bitgleich, theta = 0 beidseitig, Komplemente isometrisch.')
            neu[nk] = gewinner
        else:
            neu[nk] = e

    print(f"new registry:  {len(neu)} entries")
    print(f"collisions:     {len(kollisionen)}")
    for a, b, _ in kollisionen:
        print(f"   {a}  +  {b}")

    # ---- CHECKS
    fehler = 0
    if len(kollisionen) != 1:
        print(f"🔴 (1) expected exactly 1 collision, measured {len(kollisionen)}"); fehler += 1
    else:
        print("✅ (1) exactly one collision")
    if len(neu) != len(frozen) - 1:
        print(f"🔴 (2) expected {len(frozen)-1} entries, measured {len(neu)}"); fehler += 1
    else:
        print(f"✅ (2) {len(neu)} entries (= {len(frozen)} - 1)")
    alt_ids = {e['id'] for e in frozen.values()}
    neu_ids = {e['id'] for e in neu.values()}
    weg = alt_ids - neu_ids
    if weg != {VERSCHWINDET}:
        print(f"🔴 (3) vanished IDs: {weg} — expected exactly {{{VERSCHWINDET}}}"); fehler += 1
    else:
        print(f"✅ (3) exactly one ID vanishes: {VERSCHWINDET}")
    dazu = neu_ids - alt_ids
    if dazu:
        print(f"🔴 (3b) NEW IDs appeared: {dazu}"); fehler += 1
    else:
        print("✅ (3b) no new ID appeared")
    rest = [k for k in neu if 'x^18' in k or 'x^17' in k]
    if rest:
        print(f"🔴 (4) {len(rest)} keys still carry the old alex string"); fehler += 1
    else:
        print("✅ (4) no key carries the old alex string any more")

    if fehler:
        print(f"\n🔴 {fehler} CHECK(S) FAILED — nothing written")
        return 1
    print("\n✅ ALL CHECKS PASSED")

    if apply:
        bak = REG.with_suffix('.json.bak-pre-alexnorm-20260919')
        shutil.copy(REG, bak)
        d['frozen'] = neu
        d.setdefault('_doc', [])
        if isinstance(d['_doc'], list):
            d['_doc'].append(
                f"{date.today().isoformat()} — MIGRATED to the alex normal form (Disc #809): "
                f"the key at index 2 no longer carries the Alexander STRING but its "
                f"Normalform (fuehrender Koeffizient positiv, nach Grad sortiert). Grund: "
                f"Δ is only defined up to ±t^k, and two spellings of the same Δ "
                f"separated {BLEIBT} from {VERSCHWINDET} — with bit-identical Jones and HOMFLY, "
                f"theta = 0 on both sides and ISOMETRIC complements. "
                f"{VERSCHWINDET} ist in {BLEIBT} aufgegangen. Backup: {bak.name}")
        REG.write_text(json.dumps(d, ensure_ascii=False, indent=1))
        print(f"== WRITTEN: {REG}  (backup {bak.name}) ==")
    else:
        print("== DRY-RUN (no --apply) ==")
    return 0


if __name__ == '__main__':
    raise SystemExit(main('--apply' in sys.argv))
