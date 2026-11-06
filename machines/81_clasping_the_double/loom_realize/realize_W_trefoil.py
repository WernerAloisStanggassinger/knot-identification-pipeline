#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Realize W(trefoil) (untwisted Whitehead double, bite 6 verified Δ=1∧Jones≠1)
as .fseries via Loom-7. Braid word from spherogram (W7 construction), 7 strands, len 30.
Gates: station_crossings==30, embedding, Nyquist. Output: W_trefoil_untw.fseries."""
import sys, json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))                       # generate_loom7
sys.path.insert(0, str(HERE.parent.parent))         # maschinen_ideen (machine_common)

from generate_loom7 import (build_station, connector, station_crossings,
                            word_perm, N_STR)
from machine_common import densify, embedding_report, polyline_to_fseries

# 7-strand braid word of W(trefoil) untwisted (spherogram braid_word, bite 6d)
WORD_RAW = [-1, -2, -3, 2, 4, 5, 4, 6, -3, -5, -2, 4, -3, 2, 1, -3, -2, 3,
            -4, 3, 5, -4, -6, -3, 2, -3, -2, -4, -5, 4]


def braid_word():
    return [(abs(x), +1 if x > 0 else -1) for x in WORD_RAW]


def build_curve():
    word = braid_word()
    paths, end_lane = build_station(word)
    pts, cur, visited = [], 0, []
    for _ in range(N_STR):
        visited.append(cur)
        pts += paths[cur]
        lane = end_lane[cur]
        pts += connector(lane)
        cur = lane
    if cur != 0 or sorted(visited) != list(range(N_STR)):
        raise RuntimeError("permutation not a %d-cycle — not a knot (perm=%s)"
                           % (N_STR, word_perm(word)))
    return np.asarray(pts, float)


if __name__ == '__main__':
    print('W(trefoil) Loom-7 realization')
    print('  word length (target station_crossings):', len(WORD_RAW))
    print('  Permutation:', word_perm(braid_word()),
          '(7-cycle?)', sorted(set(word_perm(braid_word()))) == list(range(7)))
    P = build_curve()
    n_c = station_crossings(P)
    print('  station_crossings:', n_c, '/ target', len(WORD_RAW),
          '->', 'OK' if n_c == len(WORD_RAW) else 'MISMATCH')
    Pd = densify(P)
    emb = embedding_report(Pd)
    print('  embedding:', {k: emb[k] for k in emb if k in ('embedded', 'min_clearance', 'self_min')})
    if not emb.get('embedded'):
        print('  -> EMBEDDING GATE RED (collision). Loom-7 geometry needs adjustment.')
        sys.exit(0)
    folder = HERE / 'W_trefoil_untw'
    folder.mkdir(parents=True, exist_ok=True)
    rep = None
    for nt in (256, 384, 512, 768, 1024):
        try:
            rep = polyline_to_fseries(P, folder / 'W_trefoil_untw.fseries',
                                      knot_name='W_trefoil_untw',
                                      source='maschinen_ideen/81 Loom-7 (Biss6)',
                                      n_terms=nt, nyquist_max_rel=0.02)
            print('  fseries OK at n_terms=%d, dev=%.2e' % (nt, rep.get('nyquist_dev', -1)))
            break
        except RuntimeError as e:
            print('  n_terms=%d Nyquist red' % nt)
    if rep is None:
        print('  -> NYQUIST RED up to 1024')
        sys.exit(0)
    print('  -> .fseries written:', folder / 'W_trefoil_untw.fseries')
