# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 96 — swapping_pretzel_strands (classical mutant pairs).

4-band pretzel P(p,q,r,s): tangle SWAP = mutation — and with
four bands (unlike three, where all permutations are isotopic)
it generally really changes the knot: P(2,3,5,7) vs P(2,5,3,7) =
citable mutant pair. 74er mechanics generalized to n bands
(clasp-routing lessons inherited).
MAIL THREAD: Dror — classically named pairs supplement the 90er
random pairs. EXACT BETS today: component trace = 1 for both
(one even pᵢ); crossing target ≥ Σ|pᵢ|; pair embedded.
IGNITION: polynomial collision (mandatory) + type separation + θ.
REVERSAL CLAUSE: swap may be isotopic → ignition judges.
Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import (embedding_report, polyline_to_fseries,  # noqa: E402
                            count_projection_crossings)

PTS_PER_TWIST = 120
BAND_HW = 0.28
Y_TOP, Y_BOT = 1.5, -1.5


def band_strands(x0, p):
    n = max(60, PTS_PER_TWIST * max(1, abs(p)) // 2)
    t = np.linspace(0, 1, n)
    ang = np.pi * p * t
    y = Y_TOP + (Y_BOT - Y_TOP) * t
    sA = np.column_stack([x0 + BAND_HW * np.cos(ang), y,
                          BAND_HW * np.sin(ang)])
    sB = np.column_stack([x0 - BAND_HW * np.cos(ang), y,
                          -BAND_HW * np.sin(ang)])
    return sA, sB


def arc(p0, p1, lift, z_sign=1.0, n=40):
    t = np.linspace(0, 1, n)[1:-1]
    base = np.outer(1 - t, p0) + np.outer(t, p1)
    bulge = np.sin(np.pi * t)
    base[:, 1] += lift * bulge
    base[:, 2] += z_sign * 0.8 * abs(lift) * bulge
    return base


def pretzel_n(ps):
    """n-band pretzel: strands + gluing graph + component trace."""
    n = len(ps)
    xs = np.linspace(-(n - 1), (n - 1), n)
    strands = {}
    for i, (x0, p) in enumerate(zip(xs, ps)):
        sA, sB = band_strands(x0, p)
        strands[(i, "A")] = sA
        strands[(i, "B")] = sB
    def bot_side(i, s):
        return s if ps[i] % 2 == 0 else ("B" if s == "A" else "A")
    top_pairs = [((i, "B"), (i + 1, "A")) for i in range(n - 1)] \
        + [((0, "A"), (n - 1, "B"))]
    bot_pairs = list(top_pairs)
    def strand_with_botside(i, side):
        for s in ("A", "B"):
            if bot_side(i, s) == side:
                return (i, s)
    adj = {}
    for (u, v) in top_pairs:
        adj.setdefault((u, "top"), []).append((v, "top"))
        adj.setdefault((v, "top"), []).append((u, "top"))
    for (i1, sd1), (i2, sd2) in bot_pairs:
        u = strand_with_botside(i1, sd1)
        v = strand_with_botside(i2, sd2)
        adj.setdefault((u, "bot"), []).append((v, "bot"))
        adj.setdefault((v, "bot"), []).append((u, "bot"))
    seen, n_comp, order = set(), 0, []
    for start in [(i, s) for i in range(n) for s in ("A", "B")]:
        if start in seen:
            continue
        n_comp += 1
        cyc, cur, lvl = [], start, "top"
        while cur not in seen:
            seen.add(cur)
            cyc.append((cur, lvl))
            other = "bot" if lvl == "top" else "top"
            cur, lvl = adj[(cur, other)][0]
        order.append(cyc)
    return strands, order, n_comp, xs


def assemble(ps):
    strands, order, n_comp, xs = pretzel_n(ps)
    if n_comp != 1:
        return None, n_comp
    pieces = []
    cyc = order[0]
    n_x = max(abs(xs[0]), 1)
    for idx, ((i, s), lvl) in enumerate(cyc):
        S = strands[(i, s)]
        seg = S if lvl == "top" else S[::-1]
        pieces.append(seg)
        nxt_strand, nxt_lvl = cyc[(idx + 1) % len(cyc)]
        SN = strands[nxt_strand]
        target = SN[0] if nxt_lvl == "top" else SN[-1]
        sign = 1.0 if abs(seg[-1][1] - Y_TOP) < 0.5 else -1.0
        span = abs(float(target[0] - seg[-1][0]))
        mid_x = 0.5 * float(target[0] + seg[-1][0])
        if span > 1.6 * n_x:
            z_sign = -1.0                       # outer clasp to the back
        else:
            # Neighbor clasp: amplitude staggering over the x position
            z_sign = 0.45 + 0.55 * (mid_x - xs[0]) / (xs[-1] - xs[0])
        pieces.append(arc(seg[-1], target, sign * (0.45 + 0.22 * span),
                          z_sign=z_sign))
    return np.vstack(pieces), n_comp


def make_pair(ps_a, ps_b, out_root):
    name = "pz4_" + "_".join(map(str, ps_a)) + "__vs__" \
        + "_".join(map(str, ps_b))
    entry = {"name": name, "machine": "96_swapping_pretzel_strands",
             "pair": [list(ps_a), list(ps_b)],
             "predictions": {
                 "zuendung": "Polynome kollidieren (Mutation), Typ/θ "
                             "entscheiden", "schiene": "Dror"}}
    folder = Path(out_root) / name
    for tag, ps in (("K", ps_a), ("Kmut", ps_b)):
        C, n_comp = assemble(list(ps))
        if C is None:
            entry["rejected"] = f"{tag}: {n_comp} Komponenten"
            return entry, None
        nx = count_projection_crossings(C[::4], axis=2)
        entry[f"crossings_{tag}"] = int(nx)
        emb = embedding_report(C)
        if not emb["embedded"]:
            entry["rejected"] = f"embedding gate ({tag}, " \
                f"{emb['d_min_rel']:.5f})"
            return entry, None
        folder.mkdir(parents=True, exist_ok=True)
        rep = None
        for nt in (160, 256, 384):
            try:
                rep = polyline_to_fseries(
                    C, folder / f"{name}_{tag}.fseries",
                    knot_name=f"{name}_{tag}",
                    source="maschinen_ideen/96", n_terms=nt,
                    nyquist_max_rel=0.02)
                break
            except RuntimeError:
                continue
        if rep is None:
            entry["rejected"] = f"Nyquist rot ({tag})"
            return entry, None
        entry[f"fseries_{tag}"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] smoke test machine 96 (pretzel mutants)")
    # (1) component rules: P(2,3,5,7)=1; P(3,5,7)=1 (inherited from 74);
    #     P(3,3,5,7) [all odd, n even] ≠ 1
    c1 = pretzel_n([2, 3, 5, 7])[2]
    c2 = pretzel_n([2, 5, 3, 7])[2]
    c3 = pretzel_n([3, 3, 5, 7])[2]
    ok1 = c1 == 1 and c2 == 1 and c3 != 1
    print(f"  (1) components: P(2,3,5,7)={c1}, P(2,5,3,7)={c2}, "
          f"P(3,3,5,7)={c3} (all-odd 4-tuple = link)  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) THE PAIR: P(2,3,5,7) vs P(2,5,3,7) built
    e2, f2 = make_pair((2, 3, 5, 7), (2, 5, 3, 7), HERE / "demo_smoke")
    ok2 = f2 is not None
    print(f"  (2) {e2['name']}: crossings "
          f"{e2.get('crossings_K')}/{e2.get('crossings_Kmut')} "
          f"(expected twists ≥17), nyq="
          f"{e2.get('fseries_K', {}).get('nyquist_dev_rel')}/"
          f"{e2.get('fseries_Kmut', {}).get('nyquist_dev_rel')}  "
          f"{'✓ classical pair in place' if ok2 else '✗ FAIL: ' + str(e2.get('rejected'))}")
    assert ok2
    # (3) Second pair (smaller): P(2,3,5) has 3 bands → take
    #     P(2,3,5,3) vs P(2,5,3,3)
    e3, f3 = make_pair((2, 3, 5, 3), (2, 5, 3, 3), HERE / "demo_smoke")
    ok3 = f3 is not None
    print(f"  (3) {e3['name'] if e3 else '—'}: "
          f"nyq={e3.get('fseries_K', {}).get('nyquist_dev_rel')}/"
          f"{e3.get('fseries_Kmut', {}).get('nyquist_dev_rel')}  "
          f"{'✓ series' if ok3 else '✗ FAIL: ' + str(e3.get('rejected'))}")
    assert ok3
    print("[DEMO] 3/3 green — pretzel mutant machine ready to fire "
          "(citable pairs; polynomial collision = ignition duty).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
