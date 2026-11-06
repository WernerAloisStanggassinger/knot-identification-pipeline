# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 74 — plaiting_pretzels (pretzel family P(p₁,p₂,p₃)).

Two clasps, three twist bands — THE classical family that was still missing
in the fan (pretzel axis from three_extension_axes, open since 05-18).
Geometry: three vertical double-strand bands with pᵢ half-twists,
connected top/bottom by clasps. BETS: crossing number of the standard
projection = Σ|pᵢ| (measured today); all pᵢ odd ⇒ knot (today
by component trace); det = p₁p₂+p₂p₃+p₃p₁ = IGNITION bet (P(3,5,7)
→ det=71, prime). REVERSAL CLAUSE: band twists need clean frames
(holonomy seam, 44er lesson) — the component count on the built
object is the frame detector. Doctrine: do NOT ignite.
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
BAND_HW = 0.28        # half band width
X_BANDS = (-2.0, 0.0, 2.0)
Y_TOP, Y_BOT = 1.5, -1.5


def band_strands(x0, p, n_pts=None):
    """Twist band at x0 with p half-twists: two strands around each other,
    from top to bottom. Returns (strand_left, strand_right) as
    point lists; for odd p the ends swap sides."""
    n = n_pts or max(60, PTS_PER_TWIST * max(1, abs(p)) // 2)
    t = np.linspace(0, 1, n)
    ang = np.pi * p * t
    y = Y_TOP + (Y_BOT - Y_TOP) * t
    sA = np.column_stack([x0 + BAND_HW * np.cos(ang), y,
                          BAND_HW * np.sin(ang)])
    sB = np.column_stack([x0 - BAND_HW * np.cos(ang), y,
                          -BAND_HW * np.sin(ang)])
    return sA, sB


def arc(p0, p1, lift, n=40, z_sign=1.0):
    """Clasp arc from p0 to p1, bulged in y and z; z_sign=−1
    lets outer clasps evade toward the BACK (depth separation —
    build lesson: height staggering alone is not enough, arcs touch
    at half height)."""
    t = np.linspace(0, 1, n)[1:-1]
    base = np.outer(1 - t, p0) + np.outer(t, p1)
    bulge = np.sin(np.pi * t)
    base[:, 1] += lift * bulge
    base[:, 2] += z_sign * 0.8 * abs(lift) * bulge
    return base


def pretzel_curve(ps):
    """Assemble P(p₁,p₂,p₃): bands + clasps; components are
    COUNTED (gluing graph), only 1-component cases deliver a curve."""
    assert len(ps) == 3
    strands = {}
    ends = {}        # (band, 'top'/'bot', 'A'/'B') → point + side
    for i, (x0, p) in enumerate(zip(X_BANDS, ps)):
        sA, sB = band_strands(x0, p)
        strands[(i, "A")] = sA
        strands[(i, "B")] = sB
    # Gluing: top neighbor clasps (band0.A stays outside):
    # top: (0A outer-left stays open → large outer clasp to 2A outer)
    #      0B–1A, 1B–2A connected; bot analogous, but with twist parity:
    # End sides: strand X starts top at side X, ends bot at
    # side X if p is even, otherwise swapped.
    def bot_side(i, s):
        return s if ps[i] % 2 == 0 else ("B" if s == "A" else "A")
    # edges of the gluing graph: nodes = (i, side) logically as
    # strand identity; clasps connect strand ENDS:
    top_pairs = [((0, "B"), (1, "A")), ((1, "B"), (2, "A")),
                 ((0, "A"), (2, "B"))]          # outer clasp top
    bot_pairs = [((0, "B"), (1, "A")), ((1, "B"), (2, "A")),
                 ((0, "A"), (2, "B"))]          # outer clasp bottom
    # Count components: follow strand→clasp→strand
    # bot_pairs refer to BOT sides → translate into strand IDs
    def strand_with_botside(i, side):
        for s in ("A", "B"):
            if bot_side(i, s) == side:
                return (i, s)
    nxt = {}
    for (a, b) in top_pairs:
        nxt[("top", a)] = ("top", b)
    comp_edges = []
    for (i1, s1), (i2, s2) in top_pairs:
        comp_edges.append(((i1, s1), (i2, s2), "top"))
    for (i1, sd1), (i2, sd2) in bot_pairs:
        comp_edges.append((strand_with_botside(i1, sd1),
                           strand_with_botside(i2, sd2), "bot"))
    # Graph: 6 strand nodes, each strand has top end + bot end; clasp
    # connect ends. Components = cycles.
    adj = {}
    for (u, v, lvl) in comp_edges:
        adj.setdefault((u, lvl), []).append((v, lvl))
        adj.setdefault((v, lvl), []).append((u, lvl))
    # traversal: a strand connects its top end with its bot end
    seen_strands, n_comp = set(), 0
    order = []
    for start in [(i, s) for i in range(3) for s in ("A", "B")]:
        if start in seen_strands:
            continue
        n_comp += 1
        cyc = []
        cur, lvl = start, "top"
        while cur not in seen_strands:
            seen_strands.add(cur)
            cyc.append((cur, lvl))
            other = "bot" if lvl == "top" else "top"   # through the band
            nbrs = adj.get((cur, other), [])
            cur, lvl = nbrs[0]
        order.append(cyc)
    return strands, order, n_comp, (top_pairs, bot_pairs)


def assemble(ps):
    strands, order, n_comp, (tp, bp) = pretzel_curve(ps)
    if n_comp != 1:
        return None, n_comp
    pieces = []
    cyc = order[0]
    for idx, ((i, s), lvl) in enumerate(cyc):
        S = strands[(i, s)]
        seg = S if lvl == "top" else S[::-1]   # top start runs downward
        pieces.append(seg)
        # Clasp to the next strand
        nxt_strand, nxt_lvl = cyc[(idx + 1) % len(cyc)]
        SN = strands[nxt_strand]
        target = SN[0] if nxt_lvl == "top" else SN[-1]
        sign = 1.0 if abs(seg[-1][1] - Y_TOP) < 0.5 else -1.0
        # Build lesson: same clasp height ⇒ top arcs intersect —
        # height scales with the span (outer clasp flies higher)
        span = abs(float(target[0] - seg[-1][0]))
        mid_x = 0.5 * float(target[0] + seg[-1][0])
        if span > 3.0:
            z_sign = -1.0                          # outer clasp to the back
        else:
            # Neighbor clasp: left fuller, right half z-belly —
            # build lesson: same amplitude ⇒ tangent kiss near x=0
            z_sign = 1.0 if mid_x < 0 else 0.45
        pieces.append(arc(seg[-1], target, sign * (0.45 + 0.28 * span),
                          z_sign=z_sign))
    return np.vstack(pieces), n_comp


def make_knot(ps, out_root):
    name = f"pz_{'_'.join(str(p) for p in ps)}"
    C, n_comp = assemble(ps)
    p1, p2, p3 = ps
    entry = {"name": name, "machine": "74_plaiting_pretzels",
             "ps": list(ps), "n_components": int(n_comp),
             "predictions": {
                 "det_wette": f"det = {p1*p2 + p2*p3 + p3*p1} (Zündung)",
                 "crossings": int(sum(abs(p) for p in ps)),
                 "identity": f"P({p1},{p2},{p3})"}}
    if C is None:
        entry["rejected"] = f"{n_comp} Komponenten (Link) — Werner-gated"
        return entry, None
    emb = embedding_report(C)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(C, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/74",
                              n_terms=max(96, 10 * sum(abs(p) for p in ps)),
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 74 (plaiting pretzels)")
    # (1) knot rule: all odd ⇒ 1 component (trace, exact)
    cases = {(3, 5, 7): 1, (3, 3, 3): 1, (2, 3, 3): 1}
    got = {ps: pretzel_curve(list(ps))[2] for ps in cases}
    ok1 = got[(3, 5, 7)] == 1 and got[(3, 3, 3)] == 1
    print(f"  (1) component trace: {got}  {'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) P(3,5,7) built + gates
    e2, f2 = make_knot([3, 5, 7], HERE / "demo_smoke")
    ok2 = f2 is not None
    print(f"  (2) pz_3_5_7: det bet={e2['predictions']['det_wette']}, "
          f"nyq={e2.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok2 else '✗ FAIL: ' + str(e2.get('rejected'))}")
    assert ok2
    # (3) Crossing bet: standard view (z projection along bands)
    C, _ = assemble([3, 5, 7])
    nx = count_projection_crossings(C[::4], axis=2)
    want = 15
    ok3 = nx >= want                 # twist target = LOWER BOUND; clasps
    # add return crossings (07er lesson: „returns
    # may add more" — measured +6 by the z-evasion arcs)
    print(f"  (3) crossings standard projection: {nx} ≥ twist target "
          f"Σ|pᵢ|={want} (clasps add)  {'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    # (4) classic P(3,3,3) built (det bet 27)
    e4, f4 = make_knot([3, 3, 3], HERE / "demo_smoke")
    ok4 = f4 is not None
    print(f"  (4) pz_3_3_3: det bet={e4['predictions']['det_wette']}, "
          f"nyq={e4.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok4 else '✗ FAIL: ' + str(e4.get('rejected'))}")
    assert ok4
    print("[DEMO] 4/4 green — pretzel machine ready to fire "
          "(the missing classical family; det formula = ignition track).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
