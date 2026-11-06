# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 02 — plugging_tangles (Conway tangle construction kit).

Plug plan = computation: continued fraction [a1,…,an] → p/q → 2-bridge knot b(p,q).
The identity of the result is fixed BEFORE building — det = p exactly,
Schubert class {q, q⁻¹ mod p}, mirror = b(p, p−q). [Conway/Schubert,
textbook-fixed; NO speculative formulas — the pipeline measures Alexander/sig
blindly at firing.]

GEOMETRY CONVENTION (software side; steel module design = Werner):
  4-plat: four strands on x∈{0,1,2,3}, height = z. Continued-fraction entry a_i
  is realised alternately as a twist region: odd position → σ2 pair
  (slots 1,2; „horizontal"), even position → σ1 pair (slots 0,1;
  „vertical", sign flipped). |a_i| half-twists as a helix (radius 0.5,
  Δz=1 per half-twist), sign = sense of rotation. Plat closure: semicircle
  caps (0,1)+(2,3) at the bottom and top. Projection along y shows exactly
  Σ|a_i| crossings (geometry sanity gate in the demo).
  Identity bet holds UP TO MIRROR (sign conventions about
  twist handedness × Disc-#257 world are settled by the blind test at the first firing,
  after that the convention is calibrated and stands here in the header).

Bets (bit-fixed): det = p · components = 1 if p is odd, 2 if
even (even p = link → rejected for the knot corpus) · b(p,q) class.
"""
from __future__ import annotations
import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import (embedding_report, polyline_to_fseries,  # noqa: E402
                            count_projection_crossings)

# ODD point count per half-twist: the projection crossing lies at
# θ=π/2+nπ — with an even division it falls EXACTLY on a common
# vertex of both strands and the segment-intersection test (open intervals)
# rejects them (demo finding: 0/4 crossings). 17 places them between vertices.
PTS_PER_HALFTWIST = 17


# ----------------------------------------------------------- Conway computation

def cf_value(cf) -> Fraction:
    """[a1,…,an] → value = a_n + 1/(a_{n-1} + … + 1/a_1) (signed)."""
    x = Fraction(cf[0])
    for a in cf[1:]:
        if x == 0:
            raise ValueError(f"Continued fraction {cf} runs through 0 — plug plan invalid")
        x = Fraction(a) + 1 / x
    return x


def cf_to_pq(cf):
    x = cf_value(cf)
    return abs(x.numerator), abs(x.denominator)


def schubert_class(p, q):
    """Unoriented Schubert class: {q, q⁻¹ mod p}; mirror separate."""
    q = q % p
    qinv = pow(q, -1, p)
    cls = sorted({q, qinv})
    mir = sorted({(p - q) % p, pow((p - q) % p, -1, p)})
    return cls, mir


def predictions(cf):
    p, q = cf_to_pq(cf)
    cls, mir = schubert_class(p, q) if p > 1 else ([0], [0])
    return {
        "continued_fraction": list(cf),
        "p_over_q": f"{p}/{q}",
        "det": p,
        "n_components": 1 if p % 2 else 2,
        "is_knot": bool(p % 2),
        "two_bridge_class_q": cls,
        "mirror_class_q": mir,
        "identity_bet": f"b({p},{cls[0]}) bis auf Spiegel" if p % 2 else
                        f"2-Komponenten-Link b({p},{cls[0]}) — kein Knoten-Record",
        "source": "Conway/Schubert (Lehrbuch); Alexander/sig misst Pipeline blind",
    }


# ------------------------------------------------------------- 4-Plat-Bauer

def odd_cf_of(p, q):
    """Canonical ODD-length positive continued-fraction form of p/q.

    4-plat parity: only odd-length words close to the 2-bridge KNOT
    (even length → identity-permutation trap, demo finding [2,2]).
    Textbook identity [...,a] = [...,a−1,1] flips the length parity.
    The list is delivered in build order (innermost entry first,
    compatible with the right fold of cf_to_pq).
    """
    # Standard Euclid: p/q = c1 + 1/(c2 + 1/…), all c_i ≥ 1
    cs, a, b_ = [], p, q
    while b_:
        cs.append(a // b_)
        a, b_ = b_, a % b_
    if len(cs) % 2 == 0:
        if cs[-1] >= 2:
            cs = cs[:-1] + [cs[-1] - 1, 1]
        else:
            cs = cs[:-2] + [cs[-2] + 1]
    build = list(reversed(cs))      # right fold: outermost term last
    chk = cf_to_pq(build)
    if chk[0] != p:
        raise RuntimeError(f"odd_cf_of({p},{q}): parity fix lost p — {chk}")
    return build


def cf_to_braid_word(cf):
    """Entry i → (slot pair, half-twists): alternating σ2 / σ1⁻."""
    word = []
    for i, a in enumerate(cf):
        if i % 2 == 0:
            word.append((1, a))      # σ2: Slots 1,2
        else:
            word.append((0, -a))     # σ1 with convention flip: slots 0,1
    return word


def build_plat(word):
    """4-plat polyline. Returns (P closed, n_components, z_total)."""
    slots_x = [0.0, 1.0, 2.0, 3.0]
    wires = [[np.array([slots_x[s], 0.0, 0.0])] for s in range(4)]
    slot_of = [0, 1, 2, 3]      # wire index → current slot
    wire_at = [0, 1, 2, 3]      # Slot → wire-Index
    z = 0.0
    for (i, k) in word:
        if k == 0:
            continue
        dz = abs(k) * 1.0
        wa, wb = wire_at[i], wire_at[i + 1]
        xc = (slots_x[i] + slots_x[i + 1]) / 2.0
        r = 0.5
        sgn = 1.0 if k > 0 else -1.0
        th = np.linspace(0.0, abs(k) * np.pi, abs(k) * PTS_PER_HALFTWIST + 1)[1:]
        zz = z + dz * th / th[-1]
        wires[wa] += list(np.column_stack([xc - r * np.cos(th),
                                           -sgn * r * np.sin(th), zz]))
        wires[wb] += list(np.column_stack([xc + r * np.cos(th),
                                           sgn * r * np.sin(th), zz]))
        for w in range(4):
            if w not in (wa, wb):
                wires[w].append(np.array([slots_x[slot_of[w]], 0.0, z + dz]))
        if abs(k) % 2 == 1:     # odd half-twists = slot swap
            slot_of[wa], slot_of[wb] = i + 1, i
            wire_at[i], wire_at[i + 1] = wb, wa
        z += dz

    # Caps: semicircles in the xz plane, bottom −z, top +z
    def cap(xa, xb, z0, up):
        xc, r = (xa + xb) / 2.0, abs(xb - xa) / 2.0
        th = np.linspace(0.0, np.pi, 17)[1:-1]
        s = 1.0 if up else -1.0
        return np.column_stack([xc - r * np.cos(th) * np.sign(xb - xa),
                                np.zeros_like(th), z0 + s * r * np.sin(th)])

    # Plat closure: at the bottom the STARTS of slots (0,1),(2,3) pair up,
    # at the top the ENDS. Wires start on slot = own index.
    pair = {0: 1, 1: 0, 2: 3, 3: 2}

    def trace_from(w0):
        """Follow closed walk from (w0, upward). Returns step list."""
        steps, w, direction = [], w0, +1
        while True:
            steps.append((w, direction))
            if direction == +1:
                s2 = pair[slot_of[w]]          # partner end slot at the top
                w, direction = wire_at[s2], -1
            else:
                s2 = pair[w]                   # partner start slot at the bottom
                w, direction = s2, +1          # wire starts on slot = index
            if (w, direction) == (w0, +1):
                return steps

    # Count components: every walk visits each of its wires exactly once
    # (upward OR downward) — mark both directions
    seen, n_comp = set(), 0
    for w0 in range(4):
        if w0 in seen:
            continue
        n_comp += 1
        seen |= {w for (w, d) in trace_from(w0)}
    if n_comp != 1:
        return None, n_comp, z

    P = []
    for (w, direction) in trace_from(0):
        pts = np.asarray(wires[w], float)
        if direction == +1:
            P.append(pts)
            P.append(cap(slots_x[slot_of[w]], slots_x[pair[slot_of[w]]],
                         z, up=True))
        else:
            P.append(pts[::-1])
            P.append(cap(slots_x[w], slots_x[pair[w]], 0.0, up=False))
    return np.vstack(P), 1, z


# ----------------------------------------------------------------- generation

def name_of(cf):
    return "tg_" + "_".join(("m" + str(-a)) if a < 0 else str(a) for a in cf)


def make_knot(cf, out_root):
    pred = predictions(cf)
    entry = {"name": name_of(cf), "machine": "02_plugging_tangles",
             "predictions": pred}
    if not pred["is_knot"]:
        entry["rejected"] = "p gerade = 2-Komponenten-Link"
        return entry, None
    # Build form: canonical odd-length continued fraction (knot-equivalent,
    # Schubert); negative values = mirror build (all twists negated)
    x = cf_value(cf)
    p, q = abs(x.numerator), abs(x.denominator)
    build_cf = odd_cf_of(p, q)
    entry["build_cf"] = build_cf
    word = cf_to_braid_word(build_cf)
    if x < 0:
        word = [(slot, -k) for (slot, k) in word]
        entry["mirror_build"] = True
    P, ncomp, _ = build_plat(word)
    if ncomp != 1:
        raise RuntimeError(f"{cf}: build gave {ncomp} components, prediction 1 "
                           f"(p={pred['det']} odd) — construction bug!")
    nx = count_projection_crossings(P, axis=1)
    want = sum(abs(a) for a in build_cf)
    entry["proj_crossings"] = {"got": nx, "want": want, "ok": nx == want}
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / entry["name"]
    folder.mkdir(parents=True, exist_ok=True)
    n_terms = max(64, 12 * want)
    rep = polyline_to_fseries(P, folder / f"{entry['name']}.fseries",
                              knot_name=entry["name"],
                              source="maschinen_ideen/02_plugging_tangles",
                              n_terms=n_terms)
    entry["fseries"] = rep
    return entry, folder


# ----------------------------------------------------------------------- CLI

def cmd_plan(args):
    """Prediction table of the first plug plans (knots, de-duplicated)."""
    seen, rows = {}, []
    def keyof(p, cls, mir):
        return (p, tuple(min(cls, mir)))
    cfs = []
    for a in range(2, 10):
        cfs.append([a])
    for a in range(1, 5):
        for b in range(1, 5):
            cfs.append([a, b])
            for c in range(1, 4):
                cfs.append([a, b, c])
    for cf in cfs:
        try:
            pred = predictions(cf)
        except ValueError:
            continue
        if not pred["is_knot"] or pred["det"] == 1:
            continue
        k = keyof(pred["det"], pred["two_bridge_class_q"], pred["mirror_class_q"])
        if k in seen:
            continue
        seen[k] = cf
        rows.append(pred)
    rows.sort(key=lambda r: (r["det"], r["p_over_q"]))
    out = HERE / "steckplan_vorhersagen.json"
    out.write_text(json.dumps(rows, indent=1))
    print(f"[PLAN] {len(rows)} distinct 2-bridge plug plans (mirror classes "
          f"merged) → {out.name}")
    for r in rows[:20]:
        print(f"   {str(r['continued_fraction']):14s} p/q={r['p_over_q']:>7s} "
              f"det={r['det']:3d}  {r['identity_bet']}")


def cmd_make(args):
    entry, folder = make_knot(args.cf, HERE / "out_fseries")
    man_path = HERE / "MANIFEST.json"
    man = json.loads(man_path.read_text()) if man_path.exists() else {"knots": []}
    man["knots"] = [k for k in man["knots"] if k["name"] != entry["name"]] + [entry]
    man_path.write_text(json.dumps(man, indent=1))
    if folder:
        print(f"[MAKE] {entry['name']}: det bet {entry['predictions']['det']}, "
              f"crossing gate {entry['proj_crossings']}, "
              f"Nyquist {entry['fseries']['nyquist_dev_rel']} → {folder}")
    else:
        print(f"[MAKE] {entry['name']} REJECTED: {entry['rejected']}")


def cmd_demo(args):
    print("[DEMO] Smoke test machine 02")
    # 1) Conway computation against textbook anchors
    anchors = [([3], (3, 1)), ([2, 2], (5, 2)), ([3, 2], (7, 3)),
               ([2, 3, 2], (16, 7))]
    for cf, want in anchors:
        got = cf_to_pq(cf)
        print(f"  (1) {cf} → p/q = {got[0]}/{got[1]} "
              f"{'✓' if got == want else f'✗ FAIL (target {want})'}")
        assert got == want
    # 2) Link defence: p even
    e, f = make_knot([2, 3, 2], HERE / "demo_smoke")
    print(f"  (2) [2,3,2] p=16 even → rejected: {e['rejected']}  ✓")
    assert f is None
    # 3) Build [2,2] (bet b(5,2), det 5): components + crossing gate + Nyquist
    e, f = make_knot([2, 2], HERE / "demo_smoke")
    pc = e["proj_crossings"]
    print(f"  (3) [2,2]: crossings y projection {pc['got']}/{pc['want']} "
          f"{'✓' if pc['ok'] else '✗ FAIL'} · embedding "
          f"d_min_rel={e['embedding']['d_min_rel']} "
          f"{'✓' if e['embedding']['embedded'] else '✗'} · "
          f"Nyquist dev={e['fseries']['nyquist_dev_rel']} ✓")
    assert pc["ok"] and e["embedding"]["embedded"]
    # 4) Trefoil anchor [3]
    e, f = make_knot([3], HERE / "demo_smoke")
    print(f"  (4) [3] (trefoil bet, det 3): crossings "
          f"{e['proj_crossings']['got']}/3 · embedded="
          f"{e['embedding']['embedded']}  ✓")
    assert e["proj_crossings"]["ok"] and e["embedding"]["embedded"]
    print("[DEMO] 4/4 green — machine ready to fire (firing only after Monday).")


def main():
    ap = argparse.ArgumentParser(description="Machine 02 — plugging_tangles")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("plan")
    pm = sub.add_parser("make")
    pm.add_argument("cf", type=int, nargs="+", metavar="a_i")
    sub.add_parser("demo")
    args = ap.parse_args()
    {"plan": cmd_plan, "make": cmd_make, "demo": cmd_demo}[args.cmd](args)


if __name__ == "__main__":
    main()
