# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 07 — closing_differently (closure rule as an axis of its own).

Same braid, different closure = different knot. This machine is
BUILT INDEPENDENTLY in the 07 folder — the_knotfamily_machine stays untouched
(Werner gate respected, no behaviour drift for the existing stock).

THREE CLOSURES on the same N-strand braid:
  trace : trace closure — end slot s → start slot s, return
          arcs staggered in y (classical Markov world, braid index ≤ N).
  plat  : cap pairs (0,1),(2,3),… top+bottom (N even) —
          bridge number ≤ N/2.
  shift : q-TURN OFFSET — end slot s → start slot (s+q) mod N
          (the torus rail; identical word + offset = torus/cable
          strata systematically).
PARADE BET (demo): σ1³ on 2 strands — trace = trefoil, plat = unknot.
One word, two knots: the closure IS a degree of freedom.

WORD GENERATOR (sub-axis C): periodic words (w^m = combinatorial
C_m periodicity for free) + substitution rules (Fibonacci: a→ab, b→a).

Component count is predicted EXACTLY (permutation cycles) and counter-checked at the
build. At firing: Disc-#257 convention + POLYGON_PD=1 witness.
Doctrine: do NOT ignite before Monday.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import (embedding_report, polyline_to_fseries,  # noqa: E402
                            count_projection_crossings)

PTS_PER_HALFTWIST = 17        # odd — vertex trap (machine-02 lesson)
Z_MARGIN = 0.8                # straight run-out before the closures


# -------------------------------------------------------------- braid words

def parse_word(txt):
    """'s1 s2 -s1' or 's1^3' → list (generator index from 0, half-twists)."""
    word = []
    for tok in txt.split():
        if "^" in tok:
            g, p = tok.split("^")
            k = int(p)
        else:
            g, k = tok, 1
        sgn = -1 if g.startswith("-") else 1
        i = int(g.lstrip("-s")) - 1
        word.append((i, sgn * k))
    return word


def substitution_word(rules, axiom, depth, letter_map):
    """Substitution system (e.g. Fibonacci a→ab, b→a) → braid word."""
    s = axiom
    for _ in range(depth):
        s = "".join(rules.get(ch, ch) for ch in s)
    return [letter_map[ch] for ch in s], s


def braid_permutation(word, n_strands):
    perm = list(range(n_strands))
    for (i, k) in word:
        if abs(k) % 2 == 1:
            perm[i], perm[i + 1] = perm[i + 1], perm[i]
    return perm


def predict_components(word, n_strands, closure, q=1):
    """Components = cycles of closure ∘ braid permutation (exact)."""
    perm = braid_permutation(word, n_strands)   # start slot → end slot
    if closure in ("trace", "shift"):
        shift = q if closure == "shift" else 0
        nxt = {s: (perm[s] + shift) % n_strands for s in range(n_strands)}
        seen, n_comp = set(), 0
        for s0 in range(n_strands):
            if s0 in seen:
                continue
            n_comp += 1
            s = s0
            while s not in seen:
                seen.add(s)
                s = nxt[s]
        return n_comp
    # plat: knot graph over caps bottom/top + wires
    pair = {s: s + 1 if s % 2 == 0 else s - 1 for s in range(n_strands)}
    end_of = {s: perm[s] for s in range(n_strands)}
    start_of = {v: k for k, v in end_of.items()}
    seen, n_comp = set(), 0
    for w0 in range(n_strands):
        if w0 in seen:
            continue
        n_comp += 1
        w, up = w0, True
        while True:
            seen.add(w)
            if up:
                w, up = start_of[pair[end_of[w]]], False
            else:
                w, up = pair[w], True
            if w == w0 and up:
                break
    return n_comp


# ------------------------------------------------------------------ geometry

def build_braid_wires(word, n_strands):
    """Wires of the braid (without closure). Returns (wires, slot_of, z_top)."""
    wires = [[np.array([float(s), 0.0, 0.0])] for s in range(n_strands)]
    slot_of = list(range(n_strands))
    wire_at = list(range(n_strands))
    z = 0.0
    for (i, k) in word:
        if k == 0:
            continue
        dz = abs(k) * 1.0
        wa, wb = wire_at[i], wire_at[i + 1]
        xc, r = i + 0.5, 0.5
        sgn = 1.0 if k > 0 else -1.0
        th = np.linspace(0.0, abs(k) * np.pi,
                         abs(k) * PTS_PER_HALFTWIST + 1)[1:]
        zz = z + dz * th / th[-1]
        wires[wa] += list(np.column_stack([xc - r * np.cos(th),
                                           -sgn * r * np.sin(th), zz]))
        wires[wb] += list(np.column_stack([xc + r * np.cos(th),
                                           sgn * r * np.sin(th), zz]))
        for w in range(n_strands):
            if w not in (wa, wb):
                wires[w].append(np.array([float(slot_of[w]), 0.0, z + dz]))
        if abs(k) % 2 == 1:
            slot_of[wa], slot_of[wb] = i + 1, i
            wire_at[i], wire_at[i + 1] = wb, wa
        z += dz
    # straight run-out (separation closure ↔ twist zone)
    for w in range(n_strands):
        wires[w].append(np.array([float(slot_of[w]), 0.0, z + Z_MARGIN]))
        wires[w].insert(0, np.array([float(w), 0.0, -Z_MARGIN]))
    return [np.asarray(w) for w in wires], slot_of, z + Z_MARGIN


def _arc(p0, p1, n=15, bulge=None):
    """Semicircular arc p0→p1, bulging in direction bulge."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    mid = (p0 + p1) / 2.0
    r = np.linalg.norm(p1 - p0) / 2.0
    b = np.asarray(bulge, float)
    b = b / np.linalg.norm(b)
    th = np.linspace(0, np.pi, n + 2)[1:-1]
    u = (p1 - p0) / max(np.linalg.norm(p1 - p0), 1e-12)
    return np.array([mid - u * r * np.cos(t) + b * r * np.sin(t) for t in th])


def close_braid(wires, slot_of, z_top, n_strands, closure, q=1):
    """Apply closure → ONE closed polyline (or None for >1)."""
    if closure == "plat":
        if n_strands % 2:
            raise ValueError("plat needs an even strand count")
        pair = {s: s + 1 if s % 2 == 0 else s - 1 for s in range(n_strands)}
        wire_at_end = {slot_of[w]: w for w in range(n_strands)}
        steps, w, up = [], 0, True
        while True:
            steps.append((w, up))
            if up:
                w, up = wire_at_end[pair[slot_of[w]]], False
            else:
                w, up = pair[w], True       # wire with start slot s = wires[s]
            if (w, up) == (0, True):
                break
        if len({wi for (wi, _) in steps}) < n_strands:
            return None
        P = []
        for (w, up) in steps:
            if up:
                pts = wires[w]
                P.append(pts)
                partner = wires[wire_at_end[pair[slot_of[w]]]]
                P.append(_arc(pts[-1], partner[-1], bulge=[0, 0, 1]))
            else:
                pts = wires[w][::-1]
                P.append(pts)
                P.append(_arc(pts[-1], wires[pair[w]][0], bulge=[0, 0, -1]))
        return np.vstack(P)

    # trace / shift: return staggered in y (start slot = wire index)
    shift = q if closure == "shift" else 0
    P, visited, w = [], set(), 0
    while True:
        visited.add(w)
        pts = wires[w]
        P.append(pts)
        end_slot = slot_of[w]
        target = (end_slot + shift) % n_strands
        depth = 1.2 + 0.45 * end_slot
        top = pts[-1]
        bot = np.array([float(target), 0.0, -Z_MARGIN])
        # Right-angle return (two demo lessons: diagonal descents
        # form hairpins with ~0 opening → false embedding red):
        # out in y → straight down at constant (x,y) → horizontally
        # across under the floor → in. Relay depth separates the strands.
        ret = [np.array([top[0], depth, z_top + 0.5]),
               np.array([top[0], depth, -Z_MARGIN - 0.5]),
               np.array([float(target), depth, -Z_MARGIN - 0.5]),
               np.array([float(target), depth * 0.3, -Z_MARGIN - 0.3])]
        P.append(np.asarray(ret))
        w = target
        if w in visited:
            break
    if len(visited) < n_strands:
        return None
    return np.vstack(P)


# ----------------------------------------------------------------- generation

def make_knot(name, word, n_strands, closure, q, out_root):
    n_comp = predict_components(word, n_strands, closure, q)
    entry = {"name": name, "machine": "07_closing_differently",
             "word": word, "n_strands": n_strands, "closure": closure,
             "q": q if closure == "shift" else None,
             "predicted_components": n_comp,
             "predictions": {
                 "bridge_bound": n_strands // 2 if closure == "plat" else None,
                 "braid_index_bound": n_strands if closure != "plat" else None,
             }}
    if n_comp != 1:
        entry["rejected"] = f"{n_comp} Komponenten (Link) — kein Knoten-Record"
        return entry, None
    wires, slot_of, z_top = build_braid_wires(word, n_strands)
    P = close_braid(wires, slot_of, z_top, n_strands, closure, q)
    if P is None:
        raise RuntimeError(f"{name}: built components ≠ prediction — bug!")
    want = sum(abs(k) for (_, k) in word)
    nx = count_projection_crossings(P, axis=1)
    entry["proj_crossings_twistzone_want"] = want
    entry["proj_crossings_got"] = nx        # returns may add more
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(P, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/07_closing_differently",
                              n_terms=max(64, 10 * want), nyquist_max_rel=0.02)
    entry["fseries"] = rep
    return entry, folder


# ----------------------------------------------------------------------- CLI

def cmd_make(args):
    word = parse_word(" ".join(args.word))
    name = (f"br_{args.closure}{args.q if args.closure == 'shift' else ''}_"
            + "_".join(f"{'m' if k < 0 else ''}{i+1}x{abs(k)}"
                       for (i, k) in word))
    entry, folder = make_knot(name, word, args.strands, args.closure, args.q,
                              HERE / "out_fseries")
    man_path = HERE / "MANIFEST.json"
    man = json.loads(man_path.read_text()) if man_path.exists() else {"knots": []}
    man["knots"] = [k for k in man["knots"] if k["name"] != entry["name"]] + [entry]
    man_path.write_text(json.dumps(man, indent=1))
    if folder:
        print(f"[MAKE] {entry['name']} → {folder}")
    else:
        print(f"[MAKE] {entry['name']} REJECTED: {entry['rejected']}")


def cmd_demo(args):
    print("[DEMO] Smoke test machine 07 (closure axis)")
    w3 = parse_word("s1^3")
    # 1) Parade bet: σ1³ — trace 1 component (trefoil), plat 1 (unknot)
    c_tr = predict_components(w3, 2, "trace")
    c_pl = predict_components(w3, 2, "plat")
    print(f"  (1) σ1³: trace → {c_tr} comp. (bet trefoil), plat → {c_pl} "
          f"comp. (bet UNKNOT) — one word, two knots  "
          f"{'✓' if c_tr == 1 and c_pl == 1 else '✗ FAIL'}")
    assert c_tr == 1 and c_pl == 1
    # 2) Component maths: σ1² on 2 strands trace = Hopf link (2 comp.)
    c2 = predict_components(parse_word("s1^2"), 2, "trace")
    print(f"  (2) σ1² trace: {c2} components (Hopf link target 2)  "
          f"{'✓' if c2 == 2 else '✗ FAIL'}")
    assert c2 == 2
    # 3) shift closure: identical empty word, 3 strands, q=1 → 1 comp.
    c3 = predict_components([], 3, "shift", q=1)
    c3b = predict_components([], 4, "shift", q=2)
    print(f"  (3) empty braid shift: n=3,q=1 → {c3} comp. (target 1) · n=4,q=2 → "
          f"{c3b} (target gcd=2)  {'✓' if c3 == 1 and c3b == 2 else '✗ FAIL'}")
    assert c3 == 1 and c3b == 2
    # 4) Fibonacci substitution as word generator
    word, s = substitution_word({"a": "ab", "b": "a"}, "a", 5,
                                {"a": (0, 1), "b": (1, -1)})
    print(f"  (4) Fibonacci word depth 5: '{s}' → {len(word)} letters "
          f"(Fibonacci target 13)  {'✓' if len(word) == 13 else '✗ FAIL'}")
    assert len(word) == 13
    # 5) build trace σ1³ (trefoil bet) + gates
    e, f = make_knot("demo_trace_s1x3", w3, 2, "trace", 1, HERE / "demo_smoke")
    ok5 = f is not None
    print(f"  (5) build trace σ1³: embedded={e.get('embedding', {}).get('embedded')} "
          f"Nyquist={e.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok5 else '✗ FAIL: ' + str(e.get('rejected'))}")
    assert ok5
    # 6) Build plat: 4 strands σ2³ (b(?,?) world) + gates
    e6, f6 = make_knot("demo_plat_s2x3", parse_word("s2^3"), 4, "plat", 1,
                       HERE / "demo_smoke")
    ok6 = f6 is not None
    print(f"  (6) build plat 4-strand σ2³: comp. prediction "
          f"{e6['predicted_components']}, embedded="
          f"{e6.get('embedding', {}).get('embedded')}  "
          f"{'✓' if ok6 else '✗ FAIL: ' + str(e6.get('rejected'))}")
    assert ok6
    print("[DEMO] 6/6 green — closure axis ready to fire (firing after "
          "Monday; the_knotfamily_machine untouched).")


def main():
    ap = argparse.ArgumentParser(description="Machine 07 — closure axis")
    sub = ap.add_subparsers(dest="cmd", required=True)
    pm = sub.add_parser("make")
    pm.add_argument("word", nargs="+", help="e.g. s1^3 or 's1 -s2 s1'")
    pm.add_argument("--strands", type=int, required=True)
    pm.add_argument("--closure", choices=["trace", "plat", "shift"],
                    default="trace")
    pm.add_argument("--q", type=int, default=1)
    sub.add_parser("demo")
    args = ap.parse_args()
    {"make": cmd_make, "demo": cmd_demo}[args.cmd](args)


if __name__ == "__main__":
    main()
