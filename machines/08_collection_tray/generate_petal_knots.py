# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 08/Petal — petal diagram generator (one permutation = one knot).

Petal diagram: the projection has ONE single multiple point; p petals
(p odd) run star-shaped through the centre, the height PERMUTATION of the
centre passages determines the knot completely [Adams et al., petal
projections]. Radically machine-suited: the whole parameter space is S_p.

CONSTRUCTION (star construction, details in the petal_curve docstring): Passage i =
straight chord through the centre at height h_i, entry angle a_i = π(1+i)+iπ/p;
connectors = short π/p arcs on the outer cylinder in disjoint sectors.

BETS: structure bets only (p petals ⇒ petal number ≤ p; identities
measured blind at firing). Anchor UNVERIFIED (check at firing time): trefoil
has petal number 5 — a 5-permutation from the grid must hit 3_1.

Doctrine: do NOT fire before Monday — only out_fseries/ + MANIFEST.
"""
from __future__ import annotations
import argparse
import itertools
import math
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402

PTS_PER_PETAL = 160
Z_SCALE = 0.22          # height spacing of the centre passages


def petal_curve(perm):
    """Closed 3D curve from the height permutation (star construction).

    Passage i = STRAIGHT chord through the centre at height h_i, direction
    β_i = i·π/p (p tangent lines, pairwise distinct — the one
    multiple point lies vertically stacked on the z-axis). Connectors
    run OUTSIDE on the cylinder r=1 (θ monotone +π+π/p, z linear
    h_i → h_{i+1}). Teardrop lesson: petals with 216° span in the
    plane overlap inevitably — the chord+cylinder form separates
    inside (chords, heights distinct) from outside (cylinder mantle).
    The curve is DETERMINISTIC from the permutation; matching with the
    literature petal-knot assignment = firing-time question (in the manifest).
    """
    p = len(perm)
    if p % 2 == 0:
        raise ValueError("petal number must be odd")
    heights = [Z_SCALE * (h - (p - 1) / 2.0) for h in perm]
    pts = []
    n_ch = PTS_PER_PETAL // 2 * 2 + 1        # odd → midpoint exactly 0
    for i in range(p):
        # entry angle of passage i (chord direction alternates —
        # Closure computation: a_i = π(1+i) + iπ/p; connectors are then only
        # π/p wide and lie in DISJOINT outer sectors = petal picture)
        a = np.pi * (1 + i) + i * np.pi / p
        h0, h1 = heights[i], heights[(i + 1) % p]
        d = np.array([np.cos(a), np.sin(a)])
        t = np.linspace(0.0, 1.0, n_ch)
        s = 1.0 - 2.0 * t                    # +1 (entry) → −1 (exit)
        chord = np.column_stack([s * d[0], s * d[1], np.full_like(s, h0)])
        pts.append(chord[:-1])
        # Connector: θ from a+π by +π/p, z h0 → h1 (short outer arc)
        tt = np.linspace(0.0, 1.0, PTS_PER_PETAL)[:-1]
        th = a + np.pi + tt * (np.pi / p)
        z = h0 + (h1 - h0) * tt
        pts.append(np.column_stack([np.cos(th), np.sin(th), z]))
    return np.vstack(pts)


def make_knot(perm, out_root):
    name = "pt_" + "".join(str(h) for h in perm)
    entry = {"name": name, "machine": "08_collection_tray/petal",
             "permutation": list(perm),
             "predictions": {
                 "petal_number_bound": len(perm),
                 "source": "Petal-Projektion (Adams-Linie); Identität misst "
                           "die Pipeline blind",
             }}
    P = petal_curve(perm)
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(P, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/08_collection_tray",
                              n_terms=max(64, 16 * len(perm)),
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    return entry, folder


def canon_class(perm):
    """Petal permutations are cyclically + mirror equivalent (roughly):
    Normal form via rotations of the entries + reversal (dedup in the grid)."""
    p = len(perm)
    variants = []
    for r in range(p):
        rot = tuple((perm[(i + r) % p]) for i in range(p))
        variants.append(rot)
        variants.append(rot[::-1])
    return min(variants)


def cmd_plan(args):
    """Grid: all 5-permutations with h_0=0, de-duplicated (rough classes)."""
    seen, rows = set(), []
    for rest in itertools.permutations(range(1, args.petals)):
        perm = (0,) + rest
        c = canon_class(perm)
        if c in seen:
            continue
        seen.add(c)
        rows.append(list(perm))
    out = HERE / f"petal_plan_p{args.petals}.json"
    out.write_text(json.dumps({"petals": args.petals, "n_classes": len(rows),
                               "permutations": rows}, indent=1))
    print(f"[PLAN] p={args.petals}: {len(rows)} rough permutation classes "
          f"(of {math.factorial(args.petals - 1)} with h0=0) → {out.name}")


def cmd_make(args):
    entry, folder = make_knot(tuple(args.perm), HERE / "out_fseries")
    man_path = HERE / "MANIFEST.json"
    man = json.loads(man_path.read_text()) if man_path.exists() else {"knots": []}
    man["knots"] = [k for k in man["knots"] if k["name"] != entry["name"]] + [entry]
    man_path.write_text(json.dumps(man, indent=1))
    print(f"[MAKE] {entry['name']} → {folder if folder else entry['rejected']}")


def cmd_demo(args):
    print("[DEMO] smoke test petal generator (machine 08)")
    # 1) Geometry: p passages through the centre, heights = permutation
    perm = (0, 3, 1, 4, 2)
    P = petal_curve(perm)
    r = np.hypot(P[:, 0], P[:, 1])
    center_idx = np.where(r < 1e-9)[0]
    ok1 = len(center_idx) == 5
    zs = P[center_idx, 2] / Z_SCALE + 2.0
    ok1 &= np.allclose(sorted(zs.round(6)), [0, 1, 2, 3, 4])
    print(f"  (1) 5 centre passages, heights {list(np.round(zs, 2))} = "
          f"Permutation  {'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # 2) Build + gates (trefoil candidate (0,3,1,4,2), UNVERIFIED)
    e, f = make_knot(perm, HERE / "demo_smoke")
    ok2 = f is not None
    print(f"  (2) build pt_03142: embedded={e['embedding']['embedded']} "
          f"d_min_rel={e['embedding']['d_min_rel']} Nyquist="
          f"{e.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok2 else '✗ FAIL: ' + str(e.get('rejected'))}")
    assert ok2
    # 3) Identity permutation (0,1,2,3,4) = spiral → unknot bet, buildable
    e3, f3 = make_knot((0, 1, 2, 3, 4), HERE / "demo_smoke")
    print(f"  (3) pt_01234 (staircase spiral, unknot bet): embedded="
          f"{e3['embedding']['embedded']}  {'✓' if f3 else '✗ FAIL'}")
    assert f3 is not None
    print("[DEMO] 3/3 green — petal axis ready to fire (firing after Monday).")


def main():
    ap = argparse.ArgumentParser(description="Machine 08 — petal knots")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("plan"); p1.add_argument("--petals", type=int, default=5)
    p2 = sub.add_parser("make")
    p2.add_argument("perm", type=int, nargs="+")
    sub.add_parser("demo")
    args = ap.parse_args()
    {"plan": cmd_plan, "make": cmd_make, "demo": cmd_demo}[args.cmd](args)


if __name__ == "__main__":
    main()
