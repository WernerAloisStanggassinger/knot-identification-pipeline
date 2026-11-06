# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 13 — shrinking_into_itself (scale symmetry, self-similar knots).

Iteration rule (scale-invariant, reversal crack 3): K_{d+1} arises from
K_d, by placing on the SHARPEST corner of K_d a λ-scaled copy of K_d
itself (06 mount machinery). Same rule at every stage
= self-SIMILARITY, not mere iteration.

BETS (connected-sum reading, honestly labelled as composite — crack 1):
  K_{d+1} = K_d # K_d  ⇒  Δ_d = Δ_0^(2^d),  det_d = det_0^(2^d).
  Trefoil staircase: det = 3, 9, 81, 6561 — exponentially sharp pipeline tests.
Depth hard-capped d<=3 (crack 2); λ delegated to the embedding gate.
Satellite docking (a true tower) = noted extension stage.

Doctrine: do NOT ignite before Monday.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "06_eating_own_output"))
sys.path.insert(0, str(HERE.parent))
from generate_satellite_knots import (graft, predict_connected_sum,  # noqa: E402
                                      trefoil_curve, TREFOIL_DELTA)
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402

MAX_DEPTH = 3


def self_similar(base, depth, gap_frac=0.09):
    """Scale-invariant iteration: the curve eats a λ-copy of itself.

    Docking rule (fixed, manifest): sharpest corner, 1 copy per stage.
    λ follows from gap_frac (window chord) — gate decides.
    """
    if depth > MAX_DEPTH:
        raise ValueError(f"depth {depth} > cap {MAX_DEPTH} (reversal crack 2)")
    K = np.asarray(base, float)
    history = []
    for d in range(depth):
        K_next, corners, windows = graft(K, K, 1, gap_frac=gap_frac)
        history.append({"depth": d + 1, "n_points": int(len(K_next)),
                        "corner": int(corners[0])})
        K = K_next
    return K, history


def make_knot(depth, out_root, gap_frac=0.09):
    base = trefoil_curve()
    name = f"ss_trefoil_d{depth}"
    n_copies = 2 ** depth          # K#K doubling per stage
    delta, det = predict_connected_sum([TREFOIL_DELTA] * n_copies)
    entry = {"name": name, "machine": "13_shrinking_into_itself",
             "depth": depth, "lambda_rule": f"gap_frac={gap_frac}, "
             "schärfste Ecke, 1 Kopie/Stufe (skalen-invariant)",
             "predictions": {
                 "composite": True,
                 "n_summands": n_copies,
                 "det": det,
                 "alexander_degree": len(delta) - 1,
                 "mode": "Verbundsummen-Lesart (ehrlich: composite; "
                         "Satelliten-Turm = Ausbaustufe)",
             }}
    P, hist = self_similar(base, depth, gap_frac=gap_frac)
    entry["iteration_history"] = hist
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate (λ/gap_frac zu eng für d)"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(P, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/13_shrinking_into_itself",
                              n_terms=min(384, 96 * 2 ** depth),
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    return entry, folder


def cmd_make(args):
    entry, folder = make_knot(args.depth, HERE / "out_fseries",
                              gap_frac=args.gap_frac)
    man_path = HERE / "MANIFEST.json"
    man = json.loads(man_path.read_text()) if man_path.exists() else {"knots": []}
    man["knots"] = [k for k in man["knots"] if k["name"] != entry["name"]] + [entry]
    man_path.write_text(json.dumps(man, indent=1))
    print(f"[MAKE] {entry['name']}: det bet {entry['predictions']['det']} → "
          f"{folder if folder else entry['rejected']}")


def cmd_demo(args):
    print("[DEMO] Smoke test machine 13 (self-similarity)")
    # 1) det tower bets: 3, 9, 81 (exponential)
    dets = [predict_connected_sum([TREFOIL_DELTA] * (2 ** d))[1]
            for d in range(3)]
    ok1 = dets == [3, 9, 81]
    print(f"  (1) det tower d=0,1,2: {dets} (expected 3,9,81)  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # 2) depth cap engages
    try:
        self_similar(trefoil_curve(), MAX_DEPTH + 1)
        ok2 = False
    except ValueError:
        ok2 = True
    print(f"  (2) depth cap d>{MAX_DEPTH} refused  "
          f"{'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # 3) build stage 1: trefoil # trefoil, det bet 9, all gates
    e1, f1 = make_knot(1, HERE / "demo_smoke")
    ok3 = f1 is not None and e1["predictions"]["det"] == 9
    print(f"  (3) d=1: det bet {e1['predictions']['det']} · embedded="
          f"{e1['embedding']['embedded']} · Nyquist="
          f"{e1.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok3 else '✗ FAIL: ' + str(e1.get('rejected'))}")
    assert ok3
    # 4) stage 2: (K#K)#(K#K), det bet 81 — the doll inside the doll
    e2, f2 = make_knot(2, HERE / "demo_smoke", gap_frac=0.2)
    ok4 = f2 is not None and e2["predictions"]["det"] == 81
    print(f"  (4) d=2: det bet {e2['predictions']['det']} · embedded="
          f"{e2['embedding']['embedded']} · points "
          f"{e2['iteration_history'][-1]['n_points']}  "
          f"{'✓' if ok4 else '✗ FAIL: ' + str(e2.get('rejected'))}")
    assert ok4
    print("[DEMO] 4/4 green — self-similarity machine ready to fire "
          "(firing after Monday; composite honestly labelled).")


def main():
    ap = argparse.ArgumentParser(description="Machine 13 — self-similarity")
    sub = ap.add_subparsers(dest="cmd", required=True)
    pm = sub.add_parser("make")
    pm.add_argument("--depth", type=int, default=2)
    pm.add_argument("--gap-frac", type=float, default=0.09)
    sub.add_parser("demo")
    args = ap.parse_args()
    {"make": cmd_make, "demo": cmd_demo}[args.cmd](args)


if __name__ == "__main__":
    main()
