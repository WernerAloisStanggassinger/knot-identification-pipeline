# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 118 — stacking_the_commutator (C₃ ladder, Becerra axis stage 2).

3-lane station imported from 117 (word-agnostic), but OWN
width x_stn=2.4 (the 117 width tore at k≥2 at the gate — measured);
new word: W(k) = ([[A₁₂,A₂₃],A₂₃])^k · σ₁σ₂ with A₁₂=σ₁², A₂₃=σ₂²
— double commutator = γ₃ element of P₃, 20 crossings per unit,
station target 20k+2. Stanford theorem as self-test: a₂ MUST be constant
(degree <3 frozen), v₃ may climb.

BETS: pre-registered in 118_kommutator_treppe.txt (W0 unknot ·
W1 a₂ constant, Stanford · W2 v₃ ladder, honestly two-sided · W3
Polynomiality · W4 θ staircase PAIR with 117 · W5 6D). Do NOT fire.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
YARD = HERE.parent
sys.path.insert(0, str(YARD))
sys.path.insert(0, str(YARD / "117_stepping_the_delta"))
from generate_delta_ladder import (  # noqa: E402  (station word-agnostic)
    build_station, connector, station_crossings, outside_station)
from machine_common import (  # noqa: E402
    densify, embedding_report, polyline_to_fseries)

# Own station geometry: 20k+2 crossings need room — the
# 117 width (x_stn=1.0) tore at k≥2 at the embedding gate (MEASURED,
# d_min_rel < 0.003). Width constant over ALL rungs ⟹ ladder
# stays surgical. Cross-bit check to dm_k0 is dropped as a result
# (different geometry) — replaced by firing bet „both unknot".
X_STN_118 = 2.4
X_RUN_118 = 3.0

A12 = [(1, +1), (1, +1)]                 # σ₁²
A23 = [(2, +1), (2, +1)]                 # σ₂²


def inv(word):
    return [(g, -s) for g, s in reversed(word)]


def comm(x, y):
    return x + y + inv(x) + inv(y)


def braid_word(k: int):
    """([[A₁₂,A₂₃],A₂₃])^k · σ₁σ₂ — 20k+2 crossings, permutation 3-cycle."""
    unit = comm(comm(A12, A23), A23)     # 20 crossings, γ₃ element
    return unit * k + [(1, +1), (2, +1)]


def build_curve(k: int):
    """Closed waypoint curve of rung k (117 architecture)."""
    word = braid_word(k)
    paths, end_lane = build_station(word, x_stn=X_STN_118, x_run=X_RUN_118)
    pts, cur, visited = [], 0, []
    for _ in range(3):
        visited.append(cur)
        pts += paths[cur]
        lane = end_lane[cur]
        pts += connector(lane, x_run=X_RUN_118)
        cur = lane
    if cur != 0 or sorted(visited) != [0, 1, 2]:
        raise RuntimeError(f"k={k}: permutation not a 3-cycle — closure broken")
    return np.asarray(pts, float)


def make_rung(k: int, out_root: Path):
    name = f"cm_k{k}"
    P = build_curve(k)
    n_c = station_crossings(P, x_stn=X_STN_118)
    entry = {"name": name, "machine": "118_stacking_the_commutator",
             "braid_word": f"([[s1^2,s2^2],s2^2])^{k} * s1 s2",
             "station_crossings": n_c, "station_crossings_soll": 20 * k + 2,
             "predictions": {
                 "W0": "k=0: Unknot (Basis-Sprosse, wie dm_k0 als Knoten)",
                 "W1": "a2 KONSTANT ueber die Leiter (Stanford, n=3)",
                 "W2": "v3 linear in k, Sprung != 0 (ehrlich zweiseitig)",
                 "W3": "Grad-d-Invariante = Polynom in k, Grad <= d/3",
                 "W4": "theta-Antwort auf C3-Treppe (Paar mit 117)",
                 "W5": "6D trennt cm_k0/cm_k1 (Grad<3-identisches Paar)"}}
    if n_c != 20 * k + 2:
        entry["rejected"] = f"Kreuzungszahl {n_c} != Soll {20 * k + 2}"
        return entry, None
    Pd = densify(P)
    emb = embedding_report(Pd)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = out_root / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = None
    for nt in (256, 384, 512, 768, 1024):
        try:
            rep = polyline_to_fseries(
                P, folder / f"{name}.fseries", knot_name=name,
                source="maschinen_ideen/118", n_terms=nt,
                nyquist_max_rel=0.02)
            break
        except RuntimeError:
            continue
    if rep is None:
        entry["rejected"] = "Nyquist rot bis n_terms=1024"
        return entry, None
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(_args):
    print("[DEMO] Smoke test machine 118 (commutator staircase, C3 ladder)")
    out = HERE / "demo_smoke"
    ok = 0
    # S1: combinatorics
    s1 = True
    for k in range(4):
        P = build_curve(k)
        n_c = station_crossings(P, x_stn=X_STN_118)
        s1 &= (n_c == 20 * k + 2)
        print(f"  (S1) k={k}: permutation 3-cycle ✓ · station crossings "
              f"{n_c} (target {20 * k + 2}) {'✓' if n_c == 20 * k + 2 else '✗'}")
    ok += s1
    # S2: ladder surgical (bit-identical outside the 118 station)
    O0 = outside_station(build_curve(0), x_stn=X_STN_118)
    O1 = outside_station(build_curve(1), x_stn=X_STN_118)
    O2 = outside_station(build_curve(2), x_stn=X_STN_118)
    s2 = (O0.shape == O1.shape == O2.shape
          and bool(np.array_equal(O0, O1)) and bool(np.array_equal(O1, O2)))
    print(f"  (S2) cm_k0/k1/k2 outside station: {O0.shape[0]} points "
          f"bit-identical {'✓' if s2 else '✗'}")
    ok += s2
    # S3: embedding + Nyquist (k=3 MAY fail at the format limit,
    #     then it is measured — mandatory rungs are k=0..2)
    s3 = True
    for k in range(4):
        entry, folder = make_rung(k, out)
        if folder is None:
            print(f"  (S3) cm_k{k}: REJECTED — {entry['rejected']}"
                  f"{' (format limit measured, ok)' if k == 3 else ' ✗'}")
            if k < 3:
                s3 = False
        else:
            print(f"  (S3) cm_k{k}: d_min_rel={entry['embedding']['d_min_rel']}"
                  f" · n_terms={entry['fseries']['n_terms']}"
                  f" · nyq={entry['fseries']['nyquist_dev_rel']} ✓")
    ok += s3
    print(f"[DEMO] {ok}/3 green — {'staircase stands.' if ok == 3 else 'RED.'}"
          f" Bets: 118_kommutator_treppe.txt. Doctrine: do NOT fire.")
    return 0 if ok == 3 else 1


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    a = ap.parse_args()
    if a.cmd == "demo":
        sys.exit(cmd_demo(a))


if __name__ == "__main__":
    main()
