# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 122 — deepening_the_twins (twin wells, Becerra axis stage 3).

3-lane station imported from 117 (word-agnostic, 113→81→118 pattern),
own width for up to 100 crossings. Fixed base B = (σ₁σ₂)⁴
(closure = T(3,4) = 8_19), twins W_k = γ_k·B with γ_k from the
lower central series of P₃ (γ₂=[σ₁²,σ₂²], γ_{j+1}=[γ_j,σ₂²]).
Per pair (B, γ_k·B) by Stanford EVERYTHING below degree k is frozen —
the series measures separation DEPTH (a₂/v₃/det/θ/6D).

BETS: pre-registered in 122_zwillings_brunnen.txt (W0 8_19 anchor ·
W1 Stanford well · W2 separation depth two-sided · W3 det contrast ·
W4 θ well · W5 6D · W6 no 8_19 relapse). Do NOT fire.
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

# Station geometry: 100 crossings need room. Lesson from 118:
# Capacity scales ~ with width (1.0 carries ~22, 2.4 carries 62).
# 4.8 dimensioned for the 100 rung; constant over ALL
# rungs ⟹ the well stays surgical (S2).
X_STN_122 = 4.8
X_RUN_122 = 5.4

A12 = [(1, +1), (1, +1)]                 # σ₁²
A23 = [(2, +1), (2, +1)]                 # σ₂²
BASE = [(1, +1), (2, +1)] * 4            # (σ₁σ₂)⁴ -> T(3,4) = 8_19


def inv(word):
    return [(g, -s) for g, s in reversed(word)]


def comm(x, y):
    return x + y + inv(x) + inv(y)


def gamma(k: int):
    """γ_k element of P₃: γ₂=[A₁₂,A₂₃], γ_{j+1}=[γ_j,A₂₃].
    Length: |γ₂|=8, then 2·(|γ_j|+2)."""
    if k < 2:
        raise ValueError("well starts at k=2")
    g = comm(A12, A23)
    for _ in range(k - 2):
        g = comm(g, A23)
    return g


def braid_word(k: int):
    """k=0: base B. k>=2: γ_k·B — twin of depth k."""
    if k == 0:
        return list(BASE)
    return gamma(k) + list(BASE)


def word_perm(word):
    """Permutation of the word (generator i swaps strands i-1,i)."""
    p = list(range(3))
    for g, _s in word:
        p[g - 1], p[g] = p[g], p[g - 1]
    return p


def soll_crossings(k: int):
    return len(braid_word(k))


def build_curve(k: int):
    """Closed waypoint curve (117 architecture, own width)."""
    word = braid_word(k)
    paths, end_lane = build_station(word, x_stn=X_STN_122, x_run=X_RUN_122)
    pts, cur, visited = [], 0, []
    for _ in range(3):
        visited.append(cur)
        pts += paths[cur]
        lane = end_lane[cur]
        pts += connector(lane, x_run=X_RUN_122)
        cur = lane
    if cur != 0 or sorted(visited) != [0, 1, 2]:
        raise RuntimeError(f"k={k}: permutation not a 3-cycle — closure broken")
    return np.asarray(pts, float)


def make_rung(k: int, out_root: Path):
    name = "tw_base" if k == 0 else f"tw_k{k}"
    P = build_curve(k)
    n_c = station_crossings(P, x_stn=X_STN_122)
    soll = soll_crossings(k)
    entry = {"name": name, "machine": "122_deepening_the_twins",
             "braid_word": ("(s1 s2)^4" if k == 0
                            else f"gamma_{k}(P3) * (s1 s2)^4"),
             "station_crossings": n_c, "station_crossings_soll": soll,
             "predictions": {
                 "W0": "tw_base = T(3,4) = 8_19 (KnotInfo-Eichanker)",
                 "W1": "a2(tw_k)=a2(8_19) fuer k>=3; v3 gleich fuer k>=4",
                 "W2": "a2 trennt k=2, v3 trennt k=3 (ehrlich zweiseitig)",
                 "W3": "det trennt JEDES Paar (det ist kein finite type)",
                 "W4": "theta-Brunnen: tiefste Stufe k* mit theta-Trennung",
                 "W5": "6D trennt k=4/k=5-Paar (Geometrie unter Algebra-Decke)",
                 "W6": "tw_k reduziert NICHT auf 8_19 (echter Zwilling)"}}
    if n_c != soll:
        entry["rejected"] = f"Kreuzungszahl {n_c} != Soll {soll}"
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
                source="maschinen_ideen/122", n_terms=nt,
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
    print("[DEMO] Smoke test machine 122 (twin wells, γ depth)")
    out = HERE / "demo_smoke"
    ks = [0, 2, 3, 4, 5]
    ok = 0
    # S1: combinatorics + γ purity
    s1 = True
    for k in ks:
        n_c = station_crossings(build_curve(k), x_stn=X_STN_122)
        soll = soll_crossings(k)
        hit = (n_c == soll)
        pure = (k == 0) or (word_perm(gamma(k)) == [0, 1, 2])
        s1 &= hit and pure
        print(f"  (S1) k={k}: crossings {n_c} (expected {soll}) "
              f"{'✓' if hit else '✗'} · γ-Purity "
              f"{'✓' if pure else '✗'}")
    ok += s1
    # S2: well surgical (bit-identical outside the station)
    outs = [outside_station(build_curve(k), x_stn=X_STN_122) for k in ks]
    s2 = all(o.shape == outs[0].shape and bool(np.array_equal(o, outs[0]))
             for o in outs[1:])
    print(f"  (S2) {len(ks)} rungs outside station: "
          f"{outs[0].shape[0]} points bit-identical {'✓' if s2 else '✗'}")
    ok += s2
    # S3: embedding + Nyquist (mandatory base/k2/k3/k4; k5 may hit the
    #     format limit — then it is measured)
    s3 = True
    for k in ks:
        entry, folder = make_rung(k, out)
        if folder is None:
            print(f"  (S3) {entry['name']}: REJECTED — {entry['rejected']}"
                  f"{' (format limit measured, ok)' if k == 5 else ' ✗'}")
            if k != 5:
                s3 = False
        else:
            print(f"  (S3) {entry['name']}: "
                  f"d_min_rel={entry['embedding']['d_min_rel']}"
                  f" · n_terms={entry['fseries']['n_terms']}"
                  f" · nyq={entry['fseries']['nyquist_dev_rel']} ✓")
    ok += s3
    print(f"[DEMO] {ok}/3 green — {'well stands.' if ok == 3 else 'RED.'}"
          f" Bets: 122_zwillings_brunnen.txt. Doctrine: do NOT fire.")
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
