# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 124 — forging_the_fibers (fibre forge, Stallings axis).

3-lane station (117 architecture), series F_k = (σ₁³σ₂⁻¹)^k for
k ∈ {1,2,4,5} (k=3 = link, omitted). Homogeneous (σ₁ always +, σ₂
always −) ⟹ Stallings: fibred by construction, fibre genus
g = 2k−1 DOSED (1/3/7/9). fs_k1 = trefoil (Markov calibration anchor).
Two-sided self-test: theorem checks pipeline fibered field, pipeline
checks the build (Δ monic, deg Δ = 4k−2 exactly).

BETS: pre-registered in 124_faser_schmiede.txt (W0 trefoil anchor ·
W1 fibered=true for all · W2 Δ monic+degree exact · W3 genus=2k−1 ·
W4 recurrence sister · W5 θ/6D). Do NOT fire.
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

# 118 lesson: capacity ~ width (1.0 carries ~22, 2.4 carries 62);
# top rung k=5 has 20 crossings -> 2.4 is ample.
X_STN_124 = 2.4
X_RUN_124 = 3.0

KS = (1, 2, 4, 5)            # k=3: (3-cycle)^3 = id -> link, not a knot


def braid_word(k: int):
    """(σ₁³ σ₂⁻¹)^k — homogeneous, 4k crossings, permutation 3-cycle^k."""
    return ([(1, +1)] * 3 + [(2, -1)]) * k


def word_perm(word):
    p = list(range(3))
    for g, _s in word:
        p[g - 1], p[g] = p[g], p[g - 1]
    return p


def is_homogeneous(word):
    """Sign constancy per generator — MEASURED, not claimed."""
    signs = {}
    for g, s in word:
        if g in signs and signs[g] != s:
            return False
        signs[g] = s
    return True


def soll_crossings(k: int):
    return 4 * k


def build_curve(k: int):
    word = braid_word(k)
    paths, end_lane = build_station(word, x_stn=X_STN_124, x_run=X_RUN_124)
    pts, cur, visited = [], 0, []
    for _ in range(3):
        visited.append(cur)
        pts += paths[cur]
        lane = end_lane[cur]
        pts += connector(lane, x_run=X_RUN_124)
        cur = lane
    if cur != 0 or sorted(visited) != [0, 1, 2]:
        raise RuntimeError(f"k={k}: permutation is not a 3-cycle — closure broken")
    return np.asarray(pts, float)


def make_rung(k: int, out_root: Path):
    name = f"fs_k{k}"
    P = build_curve(k)
    n_c = station_crossings(P, x_stn=X_STN_124)
    soll = soll_crossings(k)
    entry = {"name": name, "machine": "124_forging_the_fibers",
             "braid_word": f"(s1^3 s2^-1)^{k}",
             "homogeneous": is_homogeneous(braid_word(k)),
             "fiber_genus_soll": 2 * k - 1,
             "station_crossings": n_c, "station_crossings_soll": soll,
             "predictions": {
                 "W0": "fs_k1 = 3_1 (Markov-Destabilisierung, Eich-Anker)",
                 "W1": "fibered = true ALLE Sprossen (Stallings; testet "
                       "Pipeline-Feld und Bau beidseitig)",
                 "W2": "Delta monisch + deg Delta = 4k-2 exakt",
                 "W3": "genus = 2k-1 exakt (1/3/7/9)",
                 "W4": "Delta-Rekurrenz: Fit aus k=1,2,4 trifft k=5 "
                       "(Burau-Transfer, ehrlich zweiseitig)",
                 "W5": "theta pro Sprosse + 6D trennt Nachbarn"}}
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
                source="maschinen_ideen/124", n_terms=nt,
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
    print("[DEMO] Smoke test machine 124 (fibre forge, Stallings)")
    out = HERE / "demo_smoke"
    ok = 0
    # S1: combinatorics + homogeneity certificate + k=3 counter-check
    s1 = True
    for k in KS:
        w = braid_word(k)
        n_c = station_crossings(build_curve(k), x_stn=X_STN_124)
        soll = soll_crossings(k)
        hit = (n_c == soll)
        hom = is_homogeneous(w)
        s1 &= hit and hom
        print(f"  (S1) k={k}: crossings {n_c} (expected {soll}) "
              f"{'✓' if hit else '✗'} · homogeneous {'✓' if hom else '✗'} · "
              f"g_expected={2 * k - 1}")
    link3 = word_perm(braid_word(3)) == [0, 1, 2]
    s1 &= link3
    print(f"  (S1) k=3 counter-check: permutation = id ⟹ link, omitted "
          f"{'✓' if link3 else '✗'}")
    ok += s1
    # S2: forge surgically — PER CLASS (build lesson from the first flight):
    # perm(F_k) = cycle^k ⟹ the closure connector sequence hangs on the
    # permutation CLASS: {1,4} (cycle) vs {2,5} (cycle²). Bit-
    # Identity outside the station holds within each class —
    # unlike 122/123 (there all words share ONE permutation).
    outs = {k: outside_station(build_curve(k), x_stn=X_STN_124)
            for k in KS}
    s2 = True
    for cls in ((1, 4), (2, 5)):
        a, b = outs[cls[0]], outs[cls[1]]
        same = a.shape == b.shape and bool(np.array_equal(a, b))
        s2 &= same
        print(f"  (S2) class k∈{cls}: outside station "
              f"{a.shape[0]} points bit-identical {'✓' if same else '✗'}")
    diff = not np.array_equal(outs[1], outs[2])
    print(f"  (S2) class separation (cycle vs cycle² = different "
          f"connector sequence, MEASURED): {'✓' if diff else '— identical?'}")
    ok += s2
    # S3: embedding + Nyquist (mandatory k=1,2,4; k=5 may hit the format limit)
    s3 = True
    for k in KS:
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
    print(f"[DEMO] {ok}/3 green — {'forge stands.' if ok == 3 else 'RED.'}"
          f" Bets: 124_faser_schmiede.txt. Doctrine: do NOT fire.")
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
