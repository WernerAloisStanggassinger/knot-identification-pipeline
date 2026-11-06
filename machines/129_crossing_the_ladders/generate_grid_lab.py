# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 129 — crossing_the_ladders (grid lab: fibration × recurrence).

W(k,m) = (σ₁³σ₂⁻¹)^k · σ₁^{2m} on the 117 station — homogeneous for
all (k,m) ⟹ Stallings across the board (deg Δ = 4k+2m−2, genus =
2k+m−1 exact), per row k a Δ recurrence in m (fit m=0..2 predicts m=3
before, three independent predictions. m=0 row ≡ 124 rungs
BIT (cross-machine anchor). First 2-parameter bet grid.

BETS: 129_gitter_labor.txt (W0 bit anchor · W1 Stallings surface ·
W2 recurrence×row · W3 det table mod 8 · W4 sig monotonicity ·
W5 θ surface/6D). Do NOT fire.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
YARD = HERE.parent
sys.path.insert(0, str(YARD))
sys.path.insert(0, str(YARD / "117_stepping_the_delta"))
sys.path.insert(0, str(YARD / "124_forging_the_fibers"))
from generate_delta_ladder import (  # noqa: E402
    build_station, connector, station_crossings, outside_station)
from generate_fiber_forge import (  # noqa: E402  (parent machine 124)
    build_curve as build_124, X_STN_124)
from machine_common import (  # noqa: E402
    densify, embedding_report, polyline_to_fseries)

X_STN_129 = X_STN_124            # = 2.4 — condition for m=0 ≡ 124 BIT
X_RUN_129 = 3.0

KS = (1, 2, 4)                   # k=3: permutation id -> link
MS = (0, 1, 2, 3)


def braid_word(k: int, m: int):
    return ([(1, +1)] * 3 + [(2, -1)]) * k + [(1, +1)] * (2 * m)


def is_homogeneous(word):
    signs = {}
    for g, s in word:
        if g in signs and signs[g] != s:
            return False
        signs[g] = s
    return True


def word_perm(word):
    p = list(range(3))
    for g, _s in word:
        p[g - 1], p[g] = p[g], p[g - 1]
    return p


def soll_crossings(k: int, m: int):
    return 4 * k + 2 * m


def build_curve(k: int, m: int):
    word = braid_word(k, m)
    paths, end_lane = build_station(word, x_stn=X_STN_129, x_run=X_RUN_129)
    pts, cur, visited = [], 0, []
    for _ in range(3):
        visited.append(cur)
        pts += paths[cur]
        lane = end_lane[cur]
        pts += connector(lane, x_run=X_RUN_129)
        cur = lane
    if cur != 0 or sorted(visited) != [0, 1, 2]:
        raise RuntimeError(f"(k={k},m={m}): no 3-cycle — closure broken")
    return np.asarray(pts, float)


def make_field(k: int, m: int, out_root: Path):
    name = f"gl_k{k}m{m}"
    P = build_curve(k, m)
    n_c = station_crossings(P, x_stn=X_STN_129)
    soll = soll_crossings(k, m)
    entry = {"name": name, "machine": "129_crossing_the_ladders",
             "braid_word": f"(s1^3 s2^-1)^{k} * s1^{2 * m}",
             "grid": [k, m], "homogeneous": is_homogeneous(braid_word(k, m)),
             "fiber_genus_soll": 2 * k + m - 1,
             "station_crossings": n_c, "station_crossings_soll": soll,
             "predictions": {
                 "W0": "m=0 == 124er-Sprosse BIT; gl_k1m0 = Trefoil",
                 "W1": "fibered=true + deg Delta = 4k+2m-2 + genus = "
                       "2k+m-1 auf ALLEN Feldern",
                 "W2": "Delta-Rekurrenz in m pro Zeile k: Fit m=0..2 "
                       "sagt m=3 bit-exakt (3 Vorhersagen)",
                 "W3": "det(k,m)-Tafel mod 8 (Muster = Befund)",
                 "W5": "theta-FLAECHE + 6D trennt Gitter-Nachbarn"}}
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
                source="maschinen_ideen/129", n_terms=nt,
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
    print("[DEMO] Smoke test machine 129 (grid lab, fibration × recurrence)")
    out = HERE / "demo_smoke"
    ok = 0
    # S1: combinatorics on all 12 cells + k=3 counter-check
    s1 = True
    for k in KS:
        line = []
        for m in MS:
            n_c = station_crossings(build_curve(k, m), x_stn=X_STN_129)
            soll = soll_crossings(k, m)
            hom = is_homogeneous(braid_word(k, m))
            s1 &= (n_c == soll) and hom
            line.append(f"m{m}:{n_c}/{soll}{'✓' if n_c == soll else '✗'}")
        print(f"  (S1) k={k}: {' · '.join(line)} · homogeneous ✓ · "
              f"g_target={[2 * k + m - 1 for m in MS]}")
    link3 = word_perm(braid_word(3, 0)) == [0, 1, 2]
    s1 &= link3
    print(f"  (S1) k=3 counter-check: id permutation ⟹ link, skipped "
          f"{'✓' if link3 else '✗'}")
    ok += s1
    # S2: class surgery + m=0 row == 124 curves BIT
    outs = {(k, m): outside_station(build_curve(k, m), x_stn=X_STN_129)
            for k in KS for m in MS}
    s2 = True
    for cls in ((1, 4), (2,)):
        ref = outs[(cls[0], 0)]
        for k in cls:
            for m in MS:
                o = outs[(k, m)]
                s2 &= o.shape == ref.shape and bool(np.array_equal(o, ref))
    cross = all(bool(np.array_equal(build_curve(k, 0), build_124(k)))
                for k in KS)
    s2 = s2 and cross
    print(f"  (S2) class surgery ✓×12 cells {'✓' if s2 else '✗'} · "
          f"m=0 row == 124 curves pointwise "
          f"{'✓ CROSS-MACHINE' if cross else '✗'}")
    ok += s2
    # S3: embedding + Nyquist on all 12 cells (mandatory)
    s3 = True
    for k in KS:
        for m in MS:
            entry, folder = make_field(k, m, out)
            if folder is None:
                print(f"  (S3) {entry['name']}: REJECTED — "
                      f"{entry['rejected']} ✗")
                s3 = False
            else:
                print(f"  (S3) {entry['name']}: "
                      f"d_min_rel={entry['embedding']['d_min_rel']}"
                      f" · nyq={entry['fseries']['nyquist_dev_rel']} ✓")
    ok += s3
    print(f"[DEMO] {ok}/3 green — {'grid stands.' if ok == 3 else 'RED.'}"
          f" Bets: 129_gitter_labor.txt. Doctrine: do NOT fire.")
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
