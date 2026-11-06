# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 123 — ratcheting_the_signature (signature ratchet, Stoimenow axis).

machine-122 station (117 architecture, x_stn=4.8) + same base B=(σ₁σ₂)⁴:
W_k = σ₁^{2k}·B = twisted torus knots T(3,4;2,k) (Dean family).
Rungs rs_k0..rs_k5 (8+2k crossings exactly). rs_k0 ≡ 122/tw_base BIT
(cross-machine check). Main bet: 2nd-order Δ recurrence fitted from
k=0,1,2 predicts k=3,4,5 bit-exactly (skein at the
twist region) — pre-registered prediction instead of re-measurement.

BETS: pre-registered in 123_signatur_ratsche.txt (W0 bit anchor ·
W1 Δ recurrence prediction · W2 det sequence/mod-8 · W3 sig ratchet ·
W4 vol convergence · W5 θ/6D). Do NOT fire.
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

# Geometry EXACTLY as 122 — condition for the rs_k0≡tw_base bit check.
X_STN_123 = 4.8
X_RUN_123 = 5.4

BASE = [(1, +1), (2, +1)] * 4            # (σ₁σ₂)⁴ -> T(3,4) = 8_19


def braid_word(k: int):
    """σ₁^{2k} · (σ₁σ₂)⁴ — twisted torus knot T(3,4;2,k), 8+2k Krz."""
    return [(1, +1)] * (2 * k) + list(BASE)


def word_perm(word):
    p = list(range(3))
    for g, _s in word:
        p[g - 1], p[g] = p[g], p[g - 1]
    return p


def soll_crossings(k: int):
    return 8 + 2 * k


def build_curve(k: int):
    word = braid_word(k)
    paths, end_lane = build_station(word, x_stn=X_STN_123, x_run=X_RUN_123)
    pts, cur, visited = [], 0, []
    for _ in range(3):
        visited.append(cur)
        pts += paths[cur]
        lane = end_lane[cur]
        pts += connector(lane, x_run=X_RUN_123)
        cur = lane
    if cur != 0 or sorted(visited) != [0, 1, 2]:
        raise RuntimeError(f"k={k}: permutation not a 3-cycle — closure broken")
    return np.asarray(pts, float)


def make_rung(k: int, out_root: Path):
    name = f"rs_k{k}"
    P = build_curve(k)
    n_c = station_crossings(P, x_stn=X_STN_123)
    soll = soll_crossings(k)
    entry = {"name": name, "machine": "123_ratcheting_the_signature",
             "braid_word": f"s1^{2 * k} * (s1 s2)^4",
             "family": f"twisted torus T(3,4;2,{k})",
             "station_crossings": n_c, "station_crossings_soll": soll,
             "predictions": {
                 "W0": "rs_k0 = 8_19, fseries BIT == 122/tw_base",
                 "W1": "Delta-Rekurrenz 2. Ordnung: fit k=0..2 sagt "
                       "k=3..5 bit-exakt vorher (Skein/Twist-Region)",
                 "W2": "det_k erfuellt Rekurrenz bei t=-1; mod-8-Muster",
                 "W3": "sig-Ratsche: Sprung in {0,-2} konstant ab k>=1",
                 "W4": "vol monoton steigend, saettigend (Dehn-Fuellung)",
                 "W5": "theta pro Sprosse + 6D trennt Nachbar-Sprossen"}}
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
                source="maschinen_ideen/123", n_terms=nt,
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
    print("[DEMO] smoke test machine 123 (signature ratchet, twist dosing)")
    out = HERE / "demo_smoke"
    ks = list(range(6))
    ok = 0
    # S1: Kombinatorik + Twist-Purity
    s1 = True
    for k in ks:
        n_c = station_crossings(build_curve(k), x_stn=X_STN_123)
        soll = soll_crossings(k)
        hit = (n_c == soll)
        pure = word_perm([(1, +1)] * (2 * k)) == [0, 1, 2]
        s1 &= hit and pure
        print(f"  (S1) k={k}: crossings {n_c} (expected {soll}) "
              f"{'✓' if hit else '✗'} · Twist-Purity {'✓' if pure else '✗'}")
    ok += s1
    # S2: surgery + cross-machine anchor (rs_k0 == 122 base curve)
    outs = [outside_station(build_curve(k), x_stn=X_STN_123) for k in ks]
    s2 = all(o.shape == outs[0].shape and bool(np.array_equal(o, outs[0]))
             for o in outs[1:])
    sys.path.insert(0, str(YARD / "122_deepening_the_twins"))
    from generate_twin_wells import build_curve as build_122  # noqa: E402
    cross = bool(np.array_equal(build_curve(0), build_122(0)))
    s2 = s2 and cross
    print(f"  (S2) {len(ks)} rungs outside station bit-identical "
          f"{'✓' if s2 else '✗'} · rs_k0 == 122/tw_base curve "
          f"{'✓ CROSS-MACHINE' if cross else '✗'}")
    ok += s2
    # S3: embedding + Nyquist (mandatory k=0..4, k=5 may hit the format limit)
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
    print(f"[DEMO] {ok}/3 green — {'ratchet stands.' if ok == 3 else 'RED.'}"
          f" Bets: 123_signatur_ratsche.txt. Doctrine: do NOT fire.")
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
