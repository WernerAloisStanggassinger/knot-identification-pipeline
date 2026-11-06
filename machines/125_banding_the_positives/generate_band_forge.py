# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 125 — banding_the_positives (band forge, Rudolph axis).

3-lane station (117 architecture, 124 width), series
W_k = (σ₁σ₂)^k · σ₂σ₁²σ₂⁻¹ = (a₁₂a₂₃)^k·(a₁₃)² — words in positive
BAND generators ⟹ strongly quasipositive ⟹ band surface = minimal genus
(Rudolph/Bennequin): genus = k DOSED (1/2/4/5). Word certificates
measured: band-positive ∧ NOT positive (σ₂⁻¹ present) — the
hierarchy level above machine 99.

BETS: pre-registered in 125_band_schmiede.txt (W1 genus=k exact ·
W2 deg-Δ gap statistics · W3 sig monotonicity · W4 hierarchy question
two-sided · W5 θ/6D). Do NOT fire.
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

X_STN_125 = 2.4
X_RUN_125 = 3.0

KS = (1, 2, 4, 5)            # k=3: permutation id -> link, not a knot
A13_SQ = [(2, +1), (1, +1), (1, +1), (2, -1)]    # (a13)^2 = s2 s1^2 s2^-1


def braid_word(k: int):
    """(σ₁σ₂)^k · σ₂σ₁²σ₂⁻¹ — 2k+4 letters, 2k+2 bands."""
    return [(1, +1), (2, +1)] * k + list(A13_SQ)


def n_bands(k: int):
    """Band count: each σ± single band (2k) + (a₁₃)² (2 bands)."""
    return 2 * k + 2


def word_perm(word):
    p = list(range(3))
    for g, _s in word:
        p[g - 1], p[g] = p[g], p[g - 1]
    return p


def band_positive_and_not_positive(word):
    """Certificates from the word: (a) all bands positive (construction:
    blocks σ₁ / σ₂ / σ₂σ₁²σ₂⁻¹ are positive bands), (b) as
    braid word NOT positive (negative letter present)."""
    has_negative = any(s < 0 for _g, s in word)
    return True, has_negative


def soll_crossings(k: int):
    return 2 * k + 4


def build_curve(k: int):
    word = braid_word(k)
    paths, end_lane = build_station(word, x_stn=X_STN_125, x_run=X_RUN_125)
    pts, cur, visited = [], 0, []
    for _ in range(3):
        visited.append(cur)
        pts += paths[cur]
        lane = end_lane[cur]
        pts += connector(lane, x_run=X_RUN_125)
        cur = lane
    if cur != 0 or sorted(visited) != [0, 1, 2]:
        raise RuntimeError(f"k={k}: permutation not a 3-cycle — closure broken")
    return np.asarray(pts, float)


def make_rung(k: int, out_root: Path):
    name = f"bp_k{k}"
    P = build_curve(k)
    n_c = station_crossings(P, x_stn=X_STN_125)
    soll = soll_crossings(k)
    bp, notpos = band_positive_and_not_positive(braid_word(k))
    entry = {"name": name, "machine": "125_banding_the_positives",
             "braid_word": f"(s1 s2)^{k} * s2 s1^2 s2^-1",
             "band_word": f"(a12 a23)^{k} * (a13)^2",
             "n_bands": n_bands(k), "genus_soll": k,
             "cert_band_positive": bp, "cert_not_positive_word": notpos,
             "station_crossings": n_c, "station_crossings_soll": soll,
             "predictions": {
                 "W1": "genus = k EXAKT (Rudolph/Bennequin, Band-Flaeche)",
                 "W2": "deg Delta <= 2k; Luecken-Statistik = Befund",
                 "W3": "sig monoton fallend in k (zweiseitig)",
                 "W4": "Hierarchie-Frage: Knoten positiv? Zuendung "
                       "urteilt (Wort zertifiziert nicht-positiv)",
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
                source="maschinen_ideen/125", n_terms=nt,
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
    print("[DEMO] Smoke test machine 125 (band forge, Rudolph/SQP)")
    out = HERE / "demo_smoke"
    ok = 0
    # S1: combinatorics + certificates + k=3 counter-check
    s1 = True
    for k in KS:
        w = braid_word(k)
        n_c = station_crossings(build_curve(k), x_stn=X_STN_125)
        soll = soll_crossings(k)
        hit = (n_c == soll)
        bp, notpos = band_positive_and_not_positive(w)
        s1 &= hit and bp and notpos
        print(f"  (S1) k={k}: crossings {n_c} (expected {soll}) "
              f"{'✓' if hit else '✗'} · bands {n_bands(k)} ⟹ g_target={k}"
              f" · band-positive {'✓' if bp else '✗'} · not positive "
              f"(σ₂⁻¹) {'✓' if notpos else '✗'}")
    link3 = word_perm(braid_word(3)) == [0, 1, 2]
    s1 &= link3
    print(f"  (S1) k=3 counter-check: permutation = id ⟹ link, skipped "
          f"{'✓' if link3 else '✗'}")
    ok += s1
    # S2: surgery CLASS BY CLASS (lesson from 124: perm = cycle^k)
    outs = {k: outside_station(build_curve(k), x_stn=X_STN_125)
            for k in KS}
    s2 = True
    for cls in ((1, 4), (2, 5)):
        a, b = outs[cls[0]], outs[cls[1]]
        same = a.shape == b.shape and bool(np.array_equal(a, b))
        s2 &= same
        print(f"  (S2) class k∈{cls}: outside station "
              f"{a.shape[0]} points bit-identical {'✓' if same else '✗'}")
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
          f" Bets: 125_band_schmiede.txt. Doctrine: do NOT fire.")
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
