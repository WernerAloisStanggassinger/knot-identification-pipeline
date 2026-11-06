# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 117 — stepping_the_delta (C₂ clasper ladder, the Becerra machine).

3 co-directional lanes through a fixed station x∈[−1,1]; in it the
braid word W(k) = (σ₁σ₂⁻¹)^{3k} · σ₁σ₂ as exact geometry
(cos lane change, ±h z-lift per crossing). (σ₁σ₂⁻¹)³ = Borromean
block (mechanism hypothesis: = Y-clasper = ONE delta unit,
Murakami–Nakanishi); σ₁σ₂ = 3-cycle tail ⟹ closure always 1
component. Return via 3 z-staggered return rails.
Station fixed ⟹ rungs dm_k0..dm_k3 outside BIT-identical
(surgical ladder).

BETS: pre-registered in 117_delta_leiter.txt (W0 unknot · W1/W2
a₂ ladder ±1 exact, Okada · W3 a₄ quadratic, Goussarov–Habiro ·
W4 θ response = Becerra question · W5 6D). Doctrine: do NOT fire.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
YARD = HERE.parent
sys.path.insert(0, str(YARD))
from machine_common import (  # noqa: E402
    densify, embedding_report, polyline_to_fseries)

S_LANE = 0.55            # lane spacing in y
H_Z    = 0.22            # z-lift per crossing (pair spacing 2h in the crossing)
X_STN  = 1.0             # station half-window (FIXED over all k)
X_RUN  = 1.6             # run-in / run-out end
BULGE  = 0.6             # x-bulge of the return-rail ramps
RAILS  = (1.35, 1.80, 2.25)   # return-rail heights z (per lane)
LANES  = (+S_LANE, 0.0, -S_LANE)
NW     = 25              # waypoints per crossing window


def braid_word(k: int):
    """(σ₁σ₂⁻¹)^{3k} · σ₁σ₂ — 6k+2 crossings, permutation = 3-cycle."""
    return [(1, +1), (2, -1)] * (3 * k) + [(1, +1), (2, +1)]


def build_station(word, x_stn=X_STN, x_run=X_RUN):
    """Three strand paths through the fixed station. Strand ID = entry lane."""
    n_c = len(word)
    w = 2 * x_stn / n_c
    occupant = [0, 1, 2]                      # occupant[gasse] = strand
    paths = {s: [] for s in range(3)}
    for s in range(3):                        # run-in
        for x in np.linspace(-x_run, -x_stn, 8, endpoint=False):
            paths[s].append((x, LANES[s], 0.0))
    x0 = -x_stn
    ts = np.linspace(0.0, 1.0, NW, endpoint=False)
    for gi, sign in word:
        a, b = gi - 1, gi                     # lanes involved
        c = 3 - a - b                         # lane not involved
        sa, sb, sc = occupant[a], occupant[b], occupant[c]
        for t in ts:
            x = x0 + w * t
            ramp = (1 - np.cos(np.pi * t)) / 2
            ya = LANES[a] + (LANES[b] - LANES[a]) * ramp
            yb = LANES[b] + (LANES[a] - LANES[b]) * ramp
            z = H_Z * np.sin(np.pi * t)
            # σ⁺: the strand migrating a→b goes OVER at the top
            paths[sa].append((x, ya, +z if sign > 0 else -z))
            paths[sb].append((x, yb, -z if sign > 0 else +z))
            paths[sc].append((x, LANES[c], 0.0))
        occupant[a], occupant[b] = sb, sa
        x0 += w
    for lane in range(3):                     # run-out
        s = occupant[lane]
        for x in np.linspace(x_stn, x_run, 8, endpoint=False):
            paths[s].append((x, LANES[lane], 0.0))
    end_lane = {occupant[lane]: lane for lane in range(3)}
    return paths, end_lane


def connector(lane: int, x_run=X_RUN):
    """Return rail lane `lane`: lift to z=R, glide back, descend."""
    y, R = LANES[lane], RAILS[lane]
    pts = []
    for t in np.linspace(0, 1, 24, endpoint=False)[1:]:        # Lift
        pts.append((x_run + BULGE * np.sin(np.pi * t), y,
                    R * (1 - np.cos(np.pi * t)) / 2))
    for x in np.linspace(x_run, -x_run, 48, endpoint=False)[1:]:  # glide
        pts.append((x, y, R))
    for t in np.linspace(0, 1, 24, endpoint=False)[1:]:        # descent
        pts.append((-x_run - BULGE * np.sin(np.pi * t), y,
                    R * (1 + np.cos(np.pi * t)) / 2))
    return pts


def build_curve(k: int):
    """Closed waypoint curve of rung k (+ permutation check)."""
    word = braid_word(k)
    paths, end_lane = build_station(word)
    pts, cur, visited = [], 0, []
    for _ in range(3):
        visited.append(cur)
        pts += paths[cur]
        lane = end_lane[cur]
        pts += connector(lane)
        cur = lane                  # return rail ends in entry lane=strand
    if cur != 0 or sorted(visited) != [0, 1, 2]:
        raise RuntimeError(f"k={k}: permutation not a 3-cycle — closure broken")
    return np.asarray(pts, float)


def station_crossings(P: np.ndarray, x_stn=X_STN) -> int:
    """Crossings of the z projection, ONLY station segments (|x|≤x_stn, |z|≤0.5)."""
    n = len(P)
    nxt = (np.arange(n) + 1) % n
    inside = (np.abs(P[:, 0]) <= x_stn + 1e-9) & (np.abs(P[:, 2]) <= 0.5)
    segs = np.where(inside & inside[nxt])[0]
    A, B = P[segs][:, :2], P[nxt[segs]][:, :2]
    cnt = 0
    for ii in range(len(segs)):
        jj = np.arange(ii + 1, len(segs))
        adj = np.abs(segs[jj] - segs[ii]) <= 1
        p, r = A[ii], B[ii] - A[ii]
        q, s = A[jj], B[jj] - A[jj]
        denom = r[0] * s[:, 1] - r[1] * s[:, 0]
        dq = q - p
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (dq[:, 0] * s[:, 1] - dq[:, 1] * s[:, 0]) / denom
            u = (dq[:, 0] * r[1] - dq[:, 1] * r[0]) / denom
        hit = (~adj) & (np.abs(denom) > 1e-12) \
            & (t > 1e-9) & (t < 1 - 1e-9) & (u > 1e-9) & (u < 1 - 1e-9)
        cnt += int(hit.sum())
    return cnt


def outside_station(P: np.ndarray, x_stn=X_STN) -> np.ndarray:
    """Points outside the station window (approaches + return rails)."""
    m = (np.abs(P[:, 0]) > x_stn + 1e-9) | (np.abs(P[:, 2]) > 0.5)
    return P[m]


def make_rung(k: int, out_root: Path):
    name = f"dm_k{k}"
    P = build_curve(k)
    n_c = station_crossings(P)
    entry = {"name": name, "machine": "117_stepping_the_delta",
             "braid_word": f"(s1 s2^-1)^{3*k} * s1 s2",
             "station_crossings": n_c, "station_crossings_soll": 6 * k + 2,
             "predictions": {
                 "W0": "k=0: Unknot (det=1)",
                 "W1_W2": "a2-Leiter linear, Steigung exakt ±1 (Okada/MN)",
                 "W3": "a4 quadratisch in k (Goussarov-Habiro)",
                 "W4": "theta-Antwort pro Sprosse = Becerra-Frage",
                 "W5": "6D trennt dm_k0/dm_k1 (Grad<2-identisches Paar)"}}
    if n_c != 6 * k + 2:
        entry["rejected"] = f"Kreuzungszahl {n_c} != Soll {6 * k + 2}"
        return entry, None, P
    Pd = densify(P)
    emb = embedding_report(Pd)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None, P
    folder = out_root / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = None
    for nt in (256, 384, 512, 768, 1024):
        try:
            rep = polyline_to_fseries(
                P, folder / f"{name}.fseries", knot_name=name,
                source="maschinen_ideen/117", n_terms=nt,
                nyquist_max_rel=0.02)
            break
        except RuntimeError:
            continue
    if rep is None:
        entry["rejected"] = "Nyquist rot bis n_terms=1024"
        return entry, None, P
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder, P


def cmd_demo(_args):
    print("[DEMO] Smoke test machine 117 (delta ladder, C2 staircase)")
    out = HERE / "demo_smoke"
    ok = 0
    # S1: combinatorics (3-cycle via build_curve assert + crossing target)
    s1 = True
    for k in range(4):
        P = build_curve(k)
        n_c = station_crossings(P)
        s1 &= (n_c == 6 * k + 2)
        print(f"  (S1) k={k}: permutation 3-cycle ✓ · station crossings "
              f"{n_c} (target {6 * k + 2}) {'✓' if n_c == 6 * k + 2 else '✗'}")
    ok += s1
    # S2: surgical ladder — bit-identical outside the station
    O0, O1 = outside_station(build_curve(0)), outside_station(build_curve(1))
    s2 = O0.shape == O1.shape and bool(np.array_equal(O0, O1))
    print(f"  (S2) dm_k0/dm_k1 outside station: {O0.shape[0]} points "
          f"bit-identical {'✓' if s2 else '✗'}")
    ok += s2
    # S3: embedding gate + Nyquist-green fseries for all 4 rungs
    s3 = True
    for k in range(4):
        entry, folder, _P = make_rung(k, out)
        if folder is None:
            s3 = False
            print(f"  (S3) dm_k{k}: REJECTED — {entry['rejected']} ✗")
        else:
            print(f"  (S3) dm_k{k}: d_min_rel={entry['embedding']['d_min_rel']}"
                  f" · n_terms={entry['fseries']['n_terms']}"
                  f" · nyq={entry['fseries']['nyquist_dev_rel']} ✓")
    ok += s3
    print(f"[DEMO] {ok}/3 green — {'ladder stands.' if ok == 3 else 'RED.'}"
          f" Bets: 117_delta_leiter.txt. Doctrine: do NOT fire.")
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
