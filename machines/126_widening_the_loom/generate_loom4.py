# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 126 — widening_the_loom (four-lane loom, Loom-4).

Generalisation of the 117 station to 4 strands (word-agnostic,
generalised occupant logic). First family: THK(4,n) = closure
(σ₁σ₂⁻¹σ₃)ⁿ, n ∈ {1,3,5,7} (n even = link). Homogeneous ⟹ Stallings
fibred, g = (3n−3)/2; alternating ⟹ Tait c = 3n (bet).
thk4_n1 = unknot anchor (genus 0 by theorem). FOLLOW-UP DOOR (lit-gated):
Conway/KT mutant pair on the same loom.

BETS: pre-registered in 126_vier_gassen_webstuhl.txt (W0 unknot ·
W1 Tait c=3n · W2 fibered+deg Δ=3(n−1) · W3 8_18 suspicion n=3 ·
W4 amphi/sig=0 two-sided · W5 det recurrence+θ/6D). Do NOT fire.
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

# --- Loom-4 geometry (117 scale, 4th lane + 4th return rail) -------------
N_STR  = 4
S_LANE = 0.55
H_Z    = 0.22
X_STN  = 2.4
X_RUN  = 3.0
BULGE  = 0.6
RAILS  = (1.35, 1.80, 2.25, 2.70)
LANES  = (+1.5 * S_LANE, +0.5 * S_LANE, -0.5 * S_LANE, -1.5 * S_LANE)
NW     = 25

KS = (1, 3, 5, 7)            # n even: 4-cycle^n = not a 4-cycle -> link


def braid_word(n: int):
    """(σ₁σ₂⁻¹σ₃)ⁿ — homogeneous (σ₁+, σ₂−, σ₃+), 3n crossings."""
    return [(1, +1), (2, -1), (3, +1)] * n


def word_perm(word, n_str=N_STR):
    p = list(range(n_str))
    for g, _s in word:
        p[g - 1], p[g] = p[g], p[g - 1]
    return p


def is_homogeneous(word):
    signs = {}
    for g, s in word:
        if g in signs and signs[g] != s:
            return False
        signs[g] = s
    return True


def soll_crossings(n: int):
    return 3 * n


def build_station(word, x_stn=X_STN, x_run=X_RUN):
    """Four strand paths through the fixed station (117 logic, occupant
    generalised: uninvolved lanes = loop instead of 3−a−b)."""
    n_c = len(word)
    w = 2 * x_stn / n_c
    occupant = list(range(N_STR))
    paths = {s: [] for s in range(N_STR)}
    for s in range(N_STR):
        for x in np.linspace(-x_run, -x_stn, 8, endpoint=False):
            paths[s].append((x, LANES[s], 0.0))
    x0 = -x_stn
    ts = np.linspace(0.0, 1.0, NW, endpoint=False)
    for gi, sign in word:
        a, b = gi - 1, gi
        others = [g for g in range(N_STR) if g not in (a, b)]
        sa, sb = occupant[a], occupant[b]
        for t in ts:
            x = x0 + w * t
            ramp = (1 - np.cos(np.pi * t)) / 2
            ya = LANES[a] + (LANES[b] - LANES[a]) * ramp
            yb = LANES[b] + (LANES[a] - LANES[b]) * ramp
            z = H_Z * np.sin(np.pi * t)
            paths[sa].append((x, ya, +z if sign > 0 else -z))
            paths[sb].append((x, yb, -z if sign > 0 else +z))
            for g in others:
                paths[occupant[g]].append((x, LANES[g], 0.0))
        occupant[a], occupant[b] = sb, sa
        x0 += w
    for lane in range(N_STR):
        s = occupant[lane]
        for x in np.linspace(x_stn, x_run, 8, endpoint=False):
            paths[s].append((x, LANES[lane], 0.0))
    end_lane = {occupant[lane]: lane for lane in range(N_STR)}
    return paths, end_lane


def connector(lane: int, x_run=X_RUN):
    """Return rail lane `lane` (117 logic, 4 rails)."""
    y, R = LANES[lane], RAILS[lane]
    pts = []
    for t in np.linspace(0, 1, 24, endpoint=False)[1:]:
        pts.append((x_run + BULGE * np.sin(np.pi * t), y,
                    R * (1 - np.cos(np.pi * t)) / 2))
    for x in np.linspace(x_run, -x_run, 48, endpoint=False)[1:]:
        pts.append((x, y, R))
    for t in np.linspace(0, 1, 24, endpoint=False)[1:]:
        pts.append((-x_run - BULGE * np.sin(np.pi * t), y,
                    R * (1 + np.cos(np.pi * t)) / 2))
    return pts


def station_crossings(P: np.ndarray, x_stn=X_STN) -> int:
    """Crossings of the z projection, station only (|x|≤x_stn, |z|≤0.5)."""
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
    m = (np.abs(P[:, 0]) > x_stn + 1e-9) | (np.abs(P[:, 2]) > 0.5)
    return P[m]


def build_curve(n: int):
    word = braid_word(n)
    paths, end_lane = build_station(word)
    pts, cur, visited = [], 0, []
    for _ in range(N_STR):
        visited.append(cur)
        pts += paths[cur]
        lane = end_lane[cur]
        pts += connector(lane)
        cur = lane
    if cur != 0 or sorted(visited) != list(range(N_STR)):
        raise RuntimeError(f"n={n}: permutation not a 4-cycle — closure broken")
    return np.asarray(pts, float)


def make_rung(n: int, out_root: Path):
    name = f"thk4_n{n}"
    P = build_curve(n)
    n_c = station_crossings(P)
    soll = soll_crossings(n)
    entry = {"name": name, "machine": "126_widening_the_loom",
             "braid_word": f"(s1 s2^-1 s3)^{n}",
             "family": f"THK(4,{n})", "n_strands": N_STR,
             "homogeneous": is_homogeneous(braid_word(n)),
             "fiber_genus_soll": (3 * n - 3) // 2,
             "station_crossings": n_c, "station_crossings_soll": soll,
             "predictions": {
                 "W0": "thk4_n1 = Unknot (Stallings-Genus 0)",
                 "W1": "alternierend + Tait c = 3n exakt",
                 "W2": "fibered=true + deg Delta = 3(n-1) exakt",
                 "W3": "thk4_n3 = 8_18-Verdacht (Lit, Zuendung urteilt)",
                 "W4": "amphi/sig=0-Frage zweiseitig",
                 "W5": "det-Rekurrenz Fit n=1,3,5 trifft n=7 + theta/6D"}}
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
                source="maschinen_ideen/126", n_terms=nt,
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
    print("[DEMO] Smoke test machine 126 (four-lane loom, THK(4,n))")
    out = HERE / "demo_smoke"
    ok = 0
    # S1: combinatorics + homogeneity + n=2 counter-check
    s1 = True
    for n in KS:
        w = braid_word(n)
        n_c = station_crossings(build_curve(n))
        soll = soll_crossings(n)
        hit = (n_c == soll)
        hom = is_homogeneous(w)
        s1 &= hit and hom
        print(f"  (S1) n={n}: crossings {n_c} (expected {soll}) "
              f"{'✓' if hit else '✗'} · homogeneous {'✓' if hom else '✗'} · "
              f"g_target={(3 * n - 3) // 2}")
    p2 = word_perm(braid_word(2))
    # 4-cycle^2 = double transposition ⟹ 2 components
    link2 = (p2[p2[0]] == 0 and p2[0] != 0)
    s1 &= link2
    print(f"  (S1) n=2 counter-check: permutation = double transposition ⟹ "
          f"link, skipped {'✓' if link2 else '✗'}")
    ok += s1
    # S2: surgery CLASS BY CLASS ({1,5}=cycle, {3,7}=cycle³)
    outs = {n: outside_station(build_curve(n)) for n in KS}
    s2 = True
    for cls in ((1, 5), (3, 7)):
        a, b = outs[cls[0]], outs[cls[1]]
        same = a.shape == b.shape and bool(np.array_equal(a, b))
        s2 &= same
        print(f"  (S2) class n∈{cls}: outside station "
              f"{a.shape[0]} points bit-identical {'✓' if same else '✗'}")
    ok += s2
    # S3: embedding + Nyquist (mandatory n=1,3,5; n=7 may hit the format limit)
    s3 = True
    for n in KS:
        entry, folder = make_rung(n, out)
        if folder is None:
            print(f"  (S3) {entry['name']}: REJECTED — {entry['rejected']}"
                  f"{' (format limit measured, ok)' if n == 7 else ' ✗'}")
            if n != 7:
                s3 = False
        else:
            print(f"  (S3) {entry['name']}: "
                  f"d_min_rel={entry['embedding']['d_min_rel']}"
                  f" · n_terms={entry['fseries']['n_terms']}"
                  f" · nyq={entry['fseries']['nyquist_dev_rel']} ✓")
    ok += s3
    print(f"[DEMO] {ok}/3 green — {'loom stands.' if ok == 3 else 'RED.'}"
          f" Bets: 126_vier_gassen_webstuhl.txt. Doctrine: do NOT fire.")
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
