# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 121 — mirroring_the_staircase (amphichiral double staircase).

V(k) = U_k · Ū_k on the FULL ring (pure closed 3-braid,
no arc sector): U_k = ([[σ₁²,σ₂²],σ₂²])^k·σ₁σ₂ (the γ₃ unit of machine 118),
Ū = sign-flipped. Window geometry depends only on the letter ⟹ the
curve is EXACTLY −id point-symmetric ⟹ strongly positive amphichiral
BY CONSTRUCTION. Theorem cross as self-tests: Stanford (a₂ const) ·
amphi⟹v₃≡0 · Hartley–Kawauchi (det = perfect square).

BETS: pre-registered in 121_amphi_doppeltreppe.txt (W0 unknot base ·
W1 a₂ const · W2 v₃≡0 · W3 det=m² · W4 degree 4 climbs, two-sided ·
W5 θ/6D). Doctrine: do NOT fire.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
YARD = HERE.parent
sys.path.insert(0, str(YARD))
sys.path.insert(0, str(YARD / "118_stacking_the_commutator"))
sys.path.insert(0, str(YARD / "117_stepping_the_delta"))
from generate_commutator_ladder import braid_word as u_word  # noqa: E402
from machine_common import (  # noqa: E402
    densify, embedding_report, polyline_to_fseries,
    count_projection_crossings)

R_RING = 2.0
S_RAD  = 0.55
H_Z    = 0.22
NW     = 25


def mirror_word(word):
    """Ū: all signs flipped (same order — point image)."""
    return [(g, -s) for g, s in word]


def staircase_word(k: int):
    """V(k) = U_k · Ū_k — 40k+4 crossings, permutation 3-cycle."""
    U = u_word(k)
    return U + mirror_word(U)


def full_ring(word) -> np.ndarray:
    """Closed 3-braid on the FULL ring (window over [0,2π]).
    −id symmetry: window j (letter g,s) lies at θ_j, window
    j+|U| (letter g,−s) at θ_j+π — lane geometry is assignment-
    independent ⟹ waypoint set point-symmetric to float epsilon."""
    n_c = len(word)
    w = 2 * np.pi / n_c
    lanes_r = (R_RING - S_RAD, R_RING, R_RING + S_RAD)
    occupant = [0, 1, 2]
    paths = {s: [] for s in range(3)}
    th0 = 0.0
    ts = np.linspace(0.0, 1.0, NW, endpoint=False)
    for gi, sign in word:
        a, b = gi - 1, gi
        c = 3 - a - b
        sa, sb, sc = occupant[a], occupant[b], occupant[c]
        for t in ts:
            th = th0 + w * t
            ramp = (1 - np.cos(np.pi * t)) / 2
            ra = lanes_r[a] + (lanes_r[b] - lanes_r[a]) * ramp
            rb = lanes_r[b] + (lanes_r[a] - lanes_r[b]) * ramp
            z = H_Z * np.sin(np.pi * t)
            ct, st = np.cos(th), np.sin(th)
            paths[sa].append((ra * ct, ra * st, +z if sign > 0 else -z))
            paths[sb].append((rb * ct, rb * st, -z if sign > 0 else +z))
            paths[sc].append((lanes_r[c] * ct, lanes_r[c] * st, 0.0))
        occupant[a], occupant[b] = sb, sa
        th0 += w
    end_lane = {occupant[lane]: lane for lane in range(3)}
    pts, cur, visited = [], 0, []
    for _ in range(3):
        visited.append(cur)
        pts += paths[cur]
        cur = end_lane[cur]                 # θ=2π ≡ 0: same lane
    if cur != 0 or sorted(visited) != [0, 1, 2]:
        raise RuntimeError("Permutation is not a 3-cycle — closure broken (link?)")
    return np.asarray(pts, float)


def point_symmetry_residual(P: np.ndarray) -> float:
    """max_i min_j |−p_i − p_j| — amphi certificate (expected ≤ 1e-12)."""
    M = -P
    worst = 0.0
    for s in range(0, len(M), 512):
        blk = M[s:s + 512]
        dd = np.linalg.norm(blk[:, None, :] - P[None, :, :], axis=2)
        worst = max(worst, float(dd.min(axis=1).max()))
    return worst


def make_rung(k: int, out_root: Path):
    name = f"ac_k{k}"
    word = staircase_word(k)
    P = full_ring(word)
    n_c = count_projection_crossings(P, axis=2)
    soll = 40 * k + 4
    res = point_symmetry_residual(P)
    entry = {"name": name, "machine": "121_mirroring_the_staircase",
             "braid_word": f"U_k * Ubar_k, U_k=([[s1^2,s2^2],s2^2])^{k}*s1s2",
             "station_crossings": n_c, "soll": soll,
             "point_symmetry_residual": res,
             "predictions": {
                 "W0": "k=0 (Kommutator-Closure) = Unknot",
                 "W1": "a2 KONSTANT (Stanford, beide Haelften gamma_3)",
                 "W2": "v3 = 0 fuer alle k (amphi toetet ungerade Grade)",
                 "W3": "det = PERFEKTES QUADRAT (Hartley-Kawauchi, spa)",
                 "W4": "Grad 4 klettert (zweiseitig: Konstanz = Befund)",
                 "W5": "theta/6D auf der spa-Serie"}}
    if n_c != soll:
        entry["rejected"] = f"Kreuzungszahl {n_c} != Soll {soll}"
        return entry, None
    if res > 1e-12:
        entry["rejected"] = f"Punktsymmetrie-Residuum {res:.2e} > 1e-12"
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
                source="maschinen_ideen/121", n_terms=nt,
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
    print("[DEMO] Smoke test machine 121 (amphichiral double staircase, U·Ū ring)",
          flush=True)
    out = HERE / "demo_smoke"
    ok = 0
    s1 = True
    for k in range(4):
        P = full_ring(staircase_word(k))
        n_c = count_projection_crossings(P, axis=2)
        good = (n_c == 40 * k + 4)
        s1 &= good
        print(f"  (S1) ac_k{k}: 1 component ✓ · crossings {n_c}/{40*k+4} "
              f"{'✓' if good else '✗'}", flush=True)
    ok += s1
    s2 = True
    for k in range(4):
        res = point_symmetry_residual(full_ring(staircase_word(k)))
        good = res <= 1e-12
        s2 &= good
        print(f"  (S2) ac_k{k}: −id residual {res:.2e} "
              f"{'✓' if good else '✗'}", flush=True)
    ok += s2
    s3 = True
    for k in range(4):
        entry, folder = make_rung(k, out)
        if folder is None:
            print(f"  (S3) ac_k{k}: REJECTED — {entry['rejected']}"
                  f"{' (ring capacity measured, ok)' if k == 3 else ' ✗'}",
                  flush=True)
            if k < 3:
                s3 = False
        else:
            print(f"  (S3) ac_k{k}: d_min_rel={entry['embedding']['d_min_rel']}"
                  f" · n_terms={entry['fseries']['n_terms']}"
                  f" · nyq={entry['fseries']['nyquist_dev_rel']} ✓", flush=True)
    ok += s3
    print(f"[DEMO] {ok}/3 green — {'double staircase stands.' if ok == 3 else 'RED.'}"
          f" Bets: 121_amphi_doppeltreppe.txt. Doctrine: do NOT fire.",
          flush=True)
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
