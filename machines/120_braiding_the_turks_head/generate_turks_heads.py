# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 120 — braiding_the_turks_head (THK(3,n), the Stoimenow series).

THK(3,n) = closure of (σ₁σ₂⁻¹)ⁿ on the 119 ring architecture
(ring_curve, word-parametric). Series n ∈ {2,4,5,7,8} → c = 2n ∈
{4,8,10,14,16}; 3∤n ⟹ 1 component. Calibration anchors: THK(3,2)=4₁,
THK(3,4)=8_18, THK(3,5)=10_123; n=7/8 beyond KnotInfo, identity =
construction. Alternating + amphichiral + periodic; Tait ⟹ c EXACT.

BETS: pre-registered in 120_turks_head_schmiede.txt (W1 Tait c=2n ·
W2 identity chain · W3 amphi/sig=0 · W4 det≡1 mod 4 + two squares
+ recurrence · W5 θ/6D). Doctrine: do NOT fire.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
YARD = HERE.parent
sys.path.insert(0, str(YARD))
sys.path.insert(0, str(YARD / "119_ghosting_the_ladder"))
sys.path.insert(0, str(YARD / "117_stepping_the_delta"))
sys.path.insert(0, str(YARD / "81_clasping_the_double"))
sys.path.insert(0, str(YARD / "10_climbing_the_torus_ladder"))
from generate_ghost_ladder import ring_curve  # noqa: E402  (ring, harvest 26)
from machine_common import (  # noqa: E402
    densify, embedding_report, polyline_to_fseries,
    count_projection_crossings)

N_SERIES = (2, 4, 5, 7, 8)        # 3∤n (n≡0 mod 3 = Borromean link)


def word_th(n: int):
    """(σ₁σ₂⁻¹)ⁿ — 2n crossings, alternating braid word."""
    return [(1, +1), (2, -1)] * n


def alternation_certificate(P: np.ndarray):
    """Geometric alternation certificate: all crossing passages of the
    z-projection along the curve (parameter s = segment
    index + fraction, flag over/under from z) — Tait premise is
    MEASURED if the sequence strictly alternates."""
    Q = np.asarray(P, float)
    n = len(Q)
    A2 = Q[:, :2]
    B2 = Q[(np.arange(n) + 1) % n][:, :2]
    Az = Q[:, 2]
    Bz = Q[(np.arange(n) + 1) % n][:, 2]
    passages = []
    for i in range(n):
        j = np.arange(i + 2, n if i > 0 else n - 1)
        if len(j) == 0:
            continue
        p, r = A2[i], B2[i] - A2[i]
        q, sv = A2[j], B2[j] - A2[j]
        denom = r[0] * sv[:, 1] - r[1] * sv[:, 0]
        dq = q - p
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (dq[:, 0] * sv[:, 1] - dq[:, 1] * sv[:, 0]) / denom
            u = (dq[:, 0] * r[1] - dq[:, 1] * r[0]) / denom
        hit = (np.abs(denom) > 1e-12) & (t > 1e-9) & (t < 1 - 1e-9) \
            & (u > 1e-9) & (u < 1 - 1e-9)
        for jj, tt, uu in zip(j[hit], t[hit], u[hit]):
            zi = Az[i] + tt * (Bz[i] - Az[i])
            zj = Az[jj] + uu * (Bz[jj] - Az[jj])
            passages.append((i + tt, bool(zi > zj)))
            passages.append((jj + uu, bool(zj > zi)))
    passages.sort()
    flags = [f for _, f in passages]
    alternating = all(flags[m] != flags[m + 1] for m in range(len(flags) - 1))
    # close cyclically: last and first passage must also alternate
    if flags:
        alternating &= (flags[-1] != flags[0])
    return alternating, len(flags) // 2


def make_rung(n: int, out_root: Path):
    name = f"th_n{n}"
    P = ring_curve(word_th(n))
    n_c = count_projection_crossings(P, axis=2)
    alt, n_alt = alternation_certificate(P)
    entry = {"name": name, "machine": "120_braiding_the_turks_head",
             "braid_word": f"(s1 s2^-1)^{n}", "family": f"THK(3,{n})",
             "station_crossings": n_c, "soll": 2 * n,
             "alternating_certificate": alt,
             "predictions": {
                 "W1": "Tait: c EXAKT 2n (kein <=) — Serie heilt Obergrenze",
                 "W2": "n=2->4_1, n=4->8_18, n=5->10_123; n=7/8 kein Match",
                 "W3": "amphi (Mirror-Match) + sig=0 durchgehend",
                 "W4": "det=1 mod 4 + Summe zweier Quadrate + Rekurrenz",
                 "W5": "theta/6D ueber amphi-Serie 4..16 Kreuzungen"}}
    if n_c != 2 * n or not alt:
        entry["rejected"] = (f"Kombinatorik: {n_c}/{2*n} Kreuzungen, "
                             f"alternierend={alt}")
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
                source="maschinen_ideen/120", n_terms=nt,
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
    print("[DEMO] Smoke test machine 120 (Turk's-head forge THK(3,n))",
          flush=True)
    out = HERE / "demo_smoke"
    ok = 0
    # S1: component (ring_curve raises on non-3-cycle) + crossing target
    s1 = True
    for n in N_SERIES:
        P = ring_curve(word_th(n))
        n_c = count_projection_crossings(P, axis=2)
        good = (n_c == 2 * n)
        s1 &= good
        print(f"  (S1) THK(3,{n}): 1 component ✓ · crossings {n_c}/{2*n} "
              f"{'✓' if good else '✗'}", flush=True)
    ok += s1
    # S2: alternation certificate (cyclic, strict)
    s2 = True
    for n in N_SERIES:
        alt, n_pass = alternation_certificate(ring_curve(word_th(n)))
        s2 &= alt
        print(f"  (S2) THK(3,{n}): {n_pass} crossings, passage sequence "
              f"alternates cyclically {'✓' if alt else '✗'}", flush=True)
    ok += s2
    # S3: embedding + Nyquist (n=8 may tear at sector capacity)
    s3 = True
    for n in N_SERIES:
        entry, folder = make_rung(n, out)
        if folder is None:
            print(f"  (S3) th_n{n}: REJECTED — {entry['rejected']}"
                  f"{' (format limit measured, ok)' if n == 8 else ' ✗'}",
                  flush=True)
            if n != 8:
                s3 = False
        else:
            print(f"  (S3) th_n{n}: d_min_rel={entry['embedding']['d_min_rel']}"
                  f" · n_terms={entry['fseries']['n_terms']}"
                  f" · nyq={entry['fseries']['nyquist_dev_rel']} ✓", flush=True)
    ok += s3
    print(f"[DEMO] {ok}/3 green — {'forge stands.' if ok == 3 else 'RED.'}"
          f" Bets: 120_turks_head_schmiede.txt. Doctrine: do NOT fire.",
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
