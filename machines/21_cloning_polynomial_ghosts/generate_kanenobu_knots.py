# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 21 — cloning_polynomial_ghosts: KANENOBU CONSTRUCTION (unlocked).

LIT CHECK 2026-06-07 (Lobb, arXiv:1105.3985 / PNAS version, Fig. 1 +
Thm 1.1 read): K(p,q) = band diagram through infinity, |p| and
|q| half-twists in the two strand pairs (top/bottom), sign by
sign; Kanenobu: P(K(p,q)) = P(K(r,s)) ⟺-sided if p+q = r+s
(HOMFLY ghosts!), K(p,q)‾ = K(−p,−q), K(p,q) = K(q,p); ribbon.
DESIGN here: circle core, band edges ±ρ·u(s), twist boxes at π/2
(p·π) and 3π/2 (q·π), junction π = CROSS hook (two cross caps,
radially staggered = the ribbon hook), junction 0 = straight caps.
TODAY'S BIT TESTS (literature anchors as geometry!): hook-cap
clearance (measured); p↔q via R_x=diag(1,−1,−1) bit-exact in the EVEN
sector (p+q even = exactly Lobb's Thm-1.1 territory; ψ matching mod π
proven, cap swap +z↔−z under R_x). MIRROR theorem K(p,q)‾ =
K(−p,−q) is rigidly IMPOSSIBLE (hook chirality flips under mirror,
compensation = isotopy, not geometry) → ignition bet. IGNITION:
K(0,0) → 4₁#4₁? (crosscheck target) + ghost bet P = f(p+q) +
mirror bet. Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402

R0 = 3.0          # circle core
RHO = 0.45        # band half-width
GAP = 0.22        # junction gap (rad)
N_ARC = 900


def smoothstep(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3 - 2 * x)


def psi_profile(s, p, q, width=0.9):
    """Band rotation angle: p·π in the box around π/2, q·π around 3π/2."""
    ps = p * np.pi * smoothstep((s - (np.pi / 2 - width / 2)) / width)
    qs = q * np.pi * smoothstep((s - (3 * np.pi / 2 - width / 2)) / width)
    return ps + qs


def edge(s, sign, p, q):
    """Band edge at ±ρ·u(s); u rotates by ψ(s) (e1=radial, e2=z)."""
    psi = psi_profile(s, p, q)
    er = np.column_stack([np.cos(s), np.sin(s), np.zeros_like(s)])
    ez = np.array([0.0, 0.0, 1.0])
    u = (np.cos(psi)[:, None] * er
         + np.sin(psi)[:, None] * ez[None, :])
    core = R0 * er
    return core + sign * RHO * u


def kanenobu_curve(p, q):
    """Traversal: A_top → cross → B_bot → straight → B_top → cross →
    A_bot → straight → closes. Cross caps radially staggered (hook)."""
    s_top = np.linspace(GAP, np.pi - GAP, N_ARC)
    s_bot = np.linspace(np.pi + GAP, 2 * np.pi - GAP, N_ARC)
    A_top = edge(s_top, +1, p, q)
    B_top = edge(s_top, -1, p, q)
    A_bot = edge(s_bot, +1, p, q)
    B_bot = edge(s_bot, -1, p, q)

    def cap(P_from, P_to, bulge_r, bulge_z=0.0, n=40):
        t = np.linspace(0, 1, n + 2)[1:-1]
        base = np.outer(1 - t, P_from) + np.outer(t, P_to)
        prof = np.sin(np.pi * t)
        # radial belly (in the xy plane outward/inward)
        mid_dir = base[:, :2] / np.linalg.norm(base[:, :2],
                                               axis=1)[:, None]
        base[:, :2] += bulge_r * prof[:, None] * mid_dir
        base[:, 2] += bulge_z * prof          # z stagger (hook!)
        return base

    # cross caps ±z-staggered (build lesson: radial stagger collides
    # at t≈0.8 — radius curves intersect in the z=0 plane).
    # Under R_x, +z↔−z swaps exactly the two hook caps.
    parts = [A_top,
             cap(A_top[-1], B_bot[0], +0.10, +0.55),  # hook top
             B_bot,
             cap(B_bot[-1], B_top[0], +0.18),         # straight outward
             B_top,
             cap(B_top[-1], A_bot[0], +0.10, -0.55),  # hook bottom
             A_bot,
             cap(A_bot[-1], A_top[0], -0.18)]         # straight inward
    C = np.vstack(parts)
    for _ in range(2):
        C = 0.5 * C + 0.25 * (np.roll(C, 1, axis=0)
                              + np.roll(C, -1, axis=0))
    return C


def set_distance(P, Q, stride=64):
    d = 0.0
    for x in Q[::stride]:
        d = max(d, float(np.linalg.norm(P - x, axis=1).min()))
    diam = float(np.linalg.norm(P.max(0) - P.min(0)))
    return d / diam


def make_knot(p, q, out_root):
    C = kanenobu_curve(p, q)
    name = f"kb_p{p}_q{q}".replace("-", "m")
    entry = {"name": name, "machine": "21_cloning_polynomial_ghosts",
             "p": p, "q": q,
             "predictions": {
                 "geister": f"P hängt nur von p+q={p+q} ab (Kanenobu, "
                            "in-house prüfbar)",
                 "ziel": "K(0,0) → 4_1 # 4_1 (Crosscheck-Ziel)",
                 "ribbon": "Ribbon-Knoten (Lobb Thm 1.4-Kontext)"}}
    emb = embedding_report(C)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = None
    for nt in (128, 192, 256, 384):
        try:
            rep = polyline_to_fseries(C, folder / f"{name}.fseries",
                                      knot_name=name,
                                      source="maschinen_ideen/21",
                                      n_terms=nt, nyquist_max_rel=0.02)
            break
        except RuntimeError:
            continue
    if rep is None:
        entry["rejected"] = "Nyquist rot"
        return entry, None
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 21 (Kanenobu, unlocked by lit check)")
    # (1) hook clearance: the two cross caps (ribbon hook) must not
    #     come close anywhere (build lesson: radial stagger did)
    C = kanenobu_curve(2, 0)
    i0 = N_ARC                       # cap1 starts after A_top
    cap1 = C[i0:i0 + 40]
    i1 = 2 * N_ARC + 80              # cap2 after B_bot + straight cap
    cap2 = C[i1 + N_ARC:i1 + N_ARC + 40]
    dmin = min(float(np.linalg.norm(cap1 - x, axis=1).min())
               for x in cap2)
    ok1 = dmin > 0.4
    print(f"  (1) cross-cap clearance (z stagger): min distance "
          f"{dmin:.3f}  {'✓ hook clean' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) p↔q via R_x=diag(1,−1,−1), EVEN sector (Lobb Thm 1.1):
    #     R_x(K(q,p)) == K(p,q) bit-exact as a point set
    Rx = np.diag([1.0, -1.0, -1.0])
    worst = 0.0
    for (p, q) in ((0, 2), (4, 0), (3, 1)):
        d = set_distance(kanenobu_curve(p, q),
                         kanenobu_curve(q, p) @ Rx.T)
        worst = max(worst, d)
    ok2 = worst < 1e-6
    print(f"  (2) R_x(K(q,p)) vs K(p,q) for (0,2),(4,0),(3,1): "
          f"max set distance {worst:.2e}  "
          f"{'✓ p↔q bit geometry' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) ghost family p+q=2 built: (0,2), (1,1), (2,0) + K(0,0)
    built = []
    for (p, q) in ((0, 0), (0, 2), (1, 1), (2, 0), (4, 0), (2, 2)):
        e, f = make_knot(p, q, HERE / "demo_smoke")
        if f is not None:
            built.append((p, q, e["fseries"]["nyquist_dev_rel"]))
        else:
            print(f"      K({p},{q}): ✗ {e.get('rejected')}")
    ok3 = len(built) >= 5
    print(f"  (3) family: {[(p, q) for p, q, _ in built]}  "
          f"{'✓ ghost ammunition in place (incl. K(0,0) target)' if ok3 else '✗ FAIL'}")
    assert ok3
    print("[DEMO] 3/3 green — Kanenobu machine BUILT (lit gate "
          "lifted; ignition: K(0,0)→4₁#4₁ + P=f(p+q)).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
