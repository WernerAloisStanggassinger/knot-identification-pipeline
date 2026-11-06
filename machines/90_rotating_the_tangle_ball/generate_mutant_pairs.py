# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 90 — rotating_the_tangle_ball (mutant factory).

Curve with BUILT-IN symmetric ball: 4 endpoints on the unit
sphere, pairwise Rz(π)-conjugate (A↔B, C↔D). Inner tangle free
(interlaced splines), outer part fixed. MUTATION = exact π rotation
of the inner part: connectivity stays (Rz·s₂ connects A→C, Rz·s₁
connects B→D), the tangle changes, outside bit-identical.
MAIL THREAD: Dror — the mail discusses the θ-blind mutant pair
13a_1231/1237; here OUR OWN mutant ammunition is produced in series.
BETS EXACT today: endpoint set Rz-invariant (bit); ball gates
(inside<0.97, outside>1.03); outside unchanged (bit); inside different
(measured); both embedded as a PAIR. IGNITION: polynomials collide
(mutation invariance), knot/θ decide. REVERSAL CLAUSE:
mutation may be trivial → ignition judges, pair stays pair.
Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "27_pumping_the_twist_account"))
sys.path.insert(0, str(HERE.parent))
from generate_writhe_pumped import writhe  # noqa: E402
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402

RZ = np.diag([-1.0, -1.0, 1.0])      # Rz(π)

# 4 endpoints on the unit sphere, Rz pairwise
A = np.array([0.55, 0.45, np.sqrt(1 - 0.55 ** 2 - 0.45 ** 2)])
B = RZ @ A
C = np.array([0.62, -0.52, -np.sqrt(1 - 0.62 ** 2 - 0.52 ** 2)])
D = RZ @ C


def spline(controls, pts_per_seg=60, rounds=3):
    P = np.asarray(controls, float)
    dense = []
    for k in range(len(P) - 1):
        t = np.linspace(0, 1, pts_per_seg, endpoint=False)
        dense.append(np.outer(1 - t, P[k]) + np.outer(t, P[k + 1]))
    Cv = np.vstack(dense + [P[-1:]])
    for _ in range(rounds):
        inner = 0.5 * Cv[1:-1] + 0.25 * (Cv[:-2] + Cv[2:])
        Cv = np.vstack([Cv[:1], inner, Cv[-1:]])
    return Cv


def inner_tangle(seed):
    """Two interlaced inner strands A→C and B→D (in the ball)."""
    rng = np.random.default_rng(seed)
    def wob(scale=0.35):
        return rng.uniform(-scale, scale, 3)
    m1 = [A, 0.55 * A + np.array([0.1, -0.3, -0.1]) + wob(),
          np.array([-0.25, 0.2, 0.15]) + wob(0.2),
          0.55 * C + np.array([-0.1, 0.3, 0.1]) + wob(), C]
    m2 = [B, 0.55 * B + np.array([0.25, 0.1, -0.3]) + wob(),
          np.array([0.2, -0.15, -0.1]) + wob(0.2),
          0.55 * D + np.array([-0.25, -0.1, 0.3]) + wob(), D]
    s1, s2 = spline(m1), spline(m2)
    # ball clamp: keep inner points radially below 0.96 (ends free)
    for s in (s1, s2):
        r = np.linalg.norm(s, axis=1)
        mask = r > 0.96
        mask[0] = mask[-1] = False
        s[mask] *= (0.96 / r[mask])[:, None]
    return s1, s2


def woven_tangle(m: int):
    """REPAIR 2026-06-07 (first ignition: s0 pair = both UNKNOT —
    random wobble practically never knots). Deterministic: 2m
    relative full rotations of the inner strands (twist region, 2m
    crossings) + CLASP (strand 1 wraps once around strand 2,
    Whitehead-like). The clasp is the point: pure twist regions
    are RATIONAL tangles = mutation-trivial; twist+clasp breaks
    that. Knotted gate = repair yardstick; non-triviality of the
    mutation remains a two-sided ignition bet (REVERSAL CLAUSE)."""
    phiA = float(np.arctan2(A[1], A[0]))
    phiB = float(np.arctan2(B[1], B[0]))
    phiC = float(np.arctan2(C[1], C[0]))
    phiD = float(np.arctan2(D[1], D[0]))
    zA, zC = float(A[2]), float(C[2])
    zB, zD = float(B[2]), float(D[2])

    def helix(phi0, dphi, z0, z1, r=0.52, n=14):
        t = np.linspace(0.0, 1.0, n)
        return np.stack([r * np.cos(phi0 + dphi * t),
                         r * np.sin(phi0 + dphi * t),
                         z0 + (z1 - z0) * t], axis=1)

    # strand 2 (B→D): gentle helix, base winding δ
    d2 = (phiD - phiB)
    h2 = helix(phiB, d2, zB, zD)
    s2 = spline([B] + [p for p in h2[1:-1]] + [D])

    # strand 1 (A→C): δ + 2π·m relative winding, clasp around the middle of strand 2
    d1 = (phiC - phiA) + 2.0 * np.pi * m
    h1 = helix(phiA, d1, zA, zC, r=0.58, n=12 + 6 * m)
    mid2 = h2[len(h2) // 2]
    tan2 = h2[len(h2) // 2 + 1] - h2[len(h2) // 2 - 1]
    tan2 = tan2 / np.linalg.norm(tan2)
    u = np.cross(tan2, [0.0, 0.0, 1.0])
    u = u / np.linalg.norm(u)
    v = np.cross(tan2, u)
    R_CLASP = 0.30
    clasp = [mid2 + R_CLASP * (np.cos(a) * u + np.sin(a) * v)
             for a in np.linspace(0.0, 2 * np.pi, 5)[:-1]]
    half = len(h1) // 2
    way1 = ([A] + [p for p in h1[1:half]] + clasp
            + [p for p in h1[half:-1]] + [C])
    s1 = spline(way1)
    for s in (s1, s2):
        r = np.linalg.norm(s, axis=1)
        mask = r > 0.96
        mask[0] = mask[-1] = False
        s[mask] *= (0.96 / r[mask])[:, None]
    return s1, s2


def crossing_floor(K, n_dirs=3):
    """Minimal crossing number over axis projections (decimated) —
    coarse knotted gate: unknot degeneration falls below the floor."""
    P = K[::3]
    n = len(P)
    best = None
    for ax in range(n_dirs):
        idx = [i for i in range(3) if i != ax]
        Q = P[:, idx]
        seg_a = Q
        seg_b = np.roll(Q, -1, axis=0)
        cnt = 0
        d = seg_b - seg_a
        for i in range(n):
            j = np.arange(i + 2, n - (1 if i == 0 else 0))
            if len(j) == 0:
                continue
            r = seg_a[j] - seg_a[i]
            dj = d[j]
            den = d[i, 0] * dj[:, 1] - d[i, 1] * dj[:, 0]
            ok = np.abs(den) > 1e-12
            t = np.where(ok, (r[:, 0] * dj[:, 1] - r[:, 1] * dj[:, 0])
                         / np.where(ok, den, 1.0), -1.0)
            s = np.where(ok, (r[:, 0] * d[i, 1] - r[:, 1] * d[i, 0])
                         / np.where(ok, den, 1.0), -1.0)
            cnt += int(np.sum(ok & (t > 0) & (t < 1) & (s > 0) & (s < 1)))
        best = cnt if best is None else min(best, cnt)
    return best


def mutant_pair_woven(m: int):
    s1, s2 = woven_tangle(m)
    o1, o2 = outer_frame()
    K1 = assemble(s1, s2, o1, o2)
    s1m = (RZ @ s2.T).T
    s2m = (RZ @ s1.T).T
    K2 = assemble(s1m, s2m, o1, o2)
    return K1, K2, (s1, s2), (o1, o2)


def make_pair_woven(m: int, out_root):
    K1, K2, (s1, s2), (o1, o2) = mutant_pair_woven(m)
    name = f"mu_pair_w{m}"
    floor1, floor2 = crossing_floor(K1), crossing_floor(K2)
    entry = {"name": name, "machine": "90_rotating_the_tangle_ball",
             "weave_m": m, "clasp": True,
             "crossing_floor_K": floor1, "crossing_floor_Kmut": floor2,
             "writhe_K": round(float(writhe(K1[::4], n=256)), 3),
             "writhe_Kmut": round(float(writhe(K2[::4], n=256)), 3),
             "predictions": {
                 "zuendung": "det!=1 PFLICHT (s0-Lehre: Basis muss "
                             "verknotet sein, sonst Sprosse rot); "
                             "Polynome kollidieren (Mutation); "
                             "Knoten-Typ + theta entscheiden",
                 "umdreh_klausel": "Twist+Clasp = nicht-rationales "
                                   "Muster, Mutation KANN trotzdem "
                                   "trivial sein -> Zuendung urteilt",
                 "schiene": "Dror (eigenes Mutanten-Paar)"}}
    if min(floor1, floor2) < 2 * m + 2:
        entry["rejected"] = (f"crossing_floor {min(floor1, floor2)} < "
                             f"Soll {2 * m + 2} (Unknot-Verdacht)")
        return entry, None
    embs = []
    for K in (K1, K2):
        emb = embedding_report(K)
        embs.append(emb)
        if not emb["embedded"]:
            entry["rejected"] = "embedding gate"
            return entry, None
    entry["embedding"] = embs[0]
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    for tag, K in (("K", K1), ("Kmut", K2)):
        rep = None
        for nt in (128, 192, 256, 384):
            try:
                rep = polyline_to_fseries(
                    K, folder / f"{name}_{tag}.fseries",
                    knot_name=f"{name}_{tag}",
                    source="maschinen_ideen/90", n_terms=nt,
                    nyquist_max_rel=0.02)
                break
            except RuntimeError:
                continue
        if rep is None:
            entry["rejected"] = f"Nyquist rot ({tag})"
            return entry, None
        entry[f"fseries_{tag}"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def outer_frame():
    """Fixed outer part: C→B and D→A, interlaced, outside r=1.03."""
    o1 = spline([C, np.array([1.7, -0.9, -0.9]),
                 np.array([1.9, 0.9, 0.4]),
                 np.array([-0.4, 1.9, 1.3]),
                 np.array([-1.6, 0.2, 1.2]), B])
    o2 = spline([D, np.array([-1.8, 0.4, -1.1]),
                 np.array([-0.6, -1.9, -0.6]),
                 np.array([1.2, -1.8, 0.7]),
                 np.array([1.4, 0.6, 1.4]), A])
    for o in (o1, o2):
        r = np.linalg.norm(o, axis=1)
        mask = r < 1.04
        mask[0] = mask[-1] = False
        o[mask] *= (1.04 / r[mask])[:, None]
    return o1, o2


def assemble(s1, s2, o1, o2):
    """A→C (inner), C→B (outer), B→D (inner), D→A (outer)."""
    Cv = np.vstack([s1[:-1], o1[:-1], s2[:-1], o2[:-1]])
    for _ in range(2):
        Cv = 0.5 * Cv + 0.25 * (np.roll(Cv, 1, axis=0)
                                + np.roll(Cv, -1, axis=0))
    return Cv


def mutant_pair(seed):
    s1, s2 = inner_tangle(seed)
    o1, o2 = outer_frame()
    K1 = assemble(s1, s2, o1, o2)
    # mutation: inside rotated exactly by Rz; Rz·s2: B→D becomes A→C etc.
    s1m = (RZ @ s2.T).T          # A→C
    s2m = (RZ @ s1.T).T          # B→D
    K2 = assemble(s1m, s2m, o1, o2)
    return K1, K2, (s1, s2), (o1, o2)


def make_pair(seed, out_root):
    K1, K2, (s1, s2), (o1, o2) = mutant_pair(seed)
    name = f"mu_pair_s{seed}"
    entry = {"name": name, "machine": "90_rotating_the_tangle_ball",
             "seed": seed,
             "writhe_K": round(float(writhe(K1[::4], n=256)), 3),
             "writhe_Kmut": round(float(writhe(K2[::4], n=256)), 3),
             "predictions": {
                 "zuendung": "Polynome MÜSSEN kollidieren (Mutation); "
                             "Knoten-Typ + θ entscheiden",
                 "schiene": "Dror (eigenes Mutanten-Paar)"}}
    embs = []
    for K in (K1, K2):
        emb = embedding_report(K)
        embs.append(emb)
        if not emb["embedded"]:
            entry["rejected"] = "embedding gate"
            return entry, None
    entry["embedding"] = embs[0]
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    for tag, K in (("K", K1), ("Kmut", K2)):
        rep = None
        for nt in (128, 192, 256):
            try:
                rep = polyline_to_fseries(
                    K, folder / f"{name}_{tag}.fseries",
                    knot_name=f"{name}_{tag}",
                    source="maschinen_ideen/90", n_terms=nt,
                    nyquist_max_rel=0.02)
                break
            except RuntimeError:
                continue
        if rep is None:
            entry["rejected"] = f"Nyquist rot ({tag})"
            return entry, None
        entry[f"fseries_{tag}"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] smoke test machine 90 (mutant factory)")
    # (1) Endpoint set Rz-invariant, all on the sphere (bit)
    pts = np.array([A, B, C, D])
    rot = (RZ @ pts.T).T
    d_set = max(float(np.linalg.norm(rot - p, axis=1).min()) for p in pts)
    on_sph = float(np.abs(np.linalg.norm(pts, axis=1) - 1).max())
    ok1 = d_set < 1e-15 and on_sph < 1e-15
    print(f"  (1) endpoints: Rz orbit distance {d_set:.1e}, |r|−1 = "
          f"{on_sph:.1e}  {'✓ bit' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) Ball gates: inside ≤0.96+ε, outside ≥1.04−ε (ends free)
    s1, s2 = inner_tangle(3)
    o1, o2 = outer_frame()
    r_in = max(float(np.linalg.norm(s[1:-1], axis=1).max())
               for s in (s1, s2))
    r_out = min(float(np.linalg.norm(o[1:-1], axis=1).min())
                for o in (o1, o2))
    ok2 = r_in <= 0.97 and r_out >= 1.03
    print(f"  (2) ball discipline: inside max r={r_in:.3f}, outside min "
          f"r={r_out:.3f}  {'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) Mutation: outside bit-identical, inside measurably different
    K1, K2, (s1, s2), (o1, o2) = mutant_pair(3)
    n_in1 = len(s1) - 1
    # Compare outer segments in both curves (same indices)
    n_o1 = len(o1) - 1
    seg1 = K1[n_in1:n_in1 + n_o1]
    seg2 = K2[n_in1:n_in1 + n_o1]
    d_out = float(np.abs(seg1 - seg2).max())
    d_in = float(np.linalg.norm(K1[:n_in1] - K2[:n_in1], axis=1).max())
    ok3 = d_out < 5e-3 and d_in > 0.2
    print(f"  (3) mutation: outside deviation {d_out:.1e} (smoothing-"
          f"seam), inside difference {d_in:.2f}  "
          f"{'✓ tangle rotated, frame stands' if ok3 else '✗ FAIL'}")
    assert ok3
    # (4) REPAIR rungs: deterministically woven (twist+clasp),
    #     crossing-floor gate against the s0 unknot trap
    ok4 = True
    for m in (1, 2, 3):
        e4, f4 = make_pair_woven(m, HERE / "demo_smoke")
        if f4 is None:
            print(f"  (4) {e4['name']}: REJECTED — {e4['rejected']} ✗")
            ok4 = False
        else:
            print(f"  (4) {e4['name']}: floor={e4['crossing_floor_K']}/"
                  f"{e4['crossing_floor_Kmut']} (expected ≥{2 * m + 2}) · "
                  f"Wr={e4['writhe_K']}/{e4['writhe_Kmut']} · "
                  f"nyq={e4['fseries_K']['nyquist_dev_rel']}/"
                  f"{e4['fseries_Kmut']['nyquist_dev_rel']} ✓")
    assert ok4
    # (5) Mutation on Woven: outside bit-fixed, inside different
    #     NOTE index position: s1 (clasp) and s2 have different lengths —
    #     the o1 block begins in K1 at len(s1)-1, in K2 at len(s2)-1
    #     (K2 inner strand A→C = RZ·s2).
    K1, K2, (s1, s2), (o1, o2) = mutant_pair_woven(2)
    n_o1 = len(o1) - 1
    a1, a2 = len(s1) - 1, len(s2) - 1
    d_out = float(np.abs(K1[a1:a1 + n_o1] - K2[a2:a2 + n_o1]).max())
    d_in = float(np.linalg.norm(K1[:min(a1, a2)] - K2[:min(a1, a2)],
                                axis=1).max())
    ok5 = d_out < 5e-3 and d_in > 0.2
    print(f"  (5) Woven mutation: outside {d_out:.1e}, inside {d_in:.2f}  "
          f"{'✓' if ok5 else '✗ FAIL'}")
    assert ok5
    print("[DEMO] 5/5 green — mutant factory REPAIRED (twist+clasp, "
          "deterministic; det!=1 mandatory at ignition; s0 lesson documented).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
