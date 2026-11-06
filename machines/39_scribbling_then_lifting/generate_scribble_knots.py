# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 39 — scribbling_then_lifting (alternating knots by construction).

Random planar Fourier curve → find self-crossings → lift ALTERNATING:
passage parity along the traverse determines over/under. GAUSS
PARITY THEOREM (VERIFIED per build): for closed planar curves
an even number of other passages lies between the two passages of every crossing ⇒
parities always opposite ⇒ the alternating lift is ALWAYS
consistent. BET (Tait track): alternating by construction;
c <= crossing number, equality if diagram reduced (ignition check).
Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402

N = 2048


def plane_curve(spec, n=N):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    x = sum(A * np.cos(k * t + p) for (k, A, p) in spec["x"])
    y = sum(A * np.sin(k * t + p) for (k, A, p) in spec["y"])
    return np.column_stack([x, y]), t


def crossings(C):
    """Self-intersections of the closed polyline: (s_i, s_j) parameters."""
    n = len(C)
    A = C
    B = C[(np.arange(n) + 1) % n]
    out = []
    for i in range(n):
        j = np.arange(i + 2, n if i > 0 else n - 1)
        if len(j) == 0:
            continue
        p, r = A[i], B[i] - A[i]
        q_, s = A[j], B[j] - A[j]
        denom = r[0] * s[:, 1] - r[1] * s[:, 0]
        dq = q_ - p
        with np.errstate(divide="ignore", invalid="ignore"):
            tt = (dq[:, 0] * s[:, 1] - dq[:, 1] * s[:, 0]) / denom
            uu = (dq[:, 0] * r[1] - dq[:, 1] * r[0]) / denom
        hit = (np.abs(denom) > 1e-12) & (tt > 1e-9) & (tt < 1 - 1e-9) \
            & (uu > 1e-9) & (uu < 1 - 1e-9)
        for k in np.nonzero(hit)[0]:
            out.append(((i + tt[k]) / n, (j[k] + uu[k]) / n))
    return out


def gauss_parity_check(cross):
    """Sort passages; every crossing must carry opposite parity."""
    passages = []
    for cid, (s1, s2) in enumerate(cross):
        passages.append((s1, cid))
        passages.append((s2, cid))
    passages.sort()
    par = {}
    ok = True
    for idx, (s, cid) in enumerate(passages):
        if cid in par:
            ok &= (par[cid] != idx % 2)
        else:
            par[cid] = idx % 2
    return ok, passages


def lift(C, cross, h_rel=0.12):
    diam = float(np.linalg.norm(C.max(0) - C.min(0)))
    h = h_rel * diam
    ok, passages = gauss_parity_check(cross)
    if not ok:
        raise RuntimeError("Gauss parity violated — not a planar curve?")
    t = np.linspace(0, 1, len(C), endpoint=False)
    z = np.zeros(len(C))
    if passages:
        gaps = np.diff([p[0] for p in passages] +
                       [passages[0][0] + 1.0])
        sigma = max(float(gaps.min()) / 3.0, 1.0 / len(C))
        for idx, (s, cid) in enumerate(passages):
            sgn = 1.0 if idx % 2 == 0 else -1.0
            d = np.minimum(np.abs(t - s), 1.0 - np.abs(t - s))
            z += sgn * h * np.exp(-0.5 * (d / sigma) ** 2)
    return np.column_stack([C[:, 0], C[:, 1], z])


def make_knot(seed, out_root):
    rng = np.random.default_rng(seed)
    spec = {"x": [(int(k), float(rng.uniform(0.3, 1.0)),
                   float(rng.uniform(0, 6.28))) for k in (1, 2, 3, 4)],
            "y": [(int(k), float(rng.uniform(0.3, 1.0)),
                   float(rng.uniform(0, 6.28))) for k in (1, 2, 3, 4)]}
    C, _ = plane_curve(spec)
    cross = crossings(C)
    name = f"sc_s{seed}"
    entry = {"name": name, "machine": "39_scribbling_then_lifting",
             "seed": seed, "n_crossings": len(cross),
             "predictions": {
                 "alternating": "per Konstruktion (Gauß-Parität verifiziert)",
                 "c_bound": len(cross),
                 "tait": "c == Kreuzungen falls reduziert (Zünd-Prüfung)"}}
    if len(cross) < 3:
        entry["rejected"] = "zu wenige Kreuzungen (Unknot-Land)"
        return entry, None
    P = lift(C, cross)
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(P, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/39",
                              n_terms=max(96, 16 * len(cross)),
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    return entry, folder


def cmd_demo(args):
    print("[DEMO] smoke test machine 39 (scribble+lift)")
    rng = np.random.default_rng(3)
    spec = {"x": [(int(k), float(rng.uniform(0.3, 1.0)),
                   float(rng.uniform(0, 6.28))) for k in (1, 2, 3, 4)],
            "y": [(int(k), float(rng.uniform(0.3, 1.0)),
                   float(rng.uniform(0, 6.28))) for k in (1, 2, 3, 4)]}
    C, _ = plane_curve(spec)
    cr = crossings(C)
    ok1 = len(cr) >= 3
    print(f"  (1) crossings found: {len(cr)}  {'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    ok2, _ = gauss_parity_check(cr)
    print(f"  (2) GAUSS PARITY (theorem, measured): all crossings "
          f"opposite  {'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    built = []
    for seed in range(20):
        e, f = make_knot(seed, HERE / "demo_smoke")
        if f is not None:
            built.append((e["name"], e["n_crossings"]))
        if len(built) >= 2:
            break
    ok3 = len(built) >= 2
    print(f"  (3) alternating factory: {built}  {'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    print("[DEMO] 3/3 green — scribble machine ready to fire "
          "(alternating knots on the assembly line).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
