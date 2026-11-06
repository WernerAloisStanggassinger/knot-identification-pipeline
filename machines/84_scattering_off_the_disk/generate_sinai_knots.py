# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 84 — scattering_off_the_disk (Sinai billiard, 3D).

Billiard in the unit cube WITH SPHERICAL SCATTERER in the middle:
dispersive ⇒ CHAOTIC — the counterpart to the integrable machine 57.
BETS EXACT today: reflection law at wall and sphere bit-exact (energy/
|v| preserved, incidence=reflection); Lyapunov > 0 MEASURED (neighbour-orbit
divergence) vs 57-λ=0 — the order/chaos pair in billiards complete.
Closure: best return in phase space + bridge (63er pattern),
declared). REVERSAL CLAUSE: grazing sphere hit = degenerate →
grazing gate (min |angle of incidence|). Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402

R_SPHERE = 0.22
CENTER = np.array([0.5, 0.5, 0.5])
GRAZE_MIN = 0.05      # grazing gate: |cos incidence| at the sphere


def fly(x, v, n_bounces, record=True):
    """Ray in cube + sphere; returns (path points, status)."""
    pts = [x.copy()]
    grazing = False
    for _ in range(n_bounces):
        # next event: wall or sphere
        t_wall = np.inf
        wall_ax, wall_side = -1, 0
        for ax in range(3):
            if v[ax] > 1e-15:
                t = (1 - x[ax]) / v[ax]
                if t < t_wall:
                    t_wall, wall_ax, wall_side = t, ax, 1
            elif v[ax] < -1e-15:
                t = -x[ax] / v[ax]
                if t < t_wall:
                    t_wall, wall_ax, wall_side = t, ax, 0
        # sphere: |x+tv−c|² = R²
        oc = x - CENTER
        b = oc @ v
        disc = b * b - (oc @ oc - R_SPHERE ** 2)
        t_sph = np.inf
        if disc > 0:
            sq = np.sqrt(disc)
            for tc in (-b - sq, -b + sq):
                if tc > 1e-12:
                    t_sph = tc
                    break
        if t_sph < t_wall:
            x = x + t_sph * v
            # Build lesson: normalize normal EXACTLY + snap onto the sphere
            # (dividing by /R is not enough — the hit residual
            # skews the reflection, |v| drift 0.26 measured)
            nrm = x - CENTER
            nrm = nrm / np.linalg.norm(nrm)
            x = CENTER + R_SPHERE * nrm
            if abs(nrm @ v) < GRAZE_MIN:
                grazing = True
            v = v - 2 * (v @ nrm) * nrm
            x = x + 1e-9 * v          # nudge against re-hit
        else:
            x = x + t_wall * v
            v = v.copy()
            v[wall_ax] = -v[wall_ax]
        if record:
            pts.append(x.copy())
    return np.asarray(pts), v, grazing


def lyapunov_estimate(x0, v0, n_bounces=40, eps=1e-9):
    """Divergence of two neighbouring orbits (position difference per bounce)."""
    p1, v1, _ = fly(x0.copy(), v0.copy(), n_bounces)
    x0b = x0 + np.array([eps, 0, 0])
    p2, v2, _ = fly(x0b, v0.copy(), n_bounces)
    m = min(len(p1), len(p2))
    d = np.linalg.norm(p1[:m] - p2[:m], axis=1)
    d = np.maximum(d, 1e-300)
    k = np.arange(m)
    lam = np.polyfit(k[1:m // 2], np.log(d[1:m // 2]), 1)[0]
    return float(lam)


def best_return_orbit(seed, n_bounces=400):
    rng = np.random.default_rng(seed)
    x0 = rng.uniform(0.1, 0.9, 3)
    if np.linalg.norm(x0 - CENTER) < R_SPHERE + 0.1:
        x0 = np.array([0.1, 0.2, 0.85])
    v0 = rng.normal(size=3)
    v0 /= np.linalg.norm(v0)
    pts, v_end, grazing = fly(x0, v0, n_bounces)
    if grazing:
        return None, None, "streif_gate"
    d = np.linalg.norm(pts[20:] - pts[0], axis=1)
    k = int(d.argmin()) + 20
    return pts[:k + 1], float(d.min()), "ok"


def rounded(P, pts=25, rounds=2):
    dense = []
    m = len(P)
    for k in range(m):
        a, b = P[k], P[(k + 1) % m]
        t = np.linspace(0, 1, pts, endpoint=False)
        dense.append(np.outer(1 - t, a) + np.outer(t, b))
    C = np.vstack(dense)
    for _ in range(rounds):
        C = 0.5 * C + 0.25 * (np.roll(C, 1, axis=0) + np.roll(C, -1, axis=0))
    return C


def make_knot(seed, out_root, ret_max=0.12):
    path, ret, status = best_return_orbit(seed)
    name = f"sn_s{seed}"
    entry = {"name": name, "machine": "84_scattering_off_the_disk",
             "seed": seed, "status": status,
             "predictions": {"identity": "Zünd-Messung",
                             "paar": "Chaos-Gegenstück zu 57"}}
    if status != "ok":
        entry["rejected"] = status
        return entry, None
    entry["return_dist"] = round(ret, 4)
    entry["n_bounces_used"] = int(len(path) - 1)
    if ret > ret_max:
        entry["rejected"] = f"Return-Gate ({ret:.3f})"
        return entry, None
    C = rounded(path)
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
                                      source="maschinen_ideen/84",
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
    global R_SPHERE
    print("[DEMO] smoke test machine 84 (Sinai billiard)")
    # (1) Reflection law: |v| bit-preserved over many bounces
    x0 = np.array([0.15, 0.3, 0.8])
    v0 = np.array([0.6, -0.5, 0.62])
    v0 /= np.linalg.norm(v0)
    pts, v_end, _ = fly(x0.copy(), v0.copy(), 200)
    dv = abs(float(np.linalg.norm(v_end)) - 1.0)
    inside = bool((pts >= -1e-9).all() and (pts <= 1 + 1e-9).all())
    d_sph = float(np.linalg.norm(pts - CENTER, axis=1).min())
    ok1 = dv < 1e-12 and inside and d_sph > R_SPHERE - 1e-9
    print(f"  (1) |v| conservation {dv:.1e}, in cube={inside}, "
          f"min sphere distance {d_sph:.4f} ≥ R={R_SPHERE}  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) THE CHAOS PAIR: λ(Sinai) > 0, λ(empty box) ≈ 0
    lam_sinai = lyapunov_estimate(x0, v0)
    R_save = R_SPHERE
    R_SPHERE = 1e-9               # sphere removed = integrable machine-57 box
    lam_box = lyapunov_estimate(x0, v0)
    R_SPHERE = R_save
    ok2 = lam_sinai > 0.2 and abs(lam_box) < 0.05
    print(f"  (2) Lyapunov per bounce: Sinai={lam_sinai:.3f} > 0, "
          f"empty box={lam_box:.4f} ≈ 0  "
          f"{'✓ order/chaos pair complete' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) chaotic orbit knot (scan seeds, gates decide)
    e3 = f3 = None
    for seed in range(60):
        e3, f3 = make_knot(seed, HERE / "demo_smoke")
        if f3 is not None:
            break
    ok3 = f3 is not None
    print(f"  (3) {e3['name']}: return={e3.get('return_dist')}, "
          f"bounces={e3.get('n_bounces_used')}, "
          f"nyq={e3.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok3 else '✗ FAIL: ' + str(e3.get('rejected'))}")
    assert ok3
    print("[DEMO] 3/3 green — Sinai machine ready to fire "
          "(chaos billiard; 57+84 = the integrable/chaotic pair).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
