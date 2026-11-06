# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 38 — twisting_the_backbone (internal coordinates / torsion words).

Chain with FIXED bend angle θ_b and TORSION WORD τ (protein principle):
Frame (t,n,b); per step first torsion (rotation about t), then bend
(rotation about b), step along t. Closure via CLOSE-RETURN
(04 method) over repetitions of the word motif + gap gate.
Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import (embedding_report, polyline_to_fseries,  # noqa: E402
                            min_nonadjacent_dist, resample_uniform)


def rot(axis, ang):
    axis = axis / np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * (K @ K)


def chain(torsion_word, theta_b, n_repeat, step=1.0):
    t = np.array([1.0, 0, 0])
    n = np.array([0, 1.0, 0])
    P = [np.zeros(3)]
    for _ in range(n_repeat):
        for tau in torsion_word:
            b = np.cross(t, n)
            Rt = rot(t, tau)
            n = Rt @ n
            b = Rt @ b
            Rb = rot(b, theta_b)
            t = Rb @ t
            n = Rb @ n
            P.append(P[-1] + step * t)
    return np.asarray(P)


def motif_transform(word, theta_b, step=1.0):
    """One motif traversal from identity frame → (R_motif, v_motif)."""
    t = np.array([1.0, 0, 0]); n = np.array([0, 1.0, 0]); p = np.zeros(3)
    F0 = np.column_stack([t, n, np.cross(t, n)])
    for tau in word:
        b = np.cross(t, n)
        Rt = rot(t, tau); n = Rt @ n
        b = np.cross(t, n)
        Rb = rot(b, theta_b); t = Rb @ t; n = Rb @ n
        p = p + step * t
    F1 = np.column_stack([t, n, np.cross(t, n)])
    return F1 @ F0.T, p


def ring_closure(word, theta_b):
    """RING RESONANCE (build lesson: ballistic chains do not return;
    C_n logic in internal coordinates): k = 2π/rotation angle of the motif,
    closing error |Σ R^j v| relative to ring size."""
    R, v = motif_transform(word, theta_b)
    cosphi = (np.trace(R) - 1) / 2
    phi = float(np.arccos(np.clip(cosphi, -1, 1)))
    if phi < 1e-6:
        return None
    k = max(3, int(round(2 * np.pi / phi)))
    E = np.zeros(3); Rj = np.eye(3); radius = 0.0
    for _ in range(k):
        E = E + Rj @ v
        radius = max(radius, float(np.linalg.norm(E)))
        Rj = R @ Rj
    return {"k": k, "phi": phi, "err": float(np.linalg.norm(E)),
            "err_rel": float(np.linalg.norm(E) / max(radius, 1e-9))}


def make_knot(word, theta_deg, out_root):
    theta = np.deg2rad(theta_deg)
    rc = ring_closure(word, theta)
    name = ("bk_" + "_".join(f"{w:+.2f}".replace(".", "p").replace("+", "")
                             .replace("-", "m") for w in word)
            + f"_th{theta_deg}")
    entry = {"name": name, "machine": "38_twisting_the_backbone",
             "torsion_word": [round(w, 3) for w in word],
             "theta_deg": theta_deg,
             "ring": rc,
             "predictions": {"identity": "Zünd-Messung",
                             "symmetry": f"C_{rc['k']}-Resonanz-Kandidat"
                             if rc else None}}
    if rc is None or rc["err_rel"] > 0.06:
        entry["rejected"] = f"Resonanz-Gate (err_rel={rc['err_rel'] if rc else 'inf'})"
        return entry, None
    P = chain(word, theta, rc["k"])
    # distribute remaining gap linearly over the chain (small by gate)
    drift = P[-1] - P[0]
    P = P - np.outer(np.linspace(0, 1, len(P)), drift)
    P = P[:-1]
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(P, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/38", n_terms=128,
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 38 (torsion backbone)")
    # 1) frame mechanics: step length exact, bend angle constant
    P = chain([0.5, -0.3], np.deg2rad(40), 30)
    seg = np.diff(P, axis=0)
    lens = np.linalg.norm(seg, axis=1)
    cosang = np.einsum("ij,ij->i", seg[:-1], seg[1:]) / (lens[:-1] * lens[1:])
    ok1 = np.allclose(lens, 1.0) and np.allclose(np.degrees(np.arccos(
        np.clip(cosang, -1, 1))), 40, atol=1e-6)
    print(f"  (1) step=1, bend=40° constant  {'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # 2) BUILD LESSON (measured): generic torsion words have helix
    #    pitch ⇒ NEVER close; BALANCED words [a,−a,…] kill the
    #    pitch (16-echo: ±-pairs!) and close with err_rel ~1e-2
    rc_bad = ring_closure([2.0, -0.8], np.deg2rad(45))
    rc_good = ring_closure([1.2, -1.2], np.deg2rad(55))
    ok2 = (rc_bad is None or rc_bad["err_rel"] > 0.3) \
        and rc_good["err_rel"] < 0.06
    print(f"  (2) pitch law: generic err="
          f"{rc_bad['err_rel'] if rc_bad else 'inf'} vs balanced err="
          f"{rc_good['err_rel']:.3f} (k={rc_good['k']})  "
          f"{'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # 3) sampler over balanced words + θ fine scan
    rng = np.random.default_rng(5)
    built, tried = [], 0
    while len(built) < 2 and tried < 80:
        tried += 1
        a = float(rng.uniform(0.8, 2.6))
        b = float(rng.uniform(0.3, 1.8))
        word = [a, -a] if tried % 2 else [a, -a, b, -b]
        best = (9, None)
        for th in np.linspace(25, 72, 48):
            rc = ring_closure(word, np.deg2rad(th))
            if rc and 4 <= rc["k"] <= 40 and rc["err_rel"] < best[0]:
                best = (rc["err_rel"], float(th))
        if best[0] > 0.05:
            continue
        e, f = make_knot(word, round(best[1], 1), HERE / "demo_smoke")
        if f is not None:
            built.append((e["name"], e["ring"]["k"]))
    ok3 = len(built) >= 1
    print(f"  (3) balanced sampler: {built} in {tried} attempts  "
          f"{'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    print("[DEMO] 3/3 green — backbone machine ready to fire "
          "(balanced resonance rings).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
