# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 06 — eating_own_output (satellite recursion, Werner's idea 05-27).

Put knots on symmetry corners, iterate outward (receipt double-3×6
→ 116 crossings). Two mathematical readings, BOTH in the prediction tool:

  CONNECTED SUM (corner graft, geometric): Δ = Δ_base · Π Δ_corner,
    det multiplicative. FREE CALIBRATION ANCHOR in our own corpus: knots-of-love =
    3_1#…#3_1 (6×), Δ = (t²−t+1)⁶ bit-exact, det = 3⁶ = 729 (Disc #142).
  SATELLITE (module with inner structure): Δ_sat(t) = Δ_pattern(t) · Δ_companion(t^w),
    w = winding number [Seifert/Litherland]. Consequence identity (built in as
    self-test): det_sat = det_pattern · |Δ_companion((-1)^w)| — for even w
    Δ_c(1) = ±1, the companion vanishes from the determinant.

Every recursion step thereby becomes BETTABLE before it is built.

CONSTRUCTOR (corner graft): find curvature peaks → base there
cut open → open module (cut-open small knot) rigidly into the
turn gap (chord alignment) → close. Iterable on its own output.

Doctrine: do NOT fire before Monday. Iteration depth > 2 = compute question
(observe the high-crossing doctrine when firing).
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import (embedding_report, polyline_to_fseries,  # noqa: E402
                            resample_uniform, densify)

TREFOIL_DELTA = np.array([1, -1, 1], int)        # t² − t + 1


# ------------------------------------------------- Prediction (polynomial core)

def poly_mul(p, q):
    return np.convolve(np.asarray(p, int), np.asarray(q, int))


def poly_comp_tw(p, w):
    """Δ(t) → Δ(t^w): spread the coefficients with w−1 zeros."""
    p = np.asarray(p, int)
    out = np.zeros((len(p) - 1) * w + 1, int)
    out[::w] = p
    return out


def poly_eval(p, x):
    return int(sum(int(c) * x ** k for k, c in enumerate(np.asarray(p, int))))


def predict_connected_sum(deltas):
    """Δ list of the summands → (Δ_total, det)."""
    acc = np.array([1], int)
    for d in deltas:
        acc = poly_mul(acc, d)
    return acc, abs(poly_eval(acc, -1))


def predict_satellite(delta_pattern, delta_companion, w):
    """Satellite formula + det identity as self-test."""
    sat = poly_mul(delta_pattern, poly_comp_tw(delta_companion, w))
    det = abs(poly_eval(sat, -1))
    det_check = abs(poly_eval(delta_pattern, -1)
                    * poly_eval(delta_companion, (-1) ** w))
    if det != det_check:
        raise RuntimeError("Satellite det identity violated — polynomial core rotten")
    return sat, det


# ----------------------------------------------------------- Corner detector

def curvature(P):
    """Discrete curvature per corner of the closed polyline."""
    A = np.roll(P, 1, axis=0) - P
    B = np.roll(P, -1, axis=0) - P
    na = np.linalg.norm(A, axis=1)
    nb = np.linalg.norm(B, axis=1)
    cosang = np.einsum("ij,ij->i", A, B) / (na * nb + 1e-30)
    return np.pi - np.arccos(np.clip(cosang, -1, 1))


def find_corners(P, n_expect=None, min_sep_frac=0.05, mask_windows=None):
    """Curvature peaks, at least min_sep_frac of the circumference apart.

    mask_windows: list (lo_frac, hi_frac) circumference windows that are NOT populated
    (old graft sites — recursion iterates OUTWARD, not in
    into its own modules).
    """
    Q = resample_uniform(P, 2048)
    k = curvature(Q)
    if mask_windows:
        n = len(Q)
        for (lo_f, hi_f) in mask_windows:
            lo, hi = int(lo_f * n) % n, int(hi_f * n) % n
            if lo <= hi:
                k[lo:hi + 1] = -np.inf
            else:
                k[lo:] = -np.inf
                k[:hi + 1] = -np.inf
    order = np.argsort(-k)
    sep = int(len(Q) * min_sep_frac)
    picked = []
    for i in order:
        if not np.isfinite(k[i]):
            break
        if all(min(abs(i - j), len(Q) - abs(i - j)) >= sep for j in picked):
            picked.append(int(i))
        if n_expect and len(picked) == n_expect:
            break
        if not n_expect and len(picked) >= 12:
            break
    return Q, sorted(picked)


# --------------------------------------------------------------- Constructor

def _rot_a_to_b(a, b):
    """Rodrigues rotation that turns direction a onto direction b."""
    u = a / np.linalg.norm(a)
    v = b / np.linalg.norm(b)
    ax = np.cross(u, v)
    s, c = np.linalg.norm(ax), float(u @ v)
    if s < 1e-12:
        return np.eye(3) if c > 0 else -np.eye(3)
    ax = ax / s
    K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
    th = np.arccos(np.clip(c, -1, 1))
    return np.eye(3) + np.sin(th) * K + (1 - c) * (K @ K)


def mount_module(M, p_in, p_out, outward, size_rel=1.0, open_frac=0.12,
                 offset_rel=0.8, n_bridge=24):
    """GRAFT: module scaled small OUTSIDE onto the corner, two bridges.

    - cut out 12% arc → open ends with a proper opening
      (NOT argmax of the point gaps: after resample_uniform all gaps are
      equal — the mini chord caused the diameter-50 bug).
    - scaled isotropically to size_rel·|gap chord| diameter (knot-type-faithful),
    - end direction ∥ chord, centroid rotated towards `outward`,
    - Seat: gap centre + outward·offset_rel·chord, bridges p_in→e1, e2→p_out.
    """
    p_in = np.asarray(p_in, float)
    p_out = np.asarray(p_out, float)
    chord_t = p_out - p_in
    Lt = float(np.linalg.norm(chord_t))
    M = resample_uniform(M, 512)
    n_cut = max(4, int(len(M) * open_frac))
    O = M[n_cut:]                            # open stretch, gap = M[0:n_cut]
    O = O - O.mean(axis=0)
    diam = float(np.linalg.norm(O.max(0) - O.min(0)))
    O = O * (size_rel * Lt / max(diam, 1e-12))
    e1, e2 = O[0], O[-1]
    # 1. Rotation: end direction onto chord direction
    R1 = _rot_a_to_b(e2 - e1, chord_t)
    O = O @ R1.T
    # 2. Rotation about the chord axis: centroid cross-arm towards outward
    u = chord_t / Lt
    mid = (O[0] + O[-1]) / 2.0
    arm = -mid                               # centroid(=0) − ends midpoint
    arm_p = arm - (arm @ u) * u
    out_p = np.asarray(outward, float) - (np.asarray(outward, float) @ u) * u
    if np.linalg.norm(arm_p) > 1e-9 and np.linalg.norm(out_p) > 1e-9:
        R2 = _rot_a_to_b(arm_p, out_p)
        O = O @ R2.T
        mid = (O[0] + O[-1]) / 2.0
    # Seat outside in front of the corner
    seat = (p_in + p_out) / 2.0 \
        + (out_p / max(np.linalg.norm(out_p), 1e-12)) * (offset_rel * Lt)
    O = O + (seat - mid)
    b1 = np.linspace(p_in, O[0], n_bridge, endpoint=False)
    b2 = np.linspace(O[-1], p_out, n_bridge, endpoint=False)[1:]
    return np.vstack([b1, O, b2])


def graft(base, module, n_corners, gap_frac=0.09, size_rel=1.0,
          mask_windows=None):
    """Place the module on the n sharpest corners of the base (connected-sum build).

    Returns (curve, corners, graft windows as circumference fractions) — the
    windows are the mask for the next recursion step.
    """
    Q, corners = find_corners(base, n_expect=n_corners,
                              mask_windows=mask_windows)
    n = len(Q)
    centroid = Q.mean(axis=0)
    h = max(2, int(n * gap_frac / 2))
    cuts = sorted(((c - h) % n, (c + h) % n) for c in corners)
    windows = [(lo / n, hi / n) for (lo, hi) in cuts]
    curve = []
    for k, (lo, hi) in enumerate(cuts):
        p_in = Q[lo]                       # last point before the gap
        p_out = Q[(hi + 1) % n]            # first point after the gap
        corner_pt = Q[(lo + h) % n]
        outward = corner_pt - centroid
        curve.append(mount_module(module, p_in, p_out, outward,
                                  size_rel=size_rel))
        nxt_lo = cuts[(k + 1) % len(cuts)][0]
        arc_idx = []
        i = (hi + 1) % n
        while i != nxt_lo:
            arc_idx.append(i)
            i = (i + 1) % n
        curve.append(Q[arc_idx])
    return np.vstack(curve), corners, windows


# ----------------------------------------------------------------- Generation

def rounded_polygon(n_sides, rounding=0.18, n=1024):
    """C_n-symmetric base with n sharp corners (test/base curve)."""
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    r = 1.0 + rounding * np.cos(n_sides * t)
    z = 0.12 * np.sin(n_sides * t)
    return np.column_stack([r * np.cos(t), r * np.sin(t), z])


def trefoil_curve(n=512):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.column_stack([np.sin(t) + 2 * np.sin(2 * t),
                            np.cos(t) - 2 * np.cos(2 * t),
                            -np.sin(3 * t)])


def make_knot(name, base, module, n_corners, deltas, out_root,
              gap_frac=0.09, mask_windows=None):
    P, corners, windows = graft(base, module, n_corners, gap_frac=gap_frac,
                                mask_windows=mask_windows)
    delta, det = predict_connected_sum(deltas)
    entry = {"name": name, "machine": "06_eating_own_output",
             "n_corners_grafted": len(corners),
             "graft_windows": [[round(a, 4), round(b, 4)] for a, b in windows],
             "predictions": {
                 "alexander_coeffs_lowfirst": [int(c) for c in delta],
                 "det": det,
                 "mode": "Verbundsumme (Ecken-Aufsatz)",
                 "source": "Δ multiplikativ; Eichanker knots-of-love Disc #142",
             }}
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate (Graft-Fenster zu eng? --gap-frac)"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    # small module blobs = high harmonics → generous + Nyquist guard
    rep = polyline_to_fseries(P, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/06_eating_own_output",
                              n_terms=192, nyquist_max_rel=0.02)
    entry["fseries"] = rep
    return entry, folder


# ----------------------------------------------------------------------- CLI

def cmd_predict(args):
    """Bet table for recursion steps (pure mathematics)."""
    rows = []
    delta, det = TREFOIL_DELTA, 3
    for k in range(1, args.depth + 1):
        deltas = [TREFOIL_DELTA] * (args.corners ** k)
        d, dt = predict_connected_sum(deltas)
        rows.append({"iteration": k, "n_modules": args.corners ** k,
                     "det": dt, "alexander_degree": len(d) - 1})
    sat, sdet = predict_satellite(TREFOIL_DELTA, TREFOIL_DELTA, 2)
    out = {"connected_sum_recursion": rows,
           "satellite_example_w2": {
               "alexander_coeffs_lowfirst": [int(c) for c in sat],
               "det": sdet,
               "note": "w=2 gerade → Begleiter-det verschwindet (Identität "
                       "geprüft im Polynom-Kern)"}}
    p = HERE / "recursion_predictions.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"[PREDICT] trefoil modules, {args.corners} corners, depth {args.depth} "
          f"→ {p.name}")
    for r in rows:
        print(f"   Iter {r['iteration']}: {r['n_modules']} modules, "
              f"det={r['det']}, deg Δ={r['alexander_degree']}")


def cmd_make(args):
    base = rounded_polygon(args.corners)
    module = trefoil_curve()
    deltas = [TREFOIL_DELTA] * args.corners
    name = f"sat_c{args.corners}_trefoil_g{args.gap_frac}".replace(".", "p")
    entry, folder = make_knot(name, base, module, args.corners, deltas,
                              HERE / "out_fseries", gap_frac=args.gap_frac)
    man_path = HERE / "MANIFEST.json"
    man = json.loads(man_path.read_text()) if man_path.exists() else {"knots": []}
    man["knots"] = [k for k in man["knots"] if k["name"] != entry["name"]] + [entry]
    man_path.write_text(json.dumps(man, indent=1))
    if folder:
        print(f"[MAKE] {name}: det bet {entry['predictions']['det']}, "
              f"Nyquist {entry['fseries']['nyquist_dev_rel']} → {folder}")
    else:
        print(f"[MAKE] {name} REJECTED: {entry['rejected']}")


def cmd_demo(args):
    print("[DEMO] Smoke test machine 06 (satellite recursion)")
    # 1) Calibration anchor knots-of-love: 6×trefoil → Δ=(t²−t+1)⁶, det=729
    d, det = predict_connected_sum([TREFOIL_DELTA] * 6)
    anchor = np.array([1], int)
    for _ in range(6):
        anchor = poly_mul(anchor, TREFOIL_DELTA)
    ok1 = det == 729 and np.array_equal(d, anchor)
    print(f"  (1) knots-of-love anchor: det={det} (expected 729), "
          f"Δ=(t²−t+1)⁶ coefficients equal  {'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # 2) Satellite det identity w even/odd
    _, dw2 = predict_satellite(TREFOIL_DELTA, TREFOIL_DELTA, 2)
    _, dw3 = predict_satellite(TREFOIL_DELTA, TREFOIL_DELTA, 3)
    ok2 = dw2 == 3 and dw3 == 9
    print(f"  (2) Satellite det: w=2 → {dw2} (companion gone), w=3 → {dw3} "
          f"(3·3)  {'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # 3) Corner detector: C_3 base → exactly 3 corners, ~120° apart
    base = rounded_polygon(3)
    Q, corners = find_corners(base, n_expect=3)
    gaps = np.diff(corners + [corners[0] + len(Q)])
    ok3 = len(corners) == 3 and np.all(np.abs(gaps - len(Q) / 3) < len(Q) * 0.1)
    print(f"  (3) Corner detector C_3: {len(corners)} corners, distances {list(gaps)} "
          f"(expected ~{len(Q)//3})  {'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    # 4) Graft: 3 trefoil modules on C_3 base, all gates
    entry, f = make_knot("demo_sat_c3", base, trefoil_curve(), 3,
                         [TREFOIL_DELTA] * 3, HERE / "demo_smoke")
    ok4 = f is not None and entry["predictions"]["det"] == 27
    print(f"  (4) Graft 3×trefoil: det bet {entry['predictions']['det']} "
          f"(expected 27) · embedded={entry['embedding']['embedded']} · Nyquist="
          f"{entry.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok4 else '✗ FAIL: ' + str(entry.get('rejected'))}")
    assert ok4
    # 5) EAT YOUR OWN OUTPUT: iteration 2 on the result of (4),
    #    old graft windows masked (iterate outward)
    out1, _, win1 = graft(base, trefoil_curve(), 3)
    entry2, f2 = make_knot("demo_sat_c3_iter2", out1, trefoil_curve(), 3,
                           [TREFOIL_DELTA] * 6, HERE / "demo_smoke",
                           mask_windows=win1)
    ok5 = f2 is not None and entry2["predictions"]["det"] == 729
    print(f"  (5) Iteration 2 (eats output of 4, windows masked): "
          f"det bet {entry2['predictions']['det']} (3⁶=729) · embedded="
          f"{entry2['embedding']['embedded']}  "
          f"{'✓' if ok5 else '✗ FAIL: ' + str(entry2.get('rejected'))}")
    assert ok5
    print("[DEMO] 5/5 green — recursion machine ready to fire (firing after "
          "Monday; depth>2 = high-crossing doctrine).")


def main():
    ap = argparse.ArgumentParser(description="Machine 06 — satellite recursion")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("predict")
    p1.add_argument("--corners", type=int, default=3)
    p1.add_argument("--depth", type=int, default=3)
    p2 = sub.add_parser("make")
    p2.add_argument("--corners", type=int, default=3)
    p2.add_argument("--gap-frac", type=float, default=0.09)
    sub.add_parser("demo")
    args = ap.parse_args()
    {"predict": cmd_predict, "make": cmd_make, "demo": cmd_demo}[args.cmd](args)


if __name__ == "__main__":
    main()
