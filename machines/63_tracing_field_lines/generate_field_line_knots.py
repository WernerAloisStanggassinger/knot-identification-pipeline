# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 63 — tracing_field_lines (Biot–Savart field lines).

Magnetic field of a conductor loop (Biot–Savart over polygon segments),
field lines = RK4 flow + close-return gate (machine-51 mechanics).
EXACT ANCHOR: every field line of a CIRCLE coil closes through
the hole and links the coil with lk = 1 (combinatorial lk from
61 as witness). Then: field of a TREFOIL coil (T(2,3) from machine 10)
— are its field lines knotted/higher linked? (helicity = Σ lk,
plasma-physics link.) SCULPTURAL: the invisible field of the
sculpture, made visible.
REVERSAL CLAUSE: lines close to the wire explode (1/r²) → wire
distance gate; ODE drift eats closure → return gate hard.
Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "61_entangling_complete_graphs"))
sys.path.insert(0, str(HERE.parent / "10_climbing_the_torus_ladder"))
sys.path.insert(0, str(HERE.parent))
from generate_conway_gordon_knots import gauss_lk  # noqa: E402
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402

WIRE_GAP_MIN = 0.15     # min distance field line↔wire (relative to R=1)
CLOSE_R = 0.05
MIN_STEPS = 200


def biot_savart(pts, wire):
    """B field at pts (n,3) of a closed polygon wire (m,3) —
    Segment sum of the Biot–Savart integral (unit current), over the
    segments VECTORISED (build lesson: the Python segment loop made
    the RK4 flow unusably slow — 61M mini-ops per line)."""
    P = np.atleast_2d(pts)                       # (n,3)
    dl = np.roll(wire, -1, axis=0) - wire        # (m,3)
    mid = (wire + np.roll(wire, -1, axis=0)) / 2
    r = P[:, None, :] - mid[None, :, :]          # (n,m,3)
    r2 = (r ** 2).sum(axis=2)                    # (n,m)
    num = np.cross(dl[None, :, :], r)            # (n,m,3)
    B = (num / np.maximum(r2, 1e-30)[:, :, None] ** 1.5).sum(axis=1)
    return B / (4 * np.pi)


def wire_distance(p, wire):
    return float(np.linalg.norm(wire - p, axis=1).min())


def trace_line(start, wire, h=0.02, n_max=20000, wire_gap=None):
    """RK4 field-line flow (normalized direction = arc-length step)
    until close-return; gates: wire distance + return."""
    if wire_gap is None:
        wire_gap = WIRE_GAP_MIN
    def f(x):
        b = biot_savart(x[None, :], wire)[0]
        n = np.linalg.norm(b)
        return b / max(n, 1e-30)
    x = np.asarray(start, float)
    path = [x]
    for it in range(n_max):
        k1 = f(x)
        k2 = f(x + 0.5 * h * k1)
        k3 = f(x + 0.5 * h * k2)
        k4 = f(x + h * k3)
        x = x + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        if wire_distance(x, wire) < wire_gap:
            return np.asarray(path), "wire_gate"
        path.append(x)
        if it > MIN_STEPS and np.linalg.norm(x - path[0]) < CLOSE_R:
            return np.asarray(path), "closed"
    return np.asarray(path), "open"


def best_return_cut(path, diam, min_steps=MIN_STEPS, max_gap_rel=0.08):
    """PHYSICS LESSON (measured, 4/4 start offsets): field lines close
    exactly only with symmetry — generic coils (trefoil) deliver
    quasi-periodic lines. Machine closure rule: at the point of best
    cut return + bridge linearly (38er drift pattern), residual gap
    honestly documented + gated."""
    if len(path) <= min_steps + 2:
        return None, np.inf
    d = np.linalg.norm(path[min_steps:] - path[0], axis=1)
    k = int(d.argmin()) + min_steps
    gap_rel = float(d.min() / diam)
    if gap_rel > max_gap_rel:
        return None, gap_rel
    cut = path[:k + 1]
    t = np.linspace(0, 1, 24, endpoint=False)[1:]
    bridge = np.outer(1 - t, cut[-1]) + np.outer(t, cut[0])
    return np.vstack([cut, bridge]), gap_rel


def circle_wire(n=256, R=1.0):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.column_stack([R * np.cos(t), R * np.sin(t), 0 * t])


def trefoil_wire(n=512):
    sys.path.insert(0, str(HERE.parent / "10_climbing_the_torus_ladder"))
    from generate_torus_ladder import torus_coeffs
    from the_knotfamily_machine.fseries_io import reconstruct_from_fourier
    a, b = torus_coeffs(2, 3)
    P = reconstruct_from_fourier(a, b, n_samples=n)
    return 0.5 * P / np.abs(P).max()      # civil scale around the origin


def make_knot(wire, start, out_root, tag, wire_name, h=0.02,
              n_max=20000, wire_gap=None):
    diam = float(np.linalg.norm(wire.max(0) - wire.min(0)))
    path, status = trace_line(start, wire, h=h, n_max=n_max,
                              wire_gap=wire_gap)
    name = f"fl_{tag}"
    entry = {"name": name, "machine": "63_tracing_field_lines",
             "wire": wire_name, "status": status,
             "n_steps": int(len(path)),
             "predictions": {
                 "helizitaet": "Σ lk(Feldlinie, Spule) = Helizitäts-"
                               "Beitrag (gemessen)",
                 "identity": "Zünd-Messung"}}
    if status == "closed":
        loop, gap_rel = path, 0.0
    else:
        loop, gap_rel = best_return_cut(path, diam)
        if loop is None:
            entry["rejected"] = (f"Return-Gate ({status}, beste Rückkehr "
                                 f"{gap_rel:.3f}·diam)")
            return entry, None
    entry["return_gap_rel"] = round(gap_rel, 5)
    path = loop
    lk = gauss_lk(path[::max(1, len(path) // 400)], wire[::4])
    entry["lk_line_wire"] = lk
    emb = embedding_report(path)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    # Escalate terms until the Nyquist guard is green (long quasi-
    # periodic lines carry high harmonics)
    rep, err = None, None
    for n_terms in (96, 160, 256, 384):
        try:
            rep = polyline_to_fseries(path, folder / f"{name}.fseries",
                                      knot_name=name,
                                      source="maschinen_ideen/63",
                                      n_terms=n_terms, nyquist_max_rel=0.02)
            break
        except RuntimeError as e:
            err = e
    if rep is None:
        entry["rejected"] = f"Nyquist auch bei 384 Termen rot ({err})"
        return entry, None
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 63 (field lines)")
    wire = circle_wire()
    # (1) Current-reversal symmetry exact: B(−I) = −B (wire backwards)
    pts = np.array([[0.4, 0.2, 0.3], [1.6, -0.1, 0.2]])
    d = np.abs(biot_savart(pts, wire) + biot_savart(pts, wire[::-1])).max()
    ok1 = d < 1e-12
    print(f"  (1) current reversal B(−I)=−B: max|Δ|={d:.2e}  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) axial field of the circle coil against the closed formula
    z = np.array([[0, 0, 0.7]])
    B = biot_savart(z, wire)[0]
    B_exact = 1.0 / (2 * (1 + 0.7 ** 2) ** 1.5)     # unit current, R=1
    ok2 = abs(B[2] - B_exact) < 1e-4 and np.abs(B[:2]).max() < 1e-12
    print(f"  (2) axial formula: Bz={B[2]:.6f} vs exact {B_exact:.6f}  "
          f"{'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) THE ANCHOR: circle-coil field line closes + lk(line,coil)=1
    e3, f3 = make_knot(wire, np.array([0.45, 0.0, 0.0]),
                       HERE / "demo_smoke", "circle_inner", "circle")
    ok3 = f3 is not None and abs(abs(e3.get("lk_line_wire", 0)) - 1) < 1e-9
    print(f"  (3) circle anchor: status={e3['status']}, "
          f"lk={e3.get('lk_line_wire')}, steps={e3['n_steps']}  "
          f"{'✓' if ok3 else '✗ FAIL: ' + str(e3.get('rejected'))}")
    assert ok3
    # (4) Trefoil coil: line traced, best-return cut (lines
    #     of generic coils are quasi-periodic — physics, measured),
    #     lk(line, coil) as helicity witness
    tw = trefoil_wire(n=256)
    diam = float(np.linalg.norm(tw.max(0) - tw.min(0)))
    e4 = f4 = None
    for off in (0.15, 0.22, 0.10, 0.30):
        p0 = tw[0]
        t0 = tw[1] - tw[0]; t0 /= np.linalg.norm(t0)
        away = p0 / np.linalg.norm(p0)
        away -= away @ t0 * t0; away /= np.linalg.norm(away)
        e4, f4 = make_knot(tw, p0 + off * diam * away, HERE / "demo_smoke",
                           f"trefoil_off{off}", "trefoil", h=0.012,
                           n_max=30000, wire_gap=0.02 * diam)
        if f4 is not None:
            break
    ok4 = f4 is not None
    print(f"  (4) trefoil coil line: status={e4['status']}, "
          f"return_gap={e4.get('return_gap_rel')}, "
          f"lk={e4.get('lk_line_wire')}, "
          f"nyq={e4.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok4 else '✗ FAIL: ' + str(e4.get('rejected'))}")
    assert ok4
    print("[DEMO] 4/4 green — field-line machine ready to fire "
          "(the invisible field of the sculpture, now measurable).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
