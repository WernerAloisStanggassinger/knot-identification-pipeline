# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 27 — pumping_the_twist_account (Wr pump, 6D calibration).

Geometric writhe Wr via Gauss SELF-integral (05 mechanics), pump =
smooth helix coil into a curve piece (m turns ≈ ±m Wr). Calibration
bet (manifest, ignition time): the SONO normal form FORGETS the pumped
Wr — same knot, Wr variants ⇒ ONE ideal Wr. Today's trefoil
calibration (6 embeddings → one point) is the special case.
Doctrine: do NOT ignite before Monday.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "06_eating_own_output"))
sys.path.insert(0, str(HERE.parent))
from generate_satellite_knots import trefoil_curve  # noqa: E402
from machine_common import (embedding_report, polyline_to_fseries,  # noqa: E402
                            resample_uniform)


def writhe(P, n=2048):
    """Gauss self-integral over segment midpoints (neighbors masked)."""
    Q = resample_uniform(P, n)
    m = (Q + np.roll(Q, -1, axis=0)) / 2.0
    t = np.roll(Q, -1, axis=0) - Q
    d = m[:, None, :] - m[None, :, :]
    r3 = np.linalg.norm(d, axis=2) ** 3
    idx = np.arange(n)
    ring = np.abs(idx[None, :] - idx[:, None])
    ring = np.minimum(ring, n - ring)
    mask = ring > 2
    cross = np.cross(t[:, None, :], t[None, :, :])
    num = np.einsum("ijk,ijk->ij", d, cross)
    with np.errstate(divide="ignore", invalid="ignore"):
        val = np.where(mask, num / r3, 0.0)
    return float(val.sum()) / (4 * np.pi)


def pump(P, m_turns, frac=0.12, radius_rel=1.5):
    """Insert helix coil with m turns into the curve piece [0, frac]."""
    Q = resample_uniform(P, 2048)
    n_seg = int(len(Q) * frac)
    seg = Q[:n_seg]
    rest = Q[n_seg:]
    p0, p1 = seg[0], seg[-1]
    axis = p1 - p0
    L = np.linalg.norm(axis)
    u = axis / L
    # orthonormal frame
    v = np.cross(u, [0, 0, 1.0])
    if np.linalg.norm(v) < 1e-9:
        v = np.cross(u, [0, 1.0, 0])
    v /= np.linalg.norm(v)
    w = np.cross(u, v)
    r = radius_rel * L / max(1, abs(m_turns))
    s = np.linspace(0, 1, max(64, 48 * abs(m_turns)))
    th = 2 * np.pi * m_turns * s
    env = np.sin(np.pi * s)            # coil swells up/down (smoothly insertable)
    coil = (p0[None, :] + np.outer(s, axis)
            + r * env[:, None] * (np.outer(np.cos(th), v)
                                  + np.outer(np.sin(th), w)))
    return np.vstack([coil, rest])


def make_knot(m_turns, out_root):
    base = trefoil_curve()
    P = pump(base, m_turns) if m_turns else resample_uniform(base, 2048)
    wr0 = writhe(base)
    wr1 = writhe(P)
    # convergence check (flaw 2): second resolution
    wr1b = writhe(P, n=3072)
    name = f"wp_trefoil_m{m_turns}".replace("-", "neg")
    entry = {"name": name, "machine": "27_pumping_the_twist_account",
             "pumped_turns": m_turns,
             "writhe_base": round(wr0, 4),
             "writhe_pumped": round(wr1, 4),
             "writhe_convergence_dev": round(abs(wr1 - wr1b), 4),
             "predictions": {
                 "same_knot_as_base": "Trefoil (Spule = isotop, Wette — "
                                      "Embedding-Gate stützt, Pipeline prüft)",
                 "normal_form_forgets": "SONO-Ideal-Wr unabhängig von "
                                        "gepumptem Wr (6D-Kalibrier-Wette)",
             }}
    if entry["writhe_convergence_dev"] > 0.05:
        entry["rejected"] = "Wr-Integral nicht konvergiert"
        return entry, None
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate (Spule fängt Stränge?)"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(P, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/27_pumping",
                              n_terms=max(96, 24 * abs(m_turns) + 48),
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 27 (twist pump)")
    # 1) calibration: flat circle Wr ≈ 0
    t = np.linspace(0, 2 * np.pi, 512, endpoint=False)
    circ = np.column_stack([np.cos(t), np.sin(t), 0 * t])
    w0 = writhe(circ)
    ok1 = abs(w0) < 0.01
    print(f"  (1) circle: Wr = {w0:+.4f} (expected ~0)  {'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # 2) trefoil curve: |Wr| in plausible order of magnitude (~3)
    w_tr = writhe(trefoil_curve())
    ok2 = 2.0 < abs(w_tr) < 4.5
    print(f"  (2) trefoil parametrisation: Wr = {w_tr:+.3f}  "
          f"{'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # 3) pump works: +4 turns raise Wr by ~4
    P4 = pump(trefoil_curve(), 4)
    w4 = writhe(P4)
    ok3 = 3.0 < (w4 - w_tr) < 7.0
    print(f"  (3) +4 turns: ΔWr = {w4 - w_tr:+.3f} (fat coil, expected +3..+7 — exact m is not the goal, SEPARATION is)  "
          f"{'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    # 4) build calibration series: m = −4, 0, +4 (same knot, three Wr)
    built = []
    for m in (-4, 0, 4):
        e, f = make_knot(m, HERE / "demo_smoke")
        if f is not None:
            built.append((m, e["writhe_pumped"]))
    ok4 = len(built) == 3
    print(f"  (4) calibration series (m, Wr): {built}  "
          f"{'✓' if ok4 else '✗ FAIL'}")
    assert ok4
    print("[DEMO] 4/4 green — twist pump ready to fire (forgetting bet = "
          "ignition program, series is ready).")


def main():
    ap = argparse.ArgumentParser(description="Machine 27 — twist pump")
    sub = ap.add_subparsers(dest="cmd", required=True)
    pm = sub.add_parser("make")
    pm.add_argument("--turns", type=int, required=True)
    sub.add_parser("demo")
    args = ap.parse_args()
    if args.cmd == "demo":
        cmd_demo(args)
    else:
        man_path = HERE / "MANIFEST.json"
        man = json.loads(man_path.read_text()) if man_path.exists() \
            else {"knots": []}
        entry, folder = make_knot(args.turns, HERE / "out_fseries")
        man["knots"] = [k for k in man["knots"]
                        if k["name"] != entry["name"]] + [entry]
        man_path.write_text(json.dumps(man, indent=1))
        print(f"[MAKE] {entry['name']}: Wr {entry['writhe_pumped']} → "
              f"{folder if folder else entry.get('rejected')}")


if __name__ == "__main__":
    main()
