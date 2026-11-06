# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 03 — drawing_on_surfaces, stage 1: BILLIARD KNOTS.

Marble run in a mirror box: every coordinate is an independent
TRIANGLE WAVE — the angular cousin of the Lissajous knots (machine 01):
    x(t) = A·tri(p·t + φx),  y(t) = B·tri(q·t + φy),  z(t) = C·tri(r·t + φz)
with tri(θ) = (2/π)·arcsin(sin θ) (exact triangle wave, period 2π).
reflection at the wall = kink of the wave. For Werner: steel-rod polygon chain
in the mirror box — a direct build instruction.

GEOMETRICALLY EXACT (no literature belief): tri(θ+π) = −tri(θ). If ALL
frequencies are odd, P(t+π) = −P(t) BIT-EXACTLY → the curve is
point-symmetric through the origin → strongly-plus-amphichiral by
construction. The demo test measures that to machine epsilon.

Shared bets with machine 01 (lit check 06-04, billiard ↔ Lissajous closely
related [Boocher et al. / Lamm]): Arf=0 ⇒ det ≡ 1 mod 8; all_odd ⇒ det =
odd square (Hartley-Kawauchi via the exact point symmetry).

Doctrine: do NOT fire before Monday — only out_fseries/ folder + MANIFEST.
"""
from __future__ import annotations
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402

N_SAMPLES = 4096


def tri(theta):
    """Exact triangle wave, period 2π, values [−1,1], tri(θ+π) = −tri(θ)."""
    return (2.0 / np.pi) * np.arcsin(np.sin(theta))


def billiard_curve(freqs, phases, box=(1.0, 1.0, 1.0), n=N_SAMPLES):
    t = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    cols = [box[d] * tri(freqs[d] * t + phases[d]) for d in range(3)]
    return np.column_stack(cols)


def parity_class(freqs):
    n_even = sum(1 for f in freqs if f % 2 == 0)
    return "all_odd" if n_even == 0 else f"{n_even}_even"


def point_symmetry_residual(freqs, phases, box, n=2048):
    """max |P(t+π) + P(t)| / diameter — exactly 0 for all_odd."""
    t = np.linspace(0.0, np.pi, n, endpoint=False)
    P1 = np.column_stack([box[d] * tri(freqs[d] * t + phases[d])
                          for d in range(3)])
    P2 = np.column_stack([box[d] * tri(freqs[d] * (t + np.pi) + phases[d])
                          for d in range(3)])
    diam = 2.0 * max(box)
    return float(np.abs(P1 + P2).max()) / diam


def predictions(freqs, phases, box):
    par = parity_class(freqs)
    p = {"parity": par,
         "shared_bets_with_machine_01":
             "Arf=0, det≡1 mod 8 (Billard↔Lissajous-Verwandtschaft, "
             "Lit-Check 2026-06-04 geteilt)"}
    if par == "all_odd":
        res = point_symmetry_residual(freqs, phases, box)
        p.update({
            "point_symmetry_residual": res,
            "point_symmetric_exact": bool(res < 1e-12),
            "strongly_plus_amphichiral": True,
            "det_is_odd_square": True,
            "basis": "P(t+π)=−P(t) GEOMETRISCH exakt (gemessen) + "
                     "Hartley-Kawauchi 1979",
        })
    else:
        p.update({"period_2": True})
    return p


def coprime_triples(max_freq):
    out = []
    for a in range(2, max_freq + 1):
        for b in range(a + 1, max_freq + 1):
            if math.gcd(a, b) != 1:
                continue
            for c in range(b + 1, max_freq + 1):
                if math.gcd(a, c) == 1 and math.gcd(b, c) == 1:
                    out.append((a, b, c))
    return out


def make_knot(freqs, phases, out_root, box=(1.0, 1.0, 1.0)):
    name = ("bl_" + "_".join(map(str, freqs))
            + f"_ph{phases[0]:.3f}_{phases[1]:.3f}_{phases[2]:.3f}"
            ).replace(".", "p")
    P = billiard_curve(freqs, phases, box)
    entry = {"name": name, "machine": "03_drawing_on_surfaces/billiard",
             "freqs": list(freqs), "phases": [round(x, 6) for x in phases],
             "box": list(box),
             "predictions": predictions(freqs, phases, box)}
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate (singuläre Phase / Wand-Treffer)"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    # Kinks (wall reflections) = slowly decaying spectrum → generous terms;
    # the Nyquist guard in machine_common trips on under-resolution.
    rep = polyline_to_fseries(P, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/03_drawing_on_surfaces",
                              n_terms=128, nyquist_max_rel=0.02)
    entry["fseries"] = rep
    return entry, folder


def phase_raster(k):
    vals = [(i + 1) * math.pi / (k + 1) for i in range(k)]
    return [(px, py, 0.0) for px in vals for py in vals]


# ----------------------------------------------------------------------- CLI

def cmd_plan(args):
    triples = coprime_triples(args.max_freq)
    rows = [{"freqs": list(t), "parity": parity_class(t)} for t in triples]
    n_odd = sum(1 for r in rows if r["parity"] == "all_odd")
    out = HERE / f"billiard_plan_maxfreq{args.max_freq}.json"
    out.write_text(json.dumps({"triples": rows}, indent=1))
    print(f"[PLAN] {len(rows)} triples ≤ {args.max_freq}: {n_odd} all_odd "
          f"(exact point-symmetry bet) → {out.name}")
    print("  First candidates (first-run sketch): (3,4,5) box analogy, "
          "(3,5,7) all_odd, (4,5,7)")


def cmd_make(args):
    out_root = HERE / "out_fseries"
    entries, made = [], 0
    for ph in phase_raster(args.phases):
        e, p = make_knot(tuple(args.freqs), ph, out_root)
        entries.append(e)
        made += p is not None
    man_path = HERE / "MANIFEST.json"
    man = json.loads(man_path.read_text()) if man_path.exists() else {"knots": []}
    known = {k["name"] for k in man["knots"]}
    man["knots"] += [e for e in entries if e["name"] not in known]
    man_path.write_text(json.dumps(man, indent=1))
    print(f"[MAKE] {tuple(args.freqs)}: {made} written, "
          f"{len(entries) - made} rejected → out_fseries/ + MANIFEST.json")


def cmd_demo(args):
    print("[DEMO] Smoke test machine 03 (billiard)")
    # 1) Triangle wave: antisymmetry + value range
    th = np.linspace(0, 2 * np.pi, 1000)
    e1 = float(np.abs(tri(th + np.pi) + tri(th)).max())
    rng = (float(tri(th).min()), float(tri(th).max()))
    print(f"  (1) tri(θ+π)=−tri(θ): max|Δ|={e1:.2e} · range=({rng[0]:.3f},"
          f"{rng[1]:.3f})  {'✓' if e1 < 1e-12 else '✗ FAIL'}")
    assert e1 < 1e-12
    # 2) all_odd ⇒ exact point symmetry of the 3D curve
    res = point_symmetry_residual((3, 5, 7), (0.61, 1.17, 0.0), (1, 1, 1))
    print(f"  (2) (3,5,7) point-symmetry residual = {res:.2e}  "
          f"{'✓' if res < 1e-12 else '✗ FAIL'}")
    assert res < 1e-12
    # 3) an even frequency breaks the symmetry (control)
    res2 = point_symmetry_residual((3, 4, 7), (0.61, 1.17, 0.0), (1, 1, 1))
    print(f"  (3) (3,4,7) residual = {res2:.3f} (must be >0)  "
          f"{'✓' if res2 > 0.1 else '✗ FAIL'}")
    assert res2 > 0.1
    # 4) Build + gates: generic phase
    e, f = make_knot((3, 4, 7), (0.61, 1.17, 0.0), HERE / "demo_smoke")
    ok = f is not None
    print(f"  (4) (3,4,7) build: embedded={e['embedding']['embedded']} "
          f"d_min_rel={e['embedding']['d_min_rel']} · Nyquist dev="
          f"{e.get('fseries', {}).get('nyquist_dev_rel', 'n/a')}  "
          f"{'✓' if ok else '✗ FAIL'}")
    assert ok
    # 5) Singularity guard: phase 0/0/0 runs through the box corner
    e0, f0 = make_knot((3, 5, 7), (0.0, 0.0, 0.0), HERE / "demo_smoke")
    print(f"  (5) (3,5,7) φ=0 (corner path): embedded="
          f"{e0['embedding']['embedded']} → "
          f"{'rejected ✓' if f0 is None else '✗ FAIL (should have been rejected)'}")
    assert f0 is None
    print("[DEMO] 5/5 green — billiard stage ready to fire (firing after Monday).")


def main():
    ap = argparse.ArgumentParser(description="Machine 03 — billiard knots")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("plan"); p1.add_argument("--max-freq", type=int, default=12)
    p2 = sub.add_parser("make")
    p2.add_argument("freqs", type=int, nargs=3, metavar=("P", "Q", "R"))
    p2.add_argument("--phases", type=int, default=3)
    sub.add_parser("demo")
    args = ap.parse_args()
    {"plan": cmd_plan, "make": cmd_make, "demo": cmd_demo}[args.cmd](args)


if __name__ == "__main__":
    main()
