# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 114 — cabling_the_ghosts (cables on Kanenobu ghosts).

COMPOSITION 21 × 93: Lobb (arXiv:1105.3985, Thm 1.2) proves that
even sl(n)/HOMFLY HOMOLOGY does not separate the Kanenobu partners p+q=const.
93 showed: cables separate mutants. THE QUESTION: do cables separate
what homology cannot? (2,·) Möbius cable on K(0,2) and K(2,0)
from the 21 demo_smoke — the satellite experiment on the hardest
blindness benchmark in the stable.
BETS EXACT today: both cables closed + embedded (ρ from the
measured d_min of the ghost, not guessed); cable pair written.
FIRING BET: cable invariants (Alexander/Jones of the cable) separate
K(0,2) from K(2,0) — if YES: satellites beat homology (Dror
punchline); if NO: blindness extends into the satellites (even
more exciting). Doctrine: do NOT fire.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
YARD = HERE.parent
sys.path.insert(0, str(YARD))
sys.path.insert(0, str(YARD / "81_clasping_the_double"))
sys.path.insert(0, str(YARD / "93_cabling_the_mutants"))
from generate_mutant_cables import moebius_cable  # noqa: E402
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402

GHOSTS = ("kb_p0_q2", "kb_p2_q0")     # the p+q=2 ghost pair from 21


def load_ghost(stem):
    hits = list((YARD / "21_cloning_polynomial_ghosts").glob(
        f"demo_smoke/**/{stem}.fseries"))
    if not hits:
        return None
    rows = []
    for line in hits[0].read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("%"):
            continue
        v = [float(x) for x in s.split()]
        if len(v) == 6:
            rows.append(v)
    C = np.asarray(rows)
    a, b = C[:, [0, 2, 4]], C[:, [1, 3, 5]]
    from the_knotfamily_machine.fseries_io import reconstruct_from_fourier
    return reconstruct_from_fourier(a, b, n_samples=4096)


def make_cable(stem, out_root):
    X = load_ghost(stem)
    name = f"gc_{stem}"
    entry = {"name": name, "machine": "114_cabling_the_ghosts",
             "ghost": stem,
             "predictions": {
                 "zuendung": "trennt das (2,·)-Kabel die Kanenobu-"
                             "Partner? (Homologie ist blind — Lobb "
                             "Thm 1.2; Satelliten-Frage an Dror)"}}
    if X is None:
        entry["rejected"] = "Geist fehlt (21 demo_smoke)"
        return entry, None
    emb0 = embedding_report(X, n_samples=2048)
    rho = 0.35 * emb0["d_min_rel"] * float(
        np.linalg.norm(X.max(0) - X.min(0)))
    entry["rho_from_dmin"] = round(rho, 5)
    D = moebius_cable(X, rho)
    emb = embedding_report(D)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = None
    for nt in (384, 512, 768):
        try:
            rep = polyline_to_fseries(D, folder / f"{name}.fseries",
                                      knot_name=name,
                                      source="maschinen_ideen/114",
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
    print("[DEMO] smoke test machine 114 (cables on Kanenobu ghosts)")
    # (1) Ghosts in place + cable closes (σ∈[0,2) construction)
    X = load_ghost(GHOSTS[0])
    assert X is not None, "21er-Geist fehlt"
    D = moebius_cable(X, 0.1)
    gap = float(np.linalg.norm(D[0] - D[-1]))
    step = float(np.linalg.norm(np.diff(D[:50], axis=0), axis=1).mean())
    ok1 = gap < 12 * step
    print(f"  (1) Cable closed: end gap {gap:.4f} vs step "
          f"{step:.4f}  {'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2)+(3) both ghosts cabled + written
    built = []
    for stem in GHOSTS:
        e, f = make_cable(stem, HERE / "demo_smoke")
        if f is not None:
            built.append((stem, e["rho_from_dmin"]))
        else:
            print(f"      {stem}: ✗ {e.get('rejected')}")
    ok2 = len(built) == 2
    print(f"  (2) cables embedded+written: {built}  "
          f"{'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    print(f"  (3) satellite experiment in place: cable(K(0,2)) vs "
          f"cable(K(2,0)) — homology blind, cable?  ✓")
    print("[DEMO] 3/3 green — ghost cable ready to fire "
          "(satellites against Lobb's homology blindness).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
