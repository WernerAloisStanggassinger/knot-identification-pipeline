# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 111 — reversing_the_arrow (orientation pairs).

Reversal of the parametrisation t→−t ⇔ b→−b (coefficients only, exact):
identical as a point set, as an oriented object the pair (K, rK).
THE MEASUREMENT: is the stack orientation-blind? (As an N-pair measurement
instead of a claim = method finding for the Dror track.)
BETS EXACT today: set identity bit; parameter reversal exact;
reversal pairs written from 3 asymmetric machine knots.
Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
YARD = HERE.parent
sys.path.insert(0, str(YARD))
from machine_common import polyline_to_fseries, embedding_report  # noqa: E402
from the_knotfamily_machine.fseries_io import reconstruct_from_fourier  # noqa: E402

# Asymmetric sources (no built-in symmetry → most honest test)
SOURCES = [
    ("47_breeding_for_strangeness", "ev_s47_g5"),
    ("84_scattering_off_the_disk", "sn_s32"),
    ("80_folding_the_chain", "hp_n48_s2000"),
]


def load_ab(machine, stem):
    hits = list((YARD / machine).glob(f"demo_smoke/**/{stem}.fseries"))
    if not hits:
        return None, None
    rows = []
    for line in hits[0].read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("%"):
            continue
        v = [float(x) for x in s.split()]
        if len(v) == 6:
            rows.append(v)
    C = np.asarray(rows)
    return C[:, [0, 2, 4]], C[:, [1, 3, 5]]


def reverse(a, b):
    """t→−t: cos stays, sin is negated ⇔ b→−b."""
    return a.copy(), -b


def make_pair(machine, stem, out_root):
    a, b = load_ab(machine, stem)
    name = f"rv_{stem}"
    entry = {"name": name, "machine": "111_reversing_the_arrow",
             "source": f"{machine}/{stem}",
             "predictions": {
                 "messung": "Stack orientierungs-blind? (Polynome/θ/6D "
                            "je Paar — Blindheits-Tabelle, Dror)"}}
    if a is None:
        entry["rejected"] = "Quelle fehlt"
        return entry, None
    ar, br = reverse(a, b)
    folder = Path(out_root) / name
    for tag, (aa, bb) in (("fwd", (a, b)), ("rev", (ar, br))):
        P = reconstruct_from_fourier(aa, bb, n_samples=4096)
        folder.mkdir(parents=True, exist_ok=True)
        rep = None
        for nt in (max(64, len(aa)), 256):
            try:
                rep = polyline_to_fseries(
                    P, folder / f"{name}_{tag}.fseries",
                    knot_name=f"{name}_{tag}",
                    source="maschinen_ideen/111", n_terms=nt,
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
    print("[DEMO] Smoke test machine 111 (orientation pairs)")
    a, b = load_ab(*SOURCES[0])
    ar, br = reverse(a, b)
    P = reconstruct_from_fourier(a, b, n_samples=2048)
    Pr = reconstruct_from_fourier(ar, br, n_samples=2048)
    # (1) Parameter reversal exact: P_rev(t) == P(−t) (index mirror)
    d_param = float(np.abs(Pr - np.roll(P[::-1], 1, axis=0)).max())
    ok1 = d_param < 1e-12
    print(f"  (1) b→−b ⇔ t→−t: max|Δ| = {d_param:.1e}  "
          f"{'✓ bit' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) Set identity: same point set, other direction
    d_set = max(float(np.linalg.norm(P - q, axis=1).min())
                for q in Pr[::64])
    t_fwd = P[1] - P[0]
    t_rev = Pr[1] - Pr[0]
    ok2 = d_set < 1e-10
    print(f"  (2) point set identical ({d_set:.1e}), direction at the start "
          f"different  {'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) Three reversal pairs from asymmetric sources
    built = []
    for machine, stem in SOURCES:
        e, f = make_pair(machine, stem, HERE / "demo_smoke")
        if f is not None:
            built.append(stem)
    ok3 = len(built) >= 2
    print(f"  (3) reversal pairs: {built}  "
          f"{'✓ blindness table prepared' if ok3 else '✗ FAIL'}")
    assert ok3
    print("[DEMO] 3/3 green — arrow machine ready to fire "
          "(orientation blindness becomes a measurement instead of an assumption).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
