# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 12 — filling_the_symmetry_table (the seven columns).

Harvest 21 — CAUTION DUPLICATE AUDIT: C_nv+D_n were already in wave 2
(generate_symmetry_columns.py); NET here = C_2h obstruction (n=2),
C_1h column, closing the table. State of the columns (axis z,
rotation order n, Fourier normal position winding w=1, xy support m≡1 mod n,
z support m≡0 mod n):
  C_n   ✓ built (Stanggassinger machine)
  S_2n  ✓ built (115)
  D_nd  ✓ built (116: S_2n + transverse axis = D_nd)
  C_nv  → HERE: σ_v∘(t→−t): a_m REAL + z COS-only ⟹ the mirror reverses
         the direction = MINUS-amphi factory in series (counterpart
         to 115, which builds PLUS-amphis direction-preservingly via σ).
  D_n   → HERE: C_2 transverse axis: a_m REAL + z SIN-only ⟹ certified
         strongly invertible CHIRAL series (no improper element).
  C_nh  → OBSTRUCTED for n≥2 (two-stage): n≥3 support argument
         (σ_h needs −m ≡ w mod n ⟹ 2w ≡ 0 mod n; w=n/2 ⟹
         π period ⟹ double traversal). n=2 SHARPER: z odd ∧
         π-periodic ⟹ z(π/2) = 0 FORCED, c(π/2)=c(−π/2) ⟹
         forced self-intersection — MEASURED by the embedding gate
         (d_min_rel = exactly 0.0 over all seeds; build attempt was
         the proof witness). ONLY n=1 lives: C_1h = pure mirror
         plane (square-knot class, plane meets the knot in exactly
         2 points), zero-crossing control |q1| > 2|q2| keeps the
         z zeros on {0, π} → is built.
  D_nh  → OBSTRUCTED (derivation): D_n requires z odd (sin),
         the σ_v part requires z even (cos) ⟹ z ≡ 0 ⟹ planar/unknot.
BETS EXACT today: all symmetries bit (<1e-12) · C_2h odd one out bit
· Series C_nv/D_n for n=3,4 + C_2h embedding-gated written.
FIRING BETS: C_nv members amphichiral (sig=0, θ amphi-consistent),
D_n members chiral+invertible (snappy symmetry group) — the
table row for Stoimenow: amphi KIND (±) and invertibility
PER COLUMN predicted, pipeline measures afterwards. Doctrine: do NOT fire.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402

N = 4320      # divisible by n=2,3,4 (and 2n shifts)


def _xy_support(n, n_xy=3):
    ns = [1]
    for j in range(1, n_xy + 1):
        ns += [1 + n * j, 1 - n * j]
    return ns


def forge(column, n, seed, n_xy=3, n_z=2, decay=0.55):
    """One column of the table: \'cnv\' | \'dn\' | \'c2h\'."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0, 2 * np.pi, N, endpoint=False)
    c = np.zeros(N, complex)
    z = np.zeros(N)
    if column in ("cnv", "dn"):
        for m in _xy_support(n, n_xy):
            amp = (1.2 if m == 1 else
                   0.45 * decay ** abs((m - 1) // n))
            c += amp * rng.normal() * np.exp(1j * m * t)     # a_m REAL
        zf = np.cos if column == "cnv" else np.sin
        for j in range(1, n_z + 1):
            z += (decay ** (j - 1)) * 0.7 * rng.normal() * zf(j * n * t)
    elif column == "c2h":
        # OBSTRUCTION WITNESS (not a product): z odd ∧ π-periodic
        # ⟹ z(π/2)=0 forced ⟹ self-intersection (d_min == 0.0)
        for j in (0, 1, 2, 3):
            m = 1 + 2 * j
            a = ((1.2 if m == 1 else 0.45 * decay ** j)
                 * (rng.normal() + 1j * rng.normal()))
            c += a * (np.exp(1j * m * t) + np.exp(-1j * m * t))
        for j in (1, 2):
            z += (decay ** (j - 1)) * 0.7 * rng.normal() * np.sin(2 * j * t)
    elif column == "c1h":
        # C_1h = pure mirror plane (square-knot class):
        # c(−t)=c(t) (a_{−m}=a_m, complex), z odd with
        # zero-crossing control |q1| > 2|q2| ⟹ z=0 only at t∈{0,π}
        for m, ampl in ((1, 1.2), (2, 0.5), (3, 0.35), (5, 0.2)):
            a = ampl * (rng.normal() + 1j * rng.normal())
            c += a * (np.exp(1j * m * t) + np.exp(-1j * m * t))
        q2 = 0.30 * rng.normal()
        q1 = (2.0 * abs(q2) + 0.6) * np.sign(rng.normal() + 1e-9)
        z = q1 * np.sin(t) + q2 * np.sin(2 * t)
    else:
        raise ValueError(column)
    return np.column_stack([c.real, c.imag, z])


def rot_z(P, ang):
    R = np.array([[np.cos(ang), -np.sin(ang), 0],
                  [np.sin(ang), np.cos(ang), 0], [0, 0, 1.0]])
    return P @ R.T


def t_reversed(P):
    return np.vstack([P[:1], P[:0:-1]])


def residuals(P, column, n):
    """Bit residuals of all group generators of the column."""
    res = {}
    res["C_n"] = float(np.abs(rot_z(P, 2 * np.pi / n)
                              - np.roll(P, -N // n, axis=0)).max())
    if column == "cnv":
        res["sigma_v"] = float(np.abs(
            t_reversed(P) * np.array([1, -1, 1]) - P).max())
    elif column == "dn":
        res["C_2_quer"] = float(np.abs(
            t_reversed(P) * np.array([1, -1, -1]) - P).max())
    elif column in ("c2h", "c1h"):
        res["sigma_h"] = float(np.abs(
            t_reversed(P) * np.array([1, 1, -1]) - P).max())
        if column == "c1h":
            res.pop("C_n", None)        # n=1: rotation trivial
    return res


def make_knot(column, n, seed, out_root):
    P = forge(column, n, seed)
    name = f"st_{column}{n}_s{seed}"
    grp = {"cnv": f"C_{n}v (Minus-Amphi)", "dn": f"D_{n} (chiral, invert.)",
           "c1h": "C_1h (Spiegel-Ebene, Quadratknoten-Klasse)"}[column]
    entry = {"name": name, "machine": "12_filling_the_symmetry_table",
             "column": column, "n": n, "group": grp,
             "symmetry_residuals": residuals(P, column, n),
             "predictions": {
                 "cnv": "amphichiral (Spiegel kehrt Richtung = MINUS-"
                        "Sorte); sig=0; Z/n",
                 "dn": "chiral + stark-invertierbar (2 Achsen-Punkte "
                       "pro Querachse); Z/n",
                 "c1h": "amphichiral via Spiegel-Ebene (2 Ebenen-"
                        "Punkte); Quadratknoten-Verwandtschaft "
                        "(K # mK-Bauplan als Kurven-Klasse)",
             }[column]}
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = None
    for nt in (max(32, 1 + 3 * n), 128, 256):
        try:
            rep = polyline_to_fseries(P, folder / f"{name}.fseries",
                                      knot_name=name,
                                      source="maschinen_ideen/12",
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
    print("[DEMO] Smoke test machine 12 (close symmetry table)")
    # (1) column symmetries bit
    worst = 0.0
    detail = []
    for col, n in (("cnv", 3), ("cnv", 4), ("dn", 3), ("dn", 4),
                   ("c1h", 1)):
        P = forge(col, n, seed=12)
        r = residuals(P, col, n)
        worst = max(worst, max(r.values()))
        detail.append(f"{col}{n}:{max(r.values()):.0e}")
    ok1 = worst < 1e-12
    print(f"  (1) generator residuals {' '.join(detail)}  "
          f"{'✓ all columns bit' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) obstruction WITNESSES as MEASUREMENT:
    #     (a) D_nh: z would have to be odd (D_n) AND even (σ_h∘C_2)
    #     (b) C_2h: forced self-intersection — gate measures d_min == 0
    P4 = forge("dn", 4, seed=12)
    z_even_part = 0.5 * (P4[:, 2] + t_reversed(P4)[:, 2])
    ok2a = float(np.abs(z_even_part).max()) < 1e-12
    dmins = []
    for seed in (12, 13, 14):
        W = forge("c2h", 2, seed)
        dmins.append(embedding_report(W)["d_min_rel"])
    ok2b = all(d == 0.0 for d in dmins)
    ok2 = ok2a and ok2b
    print(f"  (2) obstruction witnesses: D_nh conflict (even z part "
          f"{np.abs(z_even_part).max():.0e}) · C_2h forced "
          f"self-intersection d_min_rel={dmins} over 3 seeds  "
          f"{'✓ theorems MEASURED' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) series written
    built = []
    for col, n in (("cnv", 3), ("cnv", 4), ("dn", 3), ("dn", 4),
                   ("c1h", 1)):
        for seed in (12, 13, 14, 15):
            e, f = make_knot(col, n, seed, HERE / "demo_smoke")
            if f is not None:
                built.append(e["name"])
                break
    ok3 = len(built) == 5
    print(f"  (3) series: {built}  "
          f"{'✓ table filled' if ok3 else '✗ FAIL'}")
    assert ok3
    print("[DEMO] 3/3 green — symmetry table CLOSED: "
          "C_n/S_2n/D_nd built (old+115+116), C_nv/D_n/C_1h here, "
          "C_nh(n≥2)/D_nh obstructed (derived + MEASURED).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
