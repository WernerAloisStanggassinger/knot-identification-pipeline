# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 10 — climbing_the_torus_ladder ((n,n+1) ladder, Werner's design 06-04).

One more corner per rung. Three components:

1. EXACT RUNG .fseries: T(p,q) is finitely harmonic —
     x = R·cos(pt) + r/2·cos((q+p)t) + r/2·cos((q−p)t)
     y = R·sin(pt) + r/2·sin((q+p)t) − r/2·sin((q−p)t)
     z = r·sin(qt)
   Only {p, q−p, q+p} in xy + {q} in z. For (n,n+1) q−p=1 — exactly the
   MEASURED Z7 selection rule of the kf286 addendum (k≡±1 mod 7 in xy,
   k≡0 mod 7 in z). Coefficients ANALYTIC (machine-01 principle).

2. BET PROTOCOL FROM FORMULAS (not typed in): Δ_T(p,q) by EXACT
   integer polynomial division (t^{pq}−1)(t−1)/((t^p−1)(t^q−1)) →
   det=|Δ(−1)|, genus=deg/2, KFH-rank bet=#terms; c=q(p−1) [Murasugi],
   cusp=p·q, CS=−1/(2pq), |s|=2g, vol=0, L-space, fibred.
   Bit-frozen with date BEFORE every run (blind-test hygiene).

3. BRAID SECOND WITNESS: (σ1…σ_{p−1})^q as a word file (certificate method
   Disc #255 — conjugation proof is established).

Doctrine: do NOT fire before Monday. Neighbourhood sweeps = compute plan in
PREPARATION (wall doctrine, small waves). Seven-corner design (Werner,
physical) = T(7,8) blind-test candidate number one.
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
sys.path.insert(0, str(HERE.parent / "the_knotfamily_machine"))
from fseries_io import write_fseries, reconstruct_from_fourier  # noqa: E402
from machine_common import min_nonadjacent_dist, resample_uniform  # noqa: E402

R_MAJOR, R_MINOR = 1.0, 0.45


# ------------------------------------------------ Exact polynomial arithmetic

def poly_mul(a, b):
    return np.convolve(np.asarray(a, object), np.asarray(b, object))


def poly_divmod_exact(num, den):
    """Exact division of integer polynomials (remainder must be 0)."""
    num = list(np.asarray(num, object))
    den = list(np.asarray(den, object))
    while den and den[-1] == 0:
        den.pop()
    q = [0] * (len(num) - len(den) + 1)
    r = num[:]
    for k in range(len(q) - 1, -1, -1):
        c = r[k + len(den) - 1]
        if c % den[-1] != 0:
            raise RuntimeError("division not exact — formula bug")
        q[k] = c // den[-1]
        for j, d in enumerate(den):
            r[k + j] -= q[k] * d
    if any(r):
        raise RuntimeError("remainder ≠ 0 — division not exact")
    return q


def t_power_minus_1(n):
    p = [0] * (n + 1)
    p[0], p[n] = -1, 1
    return p


def alexander_torus(p, q):
    """Δ_T(p,q) = (t^{pq}−1)(t−1) / ((t^p−1)(t^q−1)) — exact."""
    num = poly_mul(t_power_minus_1(p * q), t_power_minus_1(1))
    den = poly_mul(t_power_minus_1(p), t_power_minus_1(q))
    return poly_divmod_exact(num, den)


def rung_predictions(p, q):
    delta = alexander_torus(p, q)
    det = abs(sum(c * (-1) ** k for k, c in enumerate(delta)))
    genus = (len(delta) - 1) // 2
    return {
        "knot": f"T({p},{q})",
        "crossing_number": q * (p - 1),
        "det": int(det),
        "genus": genus,
        "signature_abs": 2 * genus,
        "alexander_n_terms": int(sum(1 for c in delta if c != 0)),
        "kfh_total_rank_bet": int(sum(1 for c in delta if c != 0)),
        "cusp_volume": p * q,
        "chern_simons": f"-1/{2 * p * q}",
        "hyperbolic_volume": 0.0,
        "l_space": True, "fibered": True, "positive_braid": True,
        "braid": f"{p} Straenge x {q*(p-1)} Kreuzungen: (s1...s{p-1})^{q}",
    }


# ------------------------------------------------------- Exact rung curve

def torus_coeffs(p, q, R=R_MAJOR, r=R_MINOR, n_terms=None):
    """Analytic Fourier coefficients of T(p,q) (finitely harmonic)."""
    n_terms = max(20, q + p) if n_terms is None else n_terms
    a = np.zeros((n_terms, 3))
    b = np.zeros((n_terms, 3))
    a[p - 1, 0] = R                  # x: R·cos(pt)
    a[q + p - 1, 0] = r / 2          #    r/2·cos((q+p)t)
    a[q - p - 1, 0] += r / 2         #    r/2·cos((q−p)t)
    b[p - 1, 1] = R                  # y: R·sin(pt)
    b[q + p - 1, 1] = r / 2          #    r/2·sin((q+p)t)
    b[q - p - 1, 1] += -r / 2        #    −r/2·sin((q−p)t)
    b[q - 1, 2] = r                  # z: r·sin(qt)
    return a, b


def make_rung(p, q, out_root):
    if math.gcd(p, q) != 1:
        raise ValueError("p,q must be coprime")
    name = f"ladder_T{p}_{q}"
    a, b = torus_coeffs(p, q)
    P = reconstruct_from_fourier(a, b, n_samples=8192)
    d_min, _ = min_nonadjacent_dist(resample_uniform(P, 4096),
                                    exclude=8)
    diam = float(np.linalg.norm(P.max(0) - P.min(0)))
    entry = {"name": name, "machine": "10_climbing_the_torus_ladder",
             "predictions": rung_predictions(p, q),
             "coefficients_exact": True,
             "harmonics": {"xy": [p, q - p, q + p], "z": [q]},
             "embedding": {"d_min_rel": round(d_min / diam, 6),
                           "embedded": bool(d_min / diam >= 0.003)}}
    if not entry["embedding"]["embedded"]:
        entry["rejected"] = "embedding gate (r/R zu groß?)"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    write_fseries(folder / f"{name}.fseries", a, b, knot_name=name,
                  source="maschinen_ideen/10_climbing_the_torus_ladder")
    word = " ".join(f"s{i}" for i in range(1, p)) + f" ^{q}"
    (folder / "braid_certificate_word.txt").write_text(
        f"T({p},{q}) Braid-Zweitzeuge: ({' '.join(f's{i}' for i in range(1, p))})^{q}\n"
        f"= {p} Straenge, {q * (p - 1)} Kreuzungen, alle gleichsinnig.\n"
        f"Identitaets-Beweis beim Zuenden: Konjugations-Zertifikat (Disc #255).\n")
    return entry, folder


# ----------------------------------------------------------------------- CLI

def cmd_protokoll(args):
    """Bit-freeze the bet protocol (blind-test hygiene: date, before the run)."""
    rungs = [rung_predictions(n, n + 1) for n in range(6, args.max_n + 1)]
    dets = [r["det"] for r in rungs]
    out = {"frozen_at": args.date, "doctrine":
           "Vorhersagen UNVERÄNDERLICH ab frozen_at; Vergleich nur bit-genau "
           "gegen Pipeline-Output. Quelle: exakte Formeln (Δ-Division, "
           "Murasugi-c, cusp=pq, CS=-1/2pq), NICHT Tabellen-Lookup.",
           "rungs": rungs,
           "det_pairs_note": f"det-Folge {dets} — Paar-Struktur "
           "(Zwillinge) ist Teil der Wette."}
    p = HERE / f"WETT_PROTOKOLL_{args.date}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"[PROTOKOLL] {len(rungs)} rungs T(6,7)..T({args.max_n},"
          f"{args.max_n+1}) bit-frozen → {p.name}")
    for r in rungs:
        print(f"   {r['knot']:9s} c={r['crossing_number']:3d} det={r['det']:2d} "
              f"g={r['genus']:2d} cusp={r['cusp_volume']:3d} CS={r['chern_simons']}")


def cmd_make(args):
    entry, folder = make_rung(args.p, args.q, HERE / "out_fseries")
    man_path = HERE / "MANIFEST.json"
    man = json.loads(man_path.read_text()) if man_path.exists() else {"knots": []}
    man["knots"] = [k for k in man["knots"] if k["name"] != entry["name"]] + [entry]
    man_path.write_text(json.dumps(man, indent=1))
    print(f"[MAKE] {entry['name']} → {folder if folder else entry['rejected']}")


def cmd_demo(args):
    print("[DEMO] smoke test machine 10 (torus ladder)")
    # 1) Δ division exact: T(2,3) = t²−t+1, det 3
    d23 = alexander_torus(2, 3)
    ok1 = d23 == [1, -1, 1]
    print(f"  (1) Δ_T(2,3) = {d23} (expected [1,-1,1]), det="
          f"{abs(sum(c*(-1)**k for k, c in enumerate(d23)))}  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # 2) T(6,7) against the founding-txt line (c=35, det=7, g=15, 11 terms)
    r67 = rung_predictions(6, 7)
    ok2 = (r67["crossing_number"] == 35 and r67["det"] == 7
           and r67["genus"] == 15 and r67["alexander_n_terms"] == 11)
    print(f"  (2) T(6,7): c={r67['crossing_number']} det={r67['det']} "
          f"g={r67['genus']} Δ terms={r67['alexander_n_terms']} "
          f"(kf286 witnesses: 35/7/15/11)  {'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # 3) det twins of the ladder 6..11
    dets = [rung_predictions(n, n + 1)["det"] for n in range(6, 11)]
    ok3 = dets == [7, 7, 9, 9, 11]
    print(f"  (3) det ladder 6..10: {dets} (twin bet 7,7,9,9,11)  "
          f"{'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    # 4) Exactness of the curve coefficients (FFT re-measurement) + z selection rule
    a, b = torus_coeffs(7, 8)
    P = reconstruct_from_fourier(a, b, n_samples=4096)
    spec = np.fft.rfft(P, axis=0)
    a_fft = (2.0 / 4096) * np.real(spec[1:a.shape[0] + 1])
    err = float(np.abs(a_fft - a).max())
    amp = np.hypot(a, b)
    active_xy = sorted(set((np.nonzero(amp[:, :2].sum(axis=1) > 1e-12)[0] + 1)))
    active_z = sorted(set((np.nonzero(amp[:, 2] > 1e-12)[0] + 1)))
    ok4 = err < 1e-12 and active_xy == [1, 7, 15] and active_z == [8]
    print(f"  (4) T(7,8) exact: FFT-max|Δ|={err:.1e}, xy harmonics "
          f"{active_xy} (expected [1,7,15] = k≡±1 mod 8), z {active_z} (expected [8])"
          f"  {'✓' if ok4 else '✗ FAIL'}")
    assert ok4
    # 5) Build a rung (T(7,8), the seven-corner blind-test partner)
    e, f = make_rung(7, 8, HERE / "demo_smoke")
    ok5 = f is not None
    print(f"  (5) ladder_T7_8: d_min_rel={e['embedding']['d_min_rel']} → "
          f"fseries + braid certificate word  "
          f"{'✓' if ok5 else '✗ FAIL: ' + str(e.get('rejected'))}")
    assert ok5
    print("[DEMO] 5/5 green — ladder ready to fire (firing after Monday; "
          "neighbourhoods = compute plan/wall doctrine).")


def main():
    ap = argparse.ArgumentParser(description="Machine 10 — torus ladder")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("protokoll")
    p1.add_argument("--max-n", type=int, default=11)
    p1.add_argument("--date", required=True)
    p2 = sub.add_parser("make")
    p2.add_argument("p", type=int)
    p2.add_argument("q", type=int)
    sub.add_parser("demo")
    args = ap.parse_args()
    {"protokoll": cmd_protokoll, "make": cmd_make, "demo": cmd_demo}[args.cmd](args)


if __name__ == "__main__":
    main()
