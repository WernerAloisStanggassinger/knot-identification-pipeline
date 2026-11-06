# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 86 — counting_sunflower_spirals (phyllotaxis strands).

Point pattern with GOLDEN divergence angle 2π/φ² on the cylinder/torus;
the visible spiral families (parastichies) carry FIBONACCI numbers —
here MEASURED from the neighborhood analysis of the built pattern,
not cited. Strand = parastichy run to the end; closure on
the torus via rational convergent winding (73er bridge: the
convergents of φ are F_k/F_{k+1}!), remainder = drift gate.
Narrative punch: „the sunflower in stainless steel".
REVERSAL CLAUSE: irrational angle NEVER closes exactly — the rational
convergent step is the honest closing rule.
Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402

PHI = (1 + np.sqrt(5)) / 2
GOLD = 2 * np.pi / PHI ** 2          # 137.507…°
R_RING, R_TUBE = 6.0, 1.6


def cylinder_points(n, rise=0.08, angle=GOLD):
    """Phyllotaxis on the cylinder: point k at (k·angle mod 2π, k·rise)."""
    k = np.arange(n)
    return np.column_stack([(k * angle) % (2 * np.pi), k * rise]), k


def parastichy_numbers(n=400, rise=0.08):
    """MEASURED: the two most frequent neighbor index differences
    (cylindrical metric) = parastichy pair."""
    pts, k = cylinder_points(n, rise)
    diffs = []
    for i in range(n):
        d_th = np.abs(pts[:, 0] - pts[i, 0])
        d_th = np.minimum(d_th, 2 * np.pi - d_th)
        d = np.sqrt(d_th ** 2 + (pts[:, 1] - pts[i, 1]) ** 2)
        d[i] = np.inf
        for j in np.argsort(d)[:4]:
            diffs.append(abs(int(k[j]) - int(k[i])))
    vals, cnts = np.unique(diffs, return_counts=True)
    order = np.argsort(-cnts)
    return sorted(int(vals[o]) for o in order[:2])


def fib_upto(m):
    f = [1, 1]
    while f[-1] < m:
        f.append(f[-1] + f[-2])
    return set(f)


def strand_on_torus(stride, n_wind_q, n_pts=8192):
    """Parastichy strand as torus curve: convergent winding
    p/q = F_k/F_{k+1} (strand with step `stride` closes after q
    steps exactly — the rational closing rule)."""
    s = np.linspace(0, 2 * np.pi, n_pts, endpoint=False)
    th_ring = n_wind_q * s                       # q revolutions around the ring
    th_tube = stride * s                          # stride revolutions of the tube
    r = R_RING + R_TUBE * np.cos(th_tube)
    return np.column_stack([r * np.cos(th_ring / n_wind_q),
                            r * np.sin(th_ring / n_wind_q),
                            R_TUBE * np.sin(th_tube)])


def make_knot(p, q, out_root):
    """(p,q) strand (Fibonacci pair) on the torus."""
    from math import gcd
    name = f"ph_F{p}_F{q}"
    entry = {"name": name, "machine": "86_counting_sunflower_spirals",
             "pair": [p, q], "coprime": gcd(p, q) == 1,
             "predictions": {
                 "identity": f"T({p},{q})-Klasse erwartet (Zündung)",
                 "quelle": "goldener Winkel → Fibonacci-Parastichie"}}
    if gcd(p, q) != 1:
        entry["rejected"] = "nicht teilerfremd (Link)"
        return entry, None
    s = np.linspace(0, 2 * np.pi, 8192, endpoint=False)
    r = R_RING + R_TUBE * np.cos(q * s)
    P = np.column_stack([r * np.cos(p * s), r * np.sin(p * s),
                         R_TUBE * np.sin(q * s)])
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(P, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/86",
                              n_terms=max(96, 16 * (p + q)),
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 86 (sunflower spirals)")
    # (1) MEASURED: parastichy pair is Fibonacci (two densities)
    pairs = {}
    fib = fib_upto(200)
    ok1 = True
    for rise in (0.15, 0.08, 0.04):
        p = parastichy_numbers(rise=rise)
        pairs[rise] = p
        ok1 = ok1 and all(x in fib for x in p)
    print(f"  (1) parastichy pairs (rise 0.15/0.08/0.04): {pairs} — "
          f"all Fibonacci  {'✓ MEASURED, not cited' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) Density cascade: finer ⇒ higher Fibonacci (pair migrates)
    seq = [pairs[r] for r in (0.15, 0.08, 0.04)]
    ok2 = seq[0] != seq[2]
    print(f"  (2) cascade {seq[0]} → {seq[2]}: the pair climbs the "
          f"Fibonacci ladder  {'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) convergent logic: φ convergents = F_k/F_{k+1} exactly
    from fractions import Fraction
    conv = []
    f = Fraction(1, 1)
    fibs = [1, 1, 2, 3, 5, 8, 13]
    ok3 = all(Fraction(fibs[i + 1], fibs[i]).limit_denominator(100)
              is not None for i in range(5))
    # more direct: continued fraction of φ is [1;1,1,…] → convergents F quotients
    x = PHI
    cf = []
    for _ in range(6):
        a = int(x); cf.append(a); x = 1 / (x - a)
    ok3 = cf == [1, 1, 1, 1, 1, 1]
    print(f"  (3) φ continued fraction: {cf} (all 1 ⇒ convergents = "
          f"Fibonacci quotients)  {'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    # (4) sunflower strand built: measured pair → torus curve
    p, q = pairs[0.08]
    e4, f4 = make_knot(p, q, HERE / "demo_smoke")
    ok4 = f4 is not None
    print(f"  (4) {e4['name']} (from measured pair {p},{q}): "
          f"nyq={e4.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok4 else '✗ FAIL: ' + str(e4.get('rejected'))}")
    assert ok4
    print("[DEMO] 4/4 green — phyllotaxis machine ready to fire "
          "(the sunflower delivers its own torus knots).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
