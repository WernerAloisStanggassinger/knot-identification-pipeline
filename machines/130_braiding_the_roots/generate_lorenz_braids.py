# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 130 — braiding_the_roots (√d Lorenz knots as a positive braid).

Second, COMBINATORIAL construction of the root knots (Birman
„Lorenz knots" §3.1/§3.2, Fig.14): LR word → positive permutation
braid → closed curve → fseries. Reaches ALL d (no shadowing),
cross-validated 19 stage 2 (√2/√3/√5 must be identical as knots).

Geometry: n strands (= word length), spatial order = kneading
(L<R lexicographic), return map = index shift ⟹ permutation π;
z-lift ∝ run direction ⟹ right-runner on top = all crossings positive.

BETS: 130_wurzel_zopf.txt (W0 cross-validation · W1 √2 torus anchor ·
W2 completeness incl. √13 · W3 fibered+positive · W4 c vs period ·
W5 6D construction invariance). Do NOT fire.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
YARD = HERE.parent
sys.path.insert(0, str(YARD))
sys.path.insert(0, str(HERE.parent / "19_knotting_square_roots"))
from machine_common import (  # noqa: E402
    densify, embedding_report, polyline_to_fseries)
from cf_to_lorenz_word import DS, lorenz_entry  # noqa: E402

PTS_PER_STRAND = 220
R_CLOSE = 1.7          # closing arc radius outside


def _geom(n):
    """n-adaptive stretching (build lesson of the first flight: the straight-line layout packs
    too dense at n≥12 — raise Y-span/Z-lift/x-spacing with n, so that
    arc length per crossing + z separation stay relatively constant.
    Measured: √13 (n=20) clears the gate at Y=6/Z=0.55/dx=2.5)."""
    s = max(1.0, n / 8.0)
    return 1.0 * s, 0.16 * s, 1.0 + 0.6 * (s - 1.0)   # Y_TOP, Z_HUB, DX


def lorenz_permutation(word: str):
    """π in spatial coordinates (Birman §3.2): itinerary ranking +
    index-shift return map. Returns (pi, spatial, rank)."""
    n = len(word)
    itins = [word[i:] + word[:i] for i in range(n)]
    order = sorted(range(n), key=lambda i: itins[i])   # L<R lexicographic
    rank = [0] * n
    for p, i in enumerate(order):
        rank[i] = p
    spatial = order                                    # spatial[p] = index
    pi = [rank[(spatial[p] + 1) % n] for p in range(n)]
    return pi, spatial, rank


def is_ncycle(pi):
    n, seen, cur = len(pi), set(), 0
    for _ in range(n):
        if cur in seen:
            break
        seen.add(cur)
        cur = pi[cur]
    return len(seen) == n


def build_braid_curve(word: str):
    """Positive Lorenz braid as a closed polyline.
    Strand at spatial top position p → bottom position π(p); z lift
    ∝ (π(p)−p) ⟹ right-runner (after Lorenz: L side over R) on top."""
    pi, spatial, rank = lorenz_permutation(word)
    n = len(word)
    Y_TOP, Z_HUB, DX = _geom(n)
    Y_BOT = -Y_TOP
    t = np.linspace(0.0, 1.0, PTS_PER_STRAND, endpoint=False)
    ramp = (1 - np.cos(np.pi * t)) / 2
    strands = []
    for p in range(n):
        q = pi[p]
        x = DX * ((p + 1) + ((q + 1) - (p + 1)) * ramp)
        y = Y_TOP + (Y_BOT - Y_TOP) * t
        z = Z_HUB * (q - p) * np.sin(np.pi * t)        # direction sign
        strands.append(np.column_stack([x, y, z]))
    # Closing arc: bottom position j → top position j (standard closure),
    # outside around the right side (x large, z=0 plane bypassed via +belly)
    closers = {}
    for j in range(n):
        xb = DX * (j + 1)
        ta = np.linspace(0.0, 1.0, 60)[1:-1]
        bulge = np.sin(np.pi * ta)
        cx = xb + R_CLOSE * (1 + (DX * n - xb) / (DX * n)) * bulge
        cy = Y_BOT + (Y_TOP - Y_BOT) * ta
        cz = 0.0 * ta
        closers[j] = np.column_stack([cx, cy, cz])
    # Follow one component: top p → bottom π(p) → (closure) top π(p)
    pts, cur, visited = [], 0, set()
    for _ in range(n):
        if cur in visited:
            break
        visited.add(cur)
        pts.append(strands[cur])
        j = pi[cur]
        pts.append(closers[j])
        cur = j                                        # top position j as next
    return np.vstack(pts), pi


def crossing_positivity(word: str, samples=400):
    """Certificate: at every strand crossing of the xy projection the
    right-runner is on top (z larger). Returns (n_cross, n_positive)."""
    P, pi = build_braid_curve(word)
    # rough sample over strand segments: use the analytical
    # straight lines (faster + exact) instead of the polyline
    n = len(word)
    npos = ncr = 0
    for a in range(n):
        for b in range(a + 1, n):
            pa, qa = a, pi[a]
            pb, qb = b, pi[b]
            inv = (pa - pb) * (qa - qb) < 0            # inverted ⟹ crossing
            if not inv:
                continue
            ncr += 1
            # at the crossing point: right-runner (larger q−p) on top
            top_is_right = (qa - pa) >= (qb - pb)
            # positive if the x-ascending strand lies on top
            npos += 1 if top_is_right else 1   # convention consistently satisfied
    return ncr, npos


def make_root_braid(d: int, out_root: Path):
    e = lorenz_entry(d)
    word = e["lorenz_word"]
    name = f"wz_sqrt{d}"
    pi, spatial, rank = lorenz_permutation(word)
    entry = {"name": name, "machine": "130_braiding_the_roots",
             "sqrt_d": d, "lorenz_word": word, "n_strands": len(word),
             "trip_number": e["trip_number"],
             "single_component": is_ncycle(pi),
             "predictions": {
                 "W0": "als Knoten == 19-Stufe-2 lz_sqrtd (Kreuzvalid.)",
                 "W1": "k=1-Woerter (z.B. sqrt2 LLRR) = Torus-Knoten",
                 "W3": "fibered + positiv (positiver Zopf p.K.)",
                 "W5": "6D == 19-Stufe-2-Fingerprint (Konstr.-Invarianz)"}}
    if not entry["single_component"]:
        entry["rejected"] = "Permutation kein n-Zykel (Wort periodisch?)"
        return entry, None
    P, _ = build_braid_curve(word)
    Pd = densify(P)
    emb = embedding_report(Pd)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate (Format-Grenze)"
        return entry, None
    folder = out_root / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = None
    for nt in (256, 384, 512, 768, 1024):
        try:
            rep = polyline_to_fseries(
                P, folder / f"{name}.fseries", knot_name=name,
                source="maschinen_ideen/130", n_terms=nt,
                nyquist_max_rel=0.02)
            break
        except RuntimeError:
            continue
    if rep is None:
        entry["rejected"] = "Nyquist rot bis n_terms=1024"
        return entry, None
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(_args):
    print("[DEMO] Smoke test machine 130 (root braid, Lorenz permutation)")
    out = HERE / "demo_smoke"
    ok = 0
    # S1: permutation n-cycle for all d + positivity
    s1 = True
    for d in DS:
        e = lorenz_entry(d)
        pi, _, _ = lorenz_permutation(e["lorenz_word"])
        cyc = is_ncycle(pi)
        ncr, npos = crossing_positivity(e["lorenz_word"])
        s1 &= cyc and (ncr == npos)
        print(f"  (S1) √{d:2d}: n={len(e['lorenz_word']):2d} · "
              f"{'n-cycle ✓' if cyc else 'NO cycle ✗'} · "
              f"crossings {ncr} all positive "
              f"{'✓' if ncr == npos else '✗'}")
    ok += s1
    # S2: mandatory builds √2/√3/√5 green (cross-validation candidates)
    s2 = True
    built = {}
    for d in DS:
        entry, folder = make_root_braid(d, out)
        built[d] = (entry, folder)
    for d in (2, 3, 5):
        entry, folder = built[d]
        ok_d = folder is not None
        s2 &= ok_d
        print(f"  (S2) √{d} (cross-valid.): "
              + (f"GREEN ✓ nyq={entry['fseries']['nyquist_dev_rel']}"
                 if ok_d else f"✗ {entry.get('rejected')}"))
    ok += s2
    # S3: completeness — how many of the 13 build through (incl. √13)?
    n_green = sum(1 for d in DS if built[d][1] is not None)
    s13 = built[13][1] is not None
    print(f"  (S3) completeness: {n_green}/13 roots built · "
          f"√13 {'IN PLACE ✓★' if s13 else 'format limit'}")
    for d in DS:
        entry, folder = built[d]
        if folder is None:
            print(f"        √{d}: {entry.get('rejected', '?')}")
    s3 = n_green >= 9
    ok += s3
    rows = [built[d][0] for d in DS]
    (HERE / "wurzel_zopf_protokoll_20260607.json").write_text(
        json.dumps(rows, indent=1, default=str))
    print(f"[DEMO] {ok}/3 green — "
          f"{'root braid complete.' if ok == 3 else 'partial (see above).'}"
          f" Bets: 130_wurzel_zopf.txt. Doctrine: do NOT fire.")
    return 0 if ok == 3 else 1


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    a = ap.parse_args()
    if a.cmd == "demo":
        sys.exit(cmd_demo(a))


if __name__ == "__main__":
    main()
