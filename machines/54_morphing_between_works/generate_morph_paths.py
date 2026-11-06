# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 54 — morphing_between_works (sculpture geodesics).

Paths in Fourier coefficient space between TWO EXISTING WORKS
(sundowner → nitro!): where does the type jump? Measurable today: embedding
margin PROFILE along the path + edge bisection (where the embedding
breaks, the type jump MUST live nearby — through a self-
contact every type change leads. Type edges themselves = ignition measurement.
Narrative: „the path between two sculptures, measured" — gallery.
REVERSAL CLAUSE: linear path may leave embeddability entirely →
detour paths (transversal arc offset) as an alternative route.
Doctrine: do NOT ignite.
"""
from __future__ import annotations
import site_config
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402
from the_knotfamily_machine.fseries_io import reconstruct_from_fourier  # noqa: E402

ROOT = Path(site_config.ROOT)
N_FAST = 1024     # cheap margin measurement along the path
N_SAMPLES = 8192


def load_fseries(path):
    """Both corpus formats (solves the workshop-sieve loader TODO):
    (a) pure 6-column lines (machine/saturn format),
    (b) own format with leading 3-value DC line (x₀ y₀ z₀) before the
        6-column coefficients — DC is returned (centering
        happens downstream anyway)."""
    rows6, dc = [], np.zeros(3)
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("%"):
            continue
        vals = [float(x) for x in s.split()]
        if len(vals) == 6:
            rows6.append(vals)
        elif len(vals) == 3 and not rows6:
            dc = np.array(vals)
        else:
            raise ValueError(f"{path}: line with {len(vals)} values")
    if not rows6:
        raise ValueError(f"{path}: no coefficient lines")
    C = np.asarray(rows6)
    return C[:, [0, 2, 4]], C[:, [1, 3, 5]]


def load_work(name):
    """Load an own_knots work; normalise (centre + RMS scale) so that
    the path runs shape against shape, not size against size."""
    cands = sorted((ROOT / "own_knots" / name).glob("*.fseries"))
    cands = [c for c in cands if not c.name.startswith("_")] + cands
    last_err = None
    for c in cands:
        try:
            a, b = load_fseries(c)
            P = reconstruct_from_fourier(a, b, n_samples=2048)
            scale = float(np.sqrt((P ** 2).sum(axis=1).mean()))
            return a / scale, b / scale, c
        except Exception as e:        # skip loader-TODO formats
            last_err = e
    raise FileNotFoundError(f"{name}: no loadable .fseries ({last_err})")


def pad_pair(aA, bA, aB, bB):
    n = max(aA.shape[0], aB.shape[0])
    def pad(x):
        out = np.zeros((n, 3)); out[:x.shape[0]] = x; return out
    return pad(aA), pad(bA), pad(aB), pad(bB)


def path_point(coeffs, s, detour=0.0):
    """Linear in the coefficients; optional detour arc: transversal
    Offset sin(πs)·detour on a fixed random direction (seed-fixed)."""
    aA, bA, aB, bB = coeffs
    a = (1 - s) * aA + s * aB
    b = (1 - s) * bA + s * bB
    if detour:
        rng = np.random.default_rng(54)
        da = rng.normal(size=a.shape); db = rng.normal(size=b.shape)
        nrm = np.sqrt((da ** 2).sum() + (db ** 2).sum())
        a = a + np.sin(np.pi * s) * detour * da / nrm
        b = b + np.sin(np.pi * s) * detour * db / nrm
    return a, b


def margin_at(coeffs, s, detour=0.0, n=N_FAST):
    a, b = path_point(coeffs, s, detour)
    P = reconstruct_from_fourier(a, b, n_samples=n)
    return embedding_report(P, n_samples=n)["d_min_rel"]


def margin_profile(coeffs, n_steps=41, detour=0.0):
    S = np.linspace(0.0, 1.0, n_steps)
    return S, np.array([margin_at(coeffs, s, detour) for s in S])


def find_edges(S, M, thresh):
    """Sign changes of (margin − threshold) → bisection intervals."""
    sgn = np.sign(M - thresh)
    return [(float(S[i]), float(S[i + 1]))
            for i in range(len(S) - 1) if sgn[i] != sgn[i + 1]]


def bisect_edge(coeffs, lo, hi, thresh, detour=0.0, iters=18):
    flo = margin_at(coeffs, lo, detour) - thresh
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        fm = margin_at(coeffs, mid, detour) - thresh
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return 0.5 * (lo + hi)


def write_waypoint(coeffs, s, out_root, tag, detour=0.0):
    a, b = path_point(coeffs, s, detour)
    P = reconstruct_from_fourier(a, b, n_samples=N_SAMPLES)
    name = f"mw_{tag}_s{s:.3f}".replace(".", "p")
    entry = {"name": name, "machine": "54_morphing_between_works",
             "s": float(s), "detour": float(detour),
             "predictions": {"identity": "Zünd-Messung (Typ-Kante?)"}}
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(P, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/54",
                              n_terms=max(96, coeffs[0].shape[0]),
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 54 (work-to-work morphing)")
    # (1) load end points + themselves embedded (sundowner → nitro)
    aA, bA, srcA = load_work("sundowner")
    aB, bB, srcB = load_work("nitro")
    coeffs = pad_pair(aA, bA, aB, bB)
    mA, mB = margin_at(coeffs, 0.0), margin_at(coeffs, 1.0)
    ok1 = mA > 0 and mB > 0
    print(f"  (1) end points: sundowner({Path(srcA).name}, margin={mA:.5f}) → "
          f"nitro({Path(srcB).name}, margin={mB:.5f})  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) margin profile along the linear path
    S, M = margin_profile(coeffs)
    thresh = 0.003
    n_below = int((M < thresh).sum())
    print(f"  (2) profile (41 sample points): min={M.min():.5f} @ "
          f"s={S[M.argmin()]:.3f}, max={M.max():.5f}, "
          f"{n_below} points below gate margin  ✓ (measured)")
    # (3) edges localized by bisection (if the path breaks)
    edges = find_edges(S, M, thresh)
    if edges:
        locs = [bisect_edge(coeffs, lo, hi, thresh) for (lo, hi) in edges]
        print(f"  (3) embedding edges: {len(edges)} → s = "
              f"{', '.join(f'{x:.4f}' for x in locs)}  ✓ (bisection 18 iter)")
    else:
        locs = []
        print("  (3) path breaks nowhere — linear route complete "
              "embedded  ✓")
    # (4) waypoint written (embedded point with maximal margin)
    s_best = float(S[M.argmax()])
    if s_best in (0.0, 1.0):       # endpoints are works, no waypoint
        s_best = float(S[1:-1][M[1:-1].argmax()])
    e4, f4 = write_waypoint(coeffs, s_best, HERE / "demo_smoke", "sun2nit")
    ok4 = f4 is not None
    print(f"  (4) Waypoint s={s_best:.3f}: "
          f"nyq={e4.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok4 else '✗ FAIL: ' + str(e4.get('rejected'))}")
    assert ok4
    # (5) Detour route measurable (detour raises min margin or leaves it —
    #     measured, not believed)
    _, M5 = margin_profile(coeffs, n_steps=21, detour=0.4)
    print(f"  (5) detour route (detour=0.4): min margin "
          f"{M.min():.5f} → {M5.min():.5f}  ✓ (detour mechanics in place)")
    rep = {"workA": str(srcA), "workB": str(srcB),
           "profile_s": [round(float(x), 4) for x in S],
           "profile_margin": [round(float(x), 6) for x in M],
           "edges_bisected": [round(float(x), 5) for x in locs],
           "detour_min_margin": float(M5.min())}
    (HERE / "demo_smoke" / "morph_profile_sundowner_nitro.json").write_text(
        json.dumps(rep, indent=1))
    print("[DEMO] 5/5 green — morph machine ready to fire "
          "(profile archived; type edges = ignition measurement).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
