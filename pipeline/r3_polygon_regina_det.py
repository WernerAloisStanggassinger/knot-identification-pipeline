#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""
r3_polygon_regina_det.py — THIRD, completely independent tool witness for the
poly_cN dets (Disc #244). Builds the PD code of the sharp control polygon P
directly from the 3D geometry (same proven pd_from_polygon logic as the
Sage witness, c04 -> 2535) and computes the determinant with regina's native
Alexander routine (C++, NO timeout at ~99 crossings — exactly the wall at
which Sage's Link.determinant() died at c05/c09/s08).

RUNS IN base: conda run -n base python r3_polygon_regina_det.py
(base = pyknotid + regina; see memory reference-conda-env-tool-map).

det = |Alexander(-1)|, mirror-invariant. Validates itself first on the trefoil
(PD convention regina.fromPD vs Sage must agree -> det 3).
Pure diagnostics. Writes NOTHING to the corpus.
"""
import site_config
import sys, json, argparse
from collections import Counter
from pathlib import Path
import numpy as np
import regina

ROOT = Path(site_config.ROOT)
sys.path.insert(0, str(Path(site_config.MACHINE_YARD) / "the_knotfamily_machine"))   # authors' ROOT/the_knotfamily_machine is a symlink to this folder
from knot_machine import build_polygon_symmetric, polygon_generator  # noqa
from make_family_polygon import snap_closure  # noqa

MODULE = Path(site_config.MACHINE_YARD) / "the_knotfamily_machine" / "double3_module_extracted.npy"
SPECS = {"poly_c04": (4, False), "poly_c05": (5, False),
         "poly_c09": (9, False), "poly_s08": (8, True)}


def rot(seed):
    rng = np.random.default_rng(seed)
    q = rng.normal(size=4); q /= np.linalg.norm(q); w, x, y, z = q
    return np.array([[1-2*(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y)],
                     [2*(x*y+w*z), 1-2*(x*x+z*z), 2*(y*z-w*x)],
                     [2*(x*z-w*y), 2*(y*z+w*x), 1-2*(x*x+y*y)]])


def seg_intersect_2d(a0, a1, b0, b1):
    d1 = a1[:2]-a0[:2]; d2 = b1[:2]-b0[:2]
    den = d1[0]*d2[1]-d1[1]*d2[0]
    if abs(den) < 1e-12:
        return None
    t = ((b0[0]-a0[0])*d2[1]-(b0[1]-a0[1])*d2[0])/den
    s = ((b0[0]-a0[0])*d1[1]-(b0[1]-a0[1])*d1[0])/den
    if 1e-9 < t < 1-1e-9 and 1e-9 < s < 1-1e-9:
        return t, s
    return None


def pd_from_polygon(P):
    """PD code (standard convention, 2C edges) from a closed 3D polyline.
    Identical logic to r3_polygon_sage_det.pd_from_polygon."""
    N = len(P)
    cross = []
    for i in range(N):
        a0, a1 = P[i], P[(i+1) % N]
        for j in range(i+1, N):
            if j == i or (i+1) % N == j or (j+1) % N == i:
                continue
            b0, b1 = P[j], P[(j+1) % N]
            r = seg_intersect_2d(a0, a1, b0, b1)
            if r is None:
                continue
            t, s = r
            za = a0[2]+t*(a1[2]-a0[2]); zb = b0[2]+s*(b1[2]-b0[2])
            cross.append({"i": i, "t": t, "j": j, "s": s, "over_is_i": za > zb})
    C = len(cross)
    if C == 0:
        return None, 0
    enc = []
    for ci, c in enumerate(cross):
        enc.append((c["i"]+c["t"], ci, "i"))
        enc.append((c["j"]+c["s"], ci, "j"))
    enc.sort()
    M = len(enc)
    pos_of = {}
    for idx, (_, ci, wh) in enumerate(enc):
        pos_of[(ci, wh)] = idx

    def in_edge(idx): return ((idx-1) % M)+1
    def out_edge(idx): return (idx % M)+1

    pd = []
    for ci, c in enumerate(cross):
        under_wh = "j" if c["over_is_i"] else "i"
        over_wh = "i" if c["over_is_i"] else "j"
        iu = pos_of[(ci, under_wh)]; io = pos_of[(ci, over_wh)]
        a = in_edge(iu); cc = out_edge(iu)
        e_oi = in_edge(io); e_oo = out_edge(io)
        ui, ti = (c["i"], c["t"]) if under_wh == "i" else (c["j"], c["s"])
        oi, to = (c["i"], c["t"]) if over_wh == "i" else (c["j"], c["s"])
        ud = (P[(ui+1) % N]-P[ui])[:2]; od = (P[(oi+1) % N]-P[oi])[:2]
        dirs = {a: -ud, cc: ud, e_oi: -od, e_oo: od}
        ang0 = np.arctan2(dirs[a][1], dirs[a][0])
        others = [a, cc, e_oi, e_oo]
        order = sorted(others, key=lambda e: (np.arctan2(dirs[e][1], dirs[e][0]) - ang0) % (2*np.pi))
        pd.append([int(x) for x in order])
    return pd, C


def regina_det_deg(pd):
    """det = |Alexander(-1)|, degree = Alexander degree, via regina (C++).
    Evaluation directly via the integer polynomial coefficients (ax[i]) —
    NO string parse (regina writes '9 x' implicitly, sympify stumbles)."""
    L = regina.Link.fromPD(pd)
    # Reidemeister simplification BEFORE Alexander: det is invariant, but the
    # Alexander routine is exponentially slow on ~200-crossing projections.
    # simplify() reduces the polygon shadow diagram to near the minimal
    # crossing number -> same det, seconds instead of hours.
    L.simplify()
    ax = L.alexander()                      # regina PolynomialInt in x, from x^0
    deg = ax.degree()
    det = abs(sum(int(str(ax[i])) * ((-1) ** i) for i in range(deg + 1)))
    return int(det), int(deg)


def selftest():
    """Trefoil PD (standard convention) must give det=3 -> PD convention ok."""
    tref = [[1, 4, 2, 5], [3, 6, 4, 1], [5, 2, 6, 3]]
    det, deg = regina_det_deg(tref)
    assert det == 3, f"Trefoil-Selbsttest FEHLGESCHLAGEN: det={det} (erwartet 3)"
    print(f"[selftest] Trefoil det={det} deg={deg} -> PD convention regina ok")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stems", nargs="*", default=list(SPECS))
    ap.add_argument("-k", "--seeds", type=int, default=12)
    ap.add_argument("--out", default=str(ROOT / "claude_setup/tmp_durable/r3_polygon_regina_det.json"))
    args = ap.parse_args()
    selftest()
    stems = args.stems or list(SPECS)
    M = np.load(MODULE)
    results = {}
    for stem in stems:
        if stem not in SPECS:
            print(f"{stem}: no SPEC — skip"); continue
        n, mir = SPECS[stem]
        G = polygon_generator(n, mir)
        base = build_polygon_symmetric(snap_closure(M, G), n=n, mirror=mir,
                                       drop_seam_duplicates=True, closure_tol=1e-6)
        dets, degs, Cs, errs = [], [], [], 0
        for seed in range(args.seeds):
            P = base @ rot(5000+seed).T
            pd, C = pd_from_polygon(P)
            if pd is None:
                errs += 1; continue
            try:
                d, dg = regina_det_deg(pd)
                dets.append(d); degs.append(dg); Cs.append(C)
            except Exception as e:
                errs += 1
                if errs <= 2:
                    print(f"  {stem} seed{seed}: ERR {type(e).__name__}: {e}")
        cnt = Counter(dets)
        verdict = (f"STABIL det={cnt.most_common(1)[0][0]}" if len(cnt) == 1 and dets
                   else f"MISCH {dict(cnt)}" if dets else "alle degeneriert")
        degc = Counter(degs)
        results[stem] = {
            "n_mirror": [n, mir], "polygon_verts": int(len(base)),
            "n_valid": len(dets), "n_seeds": args.seeds, "errs": errs,
            "cross": Cs[0] if Cs else None,
            "det_counter": dict(cnt), "deg_counter": dict(degc),
            "det": cnt.most_common(1)[0][0] if len(cnt) == 1 and dets else None,
            "stable": len(cnt) == 1 and len(dets) == args.seeds,
        }
        print(f"{stem:9} (~{Cs[0] if Cs else '?'} cross, {len(base)} verts): "
              f"{len(dets)}/{args.seeds} valid, errs={errs} -> {verdict}, deg={dict(degc)}")
    Path(args.out).write_text(json.dumps(results, indent=2))
    print(f"\nwritten: {args.out}")


if __name__ == "__main__":
    main()
