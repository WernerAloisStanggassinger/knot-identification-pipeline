#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""sweep.py — control sweep over Werner's lattice. Full box."""
import site_config
import itertools, json, os, sys, traceback
import multiprocessing as mp
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, (site_config.PIPELINE_DIR))
sys.path.insert(0, (site_config.MACHINE_YARD + '/the_knotfamily_machine'))

def eine(arg):
    k, skala, spalt = arg
    try:
        from build_lattice import baue
        from embedding import min_segmentabstand
        from fseries_io import polyline_to_fourier, write_fseries
        from aggressive_simplify import parse_fseries, reconstruct_curve
        name = f"gitter_d3_k{k}_s{skala:g}_g{spalt:g}".replace(".", "p")
        P, ecken = baue(k, skala, spalt)
        terme = 1400
        fa, fb = polyline_to_fourier(P, n_terms=terme, n_samples=max(8192, 4*terme))
        pfad = site_config.MACHINE_YARD + f"/258_double_3_lattice_module/kandidaten/{name}.fseries"
        os.makedirs(os.path.dirname(pfad), exist_ok=True)
        write_fseries(pfad, fa, fb, knot_name=name, source="Maschine 258 sweep")
        off, ax, bx, ay, by, az, bz = parse_fseries(pfad)
        werte = []
        for N in (60000, 120000):
            R = reconstruct_curve(off, ax, bx, ay, by, az, bz, N).astype(np.float64)
            d, wo, h = min_segmentabstand(R)
            L = float(np.linalg.norm(np.diff(np.vstack([R, R[:1]]), axis=0), axis=1).sum())
            werte.append((d, L))
        stabil = abs(werte[0][0]-werte[1][0]) < 1e-9*max(1,werte[0][0]) or \
                 abs(werte[0][0]-werte[1][0])/max(werte[0][0],1e-12) < 1e-3
        d, L = werte[1]
        return dict(name=name, k=k, skala=skala, spalt=spalt, pfad=pfad,
                    min_abstand=d, laenge=L, verhaeltnis=d/L, stabil=bool(stabil),
                    punkte=int(len(P)))
    except Exception as e:
        return dict(name=f"k{k}_s{skala}_g{spalt}", fehler=f"{e.__class__.__name__}: {e}",
                    tb=traceback.format_exc()[-400:])

if __name__ == "__main__":
    gitter = list(itertools.product([3, 6], [1.2, 1.35, 1.5, 1.8, 2.2], [0.10, 0.25]))
    print(f"{len(gitter)} variants, {min(20, len(gitter))} workers", flush=True)
    with mp.Pool(min(20, len(gitter))) as pool:
        ergebnisse = [r for r in pool.imap_unordered(eine, gitter)]
    ok = [r for r in ergebnisse if "fehler" not in r]
    bad = [r for r in ergebnisse if "fehler" in r]
    ok.sort(key=lambda r: -r["verhaeltnis"])
    print(f"\n{'Name':<34}{'d':>11}{'Length':>10}{'d/L':>12}  stable")
    for r in ok:
        print(f"{r['name']:<34}{r['min_abstand']:11.6f}{r['laenge']:10.2f}{r['verhaeltnis']:12.3e}  {'✅' if r['stabil'] else '🔴'}")
    for r in bad:
        print(f"🔴 {r['name']}: {r['fehler']}")
    json.dump(ergebnisse, open((site_config.MACHINE_YARD + '/258_double_3_lattice_module/sweep.json'),"w",encoding="utf-8"), indent=1, ensure_ascii=False)
    print(f"\n{len(ok)} built, {len(bad)} failed → sweep.json")
