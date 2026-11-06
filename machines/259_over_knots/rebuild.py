#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Rebuild and measure the UEBER series with the repaired seam.

Werner 25.08.2026: "let it run in tangent-continuously, then rebuild."

The FOURIER curve is measured, not the polyline — the Fourier curve is
stored and identified, so it is the truth.
Reference is always double-3: quality 1.000, max kink 1.72 degrees.
"""
import site_config
import sys, os, json, numpy as np
H=(site_config.MACHINE_YARD + '/259_over_knots')
for p in (H, (site_config.MACHINE_YARD + '/258_double_3_lattice_module'),
          (site_config.PIPELINE_DIR),
          (site_config.MACHINE_YARD + '/the_knotfamily_machine')):
    sys.path.insert(0,p)
import importlib.util
sp=importlib.util.spec_from_file_location("build", f"{H}/build.py")
build=importlib.util.module_from_spec(sp); sp.loader.exec_module(build)
from embedding import min_segmentabstand
from fseries_io import polyline_to_fourier, write_fseries
from build_curve_json import reconstruct_curve

REF = 1.127e-3        # double-3 d/L
OUT = f"{H}/kurven_neu"; os.makedirs(OUT, exist_ok=True)

def profil(P):
    d=np.diff(np.vstack([P,P[:1]]),axis=0); L=np.linalg.norm(d,axis=1)
    u=d/np.maximum(L[:,None],1e-15)
    w=np.degrees(np.arccos(np.clip((u*np.roll(u,-1,axis=0)).sum(1),-1,1)))
    return float(w.max()), int((w>5).sum()), float(L.sum())

# 22400 instead of 1400: on the REPAIRED polyline the series converges only here.
# Measured 25.08. on UEBER_k3_s3 (polyline: max kink 0.96°, 0 over 5°):
#   1400 -> 10.55° / 334 kinks over 5°
#   2800 -> 14.48° / 339      (worse — the series replicates fine structure)
#   5600 ->  8.85° / 210
#  11200 ->  4.57° /   0      <- clean from here on
#  22400 ->  2.29° /   0      <- double-3 sits at 1.72°
# On the BROKEN polyline (12 zero segments, 90° kinks) more was always
# worse — hence first the seam was needed, then the terms. Individually
# neither of the two would have sufficed.
TERME = 22400
print(f"{'Curve':16s} {'maxKink':>9s} {'>5°':>5s} {'d_min':>10s} {'Quality':>8s}  (double-3: 1.72°, 0, quality 1.000)")
print("-"*74)
erg={}
for k in (3,6):
    for s in (1.5, 2.0, 3.0, 4.0):
        name=f"UEBER_k{k}_s{s:g}".replace(".","p")
        try:
            P,_ = build.baue(k, s, 0.10, nmod=40000, richtung="ueber")
            fa, fb = polyline_to_fourier(P, n_terms=TERME, n_samples=max(8192,4*TERME))
            pfad=f"{OUT}/{name}.fseries"
            write_fseries(pfad, fa, fb, knot_name=name, source="Maschine 259 tangentenstetige Naht")
            R = reconstruct_curve(np.zeros(3), np.asarray(fa), np.asarray(fb), max(8*TERME, 22400))
            mk, n5, L = profil(R)
            d,_,_ = min_segmentabstand(R, M=3000)
            g=(d/L)/REF
            erg[name]=dict(maxknick=mk, ueber5=n5, d=d, L=L, guete=g)
            print(f"{name:16s} {mk:8.2f}° {n5:>5d} {d:>10.5f} {g:>8.4f}", flush=True)
        except Exception as e:
            print(f"{name:16s} ERROR {type(e).__name__}: {str(e)[:60]}", flush=True)
json.dump(erg, open(f"{H}/rebuild_measurement.json","w"), indent=1)
print(f"\nwritten to {OUT}/ · measurement in rebuild_measurement.json")
