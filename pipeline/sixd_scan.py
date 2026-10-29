#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Knot-Identification Pipeline — Werner Alois Stanggassinger & Claude Code (Anthropic).
# ============================================================================
# sixd_scan.py — 6D geometry machine (2026-06-05, Werner's go: "build the 6D machine, full throttle")
#
# Measures, per record object (= served .fseries, decision A / Disc #251), the
# geometric axes that the PD shadow throws away:
#   - curvature kappa(s)        (2nd order, analytic from Fourier derivatives)
#   - torsion  tau(s)           (3rd order — THE 6th axis, mutant candidate)
#   - writhe   Wr               (Gauss double integral, odd under chirality)
#   - ACN                       (average crossing number, |Integrand|)
#
# All derivatives ANALYTIC (term-wise Fourier differentiation) — no
# finite-difference noise. All aggregates are invariant under:
#   reparametrisation (ds weights), similarity (scale normalisation kappa*L,
#   tau*L), rotation/translation, start-point shift + orientation reversal
#   (spectra = |FFT|). Mirroring flips the sign of tau/Wr (feature: chirality).
# NOT isotopy-invariant by construction — the cluster collapse measurement
# (sixd_analyze.py) tests empirically WHICH aggregates are isotopy-stable.
#
# Doctrine points:
#   - Nyquist/Densify (Memory feedback-resolution-nyquist-triangulation):
#     N_geo = max(NGEO_MIN, NGEO_PER_H * H), a second resolution 2x as witness.
#   - pi-periodic .fseries (only even j, Disc #70) -> period [0,pi),
#     otherwise double traversal = singular Gauss integral.
#   - No table lookup. corpus_db only for paths/stems (hard rule 1 ok).
#   - Full throttle: all cores, generous defaults (hard rule 6).
# ============================================================================
import argparse
import site_config
import hashlib
import json
import math
import os
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

ROOT = site_config.ROOT
CORPUS_DB = ROOT + "/freeze_manifest/pipeline_scripts/corpus_db.json"
OUT_DIR = ROOT + "/sixd_machine/out"

NGEO_MIN = 8192          # geometry grid lower bound
NGEO_PER_H = 64          # points per harmonic (>> Nyquist 2H)
NGEO_CAP = 65536         # upper bound (keeps j^3-amplified high harmonics clean)
NWR_DEFAULT = 16384      # writhe/ACN grid (O(N^2), generous)
NSPEC = 2048             # uniform arc-length grid for spectra
NSPEC_KEEP = 24          # |FFT| coefficients kept
QUANTILES = [0.05, 0.25, 0.50, 0.75, 0.95]


# ---------------------------------------------------------------- fseries IO
def parse_fseries(path):
    """Identical to the knot_scanner logic (lines with '%' = comment,
    first 3-column line = offset, otherwise 6 columns ax bx ay by az bz, j=1..H)."""
    coeffs = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('%'):
                continue
            coeffs.append([float(v) for v in line.split()])
    if len(coeffs[0]) == 3:
        offset = np.array(coeffs[0])
        harmonics = coeffs[1:]
    else:
        offset = np.zeros(3)
        harmonics = coeffs
    H = np.array(harmonics)                      # (H, 6)
    A = H[:, [0, 2, 4]].copy()                   # cos coefficients (H,3) x,y,z
    B = H[:, [1, 3, 5]].copy()                   # sin coefficients (H,3)
    return offset, A, B


def detect_period_pi(A, B, tol=1e-9):
    """Disc #70: only even j non-zero -> curve is pi-periodic -> period pi."""
    mag = np.abs(A) + np.abs(B)
    scale = mag.max()
    if scale == 0:
        return False
    odd = mag[0::2]                              # j=1,3,5,... (0-indexed)
    return bool(odd.max() < tol * scale)


# ------------------------------------------------------- analytic Frenet
def eval_derivatives(offset, A, B, t):
    """r, r', r'', r''' analytically. t: (N,). Returns (N,3) each."""
    Hn = A.shape[0]
    j = np.arange(1, Hn + 1, dtype=np.float64)
    th = np.outer(t, j)                          # (N,H)
    C, S = np.cos(th), np.sin(th)
    j1 = j[:, None]
    r = offset + C @ A + S @ B
    r1 = -S @ (A * j1) + C @ (B * j1)
    r2 = -C @ (A * j1**2) - S @ (B * j1**2)
    r3 = S @ (A * j1**3) - C @ (B * j1**3)
    return r, r1, r2, r3


def frenet_profile(offset, A, B, N, period):
    t = np.linspace(0.0, period, N, endpoint=False)
    dt = period / N
    r, r1, r2, r3 = eval_derivatives(offset, A, B, t)
    v = np.linalg.norm(r1, axis=1)               # speed
    cr = np.cross(r1, r2)
    ncr = np.linalg.norm(cr, axis=1)
    ds = v * dt
    L = ds.sum()
    kappa = ncr / np.maximum(v, 1e-300) ** 3
    # tau = det(r',r'',r''')/|r'xr''|^2 ; mask where kappa ~ 0 (Frenet singularity)
    eps = 1e-9 * np.median(ncr) if np.median(ncr) > 0 else 1e-300
    good = ncr > eps
    tau = np.zeros(N)
    tau[good] = np.einsum('ij,ij->i', cr[good], r3[good]) / ncr[good] ** 2
    return dict(t=t, r=r, r1=r1, v=v, ds=ds, L=L,
                kappa=kappa, tau=tau, good=good, dt=dt)


# --------------------------------------------------- Gauss double integral
def writhe_acn(r, T, ds, block=512):
    """Wr = 1/(4pi) sum_{i!=j} ((Ti x Tj).(ri-rj))/|ri-rj|^3 dsi dsj ; ACN likewise with | |."""
    N = len(ds)
    wr = 0.0
    acn = 0.0
    for i0 in range(0, N, block):
        i1 = min(i0 + block, N)
        d = r[i0:i1, None, :] - r[None, :, :]            # (b,N,3)
        dist3 = np.linalg.norm(d, axis=2) ** 3            # (b,N)
        TxT = np.cross(T[i0:i1, None, :], T[None, :, :])  # (b,N,3)
        num = np.einsum('bnk,bnk->bn', TxT, d)            # (b,N)
        w = np.outer(ds[i0:i1], ds)                       # (b,N)
        idx = np.arange(i0, i1)
        dist3[np.arange(i1 - i0), idx] = np.inf           # drop the diagonal
        integ = num / dist3 * w
        wr += integ.sum()
        acn += np.abs(integ).sum()
    return wr / (4 * math.pi), acn / (4 * math.pi)


# --------------------------------------------------------------- Aggregation
def weighted_quantiles(x, w, qs):
    o = np.argsort(x)
    x, w = x[o], w[o]
    cw = np.cumsum(w) - 0.5 * w
    cw /= w.sum()
    return [float(np.interp(q, cw, x)) for q in qs]


def spectrum(profile_x, s_cum, L, nspec, keep):
    """uniform arc-length resampling -> |rFFT|/M, invariant under start point/reversal."""
    su = np.linspace(0.0, L, nspec, endpoint=False)
    xu = np.interp(su, s_cum, profile_x, period=L)
    mag = np.abs(np.fft.rfft(xu)) / nspec
    return [float(m) for m in mag[1:keep + 1]]   # drop DC (it is in the scalars)


def aggregate(fp):
    ds, L = fp["ds"], fp["L"]
    w = ds / L
    kn = fp["kappa"] * L                          # dimensionless
    tn = fp["tau"] * L
    g = fp["good"]
    s_cum = np.cumsum(ds) - 0.5 * ds
    agg = {
        "L": float(L),
        "total_curvature_over_2pi": float((fp["kappa"] * ds).sum() / (2 * math.pi)),
        "total_torsion_signed_over_2pi": float((fp["tau"] * ds).sum() / (2 * math.pi)),
        "total_torsion_abs_over_2pi": float((np.abs(fp["tau"]) * ds).sum() / (2 * math.pi)),
        "tau_masked_fraction": float(1.0 - g.mean()),
        "kappaL_quantiles": weighted_quantiles(kn, w, QUANTILES),
        "tauL_quantiles": weighted_quantiles(tn[g], w[g], QUANTILES) if g.any() else None,
        "kappaL_spectrum": spectrum(kn, s_cum, L, NSPEC, NSPEC_KEEP),
        "tauL_spectrum": spectrum(tn, s_cum, L, NSPEC, NSPEC_KEEP),
    }
    return agg


# ------------------------------------------------------------------ Record
def scan_record(job):
    stem, fpath, nwr = job
    try:
        offset, A, B = parse_fseries(fpath)
        Hn = A.shape[0]
        period = math.pi if detect_period_pi(A, B) else 2 * math.pi
        n_geo = min(max(NGEO_MIN, NGEO_PER_H * Hn), NGEO_CAP)

        results = {}
        for label, N in (("res1", n_geo), ("res2", min(2 * n_geo, 2 * NGEO_CAP))):
            fp = frenet_profile(offset, A, B, N, period)
            results[label] = aggregate(fp)

        # writhe/ACN on its own grid (O(N^2))
        n_wr = min(nwr, NGEO_CAP)
        fpw = frenet_profile(offset, A, B, n_wr, period)
        T = fpw["r1"] / fpw["v"][:, None]
        wr, acn = writhe_acn(fpw["r"], T, fpw["ds"])
        for label in results:
            results[label]["writhe"] = float(wr)
            results[label]["acn"] = float(acn)

        # stability witness res1 vs res2 (Nyquist doctrine: >=2 resolutions)
        keys = ["total_curvature_over_2pi", "total_torsion_signed_over_2pi",
                "total_torsion_abs_over_2pi"]
        rel = {}
        for k in keys:
            a, b = results["res1"][k], results["res2"][k]
            rel[k] = abs(a - b) / max(abs(a), abs(b), 1e-12)
        payload = hashlib.sha256(
            "".join(l for l in open(fpath) if not l.strip().startswith('%'))
            .encode()).hexdigest()[:16]
        out = {
            "stem": stem,
            "fseries": fpath,
            "fseries_payload_sha16": payload,
            "n_harmonics": int(Hn),
            "period": "pi" if period == math.pi else "2pi",
            "n_geo": int(n_geo),
            "n_writhe": int(n_wr),
            "fingerprint": results["res2"],          # finer resolution = main
            "fingerprint_res1": results["res1"],
            "resolution_stability_reldiff": rel,
            "machine": "sixd_scan v1 2026-06-05",
        }
        op = os.path.join(OUT_DIR, stem + "_sixd.json")
        with open(op, "w") as f:
            json.dump(out, f, indent=1)
        flag = "OK " if max(rel.values()) < 1e-4 else "WARN-RES "
        return f"{flag}{stem} H={Hn} per={out['period']} " \
               f"Wr={wr:+.4f} ACN={acn:.2f} " \
               f"K={results['res2']['total_curvature_over_2pi']:.3f} " \
               f"Tsig={results['res2']['total_torsion_signed_over_2pi']:+.4f}"
    except Exception:
        return f"FAIL {stem}\n{traceback.format_exc()}"


def collect_jobs(nwr, only=None):
    db = json.load(open(CORPUS_DB))
    jobs, missing = [], []
    pool = list(db["knots"]) + list(db.get("validation_knots", []))
    for rec in pool:
        stem = rec["stem"]
        if only and stem not in only:
            continue
        folder = os.path.dirname(rec["pflicht_path"])
        fpath = os.path.join(folder, stem + ".fseries")
        if not os.path.isfile(fpath):
            cands = [f for f in os.listdir(folder)
                     if f.endswith(".fseries") and not f.startswith("_sculpture_")] \
                if os.path.isdir(folder) else []
            if len(cands) == 1:
                fpath = os.path.join(folder, cands[0])
            else:
                missing.append(stem)
                continue
        jobs.append((stem, fpath, nwr))
    return jobs, missing


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--nwr", type=int, default=NWR_DEFAULT)
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    jobs, missing = collect_jobs(args.nwr, set(args.only) if args.only else None)
    print(f"[6D] {len(jobs)} Records, {args.workers} Worker, "
          f"N_wr={args.nwr}; missing fseries: {len(missing)}")
    for m in missing:
        print(f"  MISSING {m}")

    n_ok = n_fail = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(scan_record, j) for j in jobs]
        for i, fu in enumerate(as_completed(futs), 1):
            line = fu.result()
            print(f"[{i}/{len(jobs)}] {line}", flush=True)
            if line.startswith("FAIL"):
                n_fail += 1
            else:
                n_ok += 1
    print(f"[6D] DONE: {n_ok} ok, {n_fail} fail, {len(missing)} missing")


if __name__ == "__main__":
    main()
