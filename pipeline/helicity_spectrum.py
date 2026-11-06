#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Part of the Knot-Identification Pipeline.
"""helicity_spectrum.py — writhe (= helicity/Gamma^2 for a twist-free vortex tube)
over the whole corpus, directly from the served .fseries centerlines.

VORTEX TRACK (Werner 2026-07-23): a knot centerline IS a vortex filament. For a
thin vortex tube WITHOUT internal twist, Moffatt/Calugareanu holds:
      H / Gamma^2 = SL = Wr + Tw ,  Tw = 0  =>  H/Gamma^2 = Wr(centerline).
Writhe is computable framing-free from the space curve (Gauss double integral), so we need
NO flow simulation — only the curve we already have.

HONESTY (also stated in the result):
  * Writhe is NOT a knot invariant, but a quantity of the REALISATION. For own_knots
    (sculptures) it is the real physical helicity of the steel curve; for machine knots
    that of the generated course (depends on the parametrisation, not on the knot type).
  * Amphichiral, symmetrically realised knots should show Wr ~ 0 (mirror antisymmetry) —
    that is the built-in plausibility check, not a circle.
  * Computed on the SERVED .fseries object (project_record_object_is_served_fseries),
    band-limited to its n_terms harmonics. At high c that is a smoothed curve;
    the writhe holds for exactly this realisation, not for a minimal diagram.

Method: writhe by midpoint discretisation of the Gauss double integral (Klenin-Langowski 1a),
row-chunked, so that the memory per knot stays small (freeze-safe next to running jobs).
"""
import site_config
import json, os, sys
import numpy as np

WURZEL = site_config.ROOT
DB = site_config.PIPELINE_DIR + "/corpus_db.json"
OUT = WURZEL + "/konversationen/helizitaet_spektrum_20260723.json"


def parse_fseries(pfad):
    a, b = [], []
    with open(pfad) as fh:
        for ln in fh:
            s = ln.strip()
            if not s or s.startswith("%"):
                continue
            parts = s.split()
            if len(parts) < 6:
                continue
            try:
                vals = [float(x) for x in parts[:6]]
            except ValueError:
                continue
            ax, bx, ay, by, az, bz = vals
            a.append([ax, ay, az])
            b.append([bx, by, bz])
    if not a:
        raise ValueError("no coefficients in %s" % pfad)
    return np.array(a), np.array(b)


def reconstruct(a, b, n_samples):
    nt = a.shape[0]
    s = np.linspace(0.0, 2 * np.pi, n_samples, endpoint=False)
    out = np.zeros((n_samples, 3))
    for j in range(1, nt + 1):
        out += np.outer(np.cos(j * s), a[j - 1]) + np.outer(np.sin(j * s), b[j - 1])
    return out


def writhe(curve, chunk=256):
    """Writhe of a closed space curve (midpoint Gauss, row-chunked)."""
    r = np.asarray(curve, float)
    n = r.shape[0]
    rn = np.roll(r, -1, axis=0)
    seg = rn - r                      # tangent * ds
    mid = 0.5 * (r + rn)
    total = 0.0
    for i0 in range(0, n, chunk):
        i1 = min(i0 + chunk, n)
        mi = mid[i0:i1]               # (k,3)
        si = seg[i0:i1]               # (k,3)
        dv = mi[:, None, :] - mid[None, :, :]         # (k,n,3)
        cr = np.cross(si[:, None, :], seg[None, :, :])# (k,n,3)
        num = np.einsum("knd,knd->kn", dv, cr)        # (k,n)
        dist = np.linalg.norm(dv, axis=2)             # (k,n)
        # mask i==j (self pair); dist there 0
        gi = np.arange(i0, i1)
        dist[np.arange(i1 - i0), gi] = np.inf
        with np.errstate(divide="ignore", invalid="ignore"):
            term = num / (dist ** 3)
        term[~np.isfinite(term)] = 0.0
        total += term.sum()
    return total / (4.0 * np.pi)


def main():
    d = json.load(open(DB))
    ks = d["knots"]
    res = []
    miss = 0
    for idx, r in enumerate(ks):
        pp = r.get("pflicht_path")
        stem = r.get("stem")
        if not pp or not stem:
            miss += 1
            continue
        fpath = os.path.join(os.path.dirname(pp), stem + ".fseries")
        if not os.path.exists(fpath):
            miss += 1
            continue
        try:
            a, b = parse_fseries(fpath)
        except Exception as e:
            miss += 1
            continue
        nt = a.shape[0]
        n_samples = int(min(4000, max(2000, 12 * nt)))
        curve = reconstruct(a, b, n_samples)
        wr = float(writhe(curve))
        res.append({
            "stem": stem,
            "category": r.get("category"),
            "distinct_knot_id": r.get("distinct_knot_id"),
            "is_distinct_representative": r.get("is_distinct_representative"),
            "cluster_folder": r.get("cluster_folder"),
            "c": r.get("c"), "det": r.get("det"),
            "amphi": r.get("amphi"), "alternating": r.get("alternating"),
            "n_terms": nt, "n_samples": n_samples,
            "writhe": wr, "abs_writhe": abs(wr),
            "helicity_over_gamma2": wr,   # = Wr for a twist-free tube
        })
        if (idx + 1) % 100 == 0:
            print("  %d/%d ..." % (idx + 1, len(ks)), flush=True)
    json.dump({"meta": {"n": len(res), "missing_fseries": miss,
                        "definition": "helicity/Gamma^2 = writhe (twist-free vortex tube); "
                                      "writhe is a property of the realization, not a knot invariant"},
               "records": res}, open(OUT, "w"), indent=2)
    print("written: %s  (%d knots, %d without fseries)" % (OUT, len(res), miss), flush=True)


if __name__ == "__main__":
    main()
