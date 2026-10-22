#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# ============================================================================
# audit_topology.py — TOPOLOGICAL GATE over ALL machine demo_smoke .fseries files.
#
# QUESTION (integrity of outgoing mail): which construction machines silently produce
# the UNKNOT? (like 81 Whitehead + 90 Wobble — the demo gates check only the embedding,
# NEVER the knot type). Per .fseries: evaluate a dense polyline -> pyknotid PD ->
# regina intelligentSimplify -> minimal crossing count. reduced == 0 = unknot.
#
# Independent of the pipeline (no kinematic scanner), fast, parallel (full throttle).
# Usage (base python with pyknotid + regina): python audit_topology.py [--workers 12]
# The release imports only its functions (unknot_zweitmeinung.py); main() is the original audit run.
# ============================================================================
import argparse, glob, json, math, os, re, sys, traceback
import site_config
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

YARD = site_config.MACHINE_YARD
OUTDIR = os.path.join(YARD, "_topo_gate_audit")


def parse_fseries(path):
    rows = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith('%'):
            continue
        rows.append([float(v) for v in line.split()])
    if len(rows[0]) == 3:
        offset = np.array(rows[0]); harm = rows[1:]
    else:
        offset = np.zeros(3); harm = rows
    Hm = np.array(harm)
    A = Hm[:, [0, 2, 4]].copy(); B = Hm[:, [1, 3, 5]].copy()
    return offset, A, B


def detect_period_pi(A, B, tol=1e-9):
    mag = np.abs(A) + np.abs(B); sc = mag.max()
    if sc == 0:
        return False
    return bool(mag[0::2].max() < tol * sc)


def eval_curve(offset, A, B, N, period):
    j = np.arange(1, A.shape[0] + 1, dtype=float)
    t = np.linspace(0.0, period, N, endpoint=False)
    th = np.outer(t, j)
    return offset + np.cos(th) @ A + np.sin(th) @ B


def reduced_crossings(P):
    """Regina minimal crossing count from a 3D polyline (pyknotid PD -> regina)."""
    from pyknotid.spacecurves import Knot
    import regina
    K = Knot(np.asarray(P, float), verbose=False)
    tu = [tuple(int(x) for x in t.split(','))
          for t in re.findall(r'X_\{([\d,]+)\}', str(K.planar_diagram()))]
    if not tu:
        return 0
    s = '[' + ','.join('[' + ','.join(map(str, c)) + ']' for c in tu) + ']'
    L = regina.Link.fromPD(s)
    L.intelligentSimplify()
    return L.size()


def worker(path):
    try:
        offset, A, B = parse_fseries(path)
        period = math.pi if detect_period_pi(A, B) else 2 * math.pi
        N = max(2000, 12 * A.shape[0])
        P = eval_curve(offset, A, B, N, period)
        rc = reduced_crossings(P)
        rel = os.path.relpath(path, YARD)
        machine = rel.split(os.sep)[0]
        stem = os.path.basename(path)[:-len(".fseries")]
        return {"machine": machine, "stem": stem, "reduced": int(rc),
                "unknot": rc == 0, "H": int(A.shape[0]), "path": rel}
    except Exception as e:
        return {"machine": os.path.relpath(path, YARD).split(os.sep)[0],
                "stem": os.path.basename(path), "error": str(e)[:120],
                "trace": traceback.format_exc()[:300]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=12)
    args = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)
    files = sorted(glob.glob(os.path.join(YARD, "*", "demo_smoke", "**", "*.fseries"),
                             recursive=True))
    print(f"[TOPO-AUDIT] {len(files)} demo_smoke .fseries, {args.workers} workers", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(worker, f): f for f in files}
        for i, fu in enumerate(as_completed(futs), 1):
            r = fu.result(); rows.append(r)
            tag = ("ERR " if "error" in r else
                   ("UNKNOT" if r["unknot"] else f"knot({r['reduced']})"))
            print(f"[{i}/{len(files)}] {r['machine']:32s} {r.get('stem',''):20s} {tag}",
                  flush=True)
    rows.sort(key=lambda r: (r["machine"], r.get("stem", "")))
    json.dump(rows, open(os.path.join(OUTDIR, "topo_audit_report.json"), "w"),
              indent=1, ensure_ascii=False)
    unknots = [r for r in rows if r.get("unknot")]
    errs = [r for r in rows if "error" in r]
    print("\n" + "=" * 64)
    print(f"TOTAL {len(rows)} · UNKNOT {len(unknots)} · ERRORS {len(errs)}")
    print("\n--- UNKNOTS (reduced==0) ---")
    for r in unknots:
        print(f"  {r['machine']:34s} {r['stem']}")
    if errs:
        print("\n--- ERRORS ---")
        for r in errs:
            print(f"  {r['machine']:34s} {r.get('stem','')}  {r.get('error')}")
    print(f"\n[REPORT] {os.path.join(OUTDIR, 'topo_audit_report.json')}")


if __name__ == "__main__":
    main()
