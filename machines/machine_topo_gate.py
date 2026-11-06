#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Part of the Knot-Identification Pipeline.
# ============================================================================
# machine_topo_gate.py — TOPOLOGICAL GATE for construction machines.
#
# ROOT (Disc #264, 2026-06-08): the machine demo gates only check EMBEDDING
# (self-intersection), NEVER the KNOT TYPE → 46% of all demo_smoke fseries
# were secretly unknots (clasp/twist does not hold; "machine built" ≠
# "delivers knots"). This module supplies the missing gate, independent of the pipeline
# (no kinematic scanner): curve/fseries → pyknotid PD → regina minimal
# crossing number. reduced==0 = unknot.
#
# INTEGRATION (one line in make_knot/cmd_demo, AFTER the embedding gate):
#     from machine_topo_gate import assert_real_knot
#     assert_real_knot(D, name)        # raises AssertionError on unknot
# or soft:
#     if reduced_crossings(D) == 0: entry["rejected"] = "topo gate (unknot)"; return ...
#
# Runs in base env (pyknotid 0.5.3 + regina). NB pyknotid alexander/determinant
# are broken with new numpy (np.complex) → use regina jones/homfly directly.
# ============================================================================
import math
import re

import numpy as np


def reduced_crossings(P):
    """regina minimal crossing number from 3D polyline P (N,3). 0 = unknot.
    Controls: trefoil→3, circle→0, fig8→4 (verified 2026-06-08)."""
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


def jones_str(P):
    """regina Jones polynomial (string) of the polyline; '1' = unknot-suspect.
    For mutant-pair equality: jones_str(A)==jones_str(B)."""
    from pyknotid.spacecurves import Knot
    import regina
    K = Knot(np.asarray(P, float), verbose=False)
    tu = [tuple(int(x) for x in t.split(','))
          for t in re.findall(r'X_\{([\d,]+)\}', str(K.planar_diagram()))]
    if not tu:
        return "1"
    s = '[' + ','.join('[' + ','.join(map(str, c)) + ']' for c in tu) + ']'
    L = regina.Link.fromPD(s)
    L.intelligentSimplify()
    return str(L.jones())


# ----------------------------------------------------- fseries → Polyline
def fseries_curve(path, N=None):
    rows = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith('%'):
            continue
        rows.append([float(v) for v in line.split()])
    off = np.array(rows[0]) if len(rows[0]) == 3 else np.zeros(3)
    Hm = np.array(rows[1:] if len(rows[0]) == 3 else rows)
    A = Hm[:, [0, 2, 4]]; B = Hm[:, [1, 3, 5]]
    mag = np.abs(A) + np.abs(B)
    period = math.pi if mag[0::2].max() < 1e-9 * mag.max() else 2 * math.pi
    if N is None:
        N = max(2000, 12 * A.shape[0])
    j = np.arange(1, A.shape[0] + 1)
    t = np.linspace(0.0, period, N, endpoint=False)
    th = np.outer(t, j)
    return off + np.cos(th) @ A + np.sin(th) @ B


def assert_real_knot(curve_or_path, name="<knot>"):
    """Gate: raises AssertionError if the object is the unknot.
    curve_or_path: (N,3) polyline OR fseries path (str)."""
    P = fseries_curve(curve_or_path) if isinstance(curve_or_path, str) else curve_or_path
    rc = reduced_crossings(P)
    assert rc > 0, (f"TOPO-GATE: {name} is the UNKNOT (reduced=0) — "
                    f"construction does not close / under-resolves. Do NOT emit.")
    return rc


if __name__ == "__main__":
    # self-test: control curves
    t = np.linspace(0, 2 * np.pi, 400, endpoint=False)
    tre = np.c_[np.sin(t) + 2 * np.sin(2 * t), np.cos(t) - 2 * np.cos(2 * t), -np.sin(3 * t)]
    cir = np.c_[np.cos(t), np.sin(t), 0 * t]
    print("Trefoil reduced =", reduced_crossings(tre), "(expected 3)")
    print("Circle  reduced =", reduced_crossings(cir), "(expected 0)")
    assert reduced_crossings(tre) == 3 and reduced_crossings(cir) == 0
    print("Selbsttest OK — Gate scharf.")
