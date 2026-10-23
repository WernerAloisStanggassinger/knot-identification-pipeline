#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""unknot_second_opinion.py — independent second opinion on "is this .fseries an unknot?"

The FIRST opinion comes from the kinematic scanner (scan_triage.py: is there a
crossing-free projection?). This one comes from a completely different corner:
curve -> pyknotid PD -> regina minimal crossing count (reduced==0 = unknot).
No scanner, no camera, no (ti,tj) grid — the two cannot make
the same mistake together.

The computing functions are taken over VERBATIM from maschinen_ideen/_topo_gate_audit/
audit_topology.py (Werner's calibrated gate, Disc #264), so that the
second opinion is not a second, differently calibrated tool.

Runs in the base env (pyknotid + regina).
"""
import sys, os, re, math, json, argparse
import site_config
import numpy as np
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, site_config.TOPO_GATE_DIR)
from audit_topology import parse_fseries, detect_period_pi, eval_curve, reduced_crossings


def one_curve(path):
    try:
        offset, A, B = parse_fseries(path)
        period = math.pi if detect_period_pi(A, B) else 2 * math.pi
        N = max(2000, 12 * A.shape[0])
        rc = reduced_crossings(eval_curve(offset, A, B, N, period))
        return dict(path=path, stem=os.path.basename(path)[:-len('.fseries')],
                    reduced=int(rc), unknot=(rc == 0))
    except Exception as e:
        return dict(path=path, stem=os.path.basename(path)[:-len('.fseries')],
                    failures=str(e)[:120])


def self_test():
    """A tool that is meant to confirm UNKNOT needs the case built in
    that it must REJECT (Footgun 16). Otherwise green is worthless.
    Note: the probe names passed to pruef() are data for the gate and stay German."""
    ok = tot = 0
    def check(name, actual, expected):
        nonlocal ok, tot
        tot += 1; good = (actual == expected); ok += good
        print(f"  {'✅' if good else '🔴'} {name}: reduced={actual} (expected {expected})")
    t = np.linspace(0, 2*np.pi, 3000, endpoint=False)
    # circle — must give 0
    circle = np.stack([np.cos(t), np.sin(t), np.zeros_like(t)], 1)
    check('circle ⟹ unknot', reduced_crossings(circle), 0)
    # COUNTER-CASE trefoil T(2,3) — must give 3, otherwise the tool says "unknot" to everything
    trefoil = np.stack([np.sin(t)+2*np.sin(2*t), np.cos(t)-2*np.cos(2*t), -np.sin(3*t)], 1)
    check('trefoil ⟹ NOT unknot (counter-case)', reduced_crossings(trefoil), 3)
    # COUNTER-CASE figure-eight — 4
    fig8 = np.stack([(2+np.cos(2*t))*np.cos(3*t), (2+np.cos(2*t))*np.sin(3*t), np.sin(4*t)], 1)
    check('figure-eight ⟹ 4 (second counter-case)', reduced_crossings(fig8), 4)
    # canary: a shifted circle stays an unknot
    check('canary: shifted circle', reduced_crossings(circle + np.array([5., -3., 2.])), 0)
    print(f"\n  {ok}/{tot} probes")
    return 0 if ok == tot else 1


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('fseries', nargs='*')
    ap.add_argument('--self-test', action='store_true')
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--json')
    a = ap.parse_args()
    if a.self_test:
        sys.exit(self_test())
    if not a.fseries:
        ap.error('fseries paths or --selbsttest')
    with ProcessPoolExecutor(a.workers) as ex:
        res = list(ex.map(one_curve, a.fseries))
    u = [r for r in res if r.get('unknot')]
    k = [r for r in res if r.get('unknot') is False]
    f = [r for r in res if 'failures' in r]
    print(f"UNKNOTS {len(u)} · KNOTS {len(k)} · ERRORS {len(f)} of {len(res)}")
    for r in sorted(k, key=lambda r: -r['reduced']):
        print(f"  🔴 KNOT    {r['stem']:<40} reduced={r['reduced']}")
    for r in f:
        print(f"  ❔ ERROR   {r['stem']:<40} {r['failures']}")
    if a.json:
        json.dump(res, open(a.json, 'w'), indent=1)
        print(f"→ {a.json}")
