#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""homfly_validation.py — the calibrated own HOMFLY against the whole corpus.

Model alex_punkt: there 557/557 determinants, 0 deviations. Here the same
for HOMFLY. Measures at the same time WHERE the own tool has its own wall —
it lies at the STRAND NUMBER (cost ~ sum of the squares of the irrep dimensions
= n!), not at the crossing number. That is the flip side of the trade: regina
dies of c, the Hecke trace dies of n.
"""
import site_config
import sys, os, json, time, argparse, tempfile
sys.path.insert(0, (site_config.TOOLS_DIR))
from homfly_calibration import parse_homfly, eval_homfly, VARIANTS, norm_K


def _mit_frist(fn, sek):
    """Footgun 13/18/20: throw-away child, WNOHANG, own clock, never alarm()."""
    fd, pf = tempfile.mkstemp(prefix='hom_'); os.close(fd)
    t0 = time.time(); pid = os.fork()
    if pid == 0:
        try:
            v = fn()
            open(pf, 'w').write(json.dumps([v.real, v.imag]))
            os._exit(0)
        except BaseException:
            os._exit(1)
    while True:
        p, st = os.waitpid(pid, os.WNOHANG)
        if p:
            try:
                re_, im_ = json.load(open(pf)); os.unlink(pf)
                return complex(re_, im_), time.time() - t0
            except BaseException:
                try: os.unlink(pf)
                except OSError: pass
                return None, time.time() - t0
        if time.time() - t0 > sek:
            os.kill(pid, 9); os.waitpid(pid, 0)
            try: os.unlink(pf)
            except OSError: pass
            return None, time.time() - t0
        time.sleep(0.2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--faelle', required=True)
    ap.add_argument('--nmax', type=int, default=10)
    ap.add_argument('--frist', type=int, default=120)
    ap.add_argument('--procs', type=int, default=12)
    ap.add_argument('--a', type=float, default=1.3)
    ap.add_argument('--q', type=float, default=1.2)
    ap.add_argument('--json')
    args = ap.parse_args()

    from hecke_homfly_fast import homfly_hecke
    faelle = json.load(open(args.faelle))
    a, q = args.a, args.q
    x, y = VARIANTS[0][1](complex(a), complex(q))
    K = norm_K(complex(a), complex(q))
    print(f"Calibration point a={a} q={q} ⟹ x={x.real:.6f} y={y.real:.6f} · K={K.real:.9f}")
    print(f"{len(faelle)} records · only n<={args.nmax} · deadline {args.frist}s per case\n")

    treffer = fehl = zu_gross = tot = 0
    schlimmster = 0.0; schlimmster_name = ''
    zeit_pro_n = {}
    abweichungen = []
    for f in faelle:
        if f['n'] > args.nmax:
            zu_gross += 1; continue
        p = parse_homfly(f['homfly'])
        if p is None:
            tot += 1; continue
        r_reg = eval_homfly(p, x, y)
        v, dt = _mit_frist(lambda: homfly_hecke(list(f['w']), f['n'],
                                                complex(a), complex(q), args.procs),
                           args.frist)
        zeit_pro_n.setdefault(f['n'], []).append(dt)
        if v is None:
            tot += 1
            print(f"  ⏱ {f['name'][:44]:<46} n={f['n']:<3} aborted after {dt:.0f}s")
            continue
        r_own = v * K
        nen = max(abs(r_reg), abs(r_own), 1e-12)
        rel = abs(r_reg - r_own) / nen
        if rel > schlimmster:
            schlimmster, schlimmster_name = rel, f['name']
        if rel < 1e-7:
            treffer += 1
        else:
            fehl += 1
            abweichungen.append({'name': f['name'], 'n': f['n'], 'c': f['c'],
                                 'rel': rel, 'regina': str(r_reg), 'eigen': str(r_own)})
            print(f"  🔴 {f['name'][:44]:<46} n={f['n']:<3} c={f['c']} rel={rel:.2e}")

    print(f"\n{'='*66}")
    print(f"  Hits               {treffer}")
    print(f"  Deviations         {fehl}")
    print(f"  not computed       {tot} (deadline/unparsable) · {zu_gross} with n>{args.nmax}")
    print(f"  largest rel. dev.  {schlimmster:.2e}  ({schlimmster_name[:40]})")
    print(f"\n  Run time per strand number (median):")
    for n in sorted(zeit_pro_n):
        v = sorted(zeit_pro_n[n])
        print(f"    n={n:<3} {len(v):>4} cases  median {v[len(v)//2]:7.2f}s  max {v[-1]:7.2f}s")
    if args.json:
        json.dump({'treffer': treffer, 'fehl': fehl, 'nicht_gerechnet': tot,
                   'zu_gross': zu_gross, 'schlimmste_rel': schlimmster,
                   'abweichungen': abweichungen,
                   'zeit_pro_n': {str(k): v for k, v in zeit_pro_n.items()}},
                  open(args.json, 'w'), indent=1)
    return 0 if fehl == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
