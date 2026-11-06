#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""SIEVE WIDE (base env). Machine 115, k=2..5 -> rotary reflection of order 4,6,8,10.
sigma acts as a parameter shift => POSITIVELY amphichiral by construction.
Filter is ALGEBRAIC and cheap: nabla from Alexander, then integral Kawauchi test.
  passes -> uninteresting (no counter-case, whatever the symmetry)
  fails -> CANDIDATE: sieved at N=1440, immediately re-checked at N=4320 (double resolution),
                  PD stored -> stage B (geometry) decides.
22 permanent workers, each writes its own JSONL (incremental, a crash costs only the rest
of its chunk). No multiprocessing pool (footgun 13).
"""
import json, os, sys, time, signal
import numpy as np
import importlib.util
sp = importlib.util.spec_from_file_location('a', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'stufe_a_pd.py'))
A = importlib.util.module_from_spec(sp); sp.loader.exec_module(A)
from sympy import Poly, symbols, sympify, expand, factor_list
import re
z, w = symbols('z w')

def nabla(alex):
    e = re.sub(r'\s+', '', re.sub(r'(\d)\s*x', r'\1*x', alex.replace('^', '**')))
    p = Poly(sympify(e), symbols('x')); c = p.all_coeffs()[::-1]; d = len(c)-1
    if d % 2 or d == 0: return None
    m = d//2
    if any(c[m+i] != c[m-i] for i in range(m+1)): return None
    pol = [Poly(2, w), Poly(w, w)]
    for i in range(2, m+1): pol.append(Poly(w*pol[-1].as_expr()-pol[-2].as_expr(), w))
    ex = Poly(c[m], w)
    for i in range(1, m+1): ex = Poly(ex.as_expr()+c[m+i]*pol[i].as_expr(), w)
    nb = Poly(expand(ex.as_expr().subs(w, z**2+2)), z)
    if nb.TC() == -1: nb = -nb
    return nb if nb.TC() == 1 else None

def kawauchi_ok(nb):
    """True = factorisation exists (uninteresting). False = COUNTER-CASE."""
    if nb.degree() == 0: return True
    lc = nb.coeff_monomial(z**nb.degree())
    d2 = nb.degree()//2
    if (-1)**d2 * lc < 0: return False                      # sign impossible
    import math
    r = math.isqrt(abs(int(lc)))
    if r*r != abs(int(lc)): return False                    # leading coefficient not a square (Fact 4.4)
    _, facs = factor_list(nb.as_expr(), z)
    def norm(p):
        p = Poly(p, z); _, pp = p.primitive()
        return -pp if pp.LC() < 0 else pp
    facs = [(norm(p), e) for p, e in facs]
    used = [False]*len(facs); fp = Poly(1, z)
    for i, (p, e) in enumerate(facs):
        if used[i]: continue
        mp = norm(Poly(p.as_expr().subs(z, -z), z))
        if mp == p:
            if e % 2: return False
            used[i] = True; fp *= p**(e//2)
        else:
            j = next((kk for kk, (q, _) in enumerate(facs) if q == mp and not used[kk]), None)
            if j is None or facs[j][1] != e: return False
            used[i] = used[j] = True; fp *= p**e
    for s in (1, -1):
        f = fp*s
        if Poly(expand(f.as_expr()*f.as_expr().subs(z, -z)), z) == nb: return True
    return False

def eval_one(k, seed, nxy, nz, dc, NS):
    import regina
    P = A.forge(k, seed, nxy, nz, dc)
    if A.sigma_ok(P, k) > 1e-11: return dict(rej='sigma')
    idx = np.linspace(0, len(P), NS, endpoint=False).astype(int)
    sg = A.signed_gauss(P[idx])
    if sg is None: return dict(rej='tangenz')
    L = regina.Link.fromSignedGauss(sg); L.simplify()
    if L.size() == 0: return dict(rej='unknoten')
    al = str(L.alexander())
    return dict(c=int(L.size()), alex=al, pd=str(L.pd()), c_raw=len(sg)//3)

def worker(wid, jobs, path):
    import regina
    f = open(path, 'w')
    for (k, seed, nxy, nz, dc) in jobs:
        try:
            r = eval_one(k, seed, nxy, nz, dc, 1440)
        except Exception as e:
            r = dict(rej=f'{type(e).__name__}')
        if 'rej' in r:
            f.write(json.dumps(dict(k=k, seed=seed, nxy=nxy, nz=nz, dc=dc, **r))+'\n'); f.flush(); continue
        try:
            nb = nabla(r['alex'])
            r['kaw'] = None if nb is None else bool(kawauchi_ok(nb))
            r['nabla'] = None if nb is None else str(nb.as_expr())
        except Exception as e:
            r['kaw'] = None; r['nabla_err'] = f'{type(e).__name__}'
        if r.get('kaw') is False:
            try:
                hi = eval_one(k, seed, nxy, nz, dc, 4320)
                r['confirm_4320'] = (hi.get('alex') == r['alex'])
                r['pd_hi'] = hi.get('pd'); r['c_hi'] = hi.get('c')
            except Exception as e:
                r['confirm_4320'] = None
        f.write(json.dumps(dict(k=k, seed=seed, nxy=nxy, nz=nz, dc=dc, **r))+'\n'); f.flush()
    f.close(); os._exit(0)

if __name__ == '__main__':
    NW = int(os.environ.get('NW', '22'))
    SEEDS = int(os.environ.get('SEEDS', '1500'))
    KS = [int(x) for x in os.environ.get('KS', '2,3,4,5').split(',')]
    jobs = [(k, s, nxy, nz, dc) for k in KS for s in range(1, SEEDS+1)
            for nxy in (2, 3, 4) for nz in (2, 3) for dc in (0.45, 0.55, 0.7)]
    print(f'=== SIEVE WIDE: {len(jobs)} candidates, k={KS}, {NW} workers ===', flush=True)
    os.makedirs('breit_out', exist_ok=True)
    t0 = time.time(); pids = []
    for i in range(NW):
        mine = jobs[i::NW]
        pid = os.fork()
        if pid == 0: worker(i, mine, f'breit_out/w{i:02d}.jsonl')
        pids.append(pid)
    print(f'{NW} workers started, ~{len(jobs)//NW} candidates each', flush=True)
    done = 0
    for p in pids:
        os.waitpid(p, 0); done += 1
        print(f'  worker done {done}/{NW}  {time.time()-t0:.0f}s', flush=True)
    print(f'ALL DONE {time.time()-t0:.0f}s', flush=True)
