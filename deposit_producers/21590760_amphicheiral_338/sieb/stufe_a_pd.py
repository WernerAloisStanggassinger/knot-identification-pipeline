#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""STAGE A (base env: regina+numpy): curve -> signed Gauss -> regina -> PD + det + Alexander.

Goal of the hunt: positively amphichiral knots with rotary reflection of order 2^a, a>=2,
WITHOUT involution and WITHOUT negative amphichirality -> there the integral Kawauchi
conjecture is false (EHW 2011, so far exactly ONE published counter-case).
Machine 115 delivers the class by construction: sigma = Rz(pi/k) o (z->-z) acts on the
curve as a parameter shift P(t)->P(t+pi/k) => knot orientation PRESERVED => POSITIVELY
amphichiral, order 2k. For k=2: order 4 = exactly EHW's class.
signed_gauss() is taken verbatim from pipeline_scripts/polygon_resolve.py (only the
module import of regina there prevents the direct import) - pipeline code stays untouched.
"""
import site_config
import sys, os, json, pickle, signal, time
import numpy as np
sys.path.insert(0, (site_config.ROOT + '/maschinen_ideen'))
sys.path.insert(0, (site_config.ROOT + '/freeze_manifest'))
from machine_common import embedding_report

N = 4320

def forge(k, seed, n_xy=3, n_z=2, decay=0.55):
    """identical to machine 115"""
    rng = np.random.default_rng(seed)
    t = np.linspace(0, 2*np.pi, N, endpoint=False)
    ns_xy = [1]
    for j in range(1, n_xy+1): ns_xy += [1+2*k*j, 1-2*k*j]
    c = np.zeros(N, complex)
    for n in ns_xy:
        amp = (1.0 if n == 1 else decay**abs((n-1)//(2*k)))
        a = amp*(rng.normal()+1j*rng.normal())*(0.45 if n != 1 else 1.2)
        c += a*np.exp(1j*n*t)
    z = np.zeros(N)
    for j in range(n_z):
        n = k+2*k*j
        z += (decay**j)*0.7*(rng.normal()*np.cos(n*t)+rng.normal()*np.sin(n*t))
    return np.column_stack([c.real, c.imag, z])

def sigma_ok(P, k, tol=1e-12):
    """Certificate of the construction: sigma(P(t)) == P(t+pi/k) bit-exactly."""
    a = np.pi/k
    R = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, -1.0]])
    return float(np.abs(P@R.T - np.roll(P, -(N//(2*k)), axis=0)).max())

# --- verbatim from polygon_resolve.signed_gauss (pipeline untouched) ---
def signed_gauss(P):
    n = len(P); S = P[:, :2]; E = np.roll(S, -1, axis=0); D = E-S
    passages = []; cid = 0; CH = 512
    for i0 in range(0, n, CH):
        i = np.arange(i0, min(i0+CH, n))
        for j0 in range(i0, n, CH):
            j = np.arange(j0, min(j0+CH, n))
            I, J = np.meshgrid(i, j, indexing="ij")
            mask = (J > I+1) & ~((I == 0) & (J == n-1))
            if not mask.any(): continue
            Ii, Jj = I[mask], J[mask]
            d = (D[Ii,0]*D[Jj,1]-D[Ii,1]*D[Jj,0]); ok = np.abs(d) > 1e-15
            Ii, Jj, d = Ii[ok], Jj[ok], d[ok]
            R = S[Jj]-S[Ii]
            t_ = (R[:,0]*D[Jj,1]-R[:,1]*D[Jj,0])/d
            u_ = (R[:,0]*D[Ii,1]-R[:,1]*D[Ii,0])/d
            hit = (t_ > 0)&(t_ < 1)&(u_ > 0)&(u_ < 1)
            for ii, jj, tt, uu in zip(Ii[hit], Jj[hit], t_[hit], u_[hit]):
                if min(tt, 1-tt, uu, 1-uu) < 1e-9: return None
                z1 = P[ii,2]+tt*(P[(ii+1)%n,2]-P[ii,2]); z2 = P[jj,2]+uu*(P[(jj+1)%n,2]-P[jj,2])
                if abs(z1-z2) < 1e-12: return None
                od, ud = (ii, jj) if z1 > z2 else (jj, ii)
                s = 1 if (D[od,0]*D[ud,1]-D[od,1]*D[ud,0]) > 0 else -1
                cid += 1
                passages.append((ii+tt, cid, z1 > z2, s)); passages.append((jj+uu, cid, z2 > z1, s))
    passages.sort()
    return "".join(f"{'O' if o else 'U'}{c}{'+' if s>0 else '-'}" for _, c, o, s in passages)

def one(k, seed, n_xy, n_z, decay):
    import regina
    P = forge(k, seed, n_xy, n_z, decay)
    r = dict(k=k, seed=seed, n_xy=n_xy, n_z=n_z, decay=decay)
    r['sigma_err'] = sigma_ok(P, k)
    if r['sigma_err'] > 1e-12: r['reject'] = 'sigma nicht bitgenau'; return r
    emb = embedding_report(P)
    if not emb.get('embedded'): r['reject'] = 'embedding gate'; return r
    sg = signed_gauss(P)
    if sg is None: r['reject'] = 'tangenz'; return r
    r['c_raw'] = len(sg)//3
    L = regina.Link.fromSignedGauss(sg)
    L.simplify()
    r['c'] = int(L.size())
    if r['c'] == 0: r['reject'] = 'unknoten'; return r
    r['pd'] = str(L.pd())
    try:
        r['alex'] = str(L.alexander())
    except Exception as e:
        r['alex'] = None; r['alex_err'] = str(e)[:60]
    try:
        h = L.homfly()
        r['det'] = None
    except Exception: pass
    try:
        import regina as rg
        tri = L.complement(); tri.simplify()
        r['tets'] = int(tri.size())
    except Exception: pass
    return r

def run(args, timeout=300):
    rd, wr = os.pipe(); pid = os.fork()
    if pid == 0:
        os.close(rd); signal.alarm(timeout)
        try: res = one(*args)
        except Exception as e: res = dict(k=args[0], seed=args[1], reject=f'{type(e).__name__}:{str(e)[:70]}')
        try:
            with os.fdopen(wr, 'wb') as f: pickle.dump(res, f)
        finally: os._exit(0)
    os.close(wr)
    try:
        with os.fdopen(rd, 'rb') as f: d = f.read()
        os.waitpid(pid, 0)
        return pickle.loads(d) if d else dict(k=args[0], seed=args[1], reject='kind gestorben')
    except Exception as e:
        return dict(k=args[0], seed=args[1], reject=f'{type(e).__name__}')

if __name__ == '__main__':
    import itertools
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    s0, s1 = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (1, 9)
    MAX = int(os.environ.get('PAR', '22'))
    jobs = [(k, s, nxy, nz, dc) for s in range(s0, s1)
            for nxy in (2, 3, 4) for nz in (2, 3) for dc in (0.45, 0.55, 0.7)]
    print(f'=== STAGE A: k={k} (rotary reflection of order {2*k}), {len(jobs)} candidates, {MAX} in parallel ===', flush=True)
    t0 = time.time(); running = {}; out = []
    def reap(block):
        while running and (block or len(running) >= MAX):
            pid, _ = os.waitpid(-1, 0 if block else os.WNOHANG)
            if pid == 0: return
            fd = running.pop(pid)
            try:
                with os.fdopen(fd, 'rb') as f: d = f.read()
                out.append(pickle.loads(d) if d else {'reject': 'leer'})
            except Exception as e: out.append({'reject': f'{type(e).__name__}'})
            if len(out) % 50 == 0: print(f'  {len(out)}/{len(jobs)} {time.time()-t0:.0f}s', flush=True)
            if not block: return
    for a in jobs:
        reap(False)
        while len(running) >= MAX: reap(False); time.sleep(0.05)
        rd, wr = os.pipe(); pid = os.fork()
        if pid == 0:
            os.close(rd); signal.alarm(600)
            try: res = one(*a)
            except Exception as e: res = dict(k=a[0], seed=a[1], n_xy=a[2], n_z=a[3], decay=a[4],
                                              reject=f'{type(e).__name__}:{str(e)[:70]}')
            try:
                with os.fdopen(wr, 'wb') as f: pickle.dump(res, f)
            finally: os._exit(0)
        os.close(wr); running[pid] = rd
    while running: reap(True)
    json.dump(out, open(f'stufe_a_k{k}.json', 'w'), indent=1, default=str)
    import collections
    c = collections.Counter(r.get('reject', 'OK') for r in out)
    print(f'\ndone {time.time()-t0:.0f}s')
    for a, b in c.most_common(): print(f'  {a}: {b}')
    ok = [r for r in out if 'pd' in r]
    print(f'  usable PDs: {len(ok)}, c distribution {sorted(collections.Counter(r["c"] for r in ok).items())}')
