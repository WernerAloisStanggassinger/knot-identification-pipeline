#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""(2) Production tool: HOMFLY via Hecke trace, PARALLEL over lambda (performance version).
WRITTEN 2026-06-16 (runs left untouched) — RUN only when the machine is free.
Method validated (see MUTATION_SEPARATION_feasibility_2026_06_16.md): q-seminormal reps (Hecke relation 1e-16),
Ocneanu weights W_λ=∏(a·qᶜ−a⁻¹q⁻ᶜ)/(qʰ−q⁻ʰ) calibrated against Sage, framing F=a⁻¹.
Speedup: the 77 irreps (n=12) are independent -> multiprocessing.Pool. Wall time ≈ largest single λ instead of the sum.
CLI:  python hecke_homfly_fast.py <braid_csv> <c_cable> <a_re> <q_re> [nprocs]   (or --selftest)
"""
import sys, numpy as np, scipy.sparse as sp
from multiprocessing import Pool

# ---- tableau / rep machinery (identical to the verified version) ----
def shape_cells(lam): return [(r,c) for r,row in enumerate(lam) for c in range(row)]
def syt(lam):
    cells=shape_cells(lam); n=len(cells); out=[]
    def rec(a):
        if len(a)==n: out.append(dict(a)); return
        v=len(a)+1
        for (r,c) in cells:
            if (r,c) in a: continue
            if (c==0 or (r,c-1) in a) and (r==0 or (r-1,c) in a):
                b=dict(a); b[(r,c)]=v; rec(b)
    rec({})
    res=[]
    for d in out:
        bo={v:b for b,v in d.items()}; res.append(tuple(bo[v] for v in range(1,n+1)))
    return [list(t) for t in dict.fromkeys(res)]
def content(b): r,c=b; return c-r
def qint(d,q): return (q**d-q**-d)/(q-q**-1)
def hecke_gen_sparse(T,idx,i,q):
    n=len(T); rows=[];cols=[];vals=[]; done=set()
    for k,t in enumerate(T):
        if k in done: continue
        d=content(t[i])-content(t[i-1])
        if d==1: rows+=[k];cols+=[k];vals+=[q]; done.add(k); continue
        if d==-1: rows+=[k];cols+=[k];vals+=[-q**-1]; done.add(k); continue
        ts=list(t); ts[i-1],ts[i]=ts[i],ts[i-1]; ks=idx[tuple(ts)]
        rows+=[k,ks,k,ks]; cols+=[k,ks,ks,k]
        vals+=[(q-q**-1)/(1-q**(-2*d)),(q-q**-1)/(1-q**(2*d)),1.0,1-1/qint(d,q)**2]
        done.add(k);done.add(ks)
    return sp.csr_matrix((vals,(rows,cols)),shape=(n,n),dtype=complex)
def cable_braid(word,c):
    def bs(i):
        base=c*(i-1); w=[]
        for st in range(c):
            for k in range(c): w.append(base+(c+st-k))
        return w
    out=[]
    for g in word:
        sw=bs(abs(g)); out+= sw if g>0 else [-x for x in reversed(sw)]
    return out
def partitions(n,m=None):
    if m is None:m=n
    if n==0: yield[];return
    for k in range(min(n,m),0,-1):
        for r in partitions(n-k,k): yield [k]+r
def hc(lam):
    out=[]
    for (r,c) in shape_cells(lam):
        arm=lam[r]-c-1; leg=sum(1 for rr in range(r+1,len(lam)) if lam[rr]>c); out.append((c-r,arm+leg+1))
    return out
def W(lam,a,q):
    P=1.0+0j
    for cont,h in hc(lam): P*=(a*q**cont-a**-1*q**-cont)/(q**h-q**-h)
    return P

_JOB=None
def _worker(lam):
    word,q,a=_JOB
    T=syt(lam); idx={tuple(t):k for k,t in enumerate(T)}; dim=len(T); sc=q-q**-1
    I=sp.identity(dim,dtype=complex,format='csr'); Gp={};Gi={}; M=np.eye(dim,dtype=complex)
    for g in word:
        i=abs(g)
        if i not in Gp: Gp[i]=hecke_gen_sparse(T,idx,i,q); Gi[i]=(Gp[i]-sc*I).tocsr()
        G=Gi[i] if g<0 else Gp[i]; M=(G.T@M.T).T
    return W(lam,a,q)*complex(np.trace(M))
def _init(job):
    global _JOB; _JOB=job
def homfly_hecke(word,n,a,q,nprocs=8):
    """assembled HOMFLY (×common norm K) of the closure of word∈B_n, parallel over λ."""
    e=sum(1 if g>0 else -1 for g in word); lams=list(partitions(n))
    with Pool(nprocs,initializer=_init,initargs=((word,q,a),)) as pool:
        parts=pool.map(_worker,lams)
    return a**(-e)*sum(parts)

if __name__=="__main__":
    if "--selftest" in sys.argv:
        a,q=complex(1.3),complex(1.2)
        # Trefoil σ1^3 + figure-8 against the validated reference (K·a^-e·ΣWχ):
        for nm,w,n in [("Trefoil",[1,1,1],2),("Figure-8",[1,-2,1,-2],3)]:
            print(f"  {nm}: ΣWχ·a^-e = {homfly_hecke(w,n,a,q,4):.6f}")
        # 2-cable KT vs Conway (must be equal):
        kt=[-1,2,2,-3,-3,2,1,-2,-2,3,-2,3,-2]; cw=[-1,2,-1,2,-1,3,-2,-2,-1,3,3]
        pk=homfly_hecke(cable_braid(kt,2),8,a,q,4); pc=homfly_hecke(cable_braid(cw,2),8,a,q,4)
        print(f"  2-cable KT/Conway |Δ|={abs(pk-pc):.2e} {'EQUAL ✓' if abs(pk-pc)<1e-6 else 'FAIL'}")
        sys.exit(0)
    bw=[int(x) for x in sys.argv[1].split(",")]; c=int(sys.argv[2])
    a=complex(float(sys.argv[3])); q=complex(float(sys.argv[4])); npr=int(sys.argv[5]) if len(sys.argv)>5 else 8
    cb=cable_braid(bw,c); n=(max(abs(g) for g in bw)+1)*c
    print(f"P = {homfly_hecke(cb,n,a,q,npr):.8f}  (n={n} strands, {len(cb)} crossings, {npr} procs)")
