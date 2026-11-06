#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Braid→.fseries bridge: twisted torus T(p,q;r,s) -> dense torus-wound centerline ->
Fourier series (.fseries format of the pipeline). With DFT roundtrip self-check (reconstructed fit error),
so that we know BEFORE the pipeline whether the harmonics count carries the curve (and thus the crossings).
Co-authors: Werner Stanggassinger + Claude Code (Anthropic). 2026-06-29."""
import numpy as np, math, sys, os
import ondemand_e2e as m

def centerline_dense(p,q,r,s, N=8192, R=1.0, rt=0.45, bump=0.16):
    """N points evenly over the WHOLE closed curve (p turns = p-cycle)."""
    word=m.braid_word(p,q,r,s); L=len(word)
    # first determine slot sequence per (loop,k), then finely sample evenly into N
    # build list of (phi_frac_within_loop, slot_before, gen) segments
    segs=[]; cur=1; start=1
    for loop in range(p):
        for k,i in enumerate(word):
            segs.append((k, i, cur))
            if   cur==i:   cur=i+1
            elif cur==i+1: cur=i
        if cur==start and loop>0: break
    nloops=len(segs)//L if L else 1
    total=len(segs)
    pts=np.empty((N,3))
    for idx in range(N):
        u=idx/N*total                # global segment parameter [0,total)
        si=int(u)%total; a=u-int(u)
        k,i,cur=segs[si]
        phi=2*math.pi*((k+a)/L)
        theta=2*math.pi*(cur-1)/p; rho=rt
        if cur==i or cur==i+1:
            if cur==i:   th0=2*math.pi*(i-1)/p; th1=2*math.pi*i/p;     rho=rt+bump*math.sin(math.pi*a)
            else:        th0=2*math.pi*i/p;     th1=2*math.pi*(i-1)/p; rho=rt-bump*math.sin(math.pi*a)
            theta=th0+(th1-th0)*a
        cph,sph=math.cos(phi),math.sin(phi)
        pts[idx]=[R*cph+rho*math.cos(theta)*cph, R*sph+rho*math.cos(theta)*sph, rho*math.sin(theta)]
    return pts

def fit_fseries(pts, n_harm):
    """Real Fourier series: offset + Σ a_j cos(j t) + b_j sin(j t), t=2π k/N."""
    N=len(pts); t=np.arange(N)*2*math.pi/N
    offset=pts.mean(axis=0)
    a=np.zeros((n_harm,3)); b=np.zeros((n_harm,3))
    for j in range(1,n_harm+1):
        cj=np.cos(j*t); sj=np.sin(j*t)
        a[j-1]=(2.0/N)*(pts*cj[:,None]).sum(axis=0)
        b[j-1]=(2.0/N)*(pts*sj[:,None]).sum(axis=0)
    return offset,a,b

def reconstruct(offset,a,b,N=4096):
    t=np.linspace(0,2*math.pi,N,endpoint=False)
    out=np.tile(offset,(N,1))
    for j in range(1,a.shape[0]+1):
        out+=np.outer(np.cos(j*t),a[j-1])+np.outer(np.sin(j*t),b[j-1])
    return out

def write_fseries(path,offset,a,b):
    with open(path,'w') as f:
        f.write("% braid->fseries (Twisted-Torus on-demand) — Stanggassinger+Claude Code 2026-06-29\n")
        f.write(f"{offset[0]:.10g} {offset[1]:.10g} {offset[2]:.10g}\n")
        for j in range(a.shape[0]):
            f.write(f"{a[j,0]:.10g} {b[j,0]:.10g} {a[j,1]:.10g} {b[j,1]:.10g} {a[j,2]:.10g} {b[j,2]:.10g}\n")

def roundtrip_err(pts,offset,a,b):
    """relative max. point error reconstruct vs original (same N)."""
    rec=reconstruct(offset,a,b,N=len(pts))
    bbox=np.ptp(pts,axis=0).max()
    return float(np.linalg.norm(rec-pts,axis=1).max()/bbox)

if __name__=='__main__':
    p,q,r,s=(int(x) for x in sys.argv[1:5])
    n_harm=int(sys.argv[5]) if len(sys.argv)>5 else max(64,4*p)
    out=sys.argv[6] if len(sys.argv)>6 else 'myknot.fseries'
    pts=centerline_dense(p,q,r,s,N=8192)
    offset,a,b=fit_fseries(pts,n_harm)
    err=roundtrip_err(pts,offset,a,b)
    write_fseries(out,offset,a,b)
    print(f"T({p},{q},{r},{s}) c={m.genus_formula(p,q,r,s)[1]} -> {out}  n_harm={n_harm}  roundtrip_relerr={err:.4f}")
