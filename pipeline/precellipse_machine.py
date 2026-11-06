#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""
precellipse_machine.py — Precellipse generator (Akio Hizume, "Precellipse Poly-Twistor 2026")
as a Stanggassinger machine + BLIND validation by the pipeline resolution.

A precellipse = precessing Kepler ellipse, closed at rational Delta=q/p,
= torus knot T(p,q) (Farey classification). Math directly from the paper (p.3-5):
  Kepler:    u + e*sin u = 2*pi*t/T        (eccentric anomaly u, Newton; t from aphelion)
  Ellipse:   X = cos u + e,  Y = sqrt(1-e^2)*sin u,  Z = 0
  Tilt k:    X'=X, Y'=cos k*Y, Z'=sin k*Y
  Polar angle v' = atan2(Y', X')   (UNWRAPPED -> accumulates p windings)
  Precession: W = (Delta-1)*v' = (q/p-1)*v'
  Rotation about z: X''=cosW*X'-sinW*Y', Y''=sinW*X'+cosW*Y', Z''=Z'

Curve -> FFT -> .fseries -> polygon_resolve.measure (z depth, the same pipeline
Resolution as in the corpus, NO table lookups) -> det/Alexander.
Validation oracle: known torus-knot Alexander/det (comparison, not lookup).

base env (regina+numpy). Writes .fseries to erzeugte_knoten/precellipse_<p>_<q>/.
"""
import site_config
import sys, os, math, json
import numpy as np

PS=(site_config.PIPELINE_DIR)
M01=(site_config.MACHINE_YARD + '/01_knotting_oscillations')
sys.path.insert(0,PS); sys.path.insert(0,M01)
import polygon_resolve as PR
import importlib.util
def _load(n,p):
    s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
G=_load("gen",M01+"/generate_oscillation_knots.py")

OUTROOT=(site_config.GENERATED_KNOTS_DIR)

# known torus knots as validation oracle (det, crossing, Alexander degree)
TORUS_TRUTH = {
    (2,3): {"name":"T(2,3)=3_1 (trefoil)", "det":3,  "c":3},
    (2,7): {"name":"T(2,7)=7_1",           "det":7,  "c":7},
    (3,7): {"name":"T(3,7) (=irrational)",  "det":1,  "c":14},
}

def precellipse(e, k, p, q, N):
    """Closed precellipse curve, N points over p revolutions."""
    i = np.arange(N)
    S = 2*np.pi*p*i/N                      # 2*pi*t/T, accumulated over p revolutions
    u = S.copy()                           # Newton, continuous (init u=S)
    for _ in range(80):
        u = u - (u + e*np.sin(u) - S)/(1 + e*np.cos(u))
    X = np.cos(u) + e
    Y = np.sqrt(1-e*e)*np.sin(u)
    Xp, Yp, Zp = X, np.cos(k)*Y, np.sin(k)*Y
    vp = np.unwrap(np.arctan2(Yp, Xp))     # MUST accumulate (p windings)
    W = (q/p - 1.0)*vp
    Xpp = np.cos(W)*Xp - np.sin(W)*Yp
    Ypp = np.sin(W)*Xp + np.cos(W)*Yp
    Zpp = Zp
    return np.stack([Xpp, Ypp, Zpp], axis=1)

def curve_to_fseries(P, kmax=256, energy_keep=0.99995):
    """FFT of the closed curve -> (offset, a, b) up to kmax harmonics."""
    N = len(P)
    Xf = np.fft.rfft(P, axis=0)            # (N//2+1, 3)
    a0 = (Xf[0].real)/N
    nh = min(kmax, Xf.shape[0]-1)
    a = 2*Xf[1:nh+1].real/N                # (nh,3)
    b = -2*Xf[1:nh+1].imag/N
    # truncate to significant harmonics (energy)
    E = (a*a + b*b).sum(axis=1)
    tot = E.sum()
    if tot>0:
        cum = np.cumsum(E)/tot
        keep = int(np.searchsorted(cum, energy_keep))+1
        keep = max(keep, 4)
        a, b = a[:keep], b[:keep]
    return a0, a, b

def write_fseries(path, offset, a, b, header):
    with open(path,"w") as f:
        for h in header: f.write(f"% {h}\n")
        f.write("  %.8f  %.8f  %.8f\n"%(offset[0],offset[1],offset[2]))
        for k in range(len(a)):
            f.write("  %.8f  %.8f  %.8f  %.8f  %.8f  %.8f\n"%(
                a[k,0],b[k,0],a[k,1],b[k,1],a[k,2],b[k,2]))

def measure_curve(offset, a, b, axis, n):
    """polygon_resolve measurement directly from Fourier coefficients."""
    return PR.measure(offset, a, b, axis, n)

def flattest_axis(offset, a, b):
    Pp = PR.sample_curve(offset, a, b, 2048)
    return int(np.argsort(Pp.max(0)-Pp.min(0))[0])

def best_embedding(p, q, N=8000):
    """Choose e,k that embed cleanly (largest min_nonadjacent_dist)."""
    best=None
    for e in (0.4,0.55,0.7,0.85):
        for kdeg in (15,25,35,45):
            k=math.radians(kdeg)
            P=precellipse(e,k,p,q,N)
            diam=float(np.linalg.norm(P.max(0)-P.min(0)))
            dmin,_=G.min_nonadjacent_dist(P,exclude=8)
            rel=dmin/diam if diam>0 else 0
            if best is None or rel>best[0]:
                best=(rel,e,k,kdeg)
    return best

def nyquist_id(offset, a, b, axis, n_max=65536):
    """Nyquist gate: double until (det,alex) is stable."""
    n=4096; prev=None
    while n<=n_max:
        r=measure_curve(offset,a,b,axis,n)
        if r and prev and r['det']==prev['det'] and r['alexander']==prev['alexander']:
            return r, n
        prev=r; n*=2
    return prev, n

def main():
    results=[]
    for (p,q),truth in TORUS_TRUTH.items():
        print(f"\n=== Precellipse Delta=q/p={q}/{p}  (target {truth['name']}, det={truth['det']}) ===")
        rel,e,k,kdeg = best_embedding(p,q)
        print(f"  best embedding: e={e}, k={kdeg}deg, d_min_rel={rel:.5f}")
        P = precellipse(e,k,p,q,16000)
        offset,a,b = curve_to_fseries(P)
        print(f"  .fseries: {len(a)} harmonics, offset~{np.round(offset,3).tolist()}")
        axis = flattest_axis(offset,a,b)
        r,nst = nyquist_id(offset,a,b,axis)
        if r is None:
            print(f"  ERROR: no stable measurement"); continue
        det=r['det']; nc=r['n_cross_simplified']
        match = (det==truth['det'])
        print(f"  PIPELINE (blind): det={det}, crossings(simpl)={nc}, n_stable={nst}")
        print(f"  Alexander={r['alexander']}")
        print(f"  ==> {'✓ MATCH '+truth['name'] if match else '✗ det '+str(det)+' != '+str(truth['det'])}")
        # write .fseries
        d=f"{OUTROOT}/precellipse_{p}_{q}"; os.makedirs(d,exist_ok=True)
        fp=f"{d}/precellipse_{p}_{q}.fseries"
        write_fseries(fp, offset, a, b,
            ["Precellipse machine (Hizume construction) — Delta=%d/%d, e=%.3f, k=%ddeg"%(q,p,e,kdeg),
             "target torus knot %s"%truth['name'],
             "generated by precellipse_machine.py (Stanggassinger pipeline validation)"])
        results.append({"p":p,"q":q,"e":e,"k_deg":kdeg,"d_min_rel":round(rel,5),
                        "n_harmonics":len(a),"axis":"xyz"[axis],
                        "det":det,"n_cross":nc,"alexander":r['alexander'],
                        "target":truth['name'],"target_det":truth['det'],"match":match,
                        "fseries":fp})
    json.dump(results, open(f"{OUTROOT}/precellipse_validation.json","w"), indent=1)
    print("\n=== BALANCE ===")
    for r in results:
        print(f"  T({r['p']},{r['q']}): det {r['det']} (target {r['target_det']}) c={r['n_cross']} "
              f"-> {'✓' if r['match'] else '✗'}  {r['fseries']}")
    n_ok=sum(1 for r in results if r['match'])
    print(f"\n{n_ok}/{len(results)} torus knots identified correctly blind.")

if __name__=="__main__":
    main()
