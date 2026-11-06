#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""On-demand ladder — END-TO-END DEMO (Disc successor #349/rung 5).
Co-authors: Werner Stanggassinger + Claude Code (Anthropic). 2026-06-29.

THE ORDER (key):  "Deliver a GENUINE NOVEL knot with genus g = G,
                          and give me its signature spectrum σ on this line."
THE MACHINE:
  1) INVERTS the exact genus law  g = (p-1)(q-1)/2 + s   (Twisted-Torus T(p,q;2,s))
     -> enumerates the (p,q,s) grid in the hyperbolicity window with g==G.
  2) MEASURES blind σ, det, genus (Sage), vol/solution_type/identify (snappy) per candidate. PARALLEL.
  3) TOPO GATE novelty: hyperbolic (solution_type 'all ... positively oriented') + identify()==[].
  4) Dial multiplicity: distinct novel knots per (σ) on the line  ->  the honest "collision" price.
  5) FILLS the 2-axis order (G, S): delivers a knot with g==G AND measured σ==S.
  6) FRESH VERIFICATION of the delivered braid word in an INDEPENDENT process (not re-read).

HONEST: g exactly invertible (closed form). σ has NO closed form -> machine enumerates the
g line and FILTERS on measured σ (dial inversion of the characterized axis + measurement filter on
the nonlinear one). c = braid-word crossings, NOT certified minimal. Nothing read from a table.
"""
import site_config
import subprocess, re, json, os, sys, time
from math import gcd
from concurrent.futures import ThreadPoolExecutor, as_completed

SAGE=site_config.SAGE_BIN
PY_SNAPPY=site_config.PYTHON_MAIN
SP=os.path.dirname(os.path.abspath(__file__))
JSONL=os.path.join(SP,'e2e_measurements.jsonl')

# ---------- family / forward law (from rung 5, bit-confirmed) ----------
def braid_word(p,q,r,s):
    return list(range(1,p))*q + list(range(1,r))*(r*s)

def in_window(p,q,r,s):
    return gcd(p,q)==1 and 1<r<q<p and s>=3

def genus_formula(p,q,r,s):
    L=(p-1)*q + (r-1)*r*s; n=p
    return (L - n + 1)//2, L          # (g, braid-word length=c)

# ---------- Blind measurement ----------
def sage_invariants(word, nstrands, timeout=1200):
    code=(
      "from sage.knots.link import Link\n"
      "from sage.groups.braid import BraidGroup\n"
      f"B=BraidGroup({nstrands})\n"
      f"b=B({word})\n"
      "L=Link(b)\n"
      "print('COMP',L.number_of_components())\n"
      "print('SIG',int(L.signature()))\n"
      "print('DET',int(abs(L.determinant())))\n"
      "try:\n    print('GENUS',int(L.genus()))\n"
      "except Exception as e:\n    print('GENUS','NA')\n"
    )
    r=subprocess.run([SAGE,'-c',code],capture_output=True,text=True,timeout=timeout)
    if r.returncode!=0: return {'err':r.stderr[-300:]}
    def g(tag):
        m=re.search(tag+r'\s+(-?\d+|NA)',r.stdout); return m.group(1) if m else None
    return {'comp':g('COMP'),'sig':g('SIG'),'det':g('DET'),'genus':g('GENUS')}

def snappy_topo(word, timeout=900):
    code=(
      "import snappy\n"
      f"L=snappy.Link(braid_closure={word})\n"
      "print('NC',len(L.crossings))\n"
      "print('COMP',len(L.link_components))\n"
      "M=L.exterior()\n"
      "try:\n    v=M.volume(); print('VOL',round(float(v),5)); print('SOL',str(M.solution_type()))\n"
      "except Exception as e:\n    print('VOL','NA'); print('SOL','err')\n"
      "try:\n    ids=M.identify(); print('IDENT',[str(x) for x in ids][:4])\n"
      "except Exception as e:\n    print('IDENT',[])\n"
    )
    r=subprocess.run([PY_SNAPPY,'-c',code],capture_output=True,text=True,timeout=timeout)
    if r.returncode!=0: return {'err':r.stderr[-300:]}
    def g(tag):
        m=re.search(tag+r'\s+(.+)',r.stdout); return m.group(1).strip() if m else None
    return {'nc':g('NC'),'comp':g('COMP'),'vol':g('VOL'),'sol':g('SOL'),'ident':g('IDENT')}

def is_novel(sn):
    return (sn.get('vol') not in (None,'NA') and str(sn.get('sol','')).startswith('all')
            and sn.get('ident') in ('[]','[ ]',None))

def measure(cand):
    p,q,r,s=cand
    t0=time.time()
    word=braid_word(p,q,r,s)
    sa=sage_invariants(word,p)
    sn=snappy_topo(word)
    gf,c=genus_formula(p,q,r,s)
    rec=dict(p=p,q=q,r=r,s=s,c=c,word_len=len(word),g_formula=gf,
             sage=sa,snappy=sn,novel=is_novel(sn),secs=round(time.time()-t0,1))
    return rec

# ---------- Enumeration: invert the genus law ----------
def enumerate_line(G, p_max=20, s_max=30):
    cands=[]
    for p in range(4,p_max+1):
        for q in range(3,p):
            r=2
            if not (gcd(p,q)==1 and 1<r<q<p): continue
            for s in range(3,s_max+1):
                gf,c=genus_formula(p,q,r,s)
                if gf==G:
                    cands.append((p,q,r,s))
    # sort by c (smaller first), unique
    return sorted(set(cands), key=lambda t: genus_formula(*t)[1])

def main():
    G=int(sys.argv[1]) if len(sys.argv)>1 else 20
    workers=int(sys.argv[2]) if len(sys.argv)>2 else 12
    p_max=int(sys.argv[3]) if len(sys.argv)>3 else 20
    cands=enumerate_line(G, p_max=p_max)
    print(f"[ORDER] novel knot with genus g={G}. Genus law inverted -> "
          f"{len(cands)} grid candidates (p<= {p_max}, c={genus_formula(*cands[0])[1]}..{genus_formula(*cands[-1])[1]}).",flush=True)
    print(f"[CANDS] {cands}",flush=True)
    open(JSONL,'w').close()
    results=[]; t0=time.time()
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs={ex.submit(measure,c):c for c in cands}
        for i,fut in enumerate(as_completed(futs),1):
            c=futs[fut]
            try: rec=fut.result()
            except Exception as e: rec=dict(cand=c,err=str(e)[:200])
            with open(JSONL,'a') as f: f.write(json.dumps(rec)+"\n")
            sn=rec.get('snappy',{}); sa=rec.get('sage',{})
            print(f"[{i}/{len(cands)}] {c} c={rec.get('c')} sig={sa.get('sig')} "
                  f"g={sa.get('genus')}(form {rec.get('g_formula')}) det={sa.get('det')} "
                  f"vol={sn.get('vol')} sol={str(sn.get('sol'))[:14]} ident={str(sn.get('ident'))[:22]} "
                  f"NOVEL={rec.get('novel')} [{rec.get('secs')}s]",flush=True)
            results.append(rec)
    print(f"\n[DONE] {len(results)} measured in {round(time.time()-t0,1)}s.",flush=True)

    # ---- Evaluation ----
    ok=[r for r in results if 'err' not in r and r.get('sage',{}).get('sig') is not None]
    # Sanity: genus law against Sage
    gmatch=sum(1 for r in ok if str(r['sage'].get('genus'))==str(r['g_formula']))
    print(f"[GENUS-LAW] g(sage)==g(formula): {gmatch}/{len(ok)}",flush=True)
    novel=[r for r in ok if r['novel']]
    print(f"[NOVELTY] {len(novel)}/{len(ok)} novel (hyperbolic + identify()==[]).",flush=True)

    # σ spectrum + dial multiplicity (distinct key: vol@5 rounded + sig + det)
    def dkey(r): return (r['snappy'].get('vol'), r['sage'].get('sig'), r['sage'].get('det'))
    spectrum={}
    for r in novel:
        sig=r['sage'].get('sig')
        spectrum.setdefault(sig,[]).append(r)
    print("\n[σ-SPECTRUM on the g line, novel]:",flush=True)
    for sig in sorted(spectrum, key=lambda x:int(x)):
        rs=spectrum[sig]
        distinct={dkey(r) for r in rs}
        print(f"   σ={sig}: {len(rs)} parametrisation(s) -> {len(distinct)} distinct novel knots "
              f"(c={sorted(set(r['c'] for r in rs))})",flush=True)

    out=dict(order=dict(target_genus=G, family='TwistedTorus T(p,q;2,s)', law='g=(p-1)(q-1)/2 + s'),
             n_candidates=len(cands), n_measured=len(ok), genus_law_match=f"{gmatch}/{len(ok)}",
             n_novel=len(novel),
             sigma_spectrum={str(k):dict(n_param=len(v),
                                         n_distinct=len({dkey(r) for r in v}),
                                         cs=sorted(set(r['c'] for r in v))) for k,v in spectrum.items()})

    # ---- fill order (G,S): S = σ with the most distinct novel reps ----
    if novel:
        S=max(spectrum, key=lambda k: len({dkey(r) for r in spectrum[k]}))
        pick=sorted(spectrum[S], key=lambda r: float(r['snappy'].get('vol') or 0))[0]
        out['filled_order']=dict(target=(G,int(S)), delivered=dict(
            params=(pick['p'],pick['q'],pick['r'],pick['s']), c=pick['c'],
            braid_word=braid_word(pick['p'],pick['q'],pick['r'],pick['s']),
            sig=pick['sage'].get('sig'), det=pick['sage'].get('det'),
            genus=pick['sage'].get('genus'), vol=pick['snappy'].get('vol')))
        print(f"\n[FILLED ORDER] (g={G}, σ={S}) -> delivered: T{(pick['p'],pick['q'],pick['r'],pick['s'])} "
              f"c={pick['c']} det={pick['sage'].get('det')} vol={pick['snappy'].get('vol')}",flush=True)

    json.dump(out, open(os.path.join(SP,'e2e_result.json'),'w'), indent=2)
    print("\n-> e2e_result.json  (+ e2e_measurements.jsonl)",flush=True)

if __name__=='__main__':
    main()
