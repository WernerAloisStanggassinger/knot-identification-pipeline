# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Shared parsing/canonicalization for corpus-vs-Maguire adjudication."""
import re, tarfile, lzma, io, cmath

# ---------- DT ----------
def decode_dt(s):
    out=[]
    for ch in s:
        if not ch.isalpha(): continue
        v=2*(ord(ch.lower())-ord('a')+1)
        if ch.isupper(): v=-v
        out.append(v)
    return out

# ---------- polynomial parsing ----------
def _parse_terms(poly, varnames):
    """Parse a polynomial string into {exptuple: coeff}. varnames e.g. ['x'] or ['x','y'].
    Handles coeffs like -, +, integers; exps like x^-4, x^2, x (=^1), bare const."""
    poly=poly.replace(' ','').replace('*','').strip()
    if not poly: return {}
    # protect exponent signs so we can split on term-level +/-
    poly=poly.replace('^-','^~')
    # ensure a leading + for uniform splitting
    if poly[0] not in '+-': poly='+'+poly
    terms=re.findall(r'[+-][^+-]*', poly)
    res={}
    for t in terms:
        sign=-1 if t[0]=='-' else 1
        body=t[1:].replace('^~','^-')
        if not body: continue
        # coefficient = leading integer (before first variable letter), else 1
        mcoef=re.match(r'^(\d+)', body)
        coeff=int(mcoef.group(1)) if mcoef else 1
        coeff*=sign
        rest=body[mcoef.end():] if mcoef else body
        exps={v:0 for v in varnames}
        for v in varnames:
            for mm in re.finditer(re.escape(v)+r'(\^-?\d+)?', rest):
                e=mm.group(1)
                exps[v]= int(e[1:]) if e else 1
        key=tuple(exps[v] for v in varnames)
        res[key]=res.get(key,0)+coeff
    return {k:v for k,v in res.items() if v!=0}

def parse_jones(poly):
    d=_parse_terms(poly,['x'])
    return {k[0]:v for k,v in d.items()}

def parse_homfly(poly):
    return _parse_terms(poly,['x','y'])  # x=a, y=z

# ---------- canonicalization (mirror-invariant) ----------
def canon_jones(d):
    if not d: return ()
    a=tuple(sorted(d.items()))
    b=tuple(sorted((-e,c) for e,c in d.items()))
    return min(a,b)

def canon_homfly(d):
    if not d: return ()
    # mirror: x-exp -> -x-exp (a->a^-1), y unchanged
    a=tuple(sorted(d.items()))
    b=tuple(sorted(((-xe,ye),c) for (xe,ye),c in d.items()))
    return min(a,b)

# ---------- determinant from Jones = |V(-1)| ; t=x^2 so x=i ----------
def det_from_jones(d):
    val=sum(c*(1j**e) for e,c in d.items())
    return round(abs(val))

# ---------- stream a maguire tar.xz of csv rows ----------
def stream_maguire(path):
    """Yield (dt_str, poly_str) for each row. tar contains one csv."""
    with tarfile.open(path,'r:xz') as tf:
        for m in tf.getmembers():
            if not m.name.endswith('.csv'): continue
            f=tf.extractfile(m)
            for raw in f:
                line=raw.decode('ascii','replace').strip()
                if not line: continue
                i=line.find(',')
                if i<0: continue
                yield line[:i], line[i+1:]
