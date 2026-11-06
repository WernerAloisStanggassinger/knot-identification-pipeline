#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""
NAKANISHI INDEX — lower bound for u(K) from the ALEXANDER MODULE.

Why exactly this: Myf Evans names in 2504.01747v2 §5 the lower bound as the open spot and
refers to the classical family — among them **Nakanishi 2005** (her reference [29]) and
Ma & Qiu 2006 ([30], „modified Wirtinger presentation …"). The Nakanishi index is this direction.

    m(K) = minimal number of generators of the Alexander module H_1(X_oo) over Lambda = Z[t,1/t]
    u(K) >= m(K)                                              (Nakanishi)
    Presentation:  H_1(X_oo) = coker( V - t V^T ),  V = Seifert matrix (n x n, n = 2g)

## How the estimate from below is made — exact, not heuristic
For every maximal ideal m = (p, f) of Lambda holds
    m(K) >= dim_{k(m)} ( M (x) k(m) ) = n - rank_{k(m)} ( V - t V^T ) ,  k(m) = F_p[t]/(f) = F_{p^d}.
The maximum over all TESTED points is a valid lower bound; more points can
only raise it. ⚠️ `m_lower <= m(K) <= u(K)` — we never claim to know m(K) itself.

## How to compute over F_{p^d} without building field arithmetic
Multiplication by t on F_p[t]/(f) is the companion matrix C_f. If one replaces every entry
a_ij(t) = V_ij - t V_ji by the d x d matrix V_ij*I - V_ji*C_f, one gets
    B = kron(V, I_d) - kron(V^T, C_f)   over F_p,   and   corank_{F_{p^d}} = corank_{F_p}(B) / d .
So Gauss over a PRIME FIELD suffices — the extension field sits in the block structure.

## Test points
- Degree 1: all alpha in F_p^* for p in KLEIN; plus **alpha = -1 for every prime divisor of det**
  (that is exactly the covering bound) and all alpha for prime divisors p <= 31.
- Degree 2 and 3: all monic irreducible f over F_2, F_3 (degree 2 additionally over F_5),
  as long as n*d does not exceed the size limit.

## Counter-checks (must bite)
  N1  m_lower >= mu(Wendt)                 — structurally enforced (alpha = -1 is included)
  N2  Corank at alpha = +1 is ALWAYS 0    — det(V - V^T) = +-1 for every Seifert matrix
  N3  m_lower <= u_obere, where witnessed
  N4  ground truth: torus knots u = (p-1)(q-1)/2, computed
  --selbsttest replaces V^T by V ⟹ N2 MUST turn red.
"""
import site_config
import ast
import itertools
import json
import os
import sys
import time

import numpy as np
from sympy import factorint

import spherogram

WURZEL = site_config.ROOT
DB = WURZEL + "/freeze_manifest/pipeline_scripts/corpus_db.json"
HIER = os.path.dirname(os.path.abspath(__file__))
ZIEL = os.path.join(HIER, "nakanishi.json")

KLEIN = [2, 3, 5, 7, 11, 13]        # degree-1 points: all alpha
ALLE_ALPHA_BIS = 31                 # prime divisors of det up to here: all alpha
GRAD23 = {2: [2, 3], 3: [2, 3], 5: [2]}   # p -> which degrees as extension
MAX_BLOCK = 620                     # n*d cap for the extension points


def rang_mod_p(M, p):
    """Gaussian elimination over F_p, vectorized."""
    A = np.mod(M, p).astype(np.int64)
    n_z, n_s = A.shape
    rang, zeile = 0, 0
    for spalte in range(n_s):
        nz = np.nonzero(A[zeile:, spalte])[0]
        if nz.size == 0:
            continue
        piv = zeile + int(nz[0])
        if piv != zeile:
            A[[zeile, piv]] = A[[piv, zeile]]
        A[zeile] = (A[zeile] * pow(int(A[zeile, spalte]), -1, p)) % p
        sp = A[:, spalte].copy()
        sp[zeile] = 0
        tr = np.nonzero(sp)[0]
        if tr.size:
            A[tr] = (A[tr] - np.outer(sp[tr], A[zeile])) % p
        zeile += 1
        rang += 1
        if zeile == n_z:
            break
    return rang


def irreduzible(p, d):
    """All monic irreducible polynomials of degree d over F_p, as a coefficient list
    [a_0, ..., a_{d-1}] for t^d + a_{d-1} t^{d-1} + ... + a_0."""
    aus = []
    for koef in itertools.product(range(p), repeat=d):
        if koef[0] == 0:          # t | f  ⟹ reducible (and t is a unit)
            continue
        # Test for roots in F_p (suffices for d = 2, 3)
        hat_nst = False
        for x in range(p):
            w = pow(x, d, p)
            for i, a in enumerate(koef):
                w = (w + a * pow(x, i, p)) % p
            if w % p == 0:
                hat_nst = True
                break
        if not hat_nst:
            aus.append(list(koef))
    return aus


def begleitmatrix(koef, p):
    d = len(koef)
    C = np.zeros((d, d), dtype=np.int64)
    for i in range(1, d):
        C[i, i - 1] = 1
    for i, a in enumerate(koef):
        C[i, d - 1] = (-a) % p
    return C


def punkte(det):
    """List of the test points: (p, coefficients_of_f_or_None, alpha_or_None, degree)."""
    aus = []
    for p in KLEIN:
        for a in range(1, p):
            aus.append((p, None, a, 1))
    for p in sorted(factorint(abs(int(det))) if det not in (0, None) else {}):
        p = int(p)
        if p in KLEIN:
            continue
        if p <= ALLE_ALPHA_BIS:
            for a in range(1, p):
                aus.append((p, None, a, 1))
        else:
            aus.append((p, None, (p - 1), 1))       # alpha = -1: the covering bound
    for p, grade in GRAD23.items():
        for d in grade:
            for f in irreduzible(p, d):
                aus.append((p, f, None, d))
    return aus


def nakanishi_lower(V, det, sabotage=False):
    V = np.array(V, dtype=np.int64)
    n = V.shape[0]
    if n == 0:
        return 0, {}, 0, 0
    VT = V.copy() if sabotage else V.T.copy()
    beste, treffer, kor1, gezaehlt = 0, {}, 0, 0
    for p, f, alpha, d in punkte(det):
        if d == 1:
            B = (V - alpha * VT) % p
            gross = n
        else:
            if n * d > MAX_BLOCK:
                continue
            C = begleitmatrix(f, p)
            I = np.eye(d, dtype=np.int64)
            B = (np.kron(V, I) - np.kron(VT, C)) % p
            gross = n * d
        kor = gross - rang_mod_p(B, p)
        gezaehlt += 1
        if d > 1:
            assert kor % d == 0, "Corank over F_{p^d} must be divisible by d"
            kor //= d
        if d == 1 and alpha % p == 1 % p:
            kor1 = max(kor1, kor)
        if kor > 0:
            treffer[f"p{p}_" + (f"a{alpha}" if d == 1 else f"f{f}")] = kor
        beste = max(beste, kor)
    return beste, treffer, kor1, gezaehlt


def eine(r, sabotage=False):
    d = json.load(open(r["pflicht_path"]))
    pd = d.get("pd_list_str")
    if not pd:
        raise ValueError("no pd_list_str")
    L = spherogram.Link([tuple(x) for x in ast.literal_eval(pd)])
    V = [[int(x) for x in z] for z in L.seifert_matrix()]
    m_low, tr, kor1, gez = nakanishi_lower(V, r.get("det"), sabotage=sabotage)
    return dict(n=len(V), m_lower=m_low, treffer=tr,
                korang_bei_alpha_1=kor1, punkte_getestet=gez)


def main():
    selbsttest = "--selbsttest" in sys.argv
    neu = "--neu" in sys.argv
    c_max = next((int(a) for a in sys.argv[1:] if a.isdigit()), None)

    db = json.load(open(DB))
    wendt = json.load(open(os.path.join(HIER, "wendt_sweep.json")))["zeilen"]
    ernte = {z["stem"]: z for z in json.load(
        open(os.path.join(HIER, "ernte_u_schranken.json")))["zeilen"]}

    stand = {}
    if os.path.exists(ZIEL) and not neu:
        alt = json.load(open(ZIEL))
        stand = {k: v for k, v in alt.get("zeilen", {}).items() if v.get("m_lower") is not None}
        print(f"↻ restart: {len(stand)} decided")

    kand = sorted([r for r in db["knots"] if r["stem"] not in stand
                   and (c_max is None or (r.get("c") or 999) <= c_max)],
                  key=lambda r: (r.get("c") or 0))
    print(f"open: {len(kand)}", flush=True)

    t0 = time.time()
    for i, r in enumerate(kand, 1):
        stem = r["stem"]
        e = dict(stem=stem, c=r.get("c"), det=r.get("det"), category=r.get("category"),
                 knotinfo_id=r.get("knotinfo_id"),
                 ist_vertreter=bool(r.get("is_distinct_representative")),
                 distinct_knot_id=r.get("distinct_knot_id"),
                 mu_wendt=wendt.get(stem, {}).get("mu"),
                 u_obere=ernte.get(stem, {}).get("u_obere"),
                 u_lower_alt=ernte.get(stem, {}).get("u_lower"),
                 m_lower=None, fehler=None)
        try:
            e.update(eine(r))
        except Exception as ex:
            e["fehler"] = f"{type(ex).__name__}: {str(ex)[:110]}"
        stand[stem] = e
        if i % 20 == 0 or i == len(kand):
            json.dump({"gebaut_am": "2026-08-11", "n": len(stand), "zeilen": stand},
                      open(ZIEL + ".tmp", "w"), indent=1)
            os.replace(ZIEL + ".tmp", ZIEL)
            print(f"[{i}/{len(kand)}] {i/(time.time()-t0):.2f}/s  {stem[:30]} "
                  f"m>={e['m_lower']} mu={e['mu_wendt']}", flush=True)

    json.dump({"gebaut_am": "2026-08-11", "n": len(stand), "zeilen": stand},
              open(ZIEL + ".tmp", "w"), indent=1)
    os.replace(ZIEL + ".tmp", ZIEL)

    ok = [v for v in stand.values() if v.get("m_lower") is not None]
    n1 = [v for v in ok if v["mu_wendt"] is not None]
    n1_bad = [v for v in n1 if v["m_lower"] < v["mu_wendt"]]
    n2_bad = [v for v in ok if v.get("korang_bei_alpha_1", 0) != 0]
    n3 = [v for v in ok if v["u_obere"] is not None]
    n3_bad = [v for v in n3 if v["m_lower"] > v["u_obere"]]
    torus = {"K3_1": (2, 3), "K5_1": (2, 5), "K7_1": (2, 7), "K9_1": (2, 9),
             "K8_19": (3, 4), "K10_124": (3, 5)}
    grund = {k: (p - 1) * (q - 1) // 2 for k, (p, q) in torus.items()}
    grund["K4_1"] = 1
    for k in list(grund):
        grund[k + "m"] = grund[k]
    n4 = [v for v in ok if (v["knotinfo_id"] or "") in grund]
    n4_bad = [v for v in n4 if v["m_lower"] > grund[v["knotinfo_id"]]]

    print()
    print(f"N1  m_lower >= mu(Wendt)      : {len(n1)-len(n1_bad)}/{len(n1)} ok, {len(n1_bad)} RED")
    for v in n1_bad[:8]:
        print(f"    RED {v['stem']} (det {v['det']}): m {v['m_lower']} < mu {v['mu_wendt']}")
    print(f"N2  corank at alpha = +1 == 0: {len(ok)-len(n2_bad)}/{len(ok)} ok, {len(n2_bad)} RED")
    print(f"N3  m_lower <= u_obere        : {len(n3)-len(n3_bad)}/{len(n3)} ok, {len(n3_bad)} RED")
    for v in n3_bad[:8]:
        print(f"    RED {v['stem']}: m {v['m_lower']} > u<= {v['u_obere']}")
    print(f"N4  ground truth              : {len(n4)-len(n4_bad)}/{len(n4)} ok, {len(n4_bad)} RED")
    for v in n4_bad[:8]:
        print(f"    RED {v['stem']} ({v['knotinfo_id']}): m {v['m_lower']} > u {grund[v['knotinfo_id']]}")

    besser = [v for v in ok if v["mu_wendt"] is not None and v["m_lower"] > v["mu_wendt"]]
    hebt = [v for v in ok if v["u_lower_alt"] is not None and v["m_lower"] > v["u_lower_alt"]]
    print()
    print(f"Nakanishi BEATS the covering bound : {len(besser)} records "
          f"({len({v['distinct_knot_id'] for v in besser})} distinct)")
    for v in sorted(besser, key=lambda v: -(v["m_lower"] - v["mu_wendt"]))[:15]:
        print(f"    +{v['m_lower']-v['mu_wendt']}  m >= {v['m_lower']} instead of mu {v['mu_wendt']}  "
              f"c={v['c']:3} det={v['det']}  {v['stem'][:34]}  {list(v['treffer'])[:3]}")
    print(f"Nakanishi raises the bound over sigma/s/tau/g4 : {len(hebt)} records")
    print(f"Errors: {len([v for v in stand.values() if v.get('fehler')])}")

    if selbsttest:
        print("\n--- SELF-TEST: V^T replaced by V ---")
        proben = [r for r in db["knots"] if (r.get("c") or 999) <= 12][:25]
        schlecht = 0
        for r in proben:
            try:
                if eine(r, sabotage=True)["korang_bei_alpha_1"] != 0:
                    schlecht += 1
            except Exception:
                pass
        print(f"    N2 red in {schlecht}/{len(proben)} cases ⟹ "
              f"{'GREEN (bites)' if schlecht else 'RED (decoration!)'}")
        if not schlecht:
            sys.exit(2)


if __name__ == "__main__":
    main()
