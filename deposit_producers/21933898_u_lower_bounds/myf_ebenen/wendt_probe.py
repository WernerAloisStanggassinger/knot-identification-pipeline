#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""
PROBE (small, profile before scaling): the Wendt bound

    u(K) >= mu( H_1(Sigma_2(K)) )        Wendt 1937

mu = minimal number of generators of the first homology of the double branched cover.
H_1(Sigma_2) is presented by  M = V + V^T  (V = Seifert matrix), |det M| = det(K).

mu is NOT computed via the Smith normal form (expensive), but exactly via
    mu = max_{p | det}  dim_{F_p} coker(M mod p) = max_p ( n - rank_{F_p}(M) )
— for a finite abelian group A, mu(A) = max_p dim(A/pA). Rank over F_p
is Gaussian elimination in a prime field: fast and exact.

WHY THIS IS INTERESTING HERE: for amphichiral knots sigma = s = tau = 0,
all previous lower bounds fall to 0. Exactly there Wendt lives on.
12 of Werner's 19 sculptures are amphichiral.

COUNTER-CHECKS (must bite):
  G1  |det(V+V^T)| == det(K) from the record   (catches PD/record discrepancy)
  G2  trefoil: H_1 = Z/3, mu = 1
  G3  knots-of-love = 3_1#3_1#3_1#m3_1#m3_1#m3_1: H_1 = (Z/3)^6, mu = 6
"""
import site_config
import ast
import json
import sys
import time

from sympy import factorint

import spherogram

WURZEL = site_config.ROOT
DB = WURZEL + "/freeze_manifest/pipeline_scripts/corpus_db.json"


def rang_mod_p(M, p):
    """Gaussian elimination over F_p. M = list of lists (int). Returns the rank."""
    A = [[x % p for x in zeile] for zeile in M]
    n_z, n_s = len(A), len(A[0]) if A else 0
    rang, zeile = 0, 0
    for spalte in range(n_s):
        piv = None
        for i in range(zeile, n_z):
            if A[i][spalte] % p:
                piv = i
                break
        if piv is None:
            continue
        A[zeile], A[piv] = A[piv], A[zeile]
        inv = pow(A[zeile][spalte], -1, p)
        A[zeile] = [(x * inv) % p for x in A[zeile]]
        for i in range(n_z):
            if i != zeile and A[i][spalte]:
                f = A[i][spalte]
                A[i] = [(A[i][j] - f * A[zeile][j]) % p for j in range(n_s)]
        zeile += 1
        rang += 1
        if zeile == n_z:
            break
    return rang


def mu_von(M, det):
    """mu(H_1) = max_p (n - rank_{F_p}(M)) over the prime divisors of det."""
    n = len(M)
    if det in (0, 1, -1):
        return 0, {}
    pro_p = {}
    for p in factorint(abs(int(det))):
        pro_p[int(p)] = n - rang_mod_p(M, int(p))
    return (max(pro_p.values()) if pro_p else 0), pro_p


def det_ganz(M):
    """Exact determinant by fraction-free Gaussian elimination (Bareiss)."""
    A = [row[:] for row in M]
    n = len(A)
    if n == 0:
        return 1
    vorz, vorheriges = 1, 1
    for k in range(n - 1):
        if A[k][k] == 0:
            for i in range(k + 1, n):
                if A[i][k] != 0:
                    A[k], A[i] = A[i], A[k]
                    vorz = -vorz
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // vorheriges
        vorheriges = A[k][k]
    return vorz * A[n - 1][n - 1]


def eine(stem, pd_str, det_record, zeit_grenze=None):
    t0 = time.time()
    pd = ast.literal_eval(pd_str)
    L = spherogram.Link([tuple(x) for x in pd])
    V = L.seifert_matrix()          # spherogram returns a plain list of lists
    V = [[int(x) for x in zeile] for zeile in V]
    n = len(V)
    M = [[V[i][j] + V[j][i] for j in range(n)] for i in range(n)]
    dM = det_ganz(M)
    mu, pro_p = mu_von(M, dM)
    return dict(stem=stem, n=n, det_record=det_record, det_VVT=abs(dM),
                det_stimmt=(abs(dM) == abs(int(det_record))),
                mu=mu, pro_primzahl=pro_p, sekunden=round(time.time() - t0, 2))


def main():
    db = json.load(open(DB))
    wunsch = sys.argv[1:] or ["celtic-trefoil", "knots-of-love", "trainstation", "double-3",
                              "nitro", "ping-pong", "mirror", "irrational", "spiral",
                              "celtic-for-playing", "sundowner", "cosmic-love", "double-6"]
    nach_stem = {x["stem"]: x for x in db["knots"]}
    for stem in wunsch:
        r = nach_stem.get(stem)
        if r is None:
            print(f"{stem:22} — not in the corpus")
            continue
        d = json.load(open(r["pflicht_path"]))
        pd = d.get("pd_list_str")
        if not pd:
            print(f"{stem:22} — no pd_list_str")
            continue
        try:
            e = eine(stem, pd, r["det"])
        except Exception as ex:
            print(f"{stem:22} c={r['c']:<4} ERROR {type(ex).__name__}: {str(ex)[:70]}")
            continue
        marke = "ok " if e["det_stimmt"] else "🔴 "
        print(f"{stem:22} c={r['c']:<4} n={e['n']:<4} det {marke}{e['det_VVT']} "
              f"(Record {e['det_record']})  mu={e['mu']}  u>={e['mu']}  "
              f"p:{e['pro_primzahl']}  {e['sekunden']}s")


if __name__ == "__main__":
    main()
