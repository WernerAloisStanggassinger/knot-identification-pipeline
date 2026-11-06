# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
import sys, json
sys.argv = ['x']
exec(open('famed_unabhaengig.py').read().split("if __name__")[0])
import snappy
def mit_stoerung(tri, art):
    N = tri.size(); F = tri.countTriangles()
    X = [matrix(ZZ, N, F) for _ in range(4)]
    for i in range(N):
        t = tri.tetrahedron(i)
        for k in range(4): X[k][i, t.triangle(k).index()] = 1
    eps = [tri.tetrahedron(i).orientation() for i in range(N)]
    if art == 'eps': eps[0] = -eps[0]
    if art == 'X': X[1], X[2] = X[2], X[1]
    cA = block_matrix([[X[0]-X[1]+X[2]], [X[2]-X[3]]])
    if cA.det() == 0: return 'detA0'
    M = snappy.Manifold(tri); a, b = M.homological_longitude(); NZ = matrix(ZZ, M.gluing_equations())
    rows = lambda k: NZ[:N-1, k::3].stack(a*NZ[N, k::3] + b*NZ[N+1, k::3])
    A = rows(0)-rows(1); B = rows(2)-rows(1)
    if art == 'long': 
        rows = lambda k: NZ[:N-1, k::3].stack(NZ[N, k::3]); A = rows(0)-rows(1); B = rows(2)-rows(1)
    if B.det() == 0: return 'detB0'
    L = B.inverse()*A
    for s in (1, -1):
        E = diagonal_matrix(QQ, [s*e for e in eps]); cB = zero_matrix(QQ, N, N).stack(E)
        P = X[0].change_ring(QQ)*cA.inverse()*cB
        if L == P + P.transpose() + (E + identity_matrix(QQ, N))/2: return 'BESTANDEN'
    return '(4) verletzt'
kand = {k['did']: k for k in json.load(open('kandidaten.json'))}
for d in ['4_1', '6_1', 'pz_3_5_7', 'T_7_5_2_6', 'sd_p3_l5']:
    src = kand[d]['isosig'] if d in kand else d
    T = FAMEDTriangulation(snappy.Manifold(src), order_management='first_hit')
    for o in T.order_triangulation()[:5]:
        tri = o.triangulation
        if eigen_famed(tri)[0] is not True: continue
        print(d, 'undisturbed:', mit_stoerung(tri, None), '| flip eps:', mit_stoerung(tri, 'eps'), '| X1<->X2:', mit_stoerung(tri, 'X'), '| meridian instead of longitude:', mit_stoerung(tri, 'long'))
        break
