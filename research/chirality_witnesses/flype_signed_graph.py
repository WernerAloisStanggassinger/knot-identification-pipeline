# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Flype test for chirality of prime alternating knots via the SIGNED interlace graph (01.11.2026).
# G±(D): vertices = crossings of a reduced alternating diagram D, edges = interlaced chords of the Gauss word, vertex label = crossing sign.
# (i) A flype keeps the labelled interlace graph (lemma in our A343358 deposit, DOI 10.5281/zenodo.23855791, file DEFINITION_CHECK_20261030.md §5).
# (ii) A flype is an ambient isotopy that carries each crossing to a crossing; the sign of a crossing is a local isotopy invariant of the
#      oriented knot and does not change when the whole knot is reversed ⟹ the labels ± are carried along too. [derivation — tested below]
# (iii) Menasco–Thistlethwaite 1993: reduced alternating diagrams of the same knot are flype-related. Mirror image flips every sign.
# ⟹ K amphichiral ⟹ G±(D) ≅ G∓(D) (sign-swapping isomorphism). Contrapositive: no sign-swapping isomorphism ⟹ K chiral.
# Two implementations: Sage canonical_label with ordered partition [P,N] vs [N,P]; networkx VF2 with node_match on the sign.
# Counter-cases (must pass): amphichiral 4_1 6_3 8_3 8_9 8_12 8_17 8_18 sign-swap compatible; chiral 3_1 5_2 not.
import sys, json, bz2, csv, itertools
import sage.all
from sage.all import Graph
import networkx as nx, snappy
dt = lambda s: [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]
def signed_graph(L):
    S = L.sage_link(); gc, signs = S.oriented_gauss_code()
    assert len(gc) == 1; w = [abs(x) for x in gc[0]]; n = len(signs)
    pos = {}
    for k, x in enumerate(w): pos.setdefault(x, []).append(k)
    assert all(len(v) == 2 for v in pos.values()) and len(pos) == n
    E = [(a, b) for a, b in itertools.combinations(sorted(pos), 2)
         if (pos[a][0] < pos[b][0] < pos[a][1]) != (pos[a][0] < pos[b][1] < pos[a][1])]
    return sorted(pos), E, {x: signs[x - 1] for x in pos}, S.is_alternating()
def test(L):
    V, E, sg, alt = signed_graph(L)
    P = [v for v in V if sg[v] > 0]; N = [v for v in V if sg[v] < 0]
    G = Graph([V, E], format='vertices_and_edges')
    if len(P) != len(N): sage_ok = False
    else: sage_ok = G.canonical_label(partition=[P, N]).edges(sort=True) == G.canonical_label(partition=[N, P]).edges(sort=True)
    H = nx.Graph(); H.add_nodes_from((v, {'s': sg[v]}) for v in V); H.add_edges_from(E)
    Hm = nx.Graph(); Hm.add_nodes_from((v, {'s': -sg[v]}) for v in V); Hm.add_edges_from(E)
    nx_ok = nx.is_isomorphic(H, Hm, node_match=lambda a, b: a['s'] == b['s'])
    assert sage_ok == nx_ok, ('implementations disagree', sage_ok, nx_ok)
    return dict(alternating_diagram=bool(alt), n=len(V), writhe=sum(sg.values()), sign_swap_compatible=bool(nx_ok))
def from_dt(code):
    L = snappy.Link('DT:' + str([tuple(dt(code))])); return L

if __name__ == '__main__':
    import os as _os
    sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
    from census_codes import main
    from census_codes import is_reduced
    def witness(L):                    # verdict wrapper (driver): the argument needs a REDUCED alternating diagram of a prime knot
        r = test(L); r['reduced_diagram'] = is_reduced(L.PD_code())
        usable = r['alternating_diagram'] and r['reduced_diagram']
        r['verdict'] = ('CHIRAL (no sign-swapping isomorphism of the signed interlace graph)' if usable and not r['sign_swap_compatible'] else
                        ('sign-swap compatible (no information)' if usable else 'not a reduced alternating diagram (no information)'))
        return r
    counter = [(nm, (lambda nm=nm: snappy.Link(nm)), (lambda r, e=e: r['alternating_diagram'] and r['reduced_diagram'] and r['sign_swap_compatible'] is e)) for nm, e in
               [(k, True) for k in ('4_1', '6_3', '8_3', '8_9', '8_12', '8_17', '8_18')] + [('3_1', False), ('5_2', False)]]
    main(witness, counter, from_dt)
