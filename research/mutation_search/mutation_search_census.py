# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Exhaustive Conway-mutation search on the census diagrams of two Burton-census knots given by NAME (non-alternating hyperbolic
# 'nh' names, any crossing numbers), mirror-blind via is_isometric_to. Functions faces/components/search copied verbatim from the
# authors' mutation_check_15.py. Calibration: run on 15nh_048139 15nh_048140, it must find their known mutation (README).
# Usage: sage python mutation_search_census.py <nameA> <nameB>   (searches A->B and B->A)
import sys, itertools, json, time, bz2, csv, signal
import snappy
from spherogram import Link
signal.alarm(3000)
RAW = __import__('os').environ.get('BURTON_CENSUS_DIR', 'burton_census/raw')   # folder with Burton's <c>{a,n}-hyp.csv.bz2
names = sys.argv[1:3]
dt = lambda s: [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]
dts = {}
for n in names:
    c = int(n.split('nh_')[0])
    for r in csv.DictReader(bz2.open(f'{RAW}/{c}n-hyp.csv.bz2', 'rt')):
        if r['name'] == n: dts[n] = dt(r['dt_code']); break
assert len(dts) == 2, dts.keys()
def faces(p):
    slots = {}
    for x, cr in enumerate(p):
        for i, e in enumerate(cr): slots.setdefault(e, []).append((x, i))
    if any(len(v) != 2 for v in slots.values()): return -1
    other = {}
    for e, (s, t) in slots.items(): other[s] = t; other[t] = s
    seen, n = set(), 0
    for s in other:
        if s in seen: continue
        n += 1; cur = s
        while cur not in seen:
            seen.add(cur); x, j = other[cur]; cur = (x, (j + 1) % 4)
    return n
def components(p):
    adj = {}
    for cr in p:
        for a, b in ((cr[0], cr[2]), (cr[1], cr[3])):
            adj.setdefault(a, []).append(b); adj.setdefault(b, []).append(a)
    seen, n = set(), 0
    for e in adj:
        if e in seen: continue
        n += 1; st = [e]
        while st:
            u = st.pop()
            if u in seen: continue
            seen.add(u); st += adj[u]
    return n
def search(src, dst_exts, src_ext):
    L = Link('DT:' + str([tuple(dts[src])])); pd = [list(x) for x in L.PD_code()]; n = len(pd)
    assert faces(pd) == n + 2
    where = {}
    for x, cr in enumerate(pd):
        for e in cr: where.setdefault(e, set()).add(x)
    pairings = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
    n_sub = n_planar = 0; t0 = time.time()
    for size in range(2, n - 1):
        for T in itertools.combinations(range(n), size):
            Ts = set(T)
            bnd = sorted(e for e, xs in where.items() if len(xs & Ts) == 1)
            if len(bnd) != 4: continue
            n_sub += 1
            for pr in pairings:
                sigma = {i: i for i in range(4)}
                for a, b in pr: sigma[a], sigma[b] = b, a
                new = [list(cr) for cr in pd]; NEW = 1000
                for i, e in enumerate(bnd):
                    for x, cr in enumerate(pd):
                        for j, f in enumerate(cr):
                            if f == e and x in Ts: new[x][j] = NEW + i
                    for x, cr in enumerate(pd):
                        for j, f in enumerate(cr):
                            if f == bnd[sigma[i]] and x not in Ts: new[x][j] = NEW + i
                if faces(new) != n + 2 or components(new) != 1: continue
                n_planar += 1
                try: M = snappy.Link(new).exterior()
                except Exception: continue
                for lab, E in dst_exts:
                    try:
                        if M.is_isometric_to(E) and not M.is_isometric_to(src_ext):
                            return dict(found=True, tangle=sorted(Ts), boundary=bnd, pairing=pr, target=lab, subsets=n_sub, planar=n_planar, seconds=round(time.time()-t0))
                    except Exception: pass
    return dict(found=False, subsets=n_sub, planar=n_planar, seconds=round(time.time()-t0))

ia, ib = names
EA = Link('DT:' + str([tuple(dts[ia])])).exterior(); EB = Link('DT:' + str([tuple(dts[ib])])).exterior()
print('vol', float(EA.volume()), float(EB.volume()), '| isometric (expect False):', EA.is_isometric_to(EB), flush=True)
r1 = search(ia, [(ib, EB)], EA); print('A->B', json.dumps(r1), flush=True)
r2 = search(ib, [(ia, EA)], EB); print('B->A', json.dumps(r2), flush=True)
print('RESULT', json.dumps(dict(A=ia, B=ib, found=bool(r1['found'] or r2['found']))), flush=True)
