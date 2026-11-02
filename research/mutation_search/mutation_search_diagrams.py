# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Conway-mutation search between two Burton-census knots on MANY NON-MINIMAL diagrams of the source knot.
# A mutation may need a diagram larger than the census diagram -> search many diagrams from backtrack() and from exterior_to_link().
# Mutation construction (faces/components/reconnection/isometry test) copied from mutation_search_census.py.
# Only change: 4-boundary subsets are found by a numpy bitmask sweep instead of itertools (same set, faster).
# A FOUND record proves: the target is a Conway mutant of the source or of its mirror image (is_isometric_to is mirror-blind).
# Usage: python mutation_search_diagrams.py SRC DST OUT.jsonl SEED MAXCROSS SECONDS [MINCROSS]
import sys, json, time, bz2, csv, random, signal
import numpy as np
import snappy
from spherogram import Link
src_name, dst_name, out_path = sys.argv[1:4]
seed, maxc, budget = int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
minc = int(sys.argv[7]) if len(sys.argv) > 7 else 0
signal.alarm(budget + 900)
random.seed(seed)
RAW = __import__('os').environ.get('BURTON_CENSUS_DIR', 'burton_census/raw')   # folder with Burton's <c>{a,n}-hyp.csv.bz2
dt = lambda s: [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]
def census_link(n):
    c = int(n.split('nh_')[0])
    for r in csv.DictReader(bz2.open(f'{RAW}/{c}n-hyp.csv.bz2', 'rt')):
        if r['name'] == n: return Link('DT:' + str([tuple(dt(r['dt_code']))]))
    raise KeyError(n)
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
def four_boundary_subsets(pd):
    n = len(pd)
    where = {}
    for x, cr in enumerate(pd):
        for e in cr: where.setdefault(e, []).append(x)
    masks = np.arange(1 << n, dtype=np.int64)
    cnt = np.zeros(1 << n, dtype=np.int8)
    for e, xs in where.items():
        a, b = xs[0], xs[1]
        if a == b: continue
        cnt += (((masks >> a) ^ (masks >> b)) & 1).astype(np.int8)
    pop = np.zeros(1 << n, dtype=np.int8)
    for i in range(n): pop += ((masks >> i) & 1).astype(np.int8)
    sel = masks[(cnt == 4) & (pop >= 2) & (pop <= n - 2)]
    return [sorted(i for i in range(n) if (m >> i) & 1) for m in sel.tolist()], where
def search_pd(pd, dst_ext, src_ext):
    n = len(pd)
    if faces(pd) != n + 2: return dict(skip='nonplanar_input')
    subsets, where = four_boundary_subsets(pd)
    pairings = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
    n_planar = 0
    for T in subsets:
        Ts = set(T)
        bnd = sorted(e for e, xs in where.items() if len(set(xs) & Ts) == 1)
        if len(bnd) != 4: continue
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
            try:
                M = snappy.Link(new).exterior()
                if M.is_isometric_to(dst_ext) and not M.is_isometric_to(src_ext):
                    return dict(found=True, tangle=T, boundary=bnd, pairing=pr, source_pd=pd, mutated_pd=new,
                                subsets=len(subsets), planar=n_planar)
            except Exception:
                pass
    return dict(found=False, subsets=len(subsets), planar=n_planar)
L0 = census_link(src_name); D0 = census_link(dst_name)
src_ext, dst_ext = L0.exterior(), D0.exterior()
out = open(out_path, 'a')
def emit(rec):
    out.write(json.dumps(rec) + '\n'); out.flush()
emit(dict(start=time.strftime('%F %T'), src=src_name, dst=dst_name, seed=seed, maxc=maxc, minc=minc, budget=budget,
          isometric_src_dst=bool(src_ext.is_isometric_to(dst_ext))))
seen = set(); t0 = time.time(); k = 0
while time.time() - t0 < budget:
    k += 1
    try:
        if k % 3 == 0:
            M = src_ext.copy(); M.randomize(); L = M.exterior_to_link(); L.simplify('global'); how = 'exterior_to_link'
        else:
            L = L0.copy(); L.backtrack(steps=random.randint(1, 12)); how = 'backtrack'
            if random.random() < 0.5: L.simplify('level')
        if len(L.link_components) != 1: continue
        pd = [list(x) for x in L.PD_code()]
    except Exception as ex:
        emit(dict(k=k, error=repr(ex)[:200])); continue
    n = len(pd)
    if n > maxc or n < minc: continue
    key = tuple(sorted(tuple(sorted(c)) for c in pd))
    if key in seen: continue
    seen.add(key)
    t1 = time.time(); r = search_pd(pd, dst_ext, src_ext)
    r.update(k=k, how=how, crossings=n, seconds=round(time.time() - t1, 2))
    if not r.get('found'):
        r.pop('source_pd', None)
    emit(r)
    if r.get('found'):
        emit(dict(RESULT='FOUND', src=src_name, dst=dst_name, at=time.strftime('%F %T'))); sys.exit(0)
emit(dict(RESULT='NOT_FOUND_IN_BUDGET', src=src_name, dst=dst_name, diagrams=len(seen), at=time.strftime('%F %T')))
