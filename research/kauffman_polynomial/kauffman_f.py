#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger. Authored in collaboration with Claude Code (Anthropic).
"""kauffman_f.py — the two-variable Kauffman polynomial F(a, z) of a knot or link from a PD code.

WHY (04.11.2026): we needed F for knots with 12 and 13 crossings computed by our own code, not only read from a table;
our earlier attempts did not get there (one stopped at braid index 4).

ALGORITHM (Kauffman 1990, the descending-diagram recursion; the same scheme KnotTheory`'s Kauffman[] uses):
  L is the regular-isotopy invariant with  L(O) = 1,  L(positive kink) = a·L,  L(X+) + L(X-) = z·(L(X0) + L(Xinf)).
  Choose a base point and a direction on every component and an order of the components. A crossing is GOOD if the
  traversal first meets it on the over strand, BAD otherwise. A diagram with no bad crossing is descending, hence a
  diagram of an unlink with writhe w and c components:  L = a^w · mu^(c-1),  mu = (a + 1/a)/z - 1.
  Otherwise take the first bad crossing in traversal order and use  L(D) = -L(D_switched) + z·(L(D_0) + L(D_inf)):
  D_switched has one bad crossing fewer (the earlier crossings are untouched), the two smoothings have one crossing fewer.
  Finally  F(a, z) = a^(-w(D)) · L(D).
  Cost ~3^n for n crossings (n = 12: ~5·10^5 leaves), memoised on the canonical PD.

CONVENTIONS: PD as KnotInfo/KnotTheory: X[i,j,k,l], i = incoming under strand, labels counterclockwise, k = outgoing
  under strand. The orientation of the over strand (l -> j or j -> l) is read off from the traversal, not from label
  arithmetic. The sign of a crossing follows the right-hand rule. Calibrated on KnotInfo's F for all prime knots with
  <= 10 crossings (`--calibrate`), which fixes the mirror convention a <-> 1/a once and for all.

USAGE:
  python3 kauffman_f.py --pd "[[1,5,2,4],[3,1,4,6],[5,3,6,2]]"      # trefoil, prints F
  python3 kauffman_f.py --selftest                                    # built-in counter-cases (no network, no KnotInfo)
  <python with database_knotinfo> kauffman_f.py --calibrate 10        # all KnotInfo knots with c <= 10, needs database_knotinfo
  <python with database_knotinfo> kauffman_f.py --knotinfo 3_1 8_19   # compute from KnotInfo PD and compare with KnotInfo F
Polynomials are dicts {(i, j): coeff} meaning coeff · a^i · z^j (pure Python, exact integers).
MEMORY: `_L` is memoised with an unbounded lru_cache. In a long-lived worker call `_L.cache_clear()` after every knot —
without it twelve sweep workers grew to 2.4 GB each on 04.11.2026 (footgun 17); `kauffman_F()` does NOT clear it itself
because repeated calls on the same diagram family benefit from the memo.
"""
import sys, ast, time, functools

# ----------------------------------------------------------------------------- Laurent polynomials in a, z
def padd(p, q):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, 0) + v
        if r[k] == 0:
            del r[k]
    return r

def pmul(p, q):
    r = {}
    for (i1, j1), v1 in p.items():
        for (i2, j2), v2 in q.items():
            k = (i1 + i2, j1 + j2)
            r[k] = r.get(k, 0) + v1 * v2
    return {k: v for k, v in r.items() if v}

def pscale(p, c):
    return {k: v * c for k, v in p.items()} if c else {}

def mono(i, j, c=1):
    return {(i, j): c}

ONE = mono(0, 0)
# mu = (a + a^-1)/z - 1
MU = {(1, -1): 1, (-1, -1): 1, (0, 0): -1}

def ppow(p, n):
    r = ONE
    for _ in range(n):
        r = pmul(r, p)
    return r

def pstr(p):
    """Readable string, sorted by z-degree then a-degree."""
    if not p:
        return "0"
    out = []
    for (i, j) in sorted(p, key=lambda k: (k[1], k[0])):
        c = p[(i, j)]
        term = ("" if abs(c) == 1 else str(abs(c)) + "*") + (f"a^({i})" if i else "") + ("*" if i and j else "") + (f"z^({j})" if j else "")
        if not term:
            term = str(abs(c))
        out.append(("-" if c < 0 else "+") + term)
    s = " ".join(out)
    return s[1:] if s.startswith("+") else s

def a_span(p):
    ex = [i for (i, j) in p]
    return (min(ex), max(ex), max(ex) - min(ex)) if ex else (0, 0, 0)

# ----------------------------------------------------------------------------- diagrams
# ORIENTED PD throughout: crossing = (a, b, c, d) counterclockwise, the under strand ENTERS at a (slot 0) and EXITS at c
# (slot 2); the over strand enters at b or d, which the traversal determines. Every arc label exits exactly one slot and
# enters exactly one slot. Base point of a component = its smallest arc label, travelling in its direction; components
# in order of their smallest label. Switching a crossing keeps labels and directions, so the traversal is identical
# along the switch chain and the switched crossing turns from bad to good (Kauffman's descending recursion).

def build_adjacency(crossings):
    adj = {}
    for ci, x in enumerate(crossings):
        for s in range(4):
            adj.setdefault(x[s], []).append((ci, s))
    return adj

def orient_slots(crossings):
    """Return enter[(ci, slot)] -> bool for every slot, following the under-strand directions."""
    adj = build_adjacency(crossings)
    enter = {}
    for ci, x in enumerate(crossings):
        enter[(ci, 0)] = True; enter[(ci, 2)] = False
    # propagate along labels: a label exiting somewhere enters at its other occurrence
    changed = True
    while changed:
        changed = False
        for ci, x in enumerate(crossings):
            for s in range(4):
                if (ci, s) in enter:
                    lab = x[s]
                    others = [e for e in adj[lab] if e != (ci, s)]
                    if not others:
                        continue
                    cj, t = others[0]
                    val = not enter[(ci, s)]
                    if (cj, t) not in enter:
                        enter[(cj, t)] = val; enter[(cj, (t + 2) % 4)] = not val; changed = True
    # components that never pass under: orientation free, fix it
    for ci, x in enumerate(crossings):
        for s in (1, 3):
            if (ci, s) not in enter:
                enter[(ci, s)] = (s == 1); enter[(ci, (s + 2) % 4)] = (s != 1)
                # propagate
                changed = True
                while changed:
                    changed = False
                    for cj, y in enumerate(crossings):
                        for t in range(4):
                            if (cj, t) in enter:
                                others = [e for e in adj[y[t]] if e != (cj, t)]
                                if others and others[0] not in enter:
                                    ck, u = others[0]; enter[(ck, u)] = not enter[(cj, t)]; enter[(ck, (u + 2) % 4)] = enter[(cj, t)]; changed = True
    return enter, adj

def traverse(crossings):
    """Returns (bad, writhe, n_components, order) with bad = crossing indices first met on the under strand."""
    enter, adj = orient_slots(crossings)
    n = len(crossings)
    # entering slot of each label
    enter_at = {}
    for ci, x in enumerate(crossings):
        for s in range(4):
            if enter[(ci, s)]:
                enter_at[x[s]] = (ci, s)
    seen = [False] * n
    first_slot = [None] * n
    comps = 0
    labels_done = set()
    for lab in sorted(enter_at):
        if lab in labels_done:
            continue
        comps += 1
        cur = lab
        while cur not in labels_done:
            labels_done.add(cur)
            ci, s = enter_at[cur]
            if first_slot[ci] is None:
                first_slot[ci] = s
            cur = crossings[ci][(s + 2) % 4]
    # signs
    w = 0
    for ci, x in enumerate(crossings):
        o_in = 1 if enter[(ci, 1)] else 3
        w += 1 if o_in == 3 else -1      # under 0->2, over 3->1 = positive (KnotTheory PositiveQ: l -> j)
    bad = [ci for ci in range(n) if first_slot[ci] in (0, 2)]
    return bad, w, comps

def relabel(crossings):
    """Deterministic relabelling by traversal order -> memo key (oriented structure preserved)."""
    enter, adj = orient_slots(crossings)
    enter_at = {}
    for ci, x in enumerate(crossings):
        for s in range(4):
            if enter[(ci, s)]:
                enter_at[x[s]] = (ci, s)
    newlab = {}
    for lab in sorted(enter_at):
        if lab in newlab:
            continue
        cur = lab
        while cur not in newlab:
            newlab[cur] = len(newlab)
            ci, s = enter_at[cur]
            cur = crossings[ci][(s + 2) % 4]
    out = sorted(tuple(newlab[v] for v in x) for x in crossings)
    return tuple(out)

def switch(crossings, ci):
    """Switch crossing ci: the over strand becomes the under strand; keep the oriented-PD form (under enters at slot 0)."""
    enter, adj = orient_slots(crossings)
    a, b, c, d = crossings[ci]
    new = list(crossings)
    new[ci] = (d, a, b, c) if enter[(ci, 3)] else (b, c, d, a)
    return new

def smooth(crossings, ci, mode):
    """Remove crossing ci, join (a,b)+(c,d) for mode 0, (a,d)+(b,c) for mode 'inf'. Returns (oriented crossings, free loops)."""
    a, b, c, d = crossings[ci]
    pairs = [(a, b), (c, d)] if mode == 0 else [(a, d), (b, c)]
    parent = {}
    def find(x):
        while parent.get(x, x) != x:
            x = parent[x]
        return x
    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry
    for x, y in pairs:
        union(x, y)
    rest = [crossings[k] for k in range(len(crossings)) if k != ci]
    new = [tuple(find(v) for v in x) for x in rest]
    used = {v for x in new for v in x}
    loops = sum(1 for r in {find(v) for v in (a, b, c, d)} if r not in used)
    return reorient(new), loops

def reorient(crossings):
    """Unoriented crossings (under at slots 0/2) -> oriented PD: pick directions by walking, rotate so under enters at 0."""
    if not crossings:
        return []
    adj = build_adjacency(crossings)
    enter = {}
    for ci in range(len(crossings)):
        for s0 in range(4):
            if (ci, s0) in enter:
                continue
            # start a walk entering (ci, s0)
            cj, t = ci, s0
            while (cj, t) not in enter:
                enter[(cj, t)] = True; enter[(cj, (t + 2) % 4)] = False
                lab = crossings[cj][(t + 2) % 4]
                others = [e for e in adj[lab] if e != (cj, (t + 2) % 4)]
                if not others:
                    break
                cj, t = others[0]
    out = []
    for ci, x in enumerate(crossings):
        if enter[(ci, 0)]:
            out.append(tuple(x))
        else:
            out.append((x[2], x[3], x[0], x[1]))
    return out

@functools.lru_cache(maxsize=None)
def _L(key, loops):
    crossings = [tuple(x) for x in key]
    if not crossings:
        return ppow(MU, loops - 1) if loops >= 1 else ONE
    bad, w, comps = traverse(crossings)
    if not bad:
        return pmul(mono(w, 0), ppow(MU, comps + loops - 1))
    ci = bad[0]
    sw = switch(crossings, ci)
    s0, l0 = smooth(crossings, ci, 0)
    si, li = smooth(crossings, ci, 'inf')
    t_sw = _L(relabel(sw), loops)
    t_0 = _L(relabel(s0), loops + l0)
    t_i = _L(relabel(si), loops + li)
    return padd(pscale(t_sw, -1), pmul(mono(0, 1), padd(t_0, t_i)))

def kauffman_F(pd):
    """F(a, z) of the link with oriented PD code pd (KnotInfo/KnotTheory convention)."""
    crossings = [tuple(x) for x in pd]
    if not crossings:
        return ONE
    bad, w, comps = traverse(crossings)
    L = _L(relabel(crossings), 0)
    return pmul(mono(-w, 0), L)

def mirror(p):
    return {(-i, j): v for (i, j), v in p.items()}

# ----------------------------------------------------------------------------- KnotInfo helpers (calibration only)
def parse_knotinfo_F(s):
    """KnotInfo string such as '(-a^(-4)-2*a^(-2))*z^(0)+ (a^(-5)+ a^(-3))*z^(1)' or '(-a^(-2)-1-a^2)*z^(0)+ …' -> dict.
    Evaluated with sympy (exact), so every spelling KnotInfo uses is covered; calibration only."""
    import sympy
    a, z = sympy.symbols('a z')
    expr = sympy.sympify(s.replace('^', '**'), locals={'a': a, 'z': z})
    expr = sympy.expand(expr)
    out = {}
    for term in sympy.Add.make_args(expr):
        c, rest = term.as_coeff_Mul()
        pw = rest.as_powers_dict()
        i = int(pw.get(a, 0)); j = int(pw.get(z, 0))
        for sym in pw:
            if sym not in (a, z) and not getattr(sym, 'is_number', False):
                raise ValueError(f"unexpected symbol {sym} in {s[:60]}")
        out[(i, j)] = out.get((i, j), 0) + int(c)
    return {k: v for k, v in out.items() if v}

def knotinfo_rows(names=None, cmax=None):
    from database_knotinfo import link_list
    rows = []
    for r in link_list():
        if r['name'] == '0_1' or not str(r['crossing_number']).strip().isdigit():
            continue   # unknot and KnotInfo's header row
        if names and r['name'] not in names:
            continue
        if cmax is not None and int(r['crossing_number']) > cmax:
            continue
        rows.append(r)
    return rows

# ----------------------------------------------------------------------------- self-test and CLI
TREFOIL = [[1, 5, 2, 4], [3, 1, 4, 6], [5, 3, 6, 2]]
FIG8 = [[4, 2, 5, 1], [8, 6, 1, 5], [6, 3, 7, 4], [2, 7, 3, 8]]
TREFOIL_F = {(-4, 0): -1, (-2, 0): -2, (-5, 1): 1, (-3, 1): 1, (-4, 2): 1, (-2, 2): 1}   # KnotInfo 3_1
FIG8_F = {(-2, 0): -1, (0, 0): -1, (2, 0): -1, (-1, 1): -1, (1, 1): -1, (-2, 2): 1, (0, 2): 2, (2, 2): 1, (-1, 3): 1, (1, 3): 1}  # KnotInfo 4_1

def selftest():
    ok = True
    def chk(name, cond):
        nonlocal ok
        print(f"  {'OK ' if cond else '🔴 FAIL'} {name}")
        ok &= cond
    chk("unknot with one positive kink: F = 1", kauffman_F([[1, 2, 2, 1]]) == ONE)
    F3 = kauffman_F(TREFOIL)
    chk("trefoil F == KnotInfo 3_1 (fixes the mirror convention)", F3 == TREFOIL_F)
    chk("trefoil F != mirror(F) (chiral)", F3 != mirror(F3))
    F4 = kauffman_F(FIG8)
    chk("figure-eight F == KnotInfo 4_1", F4 == FIG8_F)
    chk("figure-eight F == mirror(F) (amphichiral)", F4 == mirror(F4))
    # counter-case: a wrong crossing (switch one crossing of the trefoil -> unknot diagram) must NOT give the trefoil F
    Fx = kauffman_F(switch(TREFOIL, 0))
    chk("switched trefoil (an unknot diagram) F == 1", Fx == ONE)
    chk("switched trefoil F != trefoil F (counter-case)", Fx != F3)
    # parser round trip on the KnotInfo trefoil string
    s = "(-a^(-4)-2*a^(-2))*z^(0)+ (a^(-5)+ a^(-3))*z^(1)+ (a^(-4)+ a^(-2))*z^(2)"
    chk("KnotInfo string parser on 3_1", parse_knotinfo_F(s) == TREFOIL_F)
    chk("a-span of trefoil F = 3", a_span(F3)[2] == 3)
    print("SELFTEST", "GREEN" if ok else "RED")
    return ok

def main(argv):
    if '--selftest' in argv or '--selbsttest' in argv:   # --selbsttest = old spelling, kept for earlier receipts
        sys.exit(0 if selftest() else 1)
    if '--pd' in argv:
        pd = ast.literal_eval(argv[argv.index('--pd') + 1])
        t = time.time(); F = kauffman_F(pd)
        print(pstr(F)); print("a-span (min, max, span):", a_span(F), f"[{time.time()-t:.2f}s]")
        return
    if '--calibrate' in argv:
        cmax = int(argv[argv.index('--calibrate') + 1])
        rows = knotinfo_rows(cmax=cmax)
        n_ok = n_bad = 0; bad = []
        t0 = time.time()
        for r in rows:
            pd = ast.literal_eval(r['pd_notation'])
            F = kauffman_F(pd); ref = parse_knotinfo_F(r['kauffman_polynomial'])
            if F == ref:
                n_ok += 1
            else:
                n_bad += 1; bad.append(r['name'])
        print(f"calibration c<={cmax}: {len(rows)} knots, {n_ok} == KnotInfo, {n_bad} differ {bad[:10]} [{time.time()-t0:.1f}s]")
        sys.exit(0 if n_bad == 0 and n_ok > 0 else 1)
    if '--knotinfo' in argv:
        names = set(argv[argv.index('--knotinfo') + 1:])
        for r in knotinfo_rows(names=names):
            pd = ast.literal_eval(r['pd_notation'])
            t = time.time(); F = kauffman_F(pd); ref = parse_knotinfo_F(r['kauffman_polynomial'])
            print(f"{r['name']}: own F == KnotInfo F: {F == ref}; a-span own {a_span(F)} KnotInfo {a_span(ref)} [{time.time()-t:.1f}s]")
            print("   own F =", pstr(F))
        return
    print(__doc__)

if __name__ == '__main__':
    main(sys.argv[1:])
