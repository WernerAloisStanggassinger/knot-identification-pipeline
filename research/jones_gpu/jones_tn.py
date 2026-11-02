#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""jones_tn.py — Jones polynomial at a point, exactly mod p, as a tensor-network contraction of a PD code (reference, NumPy).

Kauffman bracket with LOCAL loop counting (no height function needed):
  every edge e carries an index s_e in {+1,-1} = direction in which a state loop traverses e (relative to e's reference direction);
  every crossing tensor = A * (A-smoothing) + A^-1 * (B-smoothing); a smoothing joins two pairs of ccw-adjacent ports (corners);
  a corner is consistent iff the loop enters through one port and leaves through the other; it carries x^(+1) for a left turn
  (enter port p, leave port p-1) and x^(-1) for a right turn (enter p, leave p+1), ports numbered ccw;
  every edge carries x^(s_e * tau_e) where tau_e (units of pi/2, along the reference direction) solves the face equations
  "corner turns + edge turns = +4 (bounded faces) / -4 (outer face)" of the planar 4-valent graph (integer solution via a spanning
  tree of the dual graph).
  Claim (our argument, see README; checked by calibrate_jones_tn.py, not proven here): a simple closed loop then has total
  turning +-4 units, its two orientations give x^4 + x^-4; with x^4 = -A^2 that is
  -A^2 - A^-2 = d.  So the full contraction = <K> with <O> = d; <K>_1 = <K>/d; V(t) = (-A^3)^(-w) <K>_1, t = A^-4.
Arithmetic: GF(p), p = 1 mod 4 (so i = sqrt(-1) exists); x free, A = i x^2 (A^2 = -x^4), t = A^-4 = x^-8.
Which smoothing is "A" and the orientation of t are fixed by calibration against Sage (calibrate_jones_tn.py), not assumed.
"""
import numpy as np

P = 536870909          # prime, = 1 mod 4, < 2^29 (products < 2^58, sums of 16 products < 2^62)


def inv(a, p=P): return pow(int(a) % p, p - 2, p)


def sqrt_m1(p=P):
    for g in range(2, 200):
        r = pow(g, (p - 1) // 4, p)
        if r * r % p == p - 1: return r
    raise ValueError


def edge_orientation(pd):
    """For each crossing c and slot k: +1 if the edge in slot k points INTO c (knot orientation), -1 if out of c."""
    n2 = 2 * len(pd)
    io = []
    for (a, b, c, d) in pd:
        o = [+1, 0, -1, 0]                         # slot 0 = incoming under, slot 2 = outgoing under
        if (d - b) % n2 == 1: o[1], o[3] = +1, -1   # over strand b -> d
        elif (b - d) % n2 == 1: o[1], o[3] = -1, +1
        else: raise ValueError('cannot orient over strand (edge labels not consecutive)')
        io.append(o)
    return io


def edge_turnings(pd):
    """Integer tau_e (units pi/2, along the reference direction tail->head) solving the face equations."""
    nc = len(pd); where = {}
    for c, X in enumerate(pd):
        for k, e in enumerate(X): where.setdefault(e, []).append((c, k))
    io = edge_orientation(pd)
    other = {}
    for e, ((c1, k1), (c2, k2)) in where.items():
        other[(c1, k1)] = (c2, k2); other[(c2, k2)] = (c1, k1)
    # face traversal (face on the left): arrive at crossing c through slot k, leave through slot k-1
    seen, faces = set(), []
    for c in range(nc):
        for k in range(4):
            if (c, k) in seen: continue
            face = []; cur = (c, k)          # cur = (crossing, slot) through which we LEAVE crossing
            while cur not in seen:
                seen.add(cur)
                cc, kk = cur; e = pd[cc][kk]
                sign = -io[cc][kk]            # leaving cc along e: along reference direction iff e points out of cc
                face.append((e, sign))
                c2, k2 = other[cur]          # arrive at c2 through slot k2
                cur = (c2, (k2 - 1) % 4)
            faces.append(face)
    assert len(faces) == nc + 2, (len(faces), nc)
    # equations: corners (= len(face)) + sum sign*tau_e = +4, outer face (index 0) = -4
    edges = sorted(where); idx = {e: i for i, e in enumerate(edges)}
    rhs = [(-4 if f == 0 else 4) - len(faces[f]) for f in range(len(faces))]
    assert sum(rhs) == 0
    # dual graph: faces are vertices, each edge joins the two faces on its sides; spanning tree, solve from the leaves
    sides = {}
    for f, face in enumerate(faces):
        for e, s in face: sides.setdefault(e, []).append((f, s))
    tau = {e: 0 for e in edges}
    adj = {f: [] for f in range(len(faces))}
    for e, ((f1, s1), (f2, s2)) in sides.items():
        adj[f1].append((e, f2)); adj[f2].append((e, f1))
    # BFS tree from face 0, then fix tree edges bottom-up so every non-root face equation holds; root follows from sum
    parent, order, seenf = {0: None}, [0], {0}
    for f in order:
        for e, g in adj[f]:
            if g not in seenf: seenf.add(g); parent[g] = (e, f); order.append(g)
    tree = {parent[g][0] for g in order[1:]}
    resid = list(rhs)
    for g in reversed(order[1:]):
        e, f = parent[g]
        # face g equation: sum over its edges sign*tau = resid[g]; all its other edges are already fixed (tree children or non-tree = 0)
        s_g = [s for (ff, s) in sides[e] if ff == g][0]
        known = sum(s * tau[ee] for ee, s in faces[g] if ee != e)
        tau[e] = (resid[g] - known) * s_g     # s_g = +-1
    for f, face in enumerate(faces):
        assert sum(s * tau[e] for e, s in face) + len(face) == (-4 if f == 0 else 4), 'face equation'
    return tau, io


def writhe(pd, io):
    w = 0
    for (a, b, c, d), o in zip(pd, io):
        w += 1 if o[3] == +1 else -1          # crossing sign from the over strand's direction; the global sign convention is
                                              # fixed by calibrate_jones_tn.py (it agrees with Sage for wsign = +1)
    return w


def crossing_tensors(pd, x, A_is_02=True, p=P):
    """The c crossing tensors (c, 2, 2, 2, 2) mod p at the point x (edge weights on the tail slots), A, writhe."""
    tau, io = edge_turnings(pd)
    xi = inv(x, p); I = sqrt_m1(p); A = I * x % p * x % p; Ai = inv(A, p)
    where = {}
    for c, X in enumerate(pd):
        for k, e in enumerate(X): where.setdefault(e, []).append((c, k))
    tails = {}
    for c, X in enumerate(pd):
        for k, e in enumerate(X):
            if io[c][k] == -1: tails[e] = (c, k)     # edge leaves c here = tail
    def xp(k): return pow(x, k, p) if k >= 0 else pow(xi, -k, p)
    tensors = []
    for c, X in enumerate(pd):
        T = np.zeros((2, 2, 2, 2), dtype=np.int64)
        for idx in np.ndindex(2, 2, 2, 2):
            s = [1 - 2 * i for i in idx]                 # index 0 -> +1, 1 -> -1
            enter = [s[k] * io[c][k] == 1 for k in range(4)]   # loop moves toward c along slot k
            val = 0
            for smoothing, coef in (((0, 1), (2, 3)), A if A_is_02 else Ai), (((0, 3), (1, 2)), Ai if A_is_02 else A):
                w = coef
                for (k1, k2) in smoothing:
                    if enter[k1] == enter[k2]: w = 0; break
                    ein, eout = (k1, k2) if enter[k1] else (k2, k1)
                    turn = +1 if eout == (ein - 1) % 4 else -1      # leave through p-1 = left turn
                    w = w * xp(turn) % p
                val = (val + w) % p
            # edge weights x^(s_e tau_e) on the tail slots
            for k in range(4):
                e = X[k]
                if tails[e] == (c, k): val = val * xp(s[k] * tau[e]) % p
            T[idx] = val
        tensors.append(T)
    return np.array(tensors), A, writhe(pd, io)


def bracket(pd, x, A_is_02=True, p=P):
    """Full contraction <K> (<O> = d) at the point x, mod p. Returns int."""
    tensors, A, w = crossing_tensors(pd, x, A_is_02, p)
    # contract greedily in PD order, keeping open edges as axes
    cur = np.ones((), dtype=np.int64); open_e = []
    for c, X in enumerate(pd):
        T = tensors[c]; labs = list(X)
        # contract shared edges
        shared = [e for e in labs if e in open_e]
        ax_cur = [open_e.index(e) for e in shared]; ax_T = [labs.index(e) for e in shared]
        R = np.tensordot(cur, T, axes=(ax_cur, ax_T)) % p
        open_e = [e for e in open_e if e not in shared] + [e for e in labs if e not in shared]
        # an edge whose both ends are in this crossing (a kink) would appear twice in labs: trace it
        while len(set(open_e)) < len(open_e):
            e = next(e for e in open_e if open_e.count(e) == 2)
            i, j = [k for k, ee in enumerate(open_e) if ee == e]
            R = np.trace(R, axis1=i, axis2=j) % p
            open_e = [ee for k, ee in enumerate(open_e) if k not in (i, j)]
        cur = R
    assert cur.shape == (), cur.shape
    return int(cur) % p, A, w


def jones_at(pd, x, A_is_02=True, wsign=1, p=P):
    """V at t = A^-4 = x^-8 (mod p): (-A^3)^(-w) <K>/d."""
    br, A, w = bracket(pd, x, A_is_02, p)
    d = (-A * A - inv(A * A, p)) % p
    assert d != 0, 'd = 0 mod p at this point (x^8 = -1): choose another x'
    v = br * inv(d, p) % p
    m = (-pow(A, 3, p)) % p
    w = wsign * w
    f = pow(inv(m, p), w, p) if w >= 0 else pow(m, -w, p)
    return v * f % p, pow(inv(A, p), 4, p)


# ---------------------------------------------------------------- fast path: x-independent exponent form
_S = np.array([[1 - 2 * ((a >> (3 - k)) & 1) for k in range(4)] for a in range(16)], dtype=np.int64)   # (16,4) leg signs


def _corner_table(io_c):
    """For one crossing orientation pattern: per entry a (16) and smoothing m (2): valid, corner-turn exponent."""
    valid = np.zeros((16, 2), dtype=bool); turn = np.zeros((16, 2), dtype=np.int64)
    for a in range(16):
        s = _S[a]; enter = [s[k] * io_c[k] == 1 for k in range(4)]
        for m, smoothing in enumerate((((0, 1), (2, 3)), ((0, 3), (1, 2)))):
            ok, t = True, 0
            for (k1, k2) in smoothing:
                if enter[k1] == enter[k2]: ok = False; break
                ein, eout = (k1, k2) if enter[k1] else (k2, k1)
                t += +1 if eout == (ein - 1) % 4 else -1
            valid[a, m] = ok; turn[a, m] = t
    return valid, turn


_TAB = {}


def crossing_exponents(pd):
    """x-independent form of the crossing tensors: entry (c,a) = sum_m valid * i^(+-1) x^(expo), A = i x^2 (m=0: A, m=1: A^-1).
    Returns valid (c,16,2), expo (c,16,2), isign (c,16,2) = +1 for i^1, -1 for i^-1, and the writhe."""
    tau, io = edge_turnings(pd)
    c = len(pd)
    tails = np.zeros((c, 4), dtype=np.int64)
    for v, X in enumerate(pd):
        for k, e in enumerate(X):
            if io[v][k] == -1: tails[v, k] = tau[e]
    edge = tails @ _S.T                                   # (c,16): sum over tail legs of s_k * tau_e
    valid = np.zeros((c, 16, 2), dtype=bool); expo = np.zeros((c, 16, 2), dtype=np.int64)
    for v in range(c):
        key = tuple(io[v])
        if key not in _TAB: _TAB[key] = _corner_table(io[v])
        vd, tn = _TAB[key]
        valid[v] = vd; expo[v] = tn + edge[v][:, None] + np.array([2, -2])[None, :]   # A^{+-1} = i^{+-1} x^{+-2}
    isign = np.broadcast_to(np.array([1, -1]), (c, 16, 2)).copy()
    return valid, expo, isign, writhe(pd, io)


def eval_tensors(valid, expo, isign, x, p=P):
    """Crossing tensors (c,16) mod p at the point x from the exponent form."""
    I = sqrt_m1(p); Ii = inv(I, p)
    lo, hi = int(expo.min()), int(expo.max())
    xi = inv(x, p)
    pw = np.array([pow(x, e, p) if e >= 0 else pow(xi, -e, p) for e in range(lo, hi + 1)], dtype=np.int64)
    val = pw[expo - lo] * np.where(isign == 1, I, Ii) % p
    val = np.where(valid, val, 0)
    return (val[..., 0] + val[..., 1]) % p
