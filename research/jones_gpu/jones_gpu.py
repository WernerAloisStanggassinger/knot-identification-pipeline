#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""jones_gpu.py — batched GPU contraction of the Kauffman-bracket tensor networks of jones_tn.py (many knots per kernel).

Per knot (CPU): crossing tensors at the point x (jones_tn.crossing_tensors), a crossing order that keeps the frontier of open
edges small (greedy: next = the unprocessed crossing sharing most edges with the frontier), and a slot (bit) for every open edge.
Per batch (GPU, CuPy, int64 mod p): state (B, 2^F); a crossing step sums over its 16 leg assignments a:
  out[y] += state[y_in] * T[a]   for y with (new-leg bits of y) = a_new and (freed slots of y) = 0, y_in = y without new bits | a_shared.
After the last crossing the state is a scalar per knot = <K> (with <O> = d); V(t) as in jones_tn.jones_at.
Checked separately: test_gpu_vs_ref.py compares the result with jones_tn.jones_at knot by knot (there is no check inside a run).
Limit: plan() refuses diagrams with a kink (an edge with both ends at one crossing). jones_tn.bracket has code for kinks, but it is
not covered by the shipped tests (one spelling of a one-crossing diagram fails its face-equation check).
"""
import numpy as np
from jones_tn import crossing_tensors, P, inv


def plan(pd):
    """Crossing order, per step the 4 legs' slots and whether each leg is new (opens) or shared (closes); F = slots used."""
    nc = len(pd)
    done, order, open_e = set(), [], {}
    free = []; nslots = 0
    cur = 0
    steps = []
    for _ in range(nc):
        if order:
            best, bscore = None, None
            for c in range(nc):
                if c in done: continue
                sh = sum(1 for e in pd[c] if e in open_e)
                sc = (sh, -c)
                if bscore is None or sc > bscore: best, bscore = c, sc
            cur = best
        c = cur; done.add(c); order.append(c)
        legs = []
        X = pd[c]
        counts = {e: X.count(e) for e in X}
        for k, e in enumerate(X):
            if e in open_e:
                legs.append((open_e[e], 0))                  # shared: close
            else:
                legs.append((None, 1))
        # kink: an edge appearing twice at this crossing and not yet open -> both legs new, then must be traced; not supported
        if any(v == 2 for v in counts.values()):
            raise ValueError('kink (edge with both ends at one crossing) — simplify the diagram first')
        # allocate the new legs' slots BEFORE freeing the closing ones: a slot must not be closed and opened in the same step
        for k, e in enumerate(X):
            if legs[k][1] == 1:
                if free: s = free.pop()
                else: s = nslots; nslots += 1
                legs[k] = (s, 1)
        for k, e in enumerate(X):
            if legs[k][1] == 0: free.append(open_e.pop(e))
        for k, e in enumerate(X):
            if legs[k][1] == 1: open_e[e] = legs[k][0]
        steps.append(legs)
    assert not open_e
    return order, steps, nslots


def contract_batch(items, F, p=P):
    """items: list of (tensors (c,16) int64 in plan order, steps). All knots have the same c. Returns <K> per knot (np int64)."""
    import cupy as cp                      # imported here so that plan()/prepare_fast() work without CuPy
    B = len(items); c = len(items[0][1])
    state = cp.zeros((B, 1 << F), dtype=cp.int64); state[:, 0] = 1
    y = cp.arange(1 << F, dtype=cp.int64)[None, :]
    T = cp.asarray(np.stack([it[0] for it in items]))           # (B, c, 16)
    for st in range(c):
        slots = np.array([[l[0] for l in it[1][st]] for it in items], dtype=np.int64)     # (B,4)
        new = np.array([[l[1] for l in it[1][st]] for it in items], dtype=np.int64)       # (B,4)
        bit = (1 << slots)                                                               # (B,4)
        newmask = cp.asarray((bit * new).sum(1))[:, None]
        shmask = cp.asarray((bit * (1 - new)).sum(1))[:, None]
        bitg = cp.asarray(bit); newg = cp.asarray(new)
        ok = (y & shmask) == 0                                                           # freed slots must be 0 in the output
        ybase = y & ~newmask
        out = cp.zeros_like(state)
        for a in range(16):
            ab = cp.asarray(np.array([(a >> (3 - k)) & 1 for k in range(4)], dtype=np.int64))[None, :]   # leg k value bit
            pat_new = (bitg * ab * newg).sum(1)[:, None]
            pat_sh = (bitg * ab * (1 - newg)).sum(1)[:, None]
            m = ok & ((y & newmask) == pat_new)
            yin = ybase | pat_sh
            v = cp.take_along_axis(state, cp.broadcast_to(yin, state.shape), axis=1)
            out = (out + cp.where(m, v * T[:, st, a][:, None] % p, 0)) % p
        state = out
    return cp.asnumpy(state[:, 0])


def prepare(pd, x, p=P):
    order, steps, F = plan(pd)
    T, A, w = crossing_tensors(pd, x, True, p)
    T = T[order].reshape(len(pd), 16)
    return T, steps, F, A, w


def prepare_fast(pd, xs, p=P):
    """One plan and one exponent form for all points xs: returns steps, F, w, [(T (c,16) in plan order, A) per x]."""
    from jones_tn import crossing_exponents, eval_tensors, sqrt_m1
    order, steps, F = plan(pd)
    valid, expo, isign, w = crossing_exponents(pd)
    I = sqrt_m1(p); out = []
    for x in xs:
        T = eval_tensors(valid, expo, isign, x, p)[order]
        out.append((T, I * x % p * x % p))
    return steps, F, w, out


def jones_from_bracket(br, A, w, p=P):
    d = (-A * A - inv(A * A, p)) % p
    assert d != 0, 'd = 0 mod p at this point (x^8 = -1): choose another x'
    v = int(br) * inv(d, p) % p
    m = (-pow(A, 3, p)) % p
    f = pow(inv(m, p), w, p) if w >= 0 else pow(m, -w, p)
    return v * f % p
