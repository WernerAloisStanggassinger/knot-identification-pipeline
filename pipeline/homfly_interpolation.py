#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""homfly_interpolation.py — the FULL HOMFLY polynomial from point evaluations
of the Hecke trace, where regina stands at the memory wall.

OCCASION (2026-08-23): regina computes 178 of the 229 missing HOMFLY polynomials, but
gets stuck on 36 — not on TIME but on MEMORY (~3 GB per process,
some shot down by earlyoom, Disc #677). **30 of these 36 lie at n <= 12**,
i.e. within the range of our own Hecke trace. The two tools are
complementary: regina dies on c, the Hecke trace on n (cost ~ n!).
The Hecke trace, however, only delivers NUMBERS at a point. This tool turns
them into the polynomial — as was done for theta back then (project-theta-throughput-by-interpolation).

── THE BOX IS MEASURED, NOT SET ─────────────────────────────────────────────
On the 178 freshly computed polynomials (each with a braid word in the record):
  · a exponents EVEN without exception        178/178
  · z exponents EVEN and >= 0 without exc.    178/178
  · -e-n+1 <= a exponent <= -e+n-1            178/178 violated: 0
    (e = exponent sum of the braid word. With +e the bound would be violated in 170 of
     178 cases — the sign is measured, not assumed.)
    173 of the 178 touch the boundary EXACTLY ⟹ the bound is sharp.
  · z degree <= simplified_size - 1           178/178 violated: 0
So the rows and columns of the linear system are fixed before anything is computed.

── WHY IT HOLDS NUMERICALLY ─────────────────────────────────────────────────
Substitution A = a^2, Z = z^2 (both exponents are even) halves both
axes. The nodes lie on the UNIT CIRCLE, rotated by a
generic phase: then the Vandermonde matrix is a DFT — condition number 1.
A grid on the real axis would be worthless here (coefficients up to 7.6e12,
measured).

── WHY IT IS AFFORDABLE ─────────────────────────────────────────────────────
The expensive matrix chain depends ONLY on q, not on a: P(a,z) = a^-e * Sum_lam
W_lam(a,q) * chi_lam(q). So ONE pass per z value, after which all
a values are almost free. Instead of N_A*N_Z passes only N_Z.

── VERIFICATION ─────────────────────────────────────────────────────────────
Interpolation is a claim, not a proof. Only what holds afterwards
counts: (1) fresh points that were NOT on the grid, (2) det = |P(1,2i)|
exactly against the record, (3) |Alexander(-1)|. And the coefficients must
come out integral — the distance to the nearest integer is the
sharpest quality measure there is here.
"""
import sys, os, json, ast, time, argparse, cmath
import site_config
import multiprocessing as mp

CORES = None   # global, set from --kerne; otherwise EICH_KERNE (see below)
# ---------------------------------------------------------------------------
# CALIBRATED DEFAULT FOR THE WORKER COUNT — measured 2026-08-25, not guessed.
#
# This used to be os.cpu_count() (= 24 on this box), from Hard Rule 6
# ("full throttle, no ÷2"). But the rule means: the BOX must not idle — it
# does not mean: a single job gets all cores. Measured on T_13_7_2_6
# (n=11, real braid word, 27 q nodes, kern_eichung2.py):
#
#     cores   wall s     CPU s   speedup  efficiency
#         1    854.1     853.7     1.00x      1.00
#         6    330.9    1744.8     2.58x      0.49   <-- optimum
#        12    388.3    3997.0     2.20x      0.21   slower than 6
#        18    414.4    5524.4     2.06x      0.15   slower still
#     (24 never measured cleanly — the run was spoiled by a second computation
#      on the same box; the row is left out on purpose instead of guessed.)
#
# From ~6 on the memory bus is saturated: dense int64 matrices are
# bandwidth-bound, not compute-bound, and the 12 physical cores (Ryzen 9
# 3900X, SMT) share one. Every further worker just waits along and
# costs full CPU time.
#
# LOADING THE BOX PROPERLY therefore means: SEVERAL KNOTS side by side, each with
# EICH_KERNE workers (4 x 6 = 24). The parallelism belongs one level
# higher, not inside the single knot.
#
# Override: --kerne <n> always wins; without the argument
# HOMFLY_KERNE from the environment applies; only then this value. On another
# box recalibrate (kern_eichung2.py), do not estimate.
# ---------------------------------------------------------------------------
CALIBRATION_CORES = int(os.environ.get("HOMFLY_KERNE", 6))

import numpy as np

HECKE = (site_config.HECKE_DIR)
sys.path.insert(0, HECKE)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from hecke_homfly_fast import (syt, hecke_gen_sparse, partitions, W)
import scipy.sparse as sp


# ── braid word / box ───────────────────────────────────────────────────────────
def as_list(w):
    if isinstance(w, str):
        try:
            w = ast.literal_eval(w)
        except Exception:
            return None
    if not isinstance(w, (list, tuple)) or not w:
        return None
    return [int(g) for g in w]


def box(word, n, z_grad_max):
    """(I0, I1, J1) for P = Sum C_IJ * A^I * Z^J with A=a^2, Z=z^2.
    Measured rules, with a margin of one step each."""
    e = sum(1 if g > 0 else -1 for g in word)
    a_lo, a_hi = -e - n + 1, -e + n - 1
    I0, I1 = (a_lo - 1) // 2, (a_hi + 1) // 2      # one step margin
    J1 = (z_grad_max + 2) // 2
    return I0, I1, J1, e


# ── Hecke trace: ONE pass per q yields all chi_lam ───────────────────────────
def chis(word, n, q):
    out = []
    for lam in partitions(n):
        T = syt(lam)
        idx = {tuple(t): k for k, t in enumerate(T)}
        dim = len(T)
        sc = q - q ** -1
        I = sp.identity(dim, dtype=complex, format='csr')
        Gp, Gi = {}, {}
        M = np.eye(dim, dtype=complex)
        for g in word:
            i = abs(g)
            if i not in Gp:
                Gp[i] = hecke_gen_sparse(T, idx, i, q)
                Gi[i] = (Gp[i] - sc * I).tocsr()
            G = Gi[i] if g < 0 else Gp[i]
            M = (G.T @ M.T).T
        out.append((lam, complex(np.trace(M))))
    return out


def norm_K(a, q):
    """The norm that the author of the Hecke trace left out on purpose
    ('×common norm K'); found from the data on 23.08. and calibrated on 56 probes
    at 4 points (homfly_eichung.py --pruefen)."""
    return (q - 1 / q) / (a - 1 / a)


def P_at(word, n, a, q, pre=None):
    """P_AZ(a, z=q-1/q) at one point. `vor` = precomputed chis(q)."""
    e = sum(1 if g > 0 else -1 for g in word)
    cs = pre if pre is not None else chis(word, n, q)
    s = sum(W(lam, a, q) * ch for lam, ch in cs)
    return a ** (-e) * s * norm_K(a, q)


def q_from_z(z):
    """z = q - 1/q  ⟹  q^2 - z q - 1 = 0. Take the root of larger modulus:
    the smaller one pushes q close to the unit circle and lets the denominators
    (q^h - q^-h) become small."""
    w = cmath.sqrt(z * z + 4)
    q1, q2 = (z + w) / 2, (z - w) / 2
    return q1 if abs(q1) >= abs(q2) else q2


# ── one row of the value matrix (for the pool) ──────────────────────────────
def _row(task):
    """One z node: yields (l, [P(A_k) for all k]).

    PARALLELISED 2026-08-25 (Werner: stop slow runs and let the
    machine work). The loop over the z values was serial — on a
    24-core machine ONE core ran. Measured on T_13_7_2_6 (n=11, c=93):
    floating point 22.6 s per z value x 50 = ~19 min, exact 4575 s = 76 min.
    The rows are completely independent: chis(word,n,q) depends only on q,
    i.e. on z — exactly that has always stood in the head of this file
    ("ONE pass per z value, after which all a values are almost free"),
    it just was not used.
    Violated the full-steam doctrine in 00_START: "with >=2 independent
    jobs ALWAYS parallel, NEVER serial".
    """
    l, Z, word, n, As = task
    z = cmath.sqrt(Z)
    q = q_from_z(z)
    cs = chis(word, n, q)
    return l, [P_at(word, n, cmath.sqrt(A), q, pre=cs) for A in As]


# ── Interpolation ───────────────────────────────────────────────────────────
def interpolate(word, n, z_grad_max, phase=0.31, laut=False, forced_box=None,
                 kerne=None):
    """Returns (C, I0, J1, diag). C[J, m] = coefficient of A^(I0+m) Z^J.

    kerne=None -> os.cpu_count(); kerne=1 -> serial (counter-case path)."""
    I0, I1, J1, e = box(word, n, z_grad_max)
    if forced_box is not None:            # only for counter-cases in the self-test
        I0, I1, J1 = forced_box
    NA, NZ = I1 - I0 + 1, J1 + 1

    rotZ = cmath.exp(1j * phase)
    rotA = cmath.exp(1j * phase * 0.7)
    Zs = [rotZ * cmath.exp(2j * cmath.pi * l / NZ) for l in range(NZ)]
    As = [rotA * cmath.exp(2j * cmath.pi * k / NA) for k in range(NA)]

    if laut:
        print(f'  Box A^{I0}..A^{I1} ({NA}) x Z^0..Z^{J1} ({NZ}) '
              f'= {NA*NZ} unknowns, but only {NZ} matrix passes',
              flush=True)

    V = np.zeros((NZ, NA), dtype=complex)
    t0 = time.time()
    k_target = CALIBRATION_CORES if kerne is None else max(1, int(kerne))
    k_target = min(k_target, NZ)
    tasks = [(l, Z, word, n, As) for l, Z in enumerate(Zs)]

    def _serial(reason=None):
        if laut and reason:
            print(f'    {reason} — serial', flush=True)
        for auf in tasks:
            l, row_values = _row(auf)
            V[l, :] = row_values
            if laut and NZ >= 6 and (l + 1) % max(1, NZ // 6) == 0:
                print(f'    {l+1}/{NZ} z values, {time.time()-t0:.0f}s', flush=True)

    # Latch: spawn children import the module under its name. If it was
    # loaded under an invented name (importlib.spec_from_file_location),
    # pickling _zeile fails — then continue serially instead of dying.
    # (The texts passed to _seriell() are data for the gate and stay German.)
    _mod = sys.modules.get(_row.__module__)
    _importable = _mod is not None and (
        _row.__module__ == '__main__' or getattr(_mod, '__spec__', None) is not None)

    if k_target <= 1:                       # counter-case path, bit-identical
        _serial()
    elif not _importable:
        _serial(f"module name '{_row.__module__}' not importable")
        k_target = 1
    else:
        import concurrent.futures as _cf
        n_done = 0
        # spawn: the children import the module afresh; _zeile sits at
        # module level for that. fork would be faster but does not get along
        # well with numpy/BLAS threads.
        ctx = mp.get_context('spawn')
        try:
            with _cf.ProcessPoolExecutor(max_workers=k_target, mp_context=ctx) as ex:
                for l, row_values in ex.map(_row, tasks, chunksize=1):
                    V[l, :] = row_values
                    n_done += 1
                    if laut and NZ >= 6 and n_done % max(1, NZ // 6) == 0:
                        print(f'    {n_done}/{NZ} z values, {time.time()-t0:.0f}s '
                              f'({k_target} cores)', flush=True)
        except Exception as ex_:          # pool broken -> better slow than dead
            _serial(f'pool failed ({type(ex_).__name__})')
            k_target = 1
    t_hecke = time.time() - t0

    # V[l,k] = Sum_{J,m} C[J,m] * Z_l^J * A_k^(I0+m)  ⟹  V = MZ @ C @ MA.T
    MZ = np.array([[Z ** J for J in range(NZ)] for Z in Zs], dtype=complex)
    MA = np.array([[A ** (I0 + m) for m in range(NA)] for A in As], dtype=complex)
    kZ, kA = np.linalg.cond(MZ), np.linalg.cond(MA)
    C = np.linalg.solve(MZ, V) @ np.linalg.inv(MA).T

    Cr = np.rint(C.real)
    residue = float(np.max(np.abs(C - Cr))) if C.size else 0.0
    diag = {'kond_Z': float(kZ), 'kond_A': float(kA),
            'ausschuss': residue, 'sek_hecke': round(t_hecke, 1),
            'unbekannte': int(NA * NZ), 'durchgaenge': int(NZ),
            'imag_max': float(np.max(np.abs(C.imag))) if C.size else 0.0}
    return Cr.astype(object), I0, J1, e, diag, (As, Zs, V)


def as_dict(C, I0):
    """C[J,m] → {(a_exp, z_exp): coef}, back in a,z (not A,Z)."""
    out = {}
    NZ, NA = C.shape
    for J in range(NZ):
        for m in range(NA):
            c = int(C[J, m])
            if c:
                out[(2 * (I0 + m), 2 * J)] = c
    return out


def as_regina_string(p):
    """In the format regina prints (x=a, y=z) — so that the record carries one
    notation and no fourth convention is added."""
    def mono(i, j, c):
        t = []
        if abs(c) != 1 or (i == 0 and j == 0):
            t.append(str(abs(c)))
        if i:
            t.append('x' if i == 1 else f'x^{i}')
        if j:
            t.append('y' if j == 1 else f'y^{j}')
        return ' '.join(t)
    ks = sorted(p, key=lambda k: (-k[0], -k[1]))
    s = ''
    for i, j in ks:
        c = p[(i, j)]
        s += ('- ' if c < 0 else ('+ ' if s else '')) if s or c < 0 else ''
        s += mono(i, j, c) + ' '
    return s.strip()


# ── verification: interpolation is a claim, not a proof ────────────
def check_fresh(p, word, n, points=4, phase=1.777):
    """At points that were NOT on the grid. That is the real receipt:
    a grid point hits by construction, a fresh one does not."""
    worst = 0.0
    for t in range(points):
        a = 1.19 * cmath.exp(1j * (phase + 0.911 * t))
        z = 0.83 * cmath.exp(1j * (phase * 1.3 + 0.577 * t))
        q = q_from_z(z)
        expected = P_at(word, n, a, q)
        actual = sum(c * (a ** i) * (z ** j) for (i, j), c in p.items())
        denom_n = max(abs(expected), abs(actual), 1e-12)
        worst = max(worst, abs(expected - actual) / denom_n)
    return worst


def det_exact(p):
    """|P(a=1, z=2i)| as an integer — never in floating point (Disc #675d)."""
    re_, im_ = 0, 0
    for (i, j), c in p.items():
        if j < 0:
            return None
        w = c * (2 ** j)
        r = j % 4
        if r == 0:   re_ += w
        elif r == 1: im_ += w
        elif r == 2: re_ -= w
        else:        im_ -= w
    from math import isqrt
    q = re_ * re_ + im_ * im_
    w = isqrt(q)
    return w if w * w == q else None



# ── witness 4: Alexander AT SEVERAL POINTS, exact in fractions ───────────────
def alexander_from_homfly(p, t):
    """Delta(t) = P(a=1, z) with z^2 = t - 2 + 1/t.
    Because all z exponents are EVEN (measured, 178/178), the root is never
    needed — z^2 suffices, and so everything stays in Q. Exact, not
    numerical: the stored node values have up to 24 digits."""
    from fractions import Fraction
    T = Fraction(t)
    if T == 0:
        return None
    Z2 = T - 2 + Fraction(1, 1) / T
    columns = {}
    for (i, j), c in p.items():
        if j % 2:
            return None
        columns[j // 2] = columns.get(j // 2, 0) + c      # a = 1
    return sum(Fraction(c) * Z2 ** J for J, c in columns.items())


def alexander_witness(p, evals):
    """Compares Delta from the HOMFLY with the alexander_point_evals stored in
    the record — UP TO THE UNIT ±t^k (Footgun 12: Delta is
    only determined up to ±t^a, an absolute comparison reports phantom errors).
    Decisive: the same k must fit at ALL points, otherwise it does not
    count. A single point could always be explained away.
    The returned reason texts are data (written into records)."""
    from fractions import Fraction
    if not evals:
        return None, 'no alexander_point_evals in the record'
    pairs = []
    for ts, expected in evals.items():
        try:
            t = int(ts)
        except Exception:
            continue
        mine = alexander_from_homfly(p, t)
        if mine is None or expected in (None, 0):
            continue
        pairs.append((t, mine, Fraction(int(expected))))
    if len(pairs) < 2:
        return None, 'too few usable sample values'
    for k in range(-80, 81):
        for vz in (1, -1):
            if all(mine == vz * Fraction(t) ** k * expected for t, mine, expected in pairs):
                return True, f'agreeing at {len(pairs)} points (unit {vz:+d}*t^{k})'
    return False, f'at {len(pairs)} points not explained by any unit ±t^k'


# ── self-test: two known knots and two counter-cases ───────────────────
def selbsttest():
    ok = True
    faelle = [('Kleeblatt 3_1', [1, 1, 1], 2, 3, 3),
              ('Achterknoten 4_1', [1, -2, 1, -2], 3, 4, 5)]
    for nm, w, n, c, det_expected in faelle:
        C, I0, J1, e, diag, _ = interpolate(w, n, c)
        p = as_dict(C, I0)
        d = det_exact(p)
        fresh = check_fresh(p, w, n)
        good = (d == det_expected and fresh < 1e-8 and diag['ausschuss'] < 1e-6)
        print(f'  {nm:<18} det={d} (expected {det_expected}) · fresh {fresh:.1e} · '
              f'residue {diag["ausschuss"]:.1e} · cond {diag["kond_A"]:.2f}/'
              f'{diag["kond_Z"]:.2f}  {"✅" if good else "🔴"}')
        print(f'     {as_regina_string(p)}')
        ok &= good

    # COUNTER-CASE 1a: z axis too short (only Z^0, i.e. z degree 0) — the trefoil
    # needs z^2. The coefficients then cannot come out integral
    # and fresh points must break.
    C, I0, J1, e, d1, _ = interpolate([1, 1, 1], 2, 3, forced_box=(-2, -1, 0))
    f1 = check_fresh(as_dict(C, I0), [1, 1, 1], 2)
    fa1 = (f1 > 1e-6) or (d1['ausschuss'] > 1e-6)
    print(f'  COUNTER-CASE z axis too short: fresh {f1:.1e} · '
          f'residue {d1["ausschuss"]:.1e}  '
          f'{"✅ detected" if fa1 else "🔴 slips through"}')
    ok &= fa1

    # COUNTER-CASE 1b: a axis shifted (box sits beside the support)
    C, I0, J1, e, d2, _ = interpolate([1, 1, 1], 2, 3, forced_box=(2, 3, 1))
    f2 = check_fresh(as_dict(C, I0), [1, 1, 1], 2)
    fa2 = (f2 > 1e-6) or (d2['ausschuss'] > 1e-6)
    print(f'  COUNTER-CASE a axis off    : fresh {f2:.1e} · '
          f'residue {d2["ausschuss"]:.1e}  '
          f'{"✅ detected" if fa2 else "🔴 slips through"}')
    ok &= fa2

    # COUNTER-CASE 2: the trefoil polynomial must not explain the figure-eight
    C, I0, J1, e, diag, _ = interpolate([1, 1, 1], 2, 3)
    p3 = as_dict(C, I0)
    fresh = check_fresh(p3, [1, -2, 1, -2], 3)
    print(f'  COUNTER-CASE 3_1 polynomial vs 4_1 braid: fresh {fresh:.1e}  '
          f'{"✅ detected" if fresh > 1e-6 else "🔴 slips through"}')
    ok &= fresh > 1e-6

    # COUNTER-CASE 3: ONE coefficient rounded wrongly by 1. That is the error
    # this tool can make at all — the interpolation yields
    # floating point, rounding is to integers. If the gate does not catch THAT,
    # it catches nothing (Footgun 16).
    C, I0, J1, e, diag, _ = interpolate([1, -2, 1, -2], 3, 4)   # 4_1
    p4 = as_dict(C, I0)
    d_gut = det_exact(p4)
    # take the Alexander node values from the CORRECT polynomial, then falsify
    evals = {str(t): int(alexander_from_homfly(p4, t)) for t in (2, 3, -2)}
    hits = 0
    for key_name in list(p4)[:6]:
        broken = dict(p4)
        broken[key_name] = broken[key_name] + 1
        caught = (det_exact(broken) != d_gut)
        az, _txt = alexander_witness(broken, evals)
        caught = caught or (az is False)
        fr = check_fresh(broken, [1, -2, 1, -2], 3)
        caught = caught or (fr > 1e-8)
        hits += bool(caught)
    n_all = min(6, len(p4))
    print(f'  COUNTER-CASE coefficient +1: {hits}/{n_all} caught  '
          f'{"✅" if hits == n_all else "🔴 one slips through"}')
    ok &= (hits == n_all)

    # ── COUNTER-CASES for the MANTISSA LATCH ──────────────────────────────
    # A latch that rejects everything looks like one that works.
    # These six cases separate the two.
    import os as _os, re as _re
    rg = True

    # (1) the four MEASURED GREEN ones must NEVER be rejected
    for c_green in (48, 54, 60, 63):
        g = mantissa_latch(c_green, 1e-3, exakt=False)
        good = g is None
        rg &= good
        print(f'  LATCH c={c_green:>3} (measured green) : '
              f'{"lets through ✅" if good else "🔴 REJECTS — false alarm"}')

    # (2) the two aborted targets MUST be rejected
    for c_rot in (149, 151):
        g = mantissa_latch(c_rot, 1e-3, exakt=False)
        good = g is not None
        rg &= good
        print(f'  LATCH c={c_rot} (17 decades above): '
              f'{"rejects ✅" if good else "🔴 LETS THROUGH"}')

    # (3) THE LOAD-BEARING COUNTER-CASE: in --exakt there is no cancellation.
    #     If the latch fires there, it chokes exactly the path that can do it.
    g = mantissa_latch(151, 1e-3, exakt=True)
    good = g is None
    rg &= good
    print(f'  LATCH c=151 but --exakt        : '
          f'{"lets through ✅ (no floating point, no wall)" if good else "🔴 CHOKES THE EXACT PATH"}')

    # (4) the wall must MOVE with the bound, not be nailed down
    high = mantissa_latch(75, 1e-0, exakt=False)
    tight = mantissa_latch(75, 1e-6, exakt=False)
    good = (high is None) and (tight is not None)
    rg &= good
    print(f'  LATCH c=75 at bound 1e0 / 1e-6: '
          f'{"moves along ✅" if good else "🔴 hangs on a fixed number"}')

    # (5) against the REAL 40 runs: not a single green run may fall
    log = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..',
                        'tmp_durable', 'homfly_wand_20260823',
                        'interp_ziele_ABGEBROCHEN_20260824.log')
    if _os.path.exists(log):
        false_pos = false_neg = green = rot = 0
        for table_row in open(log):
            m = _re.match(r'^\[\d+/42\]\s*(✅|🔴)\s*\S+\s+n=\s*\d+ c=(\d+)', table_row)
            if not m:
                continue
            was_green = m.group(1) == '✅'
            g = mantissa_latch(int(m.group(2)), 1e-3, exakt=False)
            green += was_green; rot += (not was_green)
            if was_green and g:
                false_pos += 1
            if (not was_green) and g:
                false_neg += 1
        good = (false_pos == 0)
        rg &= good
        print(f'  LATCH against {green + rot} real runs   : '
              f'{false_pos} of {green} green ones wrongly rejected '
              f'{"✅" if good else "🔴"} · {false_neg} of {rot} red ones caught '
              f'in advance (= saved compute time)')
    else:
        print('  LATCH against real runs        : ⚠️  measurement log missing, case not run')

    # (6) COUNTER-CASE for the double-computation guard: it may only fire when
    #     the field is REALLY filled — otherwise it chokes every run.
    import tempfile as _tf
    def _mini(mit_homfly):
        d = {'pflichtliste': {'ALGEBRA': {'homfly_polynomial':
                              'x^2 - y^2' if mit_homfly else None},
                              'DIAGRAMM': {}}}
        p = _tf.NamedTemporaryFile('w', suffix='.json', delete=False)
        json.dump(d, p); p.close(); return p.name
    full, empty = _mini(True), _mini(False)
    t_full = {'name': 'X', 'record': full, 'n': 3, 'c': 10, 'zgrad': 10}
    t_empty = {'name': 'Y', 'record': empty, 'n': 3, 'c': 10, 'zgrad': 10}
    r1 = process_target(t_full)
    gut1 = r1.get('schon_da') is True
    try:
        r2 = process_target(t_empty)
        gut2 = not r2.get('schon_da')
    except Exception:
        gut2 = True          # starts computing and fails on the braid word — correct
    _os.unlink(full); _os.unlink(empty)
    rg &= (gut1 and gut2)
    print(f'  DOUBLE-COMPUTE filled/empty    : '
          f'{"skipped / computed ✅" if (gut1 and gut2) else "🔴 fires wrongly"}')

    ok &= rg
    print('  ⇒', '✅ the interpolation bites' if ok else '🔴 self-test failed')
    return 0 if ok else 1


# ── calibration against records whose HOMFLY regina has already delivered ─────────
def calibrate_against_records(entries, laut=True):
    """The decisive test: does the interpolation reproduce what regina
    computed independently — term by term? A tool that only
    LOOKS convincing is decoration (Footgun 16)."""
    from homfly_eichung import parse_homfly
    ok = misses = 0
    rows = []
    for e in entries:
        rec = json.load(open(e['record']))
        D = rec['pflichtliste']['DIAGRAMM']
        w = as_list(D.get('braid_word'))
        n = D.get('braid_n_strands')
        expected = rec['pflichtliste']['ALGEBRA'].get('homfly_polynomial')
        if not w or not n or not expected:
            continue
        zg = int(e.get('size_simplified') or D.get('crossing_number_min') or 0)
        t0 = time.time()
        C, I0, J1, ee, diag, _ = interpolate(w, int(n), zg, kerne=CORES)
        p = as_dict(C, I0)
        ps = parse_homfly(expected)
        equal = ps is not None and all(
            p.get(k, 0) == ps.get(k, 0) for k in set(p) | set(ps))
        fresh = check_fresh(p, w, int(n))
        z = {'name': e['name'], 'n': int(n), 'c': e.get('c'), 'zgrad': zg,
             'gleich': bool(equal), 'frisch': fresh,
             'ausschuss': diag['ausschuss'], 'sek': round(time.time() - t0, 1),
             'unbekannte': diag['unbekannte'], 'durchgaenge': diag['durchgaenge']}
        rows.append(z)
        if equal:
            ok += 1
        else:
            misses += 1
        if laut:
            print(f"  {'✅' if equal else '🔴'} {e['name'][:32]:<32} n={int(n):>2} "
                  f"c={e.get('c')} · residue {diag['ausschuss']:.1e} · "
                  f"fresh {fresh:.1e} · {z['sek']}s "
                  f"({diag['durchgaenge']} passes)", flush=True)
    return ok, misses, rows



# ── targets: the records where regina stands at the memory wall ───────────
def _path_index(liste):
    """name → record file. Needed because older wall runs wrote the path only
    on success — a missing field is not a missing record."""
    idx = {}
    if not liste or not os.path.exists(liste):
        return idx
    for x in json.load(open(liste)):
        p = os.path.join(site_config.ROOT, x['pfad'])
        for cand_k in (f"{x['name']}_pflicht_invarianten_complete.json",
                     f"{x['name']}_pflicht_complete.json"):
            f = os.path.join(p, cand_k)
            if os.path.exists(f):
                idx[x['name']] = f
                break
    return idx


def target_entries(wall_file, n_max=12, liste=None):
    d = json.load(open(wall_file))
    idx = _path_index(liste)
    out = []
    for e in d:
        if e.get('status') not in ('SPEICHERWAND', 'RIEGEL', 'FEHLER'):
            continue
        rec_path = e.get('record') or idx.get(e['name'])
        if not rec_path or not os.path.exists(rec_path):
            continue
        rec = json.load(open(rec_path))
        D = rec['pflichtliste']['DIAGRAMM']
        w = as_list(D.get('braid_word'))
        n = D.get('braid_n_strands')
        if not w or not n or int(n) > n_max:
            continue
        if rec['pflichtliste']['ALGEBRA'].get('homfly_polynomial'):
            continue
        zg = int(D.get('crossing_number_min') or e.get('c') or 0)
        if zg <= 0:
            continue
        out.append({'name': e['name'], 'record': rec_path, 'w': w,
                    'n': int(n), 'zgrad': zg, 'c': e.get('c')})
    out.sort(key=lambda z: (z['n'], z['zgrad']))
    return out


# ============================================================================
# MANTISSA LATCH — measured on 24.08.2026 on 40 runs (interp_ziele.log).
#
# The residue ("Ausschuss") is NOT a property of the knot but of the floating-point number:
#     residue ~ eps * sum|coefficients|,   eps = 2^-52 = 2.22e-16
# and sum|coefficients| grows exponentially in the crossing number.
# Least squares over 40 runs (c = 48..146, over 20 decades):
#     log10(residue) = 0.1968 * c - 15.613     residuals sd = 0.25 decades
# The INTERCEPT is the receipt for the cause: 10^-15.613 = 2.4e-16
# is the mantissa itself — the model finds it without having been told.
# The strand count n explains -0.019 decades/strand of the residuals,
# i.e. nothing: n drives the COST (n!), c drives the PRECISION. Separate.
#
# Measured boundary: green at c = 48/54/60/63, first red at c = 65.
# The prediction crosses 1e-3 at c = 64.1 — right in this gap.
#
# The latch only fires where even the most favourable end of the scatter band
# (-3 sd) still lies above the bound: for 1e-3 that is c >= 68. Between
# 64 and 67 we keep COMPUTING, not guessing — there the measurement
# is not unambiguous, and a run there costs seconds.
#
# It is NOT a stop before the run but a computational result: in floating point
# the result there IS not representable. The path that can do it is named in the
# label: --exakt (F_p + CRT, homfly_modular) has no cancellation, and
# for it the latch is switched off.
# The reason text returned by mantissen_riegel() is data (stored in 'gruende') and stays German.
# ============================================================================
MANTISSA_A = 0.1968          # decades of residue per crossing
MANTISSA_B = -15.613         # intercept = log10(machine epsilon)
MANTISSA_SD = 0.252          # scatter of the residuals, decades
MANTISSA_N = 40              # runs this rests on


def residue_prediction(zgrad):
    """Expected residue in floating point at crossing number zgrad."""
    return 10.0 ** (MANTISSA_A * zgrad + MANTISSA_B)


def mantissa_latch(zgrad, ausschuss_max, exakt=False):
    """Returns a REASON if floating point no longer carries this crossing
    number — otherwise None. In --exakt there is no cancellation, hence
    no latch either (the counter-case without which this latch would choke
    exactly the path that works)."""
    if exakt or not zgrad:
        return None
    import math
    c_wand = (math.log10(ausschuss_max) - MANTISSA_B + 3 * MANTISSA_SD) / MANTISSA_A
    if zgrad < c_wand:
        return None
    expected_residue = residue_prediction(zgrad)
    decades = math.log10(expected_residue / ausschuss_max)
    return (f'mantissa latch: c={zgrad} > {c_wand:.0f}; expected residue '
            f'~1e{math.log10(expected_residue):+.0f}, that is {decades:.0f} decades above '
            f'{ausschuss_max:.0e} (measured over {MANTISSA_N} runs). '
            f'Floating point cannot represent this — way out: --exakt')


def process_target(t, ausschuss_max=1e-3, fresh_max=1e-8, exakt=False,
                   no_latch=False):
    # ausschuss_max deliberately WITHIN the measurement: term-by-term equality
    # with regina is receipted up to 2.1e-3 (34 calibration records, 0 off).
    # Going higher would be extrapolation and must be a deliberate act.
    """Interpolate and GATE. Four hurdles, all must hold:
       1. residue — the coefficients must come out integral
       2. fresh points that were not on the grid
       3. det = |P(1,2i)| exactly against the determinant in the record
       4. |Alexander(-1)| from the stored Alexander polynomial
    If one breaks, NOTHING is written and the reason is given.
    The reason texts in 'gruende' and 'frisch_text' are data (returned/stored) and stay German."""
    from homfly_gegenprobe import record_det, record_alexander, det_from_alexander
    t0 = time.time()
    # The target list is frozen at START. If a second path runs alongside
    # (regina, --exakt), it fills records in the meantime — and this
    # run computes them a second time. Happened exactly so on 23.08.: the four
    # green results of the floating-point run were repetitions of
    # records that the exact path had already written at 12:00.
    # Costs one file access, saves hours.
    try:
        _r = json.load(open(t['record']))
        if _r['pflichtliste']['ALGEBRA'].get('homfly_polynomial'):
            return {'name': t['name'], 'record': t['record'], 'n': t['n'],
                    'c': t.get('c'), 'zgrad': t.get('zgrad'), 'ok': False,
                    'gruende': ['record meanwhile carries a HOMFLY (another route '
                                'was faster) — not computed twice'],
                    'schon_da': True, 'riegel': True,
                    'det': None, 'det_record': None, 'det_alexander': None,
                    'alex_punkte': None, 'alex_punkte_text': '',
                    'frisch': None, 'frisch_text': 'not computed',
                    'exakt': bool(exakt), 'sumabs': None, 'zeugen': 0,
                    'sek': round(time.time() - t0, 2), 'homfly': None,
                    'ausschuss': None, 'unbekannte': None, 'durchgaenge': None}
    except Exception:
        pass
    if not no_latch:
        reason = mantissa_latch(t.get('zgrad'), ausschuss_max, exakt)
        if reason:
            return {'name': t['name'], 'record': t['record'], 'n': t['n'],
                    'c': t.get('c'), 'zgrad': t.get('zgrad'), 'ok': False,
                    'gruende': [reason], 'riegel': True,
                    'det': None, 'det_record': None, 'det_alexander': None,
                    'alex_punkte': None, 'alex_punkte_text': '',
                    'frisch': None, 'frisch_text': 'not computed',
                    'exakt': bool(exakt), 'sumabs': None, 'zeugen': 0,
                    'sek': round(time.time() - t0, 2), 'homfly': None,
                    'ausschuss': residue_prediction(t.get('zgrad') or 0),
                    'unbekannte': None, 'durchgaenge': None}
    if exakt:
        import homfly_modular as _hm
        from homfly_modular import interpolate_exact, check_fresh_mod
        # pass --kerne through into the EXACT path. Without this line the
        # expensive path stayed serial while the floating-point path ran in parallel — exactly the
        # half that matters (25.08.).
        _hm.CORES = CORES
        p, I0, J1, e, md = interpolate_exact(t['w'], t['n'], t['zgrad'], laut=t.get('laut', False))
        diag = {'ausschuss': 0.0, 'unbekannte': len(p),
                'durchgaenge': md['primzahlen'], 'stabil': md['stabil']}
        tr, total_pts = check_fresh_mod(p, t['w'], t['n'])
        fresh = 0.0 if tr == total_pts else 1.0
        fresh_text = f'{tr}/{total_pts} fresh points exactly equal in F_p'
        if not md['stabil']:
            pass
    else:
        C, I0, J1, e, diag, _ = interpolate(t['w'], t['n'], t['zgrad'])
        p = as_dict(C, I0)
        fresh = check_fresh(p, t['w'], t['n'])
        fresh_text = f'{fresh:.1e} (floating point)'
    rec = json.load(open(t['record']))
    d_rec = record_det(rec)
    d_int = det_exact(p)
    alex = record_alexander(rec)
    d_alex = det_from_alexander(alex) if alex else None
    sumabs = sum(abs(v) for v in p.values())
    reasons = []
    if exakt and not diag.get('stabil'):
        reasons.append('CRT did not become stable (more primes needed)')
    if not exakt and diag['ausschuss'] > ausschuss_max:
        reasons.append(f"residue {diag['ausschuss']:.1e} > {ausschuss_max:.0e}")
    if (fresh > fresh_max) if not exakt else (fresh != 0.0):
        reasons.append(f"fresh points: {fresh_text}")
    if d_rec is None:
        reasons.append('no determinant in the record')
    elif d_int != d_rec:
        reasons.append(f'det {d_int} != record {d_rec}')
    if d_alex is not None and d_rec is not None and d_alex != d_rec:
        reasons.append(f'Alexander witness {d_alex} != record {d_rec}')
    # witness 4: Alexander at several points, exact.
    evals = rec['pflichtliste']['ALGEBRA'].get('alexander_point_evals')
    az, az_text = alexander_witness(p, evals)
    if az is False:
        reasons.append(f'Alexander points: {az_text}')
    witnesses = 2 + (1 if (d_alex is not None and d_alex == d_rec) else 0) \
               + (1 if az else 0)
    return {'name': t['name'], 'record': t['record'], 'n': t['n'],
            'c': t['c'], 'zgrad': t['zgrad'], 'ok': not reasons,
            'gruende': reasons, 'det': d_int, 'det_record': d_rec,
            'det_alexander': d_alex, 'alex_punkte': az, 'alex_punkte_text': az_text,
            'frisch': fresh, 'frisch_text': fresh_text,
            'exakt': bool(exakt), 'sumabs': sumabs, 'zeugen': witnesses,
            'sek': round(time.time() - t0, 1), 'homfly': as_regina_string(p),
            **{k: diag[k] for k in ('ausschuss', 'unbekannte', 'durchgaenge')}}


def write_record(r, datum):
    f = r['record']
    rec = json.load(open(f))
    rec['pflichtliste']['ALGEBRA']['homfly_polynomial'] = r['homfly']
    rec['homfly_herkunft'] = {
        'datum': datum,
        'weg': ('Hecke trace at a point (hecke_homfly_fast) on a grid, '
                'then 2-D interpolation to the full polynomial '
                '(tools/homfly_interpolation.py)'),
        'warum_nicht_regina': ('regina stops at this diagram at the '
                               'MEMORY WALL (~3 GB per process, partly '
                               'killed by earlyoom, Disc #677) — not '
                               'at the time.'),
        'konvention': 'regina AZ, variables printed as x=a, y=z',
        'box_gemessen': ('-e-n+1 <= a exponent <= -e+n-1, both axes even, '
                         'z degree <= c-1; checked on 178 freshly computed polynomials, '
                         '0 violations, 173 touch the boundary exactly'),
        'eichung': ('20 records with an existing regina HOMFLY reproduced '
                    'term by term (c=21..54), 0 off'),
        'ausschuss': r['ausschuss'],
        'ausschuss_deutung': ('largest distance of a coefficient from the '
                              'nearest integer; measured law '
                              'residue ~ 1.2e-16 * sum|coefficients|'),
        'zeuge_1': f"fresh points outside the grid: {r['frisch']:.2e}",
        'zeuge_2': f"det from HOMFLY |P(1,2i)| = {r['det']} = record {r['det_record']}",
        'zeuge_3': (f"|Alexander(-1)| = {r['det_alexander']}"
                    if r['det_alexander'] is not None else None),
        'zeuge_4': (f"Alexander at several points (exact, fractions): "
                    f"{r['alex_punkte_text']}" if r.get('alex_punkte') else None),
        'zeugen_gesamt': r['zeugen'],
        'unbekannte': r['unbekannte'], 'gitter_durchgaenge': r['durchgaenge'],
        'sekunden': r['sek'],
    }
    tmp = f + '.tmp'
    json.dump(rec, open(tmp, 'w'), indent=1)
    os.replace(tmp, f)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--eichen', help='JSON list (from homfly_wand.py) with record/name')
    ap.add_argument('--n-max', type=int, default=8)
    ap.add_argument('--limit', type=int, default=12)
    ap.add_argument('--out')
    ap.add_argument('--liste', help='n_gross.json, to look up record paths')
    ap.add_argument('--ziele', help='wand_alle.json — takes the SPEICHERWAND (memory wall) cases')
    ap.add_argument('--ziel-n-max', type=int, default=12)
    ap.add_argument('--ziel-limit', type=int)
    ap.add_argument('--exakt', action='store_true',
                    help='Hecke trace over F_p + CRT instead of floating point')
    ap.add_argument('--ausschuss-max', type=float, default=1e-3)
    ap.add_argument('--kein-mantissen-riegel', action='store_true',
                    help='compute high c despite the measured floating-point wall')
    ap.add_argument('--datum', default='2026-08-23')
    ap.add_argument('--schreiben', action='store_true')
    ap.add_argument('--selbsttest', action='store_true')
    ap.add_argument('--record', help='path to a *_complete.json')
    ap.add_argument('--laut', action='store_true')
    ap.add_argument('--kerne', type=int, default=None,
                    help='workers for the z loop (default: all cores; 1 = serial)')
    a = ap.parse_args()
    globals()['CORES'] = a.kerne     # also takes effect in functions without their own parameter
    if a.selbsttest:
        sys.exit(selbsttest())
    if a.eichen:
        d = json.load(open(a.eichen))
        d = [x for x in d if x.get('status') == 'OK' and x.get('record')]
        # draw over the whole range, not only the cheap ones
        d = [x for x in d if (x.get('n') or 99) <= a.n_max]
        d.sort(key=lambda z: (z.get('c') or 0))
        step = max(1, len(d) // a.limit)
        d = d[::step][:a.limit]
        print(f'Calibration on {len(d)} records with an existing regina HOMFLY')
        ok, misses, rows = calibrate_against_records(d)
        print(f'\n{ok} equal term by term, {misses} off')
        if a.out:
            json.dump(rows, open(a.out, 'w'), indent=1)
        sys.exit(0 if misses == 0 else 1)
    if a.ziele:
        z = target_entries(a.ziele, a.ziel_n_max, a.liste)
        if a.ziel_limit:
            z = z[:a.ziel_limit]
        print(f'{len(z)} targets (memory wall in regina, n <= {a.ziel_n_max})',
              flush=True)
        erg, good = [], 0
        for i, t in enumerate(z, 1):
            r = process_target(t, a.ausschuss_max, exakt=a.exakt,
                               no_latch=a.kein_mantissen_riegel)
            erg.append(r)
            mark = '✅' if r['ok'] else ('⛔' if r.get('riegel') else '🔴')
            print(f"[{i}/{len(z)}] {mark} {t['name'][:30]:<30} n={t['n']:>2} "
                  f"c={t['c']} · {'REJECTED' if r.get('riegel') else ('exact/' + str(r['durchgaenge']) + ' primes' if r.get('exakt') else 'residue %.1e' % r['ausschuss'])} · "
                  f"fresh {r['frisch_text']} · det {r['det']}"
                  f"{'' if r['ok'] else '  ⟵ ' + '; '.join(r['gruende'])}"
                  f" · {r['sek']}s", flush=True)
            if r['ok']:
                good += 1
                if a.schreiben:
                    write_record(r, a.datum)
            if a.out:
                json.dump(erg, open(a.out, 'w'), indent=1)
        print(f"\n{good}/{len(z)} through all four hurdles"
              f"{' — WRITTEN' if a.schreiben else ' (dry run)'}")
        sys.exit(0)
    if a.record:
        rec = json.load(open(a.record))
        D = rec['pflichtliste']['DIAGRAMM']
        w = as_list(D.get('braid_word')); n = int(D.get('braid_n_strands'))
        zg = int(D.get('crossing_number_min') or 0) + 4
        C, I0, J1, e, diag, _ = interpolate(w, n, zg, laut=a.laut, kerne=a.kerne)
        p = as_dict(C, I0)
        print(json.dumps({'det': det_exact(p),
                          'frisch': check_fresh(p, w, n), **diag}, indent=1))
        print(as_regina_string(p))
