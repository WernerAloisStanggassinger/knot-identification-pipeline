#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""homfly_modular.py — the same interpolation, but EXACT: Hecke trace over F_p,
Vandermonde over F_p, integers via CRT.

WHY (measured 2026-08-23): the floating-point version (homfly_interpolation.py)
reproduces 34 regina polynomials term by term, but on the target knots it breaks down
already at c~63 — and not through sloppiness: the twisted torus knots
have TINY determinants (11..19) with huge HOMFLY coefficients. That is
maximal cancellation, the worst conceivable case for double precision.
Measured, the error grew from 1.5e-05 (c=54) to 2.2e+04 (c=101).

In F_p there is no cancellation and no condition number. The price is a
factor in run time (one pass per prime), not an approximation.

STOPPING CRITERION INSTEAD OF AN ESTIMATE: one more prime is taken
until the CRT result no longer changes. Two consecutive
identical results = done. Whoever guesses the number of primes instead
also guesses the size of the coefficients — and that is exactly the unknown.

The matrix core stays scipy-sparse: every Hecke generator has at most TWO
entries per row, so (G @ M) with p < 2^31 fits safely in int64
(2 * (2^31)^2 < 2^63) and is reduced only afterwards.
"""
import sys, os, json, time, argparse, random
import site_config
import numpy as np
import scipy.sparse as sp

HECKE = (site_config.HECKE_DIR)
sys.path.insert(0, HECKE)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hecke_homfly_fast import syt, partitions, shape_cells, content, hc

# Primes just below 2^31 — the product of two residues stays in int64
PRIMZAHLEN = [2147483647, 2147483629, 2147483587, 2147483579, 2147483563,
              2147483549, 2147483543, 2147483497, 2147483489, 2147483477,
              2147483423, 2147483399, 2147483353, 2147483323, 2147483269,
              2147483249, 2147483237, 2147483179, 2147483171, 2147483137]


def inv(x, p):
    x %= p
    if x == 0:
        raise ZeroDivisionError('not invertible in F_p')
    return pow(x, p - 2, p)


def qint_mod(d, q, qi, p):
    """[d]_q = (q^d - q^-d)/(q - q^-1)."""
    zaehler = (pow(q, d, p) - pow(qi, d, p)) % p if d >= 0 else \
              (pow(qi, -d, p) - pow(q, -d, p)) % p
    return zaehler * inv((q - qi) % p, p) % p


def gen_mod(T, idx, i, q, qi, p):
    """Hecke generator sigma_i in the q-seminormal form, over F_p.
    Same as hecke_gen_sparse word for word, only with modular inverses."""
    n = len(T)
    rows, cols, vals = [], [], []
    done = set()
    sc = (q - qi) % p
    for k, t in enumerate(T):
        if k in done:
            continue
        d = content(t[i]) - content(t[i - 1])
        if d == 1:
            rows.append(k); cols.append(k); vals.append(q % p); done.add(k); continue
        if d == -1:
            rows.append(k); cols.append(k); vals.append((-qi) % p); done.add(k); continue
        ts = list(t); ts[i - 1], ts[i] = ts[i], ts[i - 1]
        ks = idx[tuple(ts)]
        qm2d = pow(q, -2 * d, p) if d < 0 else inv(pow(q, 2 * d, p), p)
        qp2d = pow(q, 2 * d, p) if d > 0 else inv(pow(q, -2 * d, p), p)
        a11 = sc * inv((1 - qm2d) % p, p) % p
        a22 = sc * inv((1 - qp2d) % p, p) % p
        qd = qint_mod(d, q, qi, p)
        a21 = (1 - inv(qd * qd % p, p)) % p
        rows += [k, ks, k, ks]
        cols += [k, ks, ks, k]
        vals += [a11, a22, 1, a21]
        done.add(k); done.add(ks)
    return sp.csr_matrix((np.array(vals, dtype=np.int64),
                          (np.array(rows), np.array(cols))),
                         shape=(n, n), dtype=np.int64)


def W_mod(lam, a, ai, q, qi, p):
    P = 1
    for cont, h in hc(lam):
        if cont >= 0:
            qc, qmc = pow(q, cont, p), pow(qi, cont, p)
        else:
            qc, qmc = pow(qi, -cont, p), pow(q, -cont, p)
        zaehler = (a * qc - ai * qmc) % p
        nenner = (pow(q, h, p) - pow(qi, h, p)) % p
        P = P * zaehler % p * inv(nenner, p) % p
    return P


def chis_mod(word, n, q, qi, p):
    """One pass per q: trace of all irreps. The expensive part depends ONLY on q."""
    out = []
    for lam in partitions(n):
        T = syt(lam)
        idx = {tuple(t): k for k, t in enumerate(T)}
        dim = len(T)
        sc = (q - qi) % p
        Gp, Gi = {}, {}
        M = np.eye(dim, dtype=np.int64)
        for g in word:
            i = abs(g)
            if i not in Gp:
                Gp[i] = gen_mod(T, idx, i, q, qi, p)
                Gi[i] = (Gp[i] - sc * sp.identity(dim, dtype=np.int64,
                                                  format='csr'))
                Gi[i].data %= p
                Gi[i] = Gi[i].tocsr()
            G = Gi[i] if g < 0 else Gp[i]
            M = (G @ M) % p          # at most 2 nnz per row ⟹ no overflow
        out.append((lam, int(np.trace(M) % p)))
    return out


def P_mod(word, n, a, q, p, vor=None):
    """P_AZ(a, z=q-1/q) in F_p, including the normalisation K=(q-1/q)/(a-1/a)."""
    e = sum(1 if g > 0 else -1 for g in word)
    ai, qi = inv(a, p), inv(q, p)
    cs = vor if vor is not None else chis_mod(word, n, q, qi, p)
    s = 0
    for lam, ch in cs:
        s = (s + W_mod(lam, a, ai, q, qi, p) * ch) % p
    K = (q - qi) % p * inv((a - ai) % p, p) % p
    ae = pow(ai, e, p) if e >= 0 else pow(a, -e, p)
    return s * K % p * ae % p


# ── Linear algebra over F_p (sizes <= ~100, Gauss is enough) ─────────────
def loese_mod(A, b, p):
    """A x = b over F_p, A square. b may have several columns."""
    m = len(A)
    M = [list(map(lambda v: v % p, A[i])) + list(map(lambda v: v % p, b[i]))
         for i in range(m)]
    br = len(b[0])
    for c in range(m):
        piv = next((r for r in range(c, m) if M[r][c] % p), None)
        if piv is None:
            raise ZeroDivisionError('Vandermonde singular in F_p')
        M[c], M[piv] = M[piv], M[c]
        iv = inv(M[c][c], p)
        M[c] = [v * iv % p for v in M[c]]
        for r in range(m):
            if r != c and M[r][c]:
                f = M[r][c]
                M[r] = [(M[r][j] - f * M[c][j]) % p for j in range(m + br)]
    return [row[m:] for row in M]


def crt(reste, primzahlen):
    """Chinese remainder theorem, result as the SYMMETRIC representative — otherwise
    every negative coefficient would come out as a huge positive number."""
    M = 1
    for p in primzahlen:
        M *= p
    x = 0
    for r, p in zip(reste, primzahlen):
        Mi = M // p
        x = (x + r * Mi * pow(Mi, -1, p)) % M
    return x - M if x > M // 2 else x


# ── Interpolation over F_p, then CRT until stable ──────────────────
def _box(word, n, z_grad_max):
    e = sum(1 if g > 0 else -1 for g in word)
    a_lo, a_hi = -e - n + 1, -e + n - 1
    return (a_lo - 1) // 2, (a_hi + 1) // 2, (z_grad_max + 2) // 2, e


KERNE = None    # global; set from --kerne, otherwise EICH_KERNE (see below)
# ---------------------------------------------------------------------------
# CALIBRATED DEFAULT FOR THE WORKER COUNT — measured 2026-08-25, not guessed.
#
# This used to be os.cpu_count() (= 24 on this box), from Hard Rule 6
# ("full throttle, no ÷2"). But the rule means: the BOX must not idle — it
# does not mean: a single job gets all cores. Measured on T_13_7_2_6
# (n=11, real braid word, 27 q sample points, kern_eichung2.py):
#
#     cores   wall s     CPU s   speedup  efficiency
#         1    854.1     853.7     1.00x      1.00
#         6    330.9    1744.8     2.58x      0.49   <-- optimum
#        12    388.3    3997.0     2.20x      0.21   slower than 6
#        18    414.4    5524.4     2.06x      0.15   slower still
#     (24 never measured cleanly — the run was spoiled by a second computation
#      on the same box; the row is left out deliberately instead of guessed.)
#
# From ~6 on the memory bus is saturated: dense int64 matrices are
# bandwidth-bound, not compute-bound, and the 12 physical cores (Ryzen 9
# 3900X, SMT) share one bus. Every further worker just waits along and
# costs full CPU time.
#
# PROPER LOADING therefore means: SEVERAL KNOTS side by side, each with
# EICH_KERNE workers (4 x 6 = 24). Parallelism belongs one level
# higher, not inside the single knot.
#
# Override: --kerne <n> always wins; without the argument
# HOMFLY_KERNE from the environment applies; only then this value. On another
# box, re-calibrate (kern_eichung2.py), do not estimate.
# ---------------------------------------------------------------------------
EICH_KERNE = int(os.environ.get("HOMFLY_KERNE", 6))



def _zeile_mod(aufgabe):
    """One q sample point in F_p: returns (l, [P_mod(a_k) for all k]).

    PARALLELISED 2026-08-25. The exact path is the expensive one — the 30-hour run
    of the n=11 chain ran here, serially, on ONE of 24 cores. The sample points
    are drawn sequentially beforehand (rng stays untouched, the result depends
    on the draw, not on the order of evaluation); only the computation is
    distributed, and it is independent per l — chis_mod depends on q alone.
    """
    l, q, word, n, as_, p = aufgabe
    qi = inv(q, p)
    cs = chis_mod(word, n, q, qi, p)
    return l, [P_mod(word, n, a, q, p, vor=cs) for a in as_]


def eine_primzahl(word, n, z_grad_max, p, rng, laut=False):
    """Coefficients C[J][m] in F_p. A=a^2, Z=z^2 (both exponents even)."""
    I0, I1, J1, e = _box(word, n, z_grad_max)
    NA, NZ = I1 - I0 + 1, J1 + 1

    # Draw sample points: q determines Z, a determines A. Both must be PAIRWISE
    # distinct and must not kill any denominator — otherwise draw again.
    qs, Zs = [], []
    versuche = 0
    while len(qs) < NZ:
        versuche += 1
        if versuche > 400 * NZ:
            raise RuntimeError('no usable q sample points found')
        q = rng.randrange(2, p - 1)
        try:
            qi = inv(q, p)
            Z = (q - qi) % p
            Z = Z * Z % p
            if Z in Zs or Z == 0:
                continue
            for h in range(1, n + 1):          # denominators in W
                if (pow(q, h, p) - pow(qi, h, p)) % p == 0:
                    raise ZeroDivisionError
            for d in range(-n, n + 1):         # denominators in the generators
                if abs(d) > 1 and d != 0:
                    if (1 - pow(q, 2 * d, p) % p) % p == 0:
                        raise ZeroDivisionError
                    if qint_mod(d, q, qi, p) == 0:
                        raise ZeroDivisionError
        except (ZeroDivisionError, ValueError):
            continue
        qs.append(q); Zs.append(Z)

    as_, As = [], []
    while len(as_) < NA:
        a = rng.randrange(2, p - 1)
        try:
            ai = inv(a, p)
        except ZeroDivisionError:
            continue
        if (a - ai) % p == 0:                  # K would have denominator 0
            continue
        A = a * a % p
        if A in As or A == 0:
            continue
        as_.append(a); As.append(A)

    V = [[0] * NA for _ in range(NZ)]
    aufgaben = [(l, q, word, n, as_, p) for l, q in enumerate(qs)]
    k_soll = min(KERNE if KERNE else EICH_KERNE, NZ)

    def _seriell(grund=None):
        if laut and grund:
            print(f'      {grund} — serial', flush=True)
        for auf in aufgaben:
            l, werte = _zeile_mod(auf)
            V[l] = werte
            if laut and NZ >= 6 and (l + 1) % max(1, NZ // 4) == 0:
                print(f'      {l+1}/{NZ} q sample points', flush=True)

    _mod = sys.modules.get(_zeile_mod.__module__)
    _importierbar = _mod is not None and (
        _zeile_mod.__module__ == '__main__' or getattr(_mod, '__spec__', None) is not None)

    # SELF-CALIBRATION instead of a guessed threshold (measured 25.08.): a
    # spawn pool costs ~0.9 s to set up, and it is rebuilt FOR EACH PRIME.
    # For small knots parallel was therefore SLOWER than serial
    # (3_1: 0.01 s serial vs 0.88 s parallel). So: compute the first sample point
    # serially, estimate the remaining run time from it and only then
    # decide. No constant over n or NZ — that would be guessed.
    POOL_KOSTEN = 1.2          # s, set-up + tear-down of a spawn pool, generous
    if k_soll > 1 and _importierbar and NZ > 1:
        _t = time.time()
        _l0, _w0 = _zeile_mod(aufgaben[0])
        _je = time.time() - _t
        V[_l0] = _w0
        rest = aufgaben[1:]
        seriell_rest = _je * len(rest)
        if seriell_rest < POOL_KOSTEN * 2:
            if laut:
                print(f'      {_je*1000:.0f} ms per sample point, remaining {seriell_rest:.1f}s '
                      f'— pool not worth it', flush=True)
            k_soll = 1
            aufgaben = rest
        else:
            aufgaben = rest
            if laut:
                print(f'      {_je*1000:.0f} ms per sample point, remaining serial would be '
                      f'{seriell_rest:.0f}s — pool with {k_soll} cores', flush=True)

    if k_soll <= 1:
        _seriell()
    elif not _importierbar:
        _seriell(f"Modulname '{_zeile_mod.__module__}' nicht importierbar")
    else:
        import concurrent.futures as _cf
        import multiprocessing as _mp
        try:
            fertig = 0
            with _cf.ProcessPoolExecutor(max_workers=k_soll,
                                         mp_context=_mp.get_context('spawn')) as ex:
                for l, werte in ex.map(_zeile_mod, aufgaben, chunksize=1):
                    V[l] = werte
                    fertig += 1
                    if laut and NZ >= 6 and fertig % max(1, NZ // 4) == 0:
                        print(f'      {fertig}/{NZ} q sample points ({k_soll} cores)', flush=True)
        except Exception as ex_:
            _seriell(f'Pool gescheitert ({type(ex_).__name__})')

    MZ = [[pow(Z, J, p) for J in range(NZ)] for Z in Zs]
    X = loese_mod(MZ, V, p)                              # over J
    MA = [[pow(A, I0 + m, p) if I0 + m >= 0
           else inv(pow(A, -(I0 + m), p), p) for m in range(NA)] for A in As]
    XT = [[X[J][k] for J in range(NZ)] for k in range(NA)]
    CT = loese_mod(MA, XT, p)                            # over m
    C = [[CT[m][J] for m in range(NA)] for J in range(NZ)]
    return C, I0, J1, e


def interpoliere_exakt(word, n, z_grad_max, max_prim=12, seed=23, laut=False):
    """Add primes until the CRT result stays fixed."""
    rng = random.Random(seed)
    reste, benutzt = [], []
    vorher = None
    t0 = time.time()
    for pi, p in enumerate(PRIMZAHLEN[:max_prim]):
        C, I0, J1, e = eine_primzahl(word, n, z_grad_max, p, rng, laut)
        reste.append(C); benutzt.append(p)
        if len(benutzt) < 2:
            continue
        jetzt = {}
        NZ, NA = len(C), len(C[0])
        for J in range(NZ):
            for m in range(NA):
                v = crt([r[J][m] for r in reste], benutzt)
                if v:
                    jetzt[(2 * (I0 + m), 2 * J)] = v
        if laut:
            print(f'    {len(benutzt)} primes, {len(jetzt)} terms, '
                  f'{time.time()-t0:.0f}s', flush=True)
        if vorher is not None and jetzt == vorher:
            return jetzt, I0, J1, e, {'primzahlen': len(benutzt),
                                      'sek': round(time.time() - t0, 1),
                                      'stabil': True}
        vorher = jetzt
    return (vorher or {}), I0, J1, e, {'primzahlen': len(benutzt),
                                       'sek': round(time.time() - t0, 1),
                                       'stabil': False}


def selbsttest():
    """Two known knots exactly, plus three counter-cases."""
    from homfly_interpolation import det_exakt, als_regina_string
    ok = True
    for nm, w, n, zg, det_soll, soll_str in [
            ('Kleeblatt 3_1', [1, 1, 1], 2, 3, 3, 'x^-2 y^2 + 2 x^-2 - x^-4'),
            ('Achterknoten 4_1', [1, -2, 1, -2], 3, 4, 5, None)]:
        p, I0, J1, e, d = interpoliere_exakt(w, n, zg)
        det = det_exakt(p)
        gut = (det == det_soll and d['stabil'])
        print(f'  {nm:<18} det={det} (expected {det_soll}) · {d["primzahlen"]} '
              f'primes · stable={d["stabil"]} · {d["sek"]}s '
              f'{"✅" if gut else "🔴"}')
        print(f'     {als_regina_string(p)}')
        if soll_str:
            gleich = als_regina_string(p) == soll_str
            print(f'     against the floating-point version: '
                  f'{"✅ identical" if gleich else "🔴 different"}')
            gut &= gleich
        ok &= gut

    # COUNTER-CASE 1: box too short on the z axis ⟹ det must come out wrong
    p, I0, J1, e, d = interpoliere_exakt([1, 1, 1], 2, 0)
    # z_grad_max=0 ⟹ J1=1, that is still enough; so harder: check directly
    # that a DELIBERATELY wrong box does not give the same determinant.
    from homfly_interpolation import interpoliere, als_dict
    C2, I02, J12, e2, dg, _ = interpoliere([1, 1, 1], 2, 3, zwang_box=(2, 3, 1))
    det_falsch = det_exakt(als_dict(C2, I02))
    print(f'  COUNTER-CASE wrong box    : det={det_falsch} (must not be 3) '
          f'{"✅" if det_falsch != 3 else "🔴"}')
    ok &= (det_falsch != 3)

    # COUNTER-CASE 2: CRT with ONE prime must not count as "stable"
    _, _, _, _, d1 = interpoliere_exakt([1, 1, 1], 2, 3, max_prim=1)
    print(f'  COUNTER-CASE 1 prime      : stable={d1["stabil"]} '
          f'(must be False) {"✅" if not d1["stabil"] else "🔴"}')
    ok &= not d1['stabil']

    # COUNTER-CASE 3: the symmetric CRT representative must give negative numbers
    neg = crt([(-7) % p for p in PRIMZAHLEN[:3]], PRIMZAHLEN[:3])
    print(f'  COUNTER-CASE CRT of -7    : {neg} {"✅" if neg == -7 else "🔴"}')
    ok &= (neg == -7)

    print('  ⇒', '✅ the exact version bites' if ok else '🔴 self-test failed')
    return 0 if ok else 1


def als_json(word, n, z_grad_max, max_prim=12):
    """ONE line of JSON — so that the pipeline can call this as its own process with
    its own clock (Footgun 13: isolation belongs at the process level)."""
    from homfly_interpolation import det_exakt, als_regina_string
    t0 = time.time()
    p, I0, J1, e, md = interpoliere_exakt(word, n, z_grad_max, max_prim=max_prim)
    tr, ges = pruefe_frisch_mod(p, word, n)
    return {'homfly': als_regina_string(p), 'det': det_exakt(p),
            'primzahlen': md['primzahlen'], 'stabil': md['stabil'],
            'frisch_treffer': tr, 'frisch_gesamt': ges,
            'terme': len(p), 'sek': round(time.time() - t0, 1),
            'alexander_punkte': None}


def pruefe_frisch_mod(p_koef, word, n, anzahl=3, seed=777):
    """Fresh points, but EXACT: in F_q at points that were not on the grid.
    The floating-point test is no longer usable here — with coefficients around 1e12 and
    massive cancellation it measures its own noise. A comparison in F_p
    is either equal or not; there is nothing in between."""
    rng = random.Random(seed)
    treffer = 0
    for pr in PRIMZAHLEN[-anzahl:]:
        for _ in range(500):
            q = rng.randrange(2, pr - 1)
            a = rng.randrange(2, pr - 1)
            try:
                qi, ai = inv(q, pr), inv(a, pr)
                if (a - ai) % pr == 0:
                    continue
                for h in range(1, n + 1):
                    if (pow(q, h, pr) - pow(qi, h, pr)) % pr == 0:
                        raise ZeroDivisionError
                for d in range(-n, n + 1):
                    if abs(d) > 1 and d != 0:
                        if (1 - pow(q, 2 * d, pr)) % pr == 0:
                            raise ZeroDivisionError
                        if qint_mod(d, q, qi, pr) == 0:
                            raise ZeroDivisionError
                soll = P_mod(word, n, a, q, pr)
            except (ZeroDivisionError, ValueError):
                continue
            z = (q - qi) % pr
            ist = 0
            for (i, j), c in p_koef.items():
                ai_pot = pow(a, i, pr) if i >= 0 else pow(ai, -i, pr)
                ist = (ist + c * ai_pot % pr * pow(z, j, pr)) % pr
            treffer += int(ist % pr == soll % pr)
            break
    return treffer, anzahl


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--selbsttest', action='store_true')
    ap.add_argument('--wort', help='braid word as a JSON list')
    ap.add_argument('--n', type=int)
    ap.add_argument('--zgrad', type=int)
    ap.add_argument('--max-prim', type=int, default=12)
    ap.add_argument('--alex-evals', help='alexander_point_evals as JSON')
    a = ap.parse_args()
    if a.selbsttest:
        sys.exit(selbsttest())
    if a.wort:
        from homfly_interpolation import alexander_zeuge
        try:
            w = json.loads(a.wort)
            out = als_json([int(g) for g in w], a.n, a.zgrad, a.max_prim)
            if a.alex_evals:
                from homfly_eichung import parse_homfly
                pp = parse_homfly(out['homfly'])
                az, txt = alexander_zeuge(pp, json.loads(a.alex_evals))
                out['alexander_punkte'] = az
                out['alexander_punkte_text'] = txt
        except Exception as ex:
            out = {'error': repr(ex)[:400]}
        print(json.dumps(out))
        sys.exit(0)
    ap.print_help()
