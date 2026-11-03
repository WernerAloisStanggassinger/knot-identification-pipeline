#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Dunfield–Friedl–Jackson Conjecture 1.7: the comparison of the hyperbolic torsion polynomial T_K with genus and fiberedness.
# Conjecture (DFJ arXiv:1108.3045v3 p.5; precise form Agol-Dunfield 1501.02136v3 Conj 1.1):
#   deg T_K = 4 g(K) - 2   (deg = span of the symmetric Laurent polynomial = 2 x(K), x = 2g-1)
#   K fibered  <=>  T_K monic, i.e. leading coefficient == 1 in C (Kitayama 2108.08645v2 Conj 9.3). Proven is only
#   fibered ==> monic (and deg T_K = 4g - 2), DFJ Thm 1.5 (1108.3045v3 p.4); the converse is part of Conjecture 1.7.
# T_K from SnapPy Manifold.hyperbolic_torsion (nsagetools.compute_torsion, DFJ conventions, wada_conventions=False).
# SnapPy cleans coefficients with |c| < 2^(-bits/3) to zero; returned polynomial starts at exponent 0.
import time, resource
import sage.all  # noqa: F401  (hyperbolic_torsion needs Sage)
import snappy, spherogram

TIGHT_EXP = 100   # |x| < 2^-100  => numerically zero  (at >= 500 bits)
LOOSE_EXP = 30    # |x| > 2^-30   => numerically non-zero; in between => 'uncertain' -> recompute at higher precision


def classify(x):
    """x: non-negative real (sage RR/float). Returns 'zero', 'nonzero', or 'uncertain'."""
    x = float(x) if x < 1 else 1.0
    if x < 2.0 ** -TIGHT_EXP:
        return 'zero'
    if x > 2.0 ** -LOOSE_EXP:
        return 'nonzero'
    return 'uncertain'


def manifold_from_pd(pd):
    return spherogram.Link(pd).exterior()


def torsion_record(M, bits, tries=6):
    """Compute T_K with a geometric triangulation; retries with randomize()."""
    last = None
    for k in range(tries):
        st = M.solution_type()
        if st != 'all tetrahedra positively oriented':
            M.randomize()
            last = 'solution_type ' + st
            continue
        try:
            T = M.hyperbolic_torsion(bits_prec=bits)
            return T, k, st
        except Exception as e:  # noqa: BLE001
            last = repr(e)[:200]
            M.randomize()
    raise RuntimeError('torsion failed after %d tries: %s' % (tries, last))


def analyse(T):
    coeffs = T.list()  # dense, exponent 0..deg
    deg = T.degree()
    lc = T.leading_coefficient()
    low = coeffs[0]
    d_lc1 = (lc - 1).abs()
    d_abs1 = abs(lc.abs() - 1)
    max_im = max(c.imag().abs() for c in coeffs)
    max_c = max(c.abs() for c in coeffs)
    return dict(deg=int(deg), lc_re=float(lc.real()), lc_im=float(lc.imag()),
                lc_abs=float(lc.abs()), low_abs=float(low.abs()),
                monic=classify(d_lc1), abs_lc_is_1=classify(d_abs1),
                real=classify(max_im / max(1, max_c)),
                d_lc1_log2=float(d_lc1.log(2)) if d_lc1 > 0 else -1e9,
                max_im_rel_log2=float((max_im / max(1, max_c)).log(2)) if max_im > 0 else -1e9)


def verdict(a, g, fibered):
    """Compare analysis a with HFK genus g / fibered flag. Returns (status, notes)."""
    notes = []
    deg_ok = a['deg'] == 4 * g - 2
    if not deg_ok:
        notes.append('deg %d != 4g-2 = %d' % (a['deg'], 4 * g - 2))
    if a['monic'] == 'uncertain':
        notes.append('monic uncertain')
        mon_ok = None
    else:
        is_monic = a['monic'] == 'zero'
        mon_ok = (is_monic == bool(fibered))
        if not mon_ok:
            notes.append('monic=%s fibered=%s' % (is_monic, fibered))
    if deg_ok and mon_ok:
        return 'agree', notes
    if mon_ok is None and deg_ok:
        return 'uncertain', notes
    return 'MISMATCH', notes


def check(M, g, fibered, bits=500, escalate=(1000, 2000)):
    t0 = time.time()
    T, ntries, st = torsion_record(M, bits)
    a = analyse(T)
    status, notes = verdict(a, g, fibered)
    used = bits
    history = []
    # escalate on anything that is not a clean agreement, or where a classification was uncertain
    for b in escalate:
        if status == 'agree' and a['real'] != 'uncertain':
            break
        history.append(dict(bits=used, deg=a['deg'], status=status, notes=notes))
        M.randomize()
        T, ntries2, st = torsion_record(M, b)
        a = analyse(T)
        status, notes = verdict(a, g, fibered)
        used = b
    rec = dict(a, g=g, fibered=bool(fibered), four_g_minus_2=4 * g - 2, status=status, notes=notes,
               bits=used, retries=ntries, sec=round(time.time() - t0, 3))
    if history:
        rec['history'] = history
    return rec


def set_memlimit(gb=6):
    resource.setrlimit(resource.RLIMIT_AS, (gb * 2 ** 30, gb * 2 ** 30))
