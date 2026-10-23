#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
mandatory_invariants.py
======================
signed Gauss code → complete mandatory list ("Pflichtliste") of knot invariants

Multi-env orchestrator: calls regina (base), snappy/spherogram (snappy_env),
Sage (sage_stable) and KnotJob (java), each via subprocess, and collects
the results in one JSON report.

Input:  --signed-gauss "U1+O2+..." or --kinematic-result <result.json>
Output: JSON file with all mandatory invariants

Blind-test compliant: sealing is activated in every subprocess (see
pipeline_seal.py); no knot-table lookups.
"""

import argparse
import site_config
import ast, json
import math
import os
import re
import subprocess
import sys
import time
import textwrap

# Disc #222 (2026-05-21) / Disc #586 (2026-07-27): subprocess.run + conda run -> orphans on
# timeout. The kill only hits `conda`; the grandchildren (`python`, `sage`, `sage-eval`) keep running.
# pipeline_subprocess.run_capture_with_group_kill kills the whole process group.
# mandatory_invariants_extra.py and enrich_optional_invariants.py have used it since 21.05. —
# THIS file never got it, and exactly that caused the memory drain on 27.07.
# (MemAvailable 19894 -> 1236 MB, earlyoom, reboot). Approved by Werner 27.07. 22:2x.
from pipeline_subprocess import (
    run_capture_with_group_kill as _run_with_pgkill,
    install_signal_killer,
)

SEAL_PRELUDE = (
    "import sys\n"
    "sys.path.insert(0, '@FREEZE_DIR@')\n"
    "from pipeline_seal import seal_all\n"
    "seal_all(verbose=False)\n".replace('@FREEZE_DIR@', site_config.FREEZE_DIR)
)


def run_subprocess(env, code, timeout=600):
    """Run Python code in a conda env and parse the last JSON line of its output."""
    cmd = ['conda', 'run', '-n', env, '--no-capture-output', 'python', '-c', code]
    try:
        res = _run_with_pgkill(cmd, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {'error': f'timeout ({timeout}s)'}
    if res.returncode != 0:
        return {'error': f'rc={res.returncode}', 'stderr': res.stderr[-2000:]}
    for line in reversed(res.stdout.strip().splitlines()):
        line = line.strip()
        if line.startswith('{') and line.endswith('}'):
            try:
                return json.loads(line)
            except json.JSONDecodeError:
                continue
    return {'error': 'no json output', 'stdout': res.stdout[-2000:],
            'stderr': res.stderr[-2000:]}


def run_sage(code, timeout=600):
    """Run Sage code via inline `sage -c ...` (Discovery #58 — `sage script.py`
    triggers circular imports)."""
    cmd = ['conda', 'run', '-n', 'sage_stable', '--no-capture-output',
           'sage', '-c', code]
    try:
        res = _run_with_pgkill(cmd, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {'error': f'timeout ({timeout}s)'}
    if res.returncode != 0:
        return {'error': f'rc={res.returncode}', 'stderr': res.stderr[-2000:]}
    for line in reversed(res.stdout.strip().splitlines()):
        line = line.strip()
        if line.startswith('{') and line.endswith('}'):
            try:
                return json.loads(line)
            except json.JSONDecodeError:
                continue
    return {'error': 'no json output', 'stdout': res.stdout[-2000:],
            'stderr': res.stderr[-2000:]}


# ──────────────────────────────────────────────
# REGINA (base): classical algebra + diagram codes
# ──────────────────────────────────────────────

_REGINA_HELPER = '\n        import os as _os, time as _time, tempfile as _tf\n        _T_START = _time.time()\n        def _remaining_budget():\n            """How much may ONE expensive call still cost without the OUTER\n            limit firing and taking all the cheap results with it? 85 % of the outer\n            budget minus what has already been spent. Never below 10 s."""\n            return max(10, int(_OUTER * 0.85 - (_time.time() - _T_START)))\n        def _with_deadline(_fn, _tag, _sek=None):\n            """ONE expensive call in a throwaway child with its own clock.\n            Footgun 13: isolation at process level, not at exception level.\n            Footgun 18/24: waitpid ALWAYS with WNOHANG, never blocking.\n            Footgun 20: NEVER signal.alarm — by default SIGALRM terminates the\n                        watcher, not the watched.\n            File instead of pipe: an 84-crossing polynomial overflows the 64 KB\n            pipe buffer, the child blocks on write and the parent\n            waits for a child that is waiting for it."""\n            if _sek is None:\n                _sek = _remaining_budget()\n            _fd, _pf = _tf.mkstemp(prefix=\'reg_\' + _tag + \'_\')\n            _os.close(_fd)\n            _t0 = _time.time()\n            _pid = _os.fork()\n            if _pid == 0:\n                try:\n                    _v = str(_fn())\n                    with open(_pf, \'w\') as _f:\n                        _f.write(_v)\n                    _os._exit(0)\n                except BaseException as _e:\n                    try:\n                        with open(_pf, \'w\') as _f:\n                            _f.write(\'\\x00ERR \' + repr(_e)[:300])\n                    except BaseException:\n                        pass\n                    _os._exit(1)\n            while True:\n                _p, _st = _os.waitpid(_pid, _os.WNOHANG)\n                if _p:\n                    try:\n                        _s = open(_pf).read()\n                    except BaseException:\n                        _s = \'\\x00ERR result file not readable\'\n                    try:\n                        _os.unlink(_pf)\n                    except OSError:\n                        pass\n                    if _s.startswith(\'\\x00ERR\'):\n                        return None, _s[5:]\n                    return _s, None\n                if _time.time() - _t0 > _sek:\n                    _os.kill(_pid, 9)\n                    _os.waitpid(_pid, 0)\n                    try:\n                        _os.unlink(_pf)\n                    except OSError:\n                        pass\n                    return None, \'timeout (%ds single latch, outer budget %ds)\' % (_sek, _OUTER)\n                _time.sleep(0.25)\n'


def regina_full(signed_gauss, point_eval=False):
    # POINT-EVAL CAP (approved by Werner 2026-06-17, Disc #248): homfly/jones/bracket are
    # state sums 2^c -> for c>~30 timeout (1200s) / OOM (rc=137). point_eval=True omits
    # exactly these three; Alexander (polynomial time, Δ(-1)=det by point evaluation) +
    # Goeritz/Seifert carry determinant and identity. alexander/pd/gauss stay.
    # 2026-08-23 (Disc #673): every expensive polynomial gets its OWN deadline,
    # and the CHEAP results survive when an expensive one dies.
    # Measured on a c=84 diagram, each on its own:
    #   gauss/dt/pd/writhe/alternating/knotSig/simplify  0.00 s
    #   virtualGenus 0.5 s · jones 6 s · bracket 6 s · alexander 8 s
    #   homfly  > 120 s  ⟵ ONE call took twelve finished results down with it,
    #                       pd_code included. Root of the 11 stage-6 deaths (22.08.).
    # ORDER: the most expensive last. homfly used to stand FIRST and ate up the
    # budget before jones/bracket even got their turn.
    _exp_lines = [] if point_eval else [
        "L.improveTreewidth()",
        "_v, _e = _with_deadline(lambda: L.jones(regina.Algorithm.Treewidth), 'jones')",
        "out['jones'] = _v",
        "if _e: out['jones_err'] = _e",
        "_v, _e = _with_deadline(lambda: L.bracket(regina.Algorithm.Treewidth), 'bracket')",
        "out['kauffman'] = _v",
        "if _e: out['kauffman_err'] = _e",
        "_v, _e = _with_deadline(lambda: L.homfly(regina.Algorithm.Treewidth), 'homfly')",
        "out['homfly'] = _v",
        "if _e: out['homfly_err'] = _e",
    ]
    _expensive = "\n".join("            " + _l for _l in _exp_lines)
    _outer = (int(os.environ.get('REGINA_TIMEOUT_POINTEVAL', 300)) if point_eval
               else int(os.environ.get('REGINA_TIMEOUT', 1200)))
    code = SEAL_PRELUDE + textwrap.dedent(f"        _OUTER = {_outer}\n") \
         + textwrap.dedent(_REGINA_HELPER) + textwrap.dedent(f"""
        import json, regina
        sg = {signed_gauss!r}
        out = {{}}
        out['point_eval_mode'] = {point_eval!r}
        try:
            L = regina.Link.fromSignedGauss(sg)
            out['parsed'] = True
            out['classical'] = bool(L.isClassical())
            out['virtualGenus_input'] = int(L.virtualGenus())
            out['size_input']  = int(L.size())
            L.simplify()
            out['simplified_size'] = int(L.size())
            out['gauss_code_reduced'] = str(L.gauss())
            out['signed_gauss_reduced'] = str(L.signedGauss())
            try: out['oriented_gauss']  = str(L.orientedGauss())
            except Exception as e: out['oriented_gauss_err'] = str(e)
            try: out['dt_code'] = str(L.dt())
            except Exception as e: out['dt_err'] = str(e)
            try:
                pd_str = str(L.pd())
                out['pd_code'] = pd_str
            except Exception as e:
                out['pd_err'] = str(e)
            _v, _e = _with_deadline(lambda: L.alexander(), 'alexander')
            if _v is not None: out['alexander'] = _v
            if _e: out['alexander_err'] = _e
{_expensive}
            try: out['alternating'] = bool(L.isAlternating())
            except Exception as e: out['alt_err'] = str(e)
            try: out['knot_sig'] = str(L.knotSig())
            except Exception as e: out['knot_sig_err'] = str(e)
            try: out['writhe'] = int(L.writhe())
            except Exception as e: out['writhe_err'] = str(e)
        except Exception as e:
            out['parsed'] = False
            out['error']  = str(e)
        for _k in [_k for _k, _v in out.items() if _v is None]:
            del out[_k]          # missing key as before, not None
        print(json.dumps(out))
    """)
    # Disc #192: 180s too tight because of subprocess init overhead (conda run +
    # imports + seal_all ≈ 150s); for c≥80 (e.g. kf121) too little budget
    # remains. A direct test showed: fromSignedGauss itself takes 0.0s,
    # only the wrapper init takes time. 1200s is a generous factor.
    # 2026-06-04 (approved by Werner, targeted kf286 run): overridable via ENV
    # REGINA_TIMEOUT — the polyc07 balloon trio needs more than 1200s on raw ~150-crossing
    # diagrams (SKIP path). Default UNCHANGED 1200.
    # point_eval: nothing exponential runs -> 300s is generous (Alexander c=78 ~4.5s).
    _to = (int(os.environ.get('REGINA_TIMEOUT_POINTEVAL', 300)) if point_eval
           else int(os.environ.get('REGINA_TIMEOUT', 1200)))
    return run_subprocess('base', code, timeout=_to)


def alexander_at_minus_one(alex_str):
    """Extract |Alexander(-1)| = determinant from a Regina Alexander string
    such as '9 - 42 a + 67 a^2 - 42 a^3 + 9 a^4'."""
    if not alex_str:
        return None
    # Regina's Alexander format uses 'a' as the variable
    poly = alex_str.replace(' ', '')
    # Substitute a = -1
    # parse: extract coefficients
    terms = re.findall(r'([+-]?\d+)(?:a(?:\^(\d+))?)?', poly)
    val = 0
    for coef_s, exp_s in terms:
        if not coef_s or coef_s in ('+', '-'):
            continue
        try:
            coef = int(coef_s)
        except ValueError:
            continue
        exp = int(exp_s) if exp_s else (1 if 'a' in poly[poly.index(coef_s):
                                                          poly.index(coef_s)+len(coef_s)+2]
                                          else 0)
        # Simplified: a bare 'a' means exponent=1, otherwise 0
        # (for verification; Sage is our source of truth)
        val += coef * ((-1) ** exp)
    return abs(val)


# ──────────────────────────────────────────────
# POINT EVALUATION (Disc #248): value instead of symbolic form
# ──────────────────────────────────────────────

def _eval_regina_poly(poly_str, t, var='x'):
    """Evaluates a Regina polynomial (e.g. Alexander 'x^72 - x^71 + ... - x + 1')
    EXACTLY at t (Fraction). Disc #248 lever: polynomial-time value instead of the
    exponential state-sum form. Returns None if empty/unparsable."""
    if not poly_str:
        return None
    from fractions import Fraction
    try:
        s = poly_str.replace(' ', '').replace('^-', '^NEG')  # protect negative exponents
        s = s.replace('-', '+-')
        tv = Fraction(t)
        val = Fraction(0)
        for tok in s.split('+'):
            if not tok:
                continue
            neg = tok.startswith('-')
            if neg:
                tok = tok[1:]
            tok = tok.replace('^NEG', '^-')
            if not tok:
                continue
            if var in tok:
                cp, _, ep = tok.partition(var)
                c = Fraction(cp) if cp not in ('', '+') else Fraction(1)
                e = int(ep[1:]) if ep.startswith('^') else 1
            else:
                c = Fraction(tok)
                e = 0
            if neg:
                c = -c
            val += c * tv ** e
        return val
    except Exception:
        return None


def alex_det(alex_str):
    """Determinant = |Δ(-1)| by point evaluation (polynomial time, no state sum).
    Third det witness next to Goeritz (sage) and Seifert (snappy)."""
    v = _eval_regina_poly(alex_str, -1)
    if v is None:
        return None
    return abs(int(v)) if v.denominator == 1 else abs(float(v))


def point_alexander_det(pd_str, variant=6):
    """FALLBACK when regina cannot manage the algebra (approved by Werner 2026-08-22).

    WHY: regina's Alexander doubles every ~16 crossings — measured on our own
    corpus: c=96 17.3 s · c=114 36.5 s · c=128 65.3 s · c=144 124.4 s, and at
    c=192 it runs into the time limit. The academic tools are built for tabulated knots
    (c<=20); our objects start where the tables stop.
    Our own point tool (Disc #570) computes the same number at sample points
    in F_p + interpolation + CRT: c=128 4.4 s · c=192 11 s. Factor >100.

    VARIANT 6 is the CALIBRATED convention (alex_eichung.py: 12/12 hits;
    variants 0-5 fail completely, so the test is able to reject).
    Checked against our own corpus: 557/557 determinants hit,
    0 deviations, c from 3 to 128.

    Returns |Delta(-1)| or None. NEVER raises — a fallback that crashes
    itself is no fallback.
    """
    try:
        import sys as _sys
        _sieve = site_config.KAWAUCHI_SIEVE_DIR
        if _sieve not in _sys.path:
            _sys.path.insert(0, _sieve)
        from alex_point import parse_pd, alexander_poly, normalise
        pd = parse_pd(pd_str)
        if not pd:
            return None
        co = normalise(alexander_poly(pd, variant))
        return abs(sum(c * ((-1) ** i) for i, c in enumerate(co)))
    except Exception as e:
        print(f"  Point Alexander not possible: {type(e).__name__}: {e}")
        return None


def pd_from_stage4(output_path):
    """Fetch a PD from the stage-4 artifact when regina delivers none.

    On regina_failed, pd_code becomes None -> without this, Snappy fails too,
    although stage 4 long ago produced a PD and stored it next to the report.
    The smallest one is preferred (spherogram regularly beats regina, Disc #562).
    """
    import glob as _glob
    folder = os.path.dirname(os.path.abspath(output_path)) or '.'
    hits = sorted(_glob.glob(os.path.join(folder, '*_aggressive_simplify.json')))
    if not hits:
        return None, None
    try:
        d = json.load(open(hits[0], encoding='utf-8'))
        best = d.get('best', {}) or {}
        cand_k = [(best.get('sphe_pd_best_n'), best.get('sphe_pd_best'), 'spherogram'),
                (best.get('size_simp_regina'), best.get('pd'), 'regina')]
        cand_k = [(n, s, q) for n, s, q in cand_k
                if isinstance(s, str) and s.startswith('PD[') and isinstance(n, int)]
        if not cand_k:
            return None, None
        cand_k.sort()
        return cand_k[0][1], f"{hits[0]} ({cand_k[0][2]}, c={cand_k[0][0]})"
    except Exception as e:
        print(f"  Stage-4 PD not readable: {type(e).__name__}: {e}")
        return None, None


def alex_fingerprint(alex_str, points=(2, 3, -2)):
    """Distinctness fingerprint: Alexander at several points (Disc #248 style).
    For high c, replaces the symbolic HOMFLY canonicalisation in the distinctness key,
    where HOMFLY is exponentially out of reach."""
    if not alex_str:
        return None
    out = {}
    for p in points:
        v = _eval_regina_poly(alex_str, p)
        if v is not None:
            out[str(p)] = int(v) if v.denominator == 1 else str(v)
    return out or None


# ──────────────────────────────────────────────
# SNAPPY/SPHEROGRAM (snappy_env): hyperbolic + HFK
# ──────────────────────────────────────────────

def snappy_full(pd_code, point_eval=False):
    """pd_code: PD code string such as 'PD[X[1,4,2,5], X[3,6,4,1], ...]' or
    list-of-tuples string '[(1,4,2,5),(3,6,4,1),...]'.
    point_eval=True: HFK (knot_floer_homology, exponential) is skipped;
    geometry (vol/sym/cusp) + Seifert det/signature (numerical) stay."""
    code = SEAL_PRELUDE + textwrap.dedent(f"""
        # --- length spectrum: ONE schema ---
        def _lspec_entry(g):
            # ONE schema for the length spectrum (finding 2026-08-24). The fallback paths used
            # to write strings, the main path floats — whoever read 'length_real'
            # saw NOTHING for the strings. 53 records were silently invisible that way.
            # SnapPy.Number: .real()/.imag() are METHODS (measured in snappy_env).
            try:
                return dict(length_real=float(g.length.real()),
                            length_imag=float(g.length.imag()),
                            multiplicity=int(g.multiplicity))
            except Exception as _e:
                # Never lose silently: keep the old form AND write the reason next to it.
                return dict(length=str(g.length), mult=int(g.multiplicity),
                            schema_err=str(_e))
        import json, snappy, spherogram
        pd_str = {pd_code!r}
        # Try to evaluate as Python list of tuples
        try:
            pd_list = eval(pd_str)
            L = spherogram.Link(pd_list)
        except Exception as e1:
            try:
                L = spherogram.Link(pd_str)
            except Exception as e2:
                print(json.dumps({{'error': f'parse: {{e1}} / {{e2}}'}}))
                raise SystemExit
        out = {{'parsed': True, 'pd_input': pd_str}}
        try: out['n_crossings'] = len(L.crossings)
        except: pass
        try:
            Lg = L.copy(); Lg.simplify(mode='global')
            out['simplified_global'] = len(Lg.crossings)
        except Exception as e: out['simp_err'] = str(e)
        try:
            M = L.exterior()
            out['volume'] = float(M.volume())
            out['chern_simons'] = float(M.chern_simons())
            out['num_tetrahedra'] = int(M.num_tetrahedra())
            out['solution_type'] = str(M.solution_type())
            out['homology'] = str(M.homology())
            try: out['isosig'] = M.triangulation_isosig()
            except Exception as e: out['isosig_err'] = str(e)
            try:
                ci = M.cusp_info()
                out['cusp_shape'] = str(ci[0]) if ci else None
            except Exception as e: out['cusp_err'] = str(e)
            try:
                sg = M.symmetry_group()
                out['symmetry_group_order'] = int(sg.order())
                out['is_amphicheiral']      = bool(sg.is_amphicheiral())
                out['is_invertible_knot']   = bool(sg.is_invertible_knot())
                out['symmetry_group_repr']  = str(sg)
            except Exception as e: out['symgroup_err'] = str(e)
            try:
                ls = M.length_spectrum(0.5)
                out['length_spectrum_short'] = [
                    _lspec_entry(g)
                    for g in ls[:10]
                ]
            except Exception as e: out['lspec_err'] = str(e)
        except Exception as e:
            out['exterior_err'] = str(e)
        try:
            if {point_eval!r}: raise RuntimeError('point_eval: HFK (exponential) skipped')
            hfk = L.knot_floer_homology()
            out['knot_floer'] = {{
                'seifert_genus': int(hfk['seifert_genus']),
                'tau':           int(hfk['tau']),
                'nu':            int(hfk['nu']),
                'epsilon':       int(hfk['epsilon']),
                'fibered':       bool(hfk['fibered']),
                'L_space_knot':  bool(hfk['L_space_knot']),
                'total_rank':    int(hfk['total_rank']),
                'modulus':       int(hfk['modulus']),
            }}
        except Exception as e: out['hfk_err'] = str(e)
        # Seifert matrix for det/signature by hand
        # det EXACT in integers (Disc #582, 2026-07-27): np.linalg.det computes in float64
        # and reported spurious conflicts with the stored det above ~1e15 —
        # 6 cases in the corpus, relative error up to 4e-13, all at det >= 3e13. V+V^T is
        # integral, so the third det witness belongs computed without fractions (Bareiss).
        # The signature stays with the eigenvalues: only signs count there.
        try:
            V = L.seifert_matrix()
            import numpy as np
            # spherogram seifert_matrix returns a SnapPy matrix or a list of lists
            try:
                rows = [[int(x) for x in row] for row in np.array(V).tolist()]
            except Exception:
                rows = [[int(x) for x in list(row)] for row in V]
            n_s = len(rows)
            S_int = [[rows[i][j] + rows[j][i] for j in range(n_s)] for i in range(n_s)]

            def _det_bareiss(M):
                A = [r[:] for r in M]; n = len(A); vz = 1; prev = 1
                for k in range(n - 1):
                    if A[k][k] == 0:
                        for i2 in range(k + 1, n):
                            if A[i2][k] != 0:
                                A[k], A[i2] = A[i2], A[k]; vz = -vz; break
                        else:
                            return 0
                    for i2 in range(k + 1, n):
                        for j2 in range(k + 1, n):
                            A[i2][j2] = (A[i2][j2] * A[k][k] - A[i2][k] * A[k][j2]) // prev
                    prev = A[k][k]
                return vz * A[n - 1][n - 1] if n else 1

            out['seifert_matrix_size'] = n_s
            out['det_via_seifert'] = abs(_det_bareiss(S_int))
            out['det_via_seifert_method'] = 'bareiss_exact_integer'
            evs = np.linalg.eigvalsh(np.array(S_int, dtype=float))
            out['signature'] = int((evs > 1e-10).sum() - (evs < -1e-10).sum())
            # CONVENTION, measured on 675 records (Disc #582): this signature is throughout
            # the NEGATIVE of Sage Link.signature() on the same PD. It is a witness for
            # |sigma|, NOT for the sign — whoever puts it unflipped next to the column
            # sees nothing but spurious contradictions.
            out['signature_konvention'] = 'spherogram V+V^T = -1 * sage Link.signature() (675/675 measured)'
        except Exception as e:
            out['seifert_err'] = str(e)
        print(json.dumps(out))
    """)
    return run_subprocess('snappy_env', code, timeout=300)


# ──────────────────────────────────────────────
# SAGE (sage_stable): Knot class for arf, signature
# ──────────────────────────────────────────────

def sage_invariants(pd_code, point_eval=False):
    code = SEAL_PRELUDE + textwrap.dedent(f"""
        import json
        from sage.all import *
        pd_str = {pd_code!r}
        out = {{}}
        try:
            pd_list = eval(pd_str)
            # Sage Link is 1-based, spherogram 0-based (Discovery #13)
            # if min(arc_id) == 0, shift all by +1
            all_ids = [aid for cr in pd_list for aid in cr]
            shift = 0
            if min(all_ids) == 0:
                shift = 1
            pd_shifted = [tuple(a + shift for a in cr) for cr in pd_list]
            L = Link(pd_shifted)
            out['parsed_sage_link'] = True
            out['determinant'] = int(L.determinant())
            out['signature']   = int(L.signature())
            out['is_alternating'] = bool(L.is_alternating())
            out['writhe']      = int(L.writhe())
            _pe = {point_eval!r}
            try:
                if _pe: raise RuntimeError('point_eval: HOMFLY (libhomfly does not scale) skipped')
                out['homfly_sage'] = str(L.homfly_polynomial())
            except Exception as e: out['homfly_err'] = str(e)
            try:
                if _pe: raise RuntimeError('point_eval: Jones skipped')
                out['jones_sage']  = str(L.jones_polynomial())
            except Exception as e: out['jones_err'] = str(e)
            try: out['alexander_sage'] = str(L.alexander_polynomial())
            except Exception as e: out['alex_err'] = str(e)
            try:
                B = L.braid()
                out['braid_n_strands'] = int(B.parent().strands())
                out['braid_word']      = str(list(B.Tietze()))
                out['braid_length']    = len(list(B.Tietze()))
            except Exception as e: out['braid_err'] = str(e)
            # Knot class via the braid for Arf
            try:
                K = Knot(B)
                out['arf_invariant'] = int(K.arf_invariant())
                try: out['knot_alexander'] = str(K.alexander_polynomial())
                except: pass
            except Exception as e: out['knot_err'] = str(e)
        except Exception as e:
            out['error'] = str(e)
        print(json.dumps(out))
    """)
    return run_sage(code, timeout=600)


# ──────────────────────────────────────────────
# KnotJob (java): Khovanov homology
# ──────────────────────────────────────────────

KNOTJOB_JAR = site_config.KNOTJOB_JAR
JAVA_BIN    = site_config.JAVA_BIN

def knotjob_khovanov(pd_code, name='knot'):
    """Compute Khovanov homology via KnotJob.jar.
    pd_code: '[(1,4,2,5),(3,6,4,1),...]' (1-based, Sage style)."""
    try:
        pd_list = eval(pd_code)
    except Exception as e:
        return {'error': f'pd parse: {e}'}
    all_ids = [aid for cr in pd_list for aid in cr]
    if not all_ids:
        return {'error': 'pd_list empty (regina returned PD[] — degenerate diagram)'}
    shift = 1 if min(all_ids) == 0 else 0
    pd_1 = [tuple(a + shift for a in cr) for cr in pd_list]
    pd_str = ','.join(f'[{a},{b},{c},{d}]' for (a, b, c, d) in pd_1)
    inp = f"{name} = PD[{pd_str}]"
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.txt', delete=False) as f:
        f.write(inp + '\n')
        tmppath = f.name
    try:
        res = subprocess.run(
            [JAVA_BIN, '-jar', KNOTJOB_JAR, tmppath, '-kb0', '-nf'],
            capture_output=True, text=True, timeout=600,
        )
    except subprocess.TimeoutExpired:
        os.unlink(tmppath)
        return {'error': 'timeout'}
    finally:
        try: os.unlink(tmppath)
        except: pass
    out = {
        'rc': res.returncode,
        'khovanov_stdout_tail': res.stdout[-2000:],
        'khovanov_stderr_tail': res.stderr[-500:],
    }
    return out


def knotjob_s(pd_code, name='knot'):
    """Compute the Rasmussen s invariant via KnotJob.jar -s0 (Lee/Khovanov,
    characteristic 0). pd_code as in knotjob_khovanov (1- or 0-based).
    Returns {'rasmussen_s': int} or {'error': ...}. Sign = KnotJob convention
    (for alternating knots s = -sigma; trefoil -> -2)."""
    try:
        pd_list = eval(pd_code)
    except Exception as e:
        return {'error': f'pd parse: {e}'}
    all_ids = [aid for cr in pd_list for aid in cr]
    if not all_ids:
        return {'error': 'pd_list empty (regina returned PD[] — degenerate diagram)'}
    shift = 1 if min(all_ids) == 0 else 0
    pd_1 = [tuple(a + shift for a in cr) for cr in pd_list]
    pd_str = ','.join(f'[{a},{b},{c},{d}]' for (a, b, c, d) in pd_1)
    inp = f"{name} = PD[{pd_str}]"
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.txt', delete=False) as f:
        f.write(inp + '\n')
        tmppath = f.name
    try:
        res = subprocess.run(
            [JAVA_BIN, '-jar', KNOTJOB_JAR, tmppath, '-s0', '-nf'],
            capture_output=True, text=True, timeout=600,
        )
    except subprocess.TimeoutExpired:
        return {'error': 'timeout'}
    finally:
        try: os.unlink(tmppath)
        except: pass
    m = re.search(r'S-Invariant mod 0\s*:\s*(-?\d+)', res.stdout)
    out = {
        'rc': res.returncode,
        's_stdout_tail': res.stdout[-800:],
        's_stderr_tail': res.stderr[-500:],
    }
    if m:
        out['rasmussen_s'] = int(m.group(1))
    return out


# ──────────────────────────────────────────────
# CONWAY from HOMFLY (Discovery #73)
# ──────────────────────────────────────────────

# ──────────────────────────────────────────────
# HOMFLY FALLBACK (approved by Werner 2026-08-23: "close the cause first")
# ──────────────────────────────────────────────
HOMFLY_FALLBACK_TOOLS = site_config.TOOLS_DIR


def homfly_fallback_level(report, point_eval):
    """Fills homfly_polynomial when regina did NOT deliver it — and
    labels the gap when that fails too.

    WHY THIS EXISTS (Disc #674): the point_eval cap does not even send regina's
    exponential lines. That is right — at high c they cost hours or the
    memory. The only mistake was that an EMPTY field remained afterwards, which
    looked like "impossible". On 23.08. it turned out that of 229 such gaps,
    178 were computable in seconds. Not failed — skipped.

    The second route (Hecke trace over F_p + CRT, tools/homfly_modular.py) is
    COMPLEMENTARY to regina: regina dies on the diagram size c, the
    Hecke trace on the strand count n (cost ~ n!). Where one stalls, the other
    usually runs.

    FOUR WITNESSES, otherwise nothing is written:
      1. CRT stable (one more prime no longer changes the result)
      2. fresh points in F_p that were not on the grid
      3. det = |P(1,2i)| against the determinant in the record
      4. Alexander at several points, up to the unit ±t^k
    And if one breaks: `homfly_gap` says WHY. A labelled gap
    is a finding, an empty one is a trap.
    """
    a = report['mandatory_list']['ALGEBRA']
    d = report['mandatory_list']['DIAGRAM']
    if a.get('homfly_polynomial'):
        return                                   # regina delivered
    if os.environ.get('HOMFLY_FALLBACK', '1') == '0':
        a['homfly_gap'] = {'cause': 'per HOMFLY_FALLBACK=0 abgeschaltet'}
        return

    w, n = d.get('braid_word'), d.get('braid_n_strands')
    if isinstance(w, str):
        try:
            w = ast.literal_eval(w)
        except Exception:
            w = None
    n_max = int(os.environ.get('HOMFLY_FALLBACK_NMAX', 10))
    if not w or not n:
        a['homfly_gap'] = {'cause': 'no braid word in the record — the Hecke trace '
                                       'needs one', 'route': 'hecke_f_p'}
        return
    if int(n) > n_max:
        a['homfly_gap'] = {
            'cause': f'strand number n={int(n)} > {n_max}: the Hecke trace costs ~n!, '
                     f'from n=11 on it takes hours (representation dimension 2310, '
                     f'7700 at n=12)',
            'route': 'hecke_f_p', 'n': int(n),
            'hint': 'can be filled later with HOMFLY_FALLBACK_NMAX=12 or '
                       'tools/homfly_interpolation.py --targets --exact'}
        return

    zg = d.get('crossing_number_min')
    if not zg:
        a['homfly_gap'] = {'cause': 'no crossing number — without it the '
                                       'z axis of the interpolation box is not fixed',
                              'route': 'hecke_f_p'}
        return

    latch = int(os.environ.get('HOMFLY_FALLBACK_TIMEOUT', 1800))
    cmd = [site_config.PYTHON_MAIN,
           os.path.join(HOMFLY_FALLBACK_TOOLS, 'homfly_modular.py'),
           '--braid-word', json.dumps([int(g) for g in w]),
           '--n', str(int(n)), '--zdeg', str(int(zg))]
    ev = a.get('alexander_point_evals')
    if ev:
        cmd += ['--alex-evals', json.dumps(ev)]
    print(f"\n[PFLICHT] (3b) HOMFLY fallback (Hecke trace over F_p, "
          f"n={int(n)}, c={int(zg)}) …")
    t0 = time.time()
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=latch)
    except subprocess.TimeoutExpired:
        a['homfly_gap'] = {'cause': f'Riegel {latch}s erreicht',
                              'route': 'hecke_f_p'}
        print(f"  Time limit {latch}s reached — gap labelled, run continues")
        return
    try:
        out = json.loads(res.stdout.strip().splitlines()[-1])
    except Exception:
        a['homfly_gap'] = {'cause': 'fallback level returned no JSON',
                              'route': 'hecke_f_p',
                              'stderr': (res.stderr or '')[-400:]}
        print("  no JSON returned — gap labelled")
        return
    if out.get('error'):
        a['homfly_gap'] = {'cause': out['error'], 'route': 'hecke_f_p'}
        print(f"  Error: {out['error'][:120]} — gap labelled")
        return

    det_rec = (a.get('determinant') if a.get('determinant') is not None
               else a.get('determinant_alexander'))
    reasons = []
    if not out.get('stable'):
        reasons.append('CRT did not become stable (more primes needed)')
    if out.get('fresh_hits') != out.get('fresh_total'):
        reasons.append(f"fresh points {out.get('fresh_hits')}/"
                       f"{out.get('fresh_total')} in F_p")
    if det_rec is not None and out.get('det') != int(det_rec):
        reasons.append(f"det {out.get('det')} != record {int(det_rec)}")
    if out.get('alexander_points') is False:
        reasons.append(f"Alexander points: {out.get('alexander_points_text')}")
    if reasons:
        a['homfly_gap'] = {'cause': '; '.join(reasons), 'route': 'hecke_f_p',
                              'sek': round(time.time() - t0, 1)}
        print(f"  Witnesses disagree ⟹ NOTHING written: {'; '.join(reasons)[:160]}")
        return

    a['homfly_polynomial'] = out['homfly']
    if not a.get('conway_polynomial'):
        a['conway_polynomial'] = conway_from_homfly(out['homfly'])
    a.pop('homfly_gap', None)
    report['homfly_provenance'] = {
        'datum': time.strftime('%Y-%m-%d'),
        'route': ('Hecke trace over F_p on a grid + CRT, then '
                '2-D interpolation (pipeline/homfly_modular.py)'),
        'why_not_regina': ('point_eval mode: regina\'s exponential lines '
                               '(homfly/jones/bracket) were not sent along'
                               if point_eval else
                               'regina returned no homfly (time or memory)'),
        'convention': 'regina AZ, variables printed as x=a, y=z',
        'primes': out.get('primes'), 'terms': out.get('terms'),
        'witness_1': 'CRT stable (one more prime changes nothing)',
        'witness_2': f"fresh points in F_p: {out.get('fresh_hits')}/"
                   f"{out.get('fresh_total')} exactly equal",
        'witness_3': f"det |P(1,2i)| = {out.get('det')} = record {det_rec}",
        'witness_4': out.get('alexander_points_text'),
        'seconds': round(time.time() - t0, 1),
    }
    print(f"  ✅ HOMFLY filled in ({out.get('terms')} terms, "
          f"{out.get('primes')} primes, det {out.get('det')}, "
          f"{time.time()-t0:.1f}s)")


def conway_from_homfly(homfly_str):
    """∇(z) = HOMFLY(x=1, y=z). Set x=1, collect by powers of y.
    Accepts Regina format with variables 'x', 'y' or Sage format."""
    if not homfly_str:
        return None
    # Regina format: 'x^6 - 3 x^4 y^2 + 2 x^2 y^4 - 3 x^2 y^2 - 3 x^2 + 5 y^4 + ...'
    # Protect negative exponents before the term split
    s = (homfly_str.replace(' ', '')
                   .replace('x^-', 'x^N').replace('y^-', 'y^M'))
    # Replace 'x^k' or 'x' (with optional sign) by 1 → leave coefficient
    # Tokenize into terms
    s2 = '+' + s if s[0] not in '+-' else s
    # Split into ±term
    raw_terms = re.split(r'(?=[+-])', s2)[1:]
    coefs_by_y = {}
    for t in raw_terms:
        sign = 1
        if t[0] == '-':
            sign = -1
            t = t[1:]
        elif t[0] == '+':
            t = t[1:]
        t = t.replace('x^N', 'x^-').replace('y^M', 'y^-')  # restore
        # parse: coef × x^? × y^?
        m = re.match(r'^(\d+)?(?:x(?:\^(-?\d+))?)?(?:y(?:\^(-?\d+))?)?$', t)
        if not m:
            return None
        c, xe, ye = m.groups()
        c = int(c) if c else 1
        ye_n = int(ye) if ye else (1 if 'y' in t else 0)
        coefs_by_y.setdefault(ye_n, 0)
        coefs_by_y[ye_n] += sign * c   # x^anything = 1 at x=1
    # Build the Conway string
    parts = []
    for ye in sorted(coefs_by_y.keys()):
        c = coefs_by_y[ye]
        if c == 0:
            continue
        if ye == 0:
            parts.append(f"{c}")
        elif ye == 1:
            parts.append(f"{c}z")
        else:
            parts.append(f"{c}z^{ye}")
    return ' + '.join(parts) if parts else '0'


# ──────────────────────────────────────────────
# CONSISTENCY QUARTET (Discovery #76)
# ──────────────────────────────────────────────

def jones_at_i(jones_str):
    """|Jones(i)| via the substitution x=i, x^k = i^k.
    Protects 'x^-' from the split regex (otherwise negative exponents get cut apart).

    🔴 INTEGER, NOT complex (finding 2026-08-24). This used to read `val = 0+0j`,
    i.e. floating point. i^k only cycles through (1, i, -1, -i) — the value IS an integer, and
    the determinant of a knot is always ODD. For gitter_d3_k6_s1p8_g0p25
    (c_min=141) the float route returned 1178895905925319098368 — an EVEN number, hence
    impossible; the exact value is 1178895905925318655335 (= |Delta(-1)|, independently confirmed).
    Above 2^53 a double loses integrality, and a wrong bit looks
    like a result. Related: Disc #688 (mantissa guard)."""
    if not jones_str:
        return None
    s = jones_str.replace(' ', '').replace('x^-', 'x^N')  # marker for negative exponents
    s2 = '+' + s if s[0] not in '+-' else s
    raw_terms = re.split(r'(?=[+-])', s2)[1:]
    re_part = 0
    im_part = 0
    for t in raw_terms:
        sign = 1
        if t[0] == '-': sign = -1; t = t[1:]
        elif t[0] == '+': t = t[1:]
        t = t.replace('x^N', 'x^-')  # restore
        m = re.match(r'^(\d+)?(?:x(?:\^(-?\d+))?)?$', t)
        if not m:
            continue
        c, xe = m.groups()
        c = int(c) if c else 1
        xe_n = int(xe) if xe else (1 if 'x' in t else 0)
        # i^k cyclic: 1, i, -1, -i  — two integers instead of one float
        r = xe_n % 4
        if r == 0:
            re_part += sign * c
        elif r == 1:
            im_part += sign * c
        elif r == 2:
            re_part -= sign * c
        else:
            im_part -= sign * c
    return math.isqrt(re_part * re_part + im_part * im_part)


# ──────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────

def main():
    # Hole B (OOM 2026-07-01): pass an external kill through to the conda/sage grandchildren, otherwise
    # they survive as RAM-eating orphans. As in mandatory_invariants_extra.py:857.
    install_signal_killer()
    parser = argparse.ArgumentParser(description='Full mandatory list of knot invariants.')
    parser.add_argument('--signed-gauss', help='Signed Gauss code string')
    parser.add_argument('--kinematic-result',
                        help='kinematic_result.json (extracts signed_gauss from chosen)')
    parser.add_argument('--name', default='knot', help='Name for reporting')
    parser.add_argument('--output', default='mandatory_invariants.json')
    parser.add_argument('--skip-sage', action='store_true', help='Skip Sage (slow)')
    parser.add_argument('--skip-snappy', action='store_true')
    parser.add_argument('--skip-knotjob', action='store_true')
    parser.add_argument('--point-eval', action='store_true',
                        help='High-c cap (Disc #248): omit homfly/jones/bracket/HFK/khovanov; '
                             'det via Alexander(-1)+Goeritz+Seifert (polynomial time). '
                             'Alternatively ENV MANDATORY_POINT_EVAL=1.')
    args = parser.parse_args()

    point_eval = args.point_eval or os.environ.get('MANDATORY_POINT_EVAL') == '1'
    if point_eval:
        print("[PFLICHT] POINT-EVAL mode (Disc #248): homfly/jones/bracket/HFK/khovanov "
              "skipped — value instead of symbol. det via Alexander(-1)+Goeritz+Seifert.")

    if args.kinematic_result:
        with open(args.kinematic_result) as f:
            kin = json.load(f)
        signed_gauss = kin['chosen']['signed_gauss']
    elif args.signed_gauss:
        signed_gauss = args.signed_gauss
    else:
        print("--signed-gauss or --kinematic-result required.", file=sys.stderr)
        sys.exit(2)

    print(f"[PFLICHT] Name: {args.name}")
    print(f"[PFLICHT] Signed-Gauss ({len(signed_gauss)//3} Tokens):")
    print(f"          {signed_gauss}")

    report = {
        'name':         args.name,
        'signed_gauss': signed_gauss,
        'mandatory_list': {
            'DIAGRAM':     {},
            'ALGEBRA':      {},
            'HYPERBOLIC': {},
        },
        'raw_subprocess_outputs': {},
    }

    print(f"\n[PFLICHT] (1) Regina (base) — algebra + diagram codes")
    t0 = time.time()
    reg = regina_full(signed_gauss, point_eval=point_eval)
    report['raw_subprocess_outputs']['regina'] = reg
    print(f"  → {time.time()-t0:.1f}s")
    if reg.get('parsed') and reg.get('classical'):
        d = report['mandatory_list']['DIAGRAM']
        a = report['mandatory_list']['ALGEBRA']
        d['gauss_code_reduced'] = reg.get('gauss_code_reduced')
        d['signed_gauss_reduced'] = reg.get('signed_gauss_reduced')
        d['oriented_gauss']     = reg.get('oriented_gauss')
        d['dt_code']            = reg.get('dt_code')
        d['pd_code']            = reg.get('pd_code')
        d['crossing_number_min'] = reg.get('simplified_size')
        d['alternating']        = reg.get('alternating')
        d['knot_sig']           = reg.get('knot_sig')
        d['writhe']             = reg.get('writhe')
        d['frequency']          = 1.0
        d['frequency_ratio']    = 1.0
        a['homfly_polynomial']  = reg.get('homfly')
        a['jones_polynomial']   = reg.get('jones')
        a['alexander_polynomial'] = reg.get('alexander')
        a['kauffman_bracket']   = reg.get('kauffman')
        a['conway_polynomial']  = conway_from_homfly(reg.get('homfly'))
        a['determinant']        = jones_at_i(reg.get('jones'))
        # Point evaluation (Disc #248): det = |Δ(-1)| + distinctness fingerprint.
        # ALWAYS set (also in full mode = independent third det witness).
        a['determinant_alexander'] = alex_det(reg.get('alexander'))
        a['alexander_point_evals'] = alex_fingerprint(reg.get('alexander'))
    else:
        report['regina_failed'] = True
        print(f"  Regina failed: {reg}")

    pd_code = reg.get('pd_code') if reg.get('parsed') else None

    # FALLBACK (approved by Werner 2026-08-22, "patch it"): if regina fails, the
    # run is not lost — stage 4 has a PD, and our own point tool
    # computes the determinant from it in seconds. Without this, stage 5 aborted with
    # "no det witness" (3 runs on 22.08.), although everything needed was there.
    if report.get('regina_failed'):
        ersatz_pd, source_name = pd_from_stage4(args.output)
        if ersatz_pd:
            if not pd_code:
                pd_code = ersatz_pd
                print(f"  ↻ PD taken over from stage 4: {source_name}")
                # 2026-08-23: the replacement PD must also go INTO THE REPORT, not only into
                # the local variable. The fallback of 22.08. rescued stage 5
                # (det witness), but stage 6 reads mandatory_list.DIAGRAMM.pd_code —
                # that field is only set in the regina branch. Result: 11 runs
                # died with "FEHLER: pd_code fehlt im Input" (pd_code missing in input) AFTER stage 5
                # had computed for 20+ minutes and the PD was lying right next to it.
                # Deliberately NOT set: crossing_number_min. This PD is a
                # simplified diagram; its size is an UPPER BOUND, not a
                # minimum (Disc #231: honest bound "<=N", never report it as c).
                report['mandatory_list']['DIAGRAM']['pd_code'] = ersatz_pd
                report['pd_code_provenance'] = {
                    'origin': source_name,
                    'cause': 'regina_failed — PD from stage 4 (aggressive_simplify)',
                    'warning': ('diagram size is an upper bound for the '
                                'crossing number, not a minimum.')}
            t_pk = time.time()
            d_pk = point_alexander_det(ersatz_pd)
            if d_pk is not None:
                report['mandatory_list']['ALGEBRA']['determinant_point'] = d_pk
                report['point_alexander'] = {
                    'det': d_pk, 'origin': source_name, 'variant': 6,
                    'seconds': round(time.time() - t_pk, 2),
                    'note': ('geeicht alex_eichung.py Variante 6; '
                             '557/557 Korpus-Determinanten getroffen, 0 Abweichungen')}
                print(f"  ✓ Point Alexander (fallback): |det| = {d_pk} "
                      f"({time.time() - t_pk:.1f}s)")
            else:
                print(f"  🔴 Point Alexander delivered nothing — fallback does not apply.")
        else:
            print(f"  🔴 No stage-4 PD found — fallback not possible.")

    if pd_code and not args.skip_snappy:
        print(f"\n[PFLICHT] (2) Snappy/Spherogram (snappy_env) — hyperbolic + HFK + Seifert")
        # PD code format conversion: regina delivers "PD[ X(...), X(...), ... ]"
        # spherogram needs either 'PD[X[...], ...]' or a list of tuples
        pd_for_snappy = pd_code.replace('PD[', '').replace(']', '').strip()
        # We parse as a list of tuples: 'X(1,2,3,4),X(...) ...'
        crossings = re.findall(r'X[\(\[](\-?\d+(?:,\s*\-?\d+){3})[\)\]]', pd_code)
        if not crossings:
            print(f"  PD code not parsable: {pd_code}")
            pd_list_str = '[]'
        else:
            tuples = []
            for cr in crossings:
                vals = [int(v.strip()) for v in cr.split(',')]
                tuples.append(tuple(vals))
            pd_list_str = repr(tuples)
        t0 = time.time()
        snp = snappy_full(pd_list_str, point_eval=point_eval)
        report['raw_subprocess_outputs']['snappy'] = snp
        report['pd_list_str'] = pd_list_str
        print(f"  → {time.time()-t0:.1f}s")
        h = report['mandatory_list']['HYPERBOLIC']
        a = report['mandatory_list']['ALGEBRA']
        d = report['mandatory_list']['DIAGRAM']
        h['hyperbolic_volume']   = snp.get('volume')
        h['chern_simons']        = snp.get('chern_simons')
        h['cusp_shape']          = snp.get('cusp_shape')
        h['symmetry_group']      = snp.get('symmetry_group_repr')
        h['symmetry_group_order'] = snp.get('symmetry_group_order')
        h['is_amphicheiral']     = snp.get('is_amphicheiral')
        h['is_invertible_knot']  = snp.get('is_invertible_knot')
        h['homology']            = snp.get('homology')
        h['canonical_retriangulation'] = snp.get('isosig')
        h['length_spectrum_short'] = snp.get('length_spectrum_short')
        if 'knot_floer' in snp:
            hfk = snp['knot_floer']
            a['knot_floer_homology'] = hfk
            a['genus']           = hfk.get('seifert_genus')
            a['fibered']         = hfk.get('fibered')
            a['tau_invariant']   = hfk.get('tau')
            # smooth_4_genus: NOT |tau|. Only g_4 >= |tau| holds (Ozsvath-Szabo), not
            # equality -- the same overclaim class as the old slice field
            # (Disc #555/#556). Signature and s are not yet known here either.
            # The verdict is made at the end with four_genus_verdict.apply_to_algebra().
            # slice: do NOT decide here. tau=0 is necessary, not sufficient (error log
            # 2026-07-23 [SLICE-FELD]: the old line declared every tau=0 knot a
            # slice knot -- 213 records, 29 provably wrong). Signature, Arf and determinant
            # are not yet known at this point; the verdict is made at the end in one
            # block with slice_verdict.apply_to_algebra().
            a['ribbon']          = 'undetermined (not algorithmically decidable in the pipeline)'
        if 'signature' in snp:
            # Store in the corpus convention (= Sage Link.signature) so that the field is directly
            # readable next to `signature`. The raw spherogram value is its negative
            # (675/675 measured, Disc #582) and stays alongside so that the provenance
            # is checkable. WITNESS FOR |sigma| — not for the sign.
            a['signature_seifert'] = -snp['signature']
            a['signature_seifert_raw_spherogram'] = snp['signature']
            a['signature_seifert_convention'] = snp.get('signature_convention')
        if 'det_via_seifert' in snp:
            a['determinant_seifert'] = snp['det_via_seifert']
            a['determinant_seifert_method'] = snp.get('det_via_seifert_method')
    else:
        if not pd_code: print(f"\n[PFLICHT] (2) Snappy skipped (no PD code)")
        else: print(f"\n[PFLICHT] (2) Snappy skipped (--skip-snappy)")

    if pd_code and not args.skip_sage:
        print(f"\n[PFLICHT] (3) Sage (sage_stable) — Signature, Arf, Braid")
        crossings = re.findall(r'X[\(\[](\-?\d+(?:,\s*\-?\d+){3})[\)\]]', pd_code)
        tuples = [tuple(int(v.strip()) for v in cr.split(',')) for cr in crossings]
        pd_list_str = repr(tuples)
        t0 = time.time()
        sage_out = sage_invariants(pd_list_str, point_eval=point_eval)
        report['raw_subprocess_outputs']['sage'] = sage_out
        print(f"  → {time.time()-t0:.1f}s")
        a = report['mandatory_list']['ALGEBRA']
        d = report['mandatory_list']['DIAGRAM']
        a['signature']   = sage_out.get('signature')
        a['arf_invariant'] = sage_out.get('arf_invariant')
        a['determinant_sage'] = sage_out.get('determinant')
        d['braid_n_strands'] = sage_out.get('braid_n_strands')
        d['braid_word']      = sage_out.get('braid_word')
        d['braid_length']    = sage_out.get('braid_length')
        d['braid_index']     = sage_out.get('braid_n_strands')  # NOT minimal
    else:
        if not pd_code: print(f"\n[PFLICHT] (3) Sage skipped (no PD code)")
        else: print(f"\n[PFLICHT] (3) Sage skipped (--skip-sage)")

    # Cause closed (Werner 2026-08-23): a skipped computation must not
    # leave an empty field. Either the second route fills it — or the
    # gap carries its reason. Placed here because only the Sage step sets the
    # braid word that the Hecke trace needs.
    try:
        homfly_fallback_level(report, point_eval)
    except Exception as _e:
        report['mandatory_list']['ALGEBRA']['homfly_gap'] = {
            'cause': 'fallback level itself failed: ' + repr(_e)[:200],
            'route': 'hecke_f_p'}
        print(f"  Fallback failed: {_e!r} — gap labelled, "
              f"run continues")

    if pd_code and not args.skip_knotjob and not point_eval:
        print(f"\n[PFLICHT] (4) KnotJob — Khovanov + s invariant")
        crossings = re.findall(r'X[\(\[](\-?\d+(?:,\s*\-?\d+){3})[\)\]]', pd_code)
        tuples = [tuple(int(v.strip()) for v in cr.split(',')) for cr in crossings]
        pd_list_str = repr(tuples)
        t0 = time.time()
        kj = knotjob_khovanov(pd_list_str, args.name)
        report['raw_subprocess_outputs']['knotjob'] = kj
        print(f"  → Khovanov {time.time()-t0:.1f}s")
        report['mandatory_list']['ALGEBRA']['khovanov_stdout_tail'] = (
            kj.get('khovanov_stdout_tail', '')[-1500:])
        t0 = time.time()
        ks = knotjob_s(pd_list_str, args.name)
        report['raw_subprocess_outputs']['knotjob_s'] = ks
        print(f"  → Rasmussen-s {time.time()-t0:.1f}s: "
              f"{ks.get('rasmussen_s', ks.get('error', '?'))}")
        if ks.get('rasmussen_s') is not None:
            report['mandatory_list']['ALGEBRA']['rasmussen_s'] = ks['rasmussen_s']

    # Consistency quartet (Discovery #76, stratified per Discovery #112)
    # ----------------------------------------------------------------------
    # Det paths (mandatory identity, ALWAYS):
    #   jones_at_i == sage_det == |det(V+V^T)|
    # HFK rank condition (tiered, Ozsváth-Szabó):
    #   always:       rank(HFK^) >= |det(K)|, difference even
    #   if thin:      rank(HFK^) == |det(K)|
    # If rank > det: the knot is NOT HFK-thin (no pipeline bug,
    # consistent with the theorem). Saturn/Mercury are alternating ⇒ thin
    # ⇒ rank=det. orbit is non-alternating (Khov width=6) ⇒ rank>det
    # with even difference is the correct observation.
    print(f"\n[PFLICHT] (5) Consistency quartet (Discovery #76, stratified #112):")
    a = report['mandatory_list']['ALGEBRA']
    diag = report['mandatory_list'].get('DIAGRAM', {})

    det_paths = {
        'jones_at_i':          a.get('determinant'),
        'alexander_det':       a.get('determinant_alexander'),
        'sage_determinant':    a.get('determinant_sage'),
        'seifert_determinant': a.get('determinant_seifert'),
        # fourth witness, only populated when regina failed (fallback 2026-08-22)
        'point_alexander':     a.get('determinant_point'),
    }
    det_paths = {k: v for k, v in det_paths.items() if v is not None}
    hfk_rank = (a.get('knot_floer_homology', {}).get('total_rank')
                if isinstance(a.get('knot_floer_homology'), dict) else None)

    print(f"  Det paths: {det_paths}")
    det_consistent = (len(set(det_paths.values())) == 1) if det_paths else False
    det_value = next(iter(det_paths.values())) if det_consistent else None

    if det_consistent:
        print(f"  ✓ Det paths agree: |det(K)| = {det_value}")
    else:
        print(f"  ✗ Det paths DIVERGE — pipeline bug "
              f"(sign / index base / convention)")

    # HFK rank stratification
    hfk_rank_check = {'rank': hfk_rank, 'det': det_value}
    hfk_rank_ok    = None
    hfk_rank_equals_det = None   # field name misleading until 18.09.2026 (Disc #797)
    hfk_parity_ok  = None

    if hfk_rank is not None and det_value is not None:
        diff = hfk_rank - det_value
        hfk_parity_ok = (diff % 2 == 0)
        hfk_rank_ok   = (hfk_rank >= det_value) and hfk_parity_ok
        hfk_rank_equals_det = (hfk_rank == det_value)
        hfk_rank_check['diff_rank_minus_det'] = diff
        hfk_rank_check['rank_geq_det']        = (hfk_rank >= det_value)
        hfk_rank_check['parity_even']         = hfk_parity_ok
        hfk_rank_check['hfk_rank_equals_det'] = hfk_rank_equals_det

        print(f"  HFK rank: {hfk_rank}, diff (rank-det) = {diff}")
        if hfk_rank_equals_det:
            print(f"  ✓ HFK rank = det (NOT equivalent to delta-thin)")
        elif hfk_rank_ok:
            print(f"  ✓ HFK NOT thin: rank>det with even diff — "
                  f"consistent with Ozsváth-Szabó (Discovery #112), no pipeline bug")
        else:
            if hfk_rank < det_value:
                print(f"  ✗ ERROR: rank < det violates the Ozsváth-Szabó lower bound")
            elif not hfk_parity_ok:
                print(f"  ✗ ERROR: rank-det = {diff} odd, violates the "
                      f"Ozsváth-Szabó parity theorem")

    # Cross-check with the alternating field (alternating ⇒ Khovanov-thin ⇒ HFK-thin)
    is_alt = diag.get('alternating')
    if is_alt is True and hfk_rank_equals_det is False:
        print(f"  ⚠ WARNING: alternating=True but HFK rank>det "
              f"— inconsistency, check")

    overall_consistent = det_consistent and (hfk_rank_ok is None or hfk_rank_ok is True)

    report['consistency_quartet']  = {
        'det_paths':          det_paths,
        'det_value':          det_value,
        'det_paths_consistent': det_consistent,
        'hfk_rank_check':     hfk_rank_check,
        'hfk_rank_consistent': hfk_rank_ok,
        'hfk_rank_equals_det': hfk_rank_equals_det,
        'hfk_delta_thin':     None,   # can only be filled by a delta measurement
        'hfk_thin_note':      'rank==det is NOT delta-thinness: it only forces equal Maslov parity per Alexander grading, delta and delta-2 are compatible with it (measured 2026-09-18 on kf304_gtw_1.75pi: rank 243 = det, deltas {-1:32, +1:211}). Real thinness is measured by claude_setup_v2/tools/pruef_hfk_thin.py; hfk_delta_thin stays None as long as it has not been computed for this knot. Disc #797.',
        'alternating':        is_alt,
    }
    report['consistency_ok'] = overall_consistent

    # === SLICE VERDICT (error log 2026-07-23 [SLICE-FELD]) ===
    # Here and not earlier: only now are tau, s, signature, Arf and determinant all available.
    # The pipeline can REFUTE slice, not prove it -> False / None, True only for the unknot.
    try:
        from slice_verdict import apply_to_algebra as _slice_apply
        _v = _slice_apply(report['mandatory_list']['ALGEBRA'], report.get('name'))
        print(f"  slice verdict: {_v} — {report['mandatory_list']['ALGEBRA']['slice_note'][:110]}")
    except Exception as _e:
        print(f"  ⚠ slice verdict failed: {type(_e).__name__}: {_e}")
        _v = None

    # === 4-GENUS VERDICT (same class, Disc #556) ===
    # g_4 >= |tau|, |s|/2, |signature|/2 and >= 1 if slice is refuted; g_4 <= g_3.
    # A number only when the bounds coincide (torus knots etc.), otherwise an interval.
    try:
        from four_genus_verdict import apply_to_algebra as _g4_apply
        _g4 = _g4_apply(report['mandatory_list']['ALGEBRA'], _v, report.get('name'))
        _iv = report['mandatory_list']['ALGEBRA'].get('smooth_4_genus_interval')
        print(f"  4-genus: {_g4 if _g4 is not None else 'undetermined'} (interval {_iv})")
    except Exception as _e:
        print(f"  ⚠ 4-genus verdict failed: {type(_e).__name__}: {_e}")

    if overall_consistent:
        print(f"  ✓ Pipeline consistency OK")
    elif not det_consistent:
        print(f"  ✗ Det path divergence is a pipeline bug — check the values")
    else:
        print(f"  ✗ HFK rank condition violated — check the values")

    # ── GUARD "NEVER RAISE c" (2026-07-25, Disc #562/#563) ────────────────────
    # Stage-4 frames are sampled, spherogram simplify is randomised -> a
    # fresh full run can find a WORSE bound than the one already recorded
    # in the record. Really happened on `sundown`: 19 (recorded even 18) -> 20.
    # c is an UPPER BOUND (Disc #231): an already verified smaller
    # bound stays true, so a re-run must not overwrite it.
    # Only LOWERING is allowed; raising is refused and logged.
    try:
        _new = report['mandatory_list']['DIAGRAM'].get('crossing_number_min')
        # Read both sources: the stage-5 report of the PREVIOUS run and the record
        # (complete.json) — bound sweeps lower c only in complete.json, so the guard
        # must look there too, otherwise a re-run silently raises it again.
        _candidates = []
        for _p in (args.output,
                   args.output.replace('_mandatory_invariants.json',
                                       '_mandatory_invariants_complete.json')):
            if _p and os.path.exists(_p):
                try:
                    with open(_p) as _f:
                        _v = ((json.load(_f).get('mandatory_list') or {}).get('DIAGRAM') or {}) \
                            .get('crossing_number_min')
                    if isinstance(_v, int):
                        _candidates.append(_v)
                except Exception:
                    pass
        if isinstance(_new, int) and _candidates:
            _alt = min(_candidates)
            if _alt < _new:
                report['mandatory_list']['DIAGRAM']['crossing_number_min'] = _alt
                report['_C_BOUND_KEPT_FROM_PREVIOUS'] = {
                    'kept': _alt, 'this_run_found': _new,
                    'reason': ('c is an upper bound (Disc #231) — the smaller bound receipted '
                               'earlier stays valid; stage-4 sampling is stochastic, '
                               'a worse run must not raise it (Disc #562).'),
                    'diagram_fields': ('pd_code/dt/writhe/knot_sig belong to the diagram of THIS '
                                        'run — only the bound is the old one.'),
                    'when': '2026-07-25',
                }
                print(f"\n[PFLICHT] 🛡  c guard: this run found {_new}, record carried {_alt} "
                      f"→ keeping {_alt} (never raise, Disc #562).")
    except Exception as _e:
        print(f"[PFLICHT] ⚠ c guard skipped: {type(_e).__name__}: {_e}")

    with open(args.output, 'w') as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\n[PFLICHT] Writing {args.output}")

    # Pretty print
    print(f"\n{'='*60}")
    print(f"PFLICHTLISTE — {args.name}")
    print(f"{'='*60}")
    for category, items in report['mandatory_list'].items():
        print(f"\n=== {category} ===")
        for k, v in items.items():
            if v is None:
                continue
            s = str(v)
            if len(s) > 200:
                s = s[:197] + '...'
            print(f"  {k:>30} : {s}")

    # FALSE-GREEN FIX (error log 2026-06-04 ~22:40, approved by Werner 2026-06-05):
    # tr15e SKIP: regina died (rc=137) -> ALL mandatory fields None, but
    # stage 5 ran through to exit 0 -> pipeline_run.sh reported "[PIPELINE OK]".
    # Minimum population = at least ONE det witness (jones_at_i / sage /
    # seifert). Without it the report is empty and all later stages are hopeless.
    # The report has already been written above (diagnosis is preserved).
    # MANDATORY_ALLOW_EMPTY=1 = deliberate diagnostic run without this gate
    # (same one-line switch mechanism as REGINA_TIMEOUT).
    if not det_paths and os.environ.get('MANDATORY_ALLOW_EMPTY') != '1':
        print(f"\n[PFLICHT] ERROR: minimum population violated — no det witness "
              f"(regina_failed={report.get('regina_failed', False)}). "
              f"Report kept for diagnosis in {args.output}. Exit 4.")
        sys.exit(4)


if __name__ == '__main__':
    main()
