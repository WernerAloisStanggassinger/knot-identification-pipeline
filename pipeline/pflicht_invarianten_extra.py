#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
pflicht_invarianten_extra.py
============================
Computes the 6 mandatory invariants that `pflicht_invarianten.py` did not
cover:
  1. khovanov_homology   — KnotJob -kb0
  2. length_spectrum_short — snappy with canonize() fallback (Discovery: the Dirichlet construction fails on the default triangulation)
  3. s_invariant         — = 0 if amphicheiral (Discovery #75)
  4. bridge_number       — combinatorial from the signed Gauss code (counts maximal O blocks; upper bound)
  5. unknotting_number   — Sage brute-force braid with multiprocessing.Pool(12) (Discovery #60)
  6. kauffman_polynomial — Sage CHA Markov trace n≤4 (Discovery #55); for n>4 the KeyError is documented

Loads an existing pflicht_invarianten.json, adds the 6 fields,
writes to <name>_complete.json.

Full hardware use:
  - unknotting_number: Sage subprocess + multiprocessing.Pool(12) internally
  - others: single subprocesses (no parallelisation potential)
"""

import argparse
import site_config
import json
import os
import re
import subprocess
import sys
import textwrap
import time

# Disc #222 (2026-05-21): subprocess.run + conda run -> orphaned inner python
# on timeout. pipeline_subprocess.run_capture_with_group_kill kills the
# whole process group and avoids orphans.
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

KNOTJOB_JAR = site_config.KNOTJOB_JAR
JAVA_BIN    = site_config.JAVA_BIN


# ──────────────────────────────────────────────
# 1. KHOVANOV via KnotJob — parsed into structure
# ──────────────────────────────────────────────

def parse_pd_string(pd_str):
    crossings = re.findall(r'X[\(\[](\-?\d+(?:,\s*\-?\d+){3})[\)\]]', pd_str)
    if not crossings:
        return None
    return [tuple(int(v.strip()) for v in cr.split(',')) for cr in crossings]


def parse_khovanov_polynomial(poly_str):
    """Parse a KnotJob Khovanov polynomial of the form 'c1 t^a1 q^b1 + c2 t^a2 q^b2 + ...'.

    Handles uniformly:
      - missing coefficient ('t^2 q^3')           → coef = 1
      - missing exponent    ('t q' = t^1 q^1)     → exp  = 1
      - missing variable    ('14 q', '12 t', '5') → missing exp = 0
      - negative exponents  ('t^-5 q^-11')        — '^-' is NOT mistaken for a term separator,
                                                    because only a whitespace-flanked binary ' - '
                                                    is normalised to the term separator ' + -'.
    Returns: dict {(t_exp, q_exp): coef} (terms summing to zero removed).
    """
    s = poly_str.strip()
    if not s:
        return {}
    # Step 1: binary ' - ' (whitespace-flanked) becomes ' + -'.
    # '^-5' stays intact (no whitespace between '^' and '-').
    s = re.sub(r'\s+-\s+', ' + -', s)
    # Step 2: normalise the leading sign
    if not s.startswith(('+', '-')):
        s = '+ ' + s
    # Step 3: split at top-level ' + ' (with whitespace)
    parts = [p.strip() for p in re.split(r'\s+\+\s+', s) if p.strip()]
    # Step 4: per part, extract sign + coef + t part + q part
    terms = {}
    for part in parts:
        sign = 1
        if part.startswith('-'):
            sign = -1
            part = part[1:].strip()
        elif part.startswith('+'):
            part = part[1:].strip()
        coef = sign
        t_exp = 0
        q_exp = 0
        m_coef = re.match(r'(\d+)\b', part)
        if m_coef:
            coef = sign * int(m_coef.group(1))
            part = part[m_coef.end():].strip()
        m_t = re.match(r't(?:\^(\-?\d+))?', part)
        if m_t:
            t_exp = int(m_t.group(1)) if m_t.group(1) is not None else 1
            part = part[m_t.end():].strip()
        m_q = re.match(r'q(?:\^(\-?\d+))?', part)
        if m_q:
            q_exp = int(m_q.group(1)) if m_q.group(1) is not None else 1
        terms[(t_exp, q_exp)] = terms.get((t_exp, q_exp), 0) + coef
    return {k: v for k, v in terms.items() if v != 0}


def khovanov_via_knotjob(pd_code, name='knot'):
    pd_list = parse_pd_string(pd_code)
    if pd_list is None:
        return {'error': f'pd parse: {pd_code[:200]}'}
    all_ids = [aid for cr in pd_list for aid in cr]
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
    out = {'rc': res.returncode, 'raw_stdout': res.stdout}
    out.update(parse_knotjob_khovanov_stdout(res.stdout))
    # s invariant: KnotJob -s delivers it too, but needs a separate call
    return out


def s_invariant_via_knotjob_sgr(pd_code, name='knot', timeout=900):
    """Calls KnotJob -sgr and parses 'Graded s-invariant : N. M' (Discovery #111).

    Returns: dict {rc, raw_stdout, s_invariant (int|None), parse_error (str|None)}.
    On a crash on a 14+-crossing PD (Discovery #75) KnotJob throws
    IndexOutOfBoundsException — then rc != 0, s_invariant None.
    """
    pd_list = parse_pd_string(pd_code)
    if pd_list is None:
        return {'rc': -1, 'raw_stdout': '', 's_invariant': None,
                'parse_error': f'pd parse: {pd_code[:200]}'}
    all_ids = [aid for cr in pd_list for aid in cr]
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
            [JAVA_BIN, '-Xmx8g', '-jar', KNOTJOB_JAR, tmppath, '-sgr'],
            capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        try: os.unlink(tmppath)
        except: pass
        return {'rc': -2, 'raw_stdout': '', 's_invariant': None,
                'parse_error': f'timeout {timeout}s'}
    finally:
        try: os.unlink(tmppath)
        except: pass
    out = {'rc': res.returncode, 'raw_stdout': res.stdout,
           'stderr': res.stderr, 's_invariant': None, 'parse_error': None}
    if res.returncode != 0:
        out['parse_error'] = f'rc={res.returncode}; stderr={res.stderr[:300]}'
        return out
    # 'Graded s-invariant : -2. 2' — the first number is the s invariant.
    m = re.search(r'Graded s-invariant\s*:\s*(-?\d+)', res.stdout)
    if m:
        out['s_invariant'] = int(m.group(1))
    else:
        out['parse_error'] = f'no "Graded s-invariant" line in stdout: {res.stdout[-300:]}'
    return out


def parse_knotjob_khovanov_stdout(raw_stdout):
    """Parses KnotJob Khovanov output into structured fields.

    Recognises three polynomial blocks:
      - "Integral unreduced Khovanov Homology" → free-part terms_t_q
      - "Torsion of order P" → list per order
      - "Integral reduced Khovanov Homology" → reduced terms_t_q

    Returns a dict with the following possible keys (depending on the KnotJob output):
      polynomial, terms_t_q, delta_diagonals, n_distinct_diagonals, width, thin,
      torsion (list of dicts), torsion_orders, has_torsion, z2_thin_only,
      reduced_polynomial, reduced_terms_t_q, reduced_delta_diagonals,
      reduced_n_distinct_diagonals, reduced_width, reduced_thin
    """
    out = {}
    # Free part of Khovanov
    m = re.search(r'Integral unreduced Khovanov Homology\s*:\s*(.+?)$',
                  raw_stdout, re.MULTILINE)
    if m:
        polynomial_str = m.group(1).strip()
        out['polynomial'] = polynomial_str
        terms = parse_khovanov_polynomial(polynomial_str)
        out['terms_t_q'] = {f'{te},{qe}': c for (te, qe), c in terms.items()}
        # Standard definition Lee/Khovanov:
        #   δ = q − 2t (delta grading); width = max(δ) − min(δ)
        #   thin ⇔ width ≤ 2 (Lee 2005: alternating knots are always thin).
        if terms:
            diags = sorted(set(qe - 2 * te for (te, qe) in terms.keys()))
            out['delta_diagonals'] = diags
            out['n_distinct_diagonals'] = len(diags)
            out['width'] = diags[-1] - diags[0]
            out['thin']  = (out['width'] <= 2)
    # Torsion (can deliver several orders, e.g. order 2, order 3, ...)
    torsion_pattern = re.compile(r'Torsion of order\s+(\d+)\s*:\s*(.+?)$',
                                  re.MULTILINE)
    torsions = []
    for tm in torsion_pattern.finditer(raw_stdout):
        order = int(tm.group(1))
        tor_str = tm.group(2).strip()
        if not tor_str or tor_str.lower() == 'none':
            torsions.append({'order': order, 'polynomial': '', 'terms_t_q': {},
                             'total_torsion_count': 0})
            continue
        tor_terms = parse_khovanov_polynomial(tor_str)
        torsions.append({
            'order':               order,
            'polynomial':          tor_str,
            'terms_t_q':           {f'{te},{qe}': c for (te, qe), c in tor_terms.items()},
            'total_torsion_count': sum(abs(c) for c in tor_terms.values()),
        })
    out['torsion'] = torsions
    out['has_torsion']    = bool(torsions and any(t['total_torsion_count'] > 0 for t in torsions))
    out['torsion_orders'] = sorted({t['order'] for t in torsions if t['total_torsion_count'] > 0})
    # Z2-thin-only: standard property of alternating knots (Shumakovitch);
    # 'pure-thin' if NO torsion.
    out['z2_thin_only']   = (out['torsion_orders'] == [2])
    # Reduced Khovanov
    rm = re.search(r'Integral reduced Khovanov Homology\s*:\s*(.+?)$',
                   raw_stdout, re.MULTILINE)
    if rm:
        red_str = rm.group(1).strip()
        red_terms = parse_khovanov_polynomial(red_str)
        out['reduced_polynomial'] = red_str
        out['reduced_terms_t_q']  = {f'{te},{qe}': c for (te, qe), c in red_terms.items()}
        if red_terms:
            r_diags = sorted(set(qe - 2 * te for (te, qe) in red_terms.keys()))
            out['reduced_delta_diagonals']      = r_diags
            out['reduced_n_distinct_diagonals'] = len(r_diags)
            out['reduced_width']                = r_diags[-1] - r_diags[0]
            # reduced is 'thin' in the strict sense when concentrated on a single
            # diagonal (typical for alternating knots).
            out['reduced_thin']                 = (out['reduced_width'] == 0)
    return out


# ──────────────────────────────────────────────
# 2. LENGTH_SPECTRUM_SHORT — with canonize() fallback
# ──────────────────────────────────────────────

def length_spectrum_via_canonize(pd_code, n_geos=10, timeout=300, isosig=None):
    """Short length spectrum over several routes, with multi-route fallback.

    Route 1 (PRIMARY, Discovery #103): snappy.Manifold.length_spectrum_alt_gen()
        — generator based; a numerical method independent of the Dirichlet
        polyhedron, works also for high-symmetry manifolds such as Saturn.
        The first geodesic can be expensive, further ones come quickly.
    Route 2 (FALLBACK, Discovery #95): classical length_spectrum(0.5) via
        Dirichlet, with canonize() / randomize-x_k fallback. Fails for
        D_6-symmetric manifolds.

    Input: pd_code (string with "X[a,b,c,d], ...") and optionally isosig.
    If isosig is given, the manifold is built directly from the isosig
    (bypasses spherogram.Link.exterior() — relevant when the diagram PD
    takes a different triangulation route than the one persisted by the
    mandatory block).

    Geodesics are grouped by numerical proximity (1e-9); multiplicity
    = number of hits.
    """
    iso_str = isosig if isosig else ''
    code = SEAL_PRELUDE + textwrap.dedent(f"""
        import json, snappy, spherogram, time
        pd_str = {pd_code!r}
        iso    = {iso_str!r}
        N_GEOS = {n_geos}

        if iso:
            M = snappy.Manifold(iso)
        else:
            import re
            crossings = re.findall(r'X[\\(\\[](\\-?\\d+(?:,\\s*\\-?\\d+){{3}})[\\)\\]]', pd_str)
            pd_list = [tuple(int(v.strip()) for v in cr.split(',')) for cr in crossings]
            L = spherogram.Link(pd_list)
            M = L.exterior()

        out = {{'volume': float(M.volume()), 'isosig_used': bool(iso)}}

        # --- Laengenspektrum: EIN Schema ---
        def _lspec_eintrag(g):
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
        def cluster_geodesics(geos, tol=1e-9):
            '''Group geodesics by complex proximity; multiplicity = cluster size.
            Conjugate pairs (+/-Im) are treated as DIFFERENT geodesics.'''
            clusters = []
            for g in geos:
                L = complex(g.length)
                placed = False
                for cl in clusters:
                    if abs(L - cl['L_repr']) < tol:
                        cl['count'] += 1
                        placed = True
                        break
                if not placed:
                    clusters.append({{'L_repr': L, 'count': 1,
                                      'real': L.real, 'imag': L.imag}})
            return clusters

        # Pfad 1: length_spectrum_alt_gen (numerisch, robuster)
        t0 = time.time()
        try:
            gen = M.length_spectrum_alt_gen()
            geos = []
            for i, g in enumerate(gen):
                geos.append(g)
                if len(geos) >= N_GEOS:
                    break
            clusters = cluster_geodesics(geos)
            out['method'] = 'length_spectrum_alt_gen'
            out['length_spectrum_short'] = [
                {{'length_real':       cl['real'],
                  'length_imag':       cl['imag'],
                  'multiplicity':      cl['count']}}
                for cl in clusters
            ]
            out['n_geodesics_raw'] = len(geos)
            out['duration_s']      = time.time() - t0
        except Exception as e_alt:
            out['err_alt'] = str(e_alt)
            # path 2: classical length_spectrum with Dirichlet
            try:
                ls = M.length_spectrum(0.5)
                out['method'] = 'direct_dirichlet'
                out['length_spectrum_short'] = [
                    _lspec_eintrag(g)
                    for g in ls[:N_GEOS]
                ]
                out['duration_s'] = time.time() - t0
            except Exception as e_dir:
                out['err_direct'] = str(e_dir)
                # path 3: after randomize() (short loop)
                for k in range(3):
                    try:
                        Mr = M.copy()
                        for _ in range(k+1):
                            Mr.randomize()
                        ls = Mr.length_spectrum(0.5)
                        out['method'] = f'randomize_x{{k+1}}'
                        out['length_spectrum_short'] = [
                            _lspec_eintrag(g)
                            for g in ls[:N_GEOS]
                        ]
                        out['duration_s'] = time.time() - t0
                        break
                    except Exception as e_rand:
                        out[f'err_rand_{{k+1}}'] = str(e_rand)
        print(json.dumps(out, default=str))
    """)
    cmd = ['conda', 'run', '-n', 'snappy_env', '--no-capture-output',
           'python', '-c', code]
    res = _run_with_pgkill(cmd, timeout=timeout)
    if res.returncode != 0:
        return {'error': f'rc={res.returncode}', 'stderr': res.stderr[-1000:]}
    for line in reversed(res.stdout.strip().splitlines()):
        line = line.strip()
        if line.startswith('{') and line.endswith('}'):
            try: return json.loads(line)
            except: continue
    return {'error': 'no json', 'stdout': res.stdout[-1000:]}


# ──────────────────────────────────────────────
# 2b. MFW LOWER BOUND for braid_index from HOMFLY
# ──────────────────────────────────────────────

def mfw_braid_index_bound(homfly_str):
    """Morton-Franks-Williams lower bound for braid_index(K) from HOMFLY P(v,z).

    MFW theorem: b(K) >= (max v-exp - min v-exp)/2 + 1, i.e. half the span
    of the exponents of the first HOMFLY variable plus 1.

    Accepts both HOMFLY conventions automatically:
      - Regina: variables 'x' (=v) and 'y' (=z), format 'x^6 - 3 x^4 y^2 + ...'
      - Sage:   variables 'L' (=v) and 'M' (=z), format '-L^6 + 3*L^4*M^2 + ...'
    The exponent spread is invariant under the sign convention v vs v^{-1}.

    Returns: dict with keys 'lower_bound', 'max_v_exp', 'min_v_exp', 'spread',
             'variable_used', or None on a parse error.
    """
    if not homfly_str:
        return None
    # Detect the variable: 'L' (Sage) takes precedence over 'x' (Regina), because 'L'
    # is unambiguous and 'x' could also occur in a word like 'exp'.
    var = 'L' if re.search(r'\bL\b', homfly_str) else 'x'
    s = homfly_str.replace('*', ' ')
    s = re.sub(r'\s+-\s+', ' + -', s.strip())
    if not s.startswith(('+', '-')):
        s = '+ ' + s
    parts = [p.strip() for p in re.split(r'\s+\+\s+', s) if p.strip()]

    var_pat = re.compile(rf'(?:^|[^a-zA-Z0-9]){var}(?:\s*\^\s*(\-?\d+))?')
    v_exps = []
    for part in parts:
        if part.startswith(('+', '-')):
            part = part[1:].strip()
        m = var_pat.search(part)
        if m:
            v_exps.append(int(m.group(1)) if m.group(1) is not None else 1)
        else:
            # Constant term or only the z variable → v exp = 0
            v_exps.append(0)
    if not v_exps:
        return None
    Emax = max(v_exps)
    Emin = min(v_exps)
    return {
        'lower_bound': (Emax - Emin) // 2 + 1,
        'max_v_exp':   Emax,
        'min_v_exp':   Emin,
        'spread':      Emax - Emin,
        'variable_used': var,
    }


# ──────────────────────────────────────────────
# 3. BRIDGE NUMBER — Schubert bounds from genus + crossing_number
# ──────────────────────────────────────────────

def bridge_number_bounds(genus, crossing_number, det=None, alternating=False):
    """Bridge number bounds via standard inequalities.

    Lower bound (CORRECTED — Discovery #109):
      Formerly here: 'Schubert b(K) >= g(K)+1'. That is NOT a general theorem
      (counterexample: (5,2) torus knot g=2, b=2). Current lower bound:
        - unknot (Seifert genus 0) ⇒ b = 1
        - det == ±1 (but genus > 0) ⇒ b >= 3, see below
        - otherwise: b >= 2 (nontrivial knot)
      Sharpenings via the Wirtinger meridional_rank come in
      enrich_optional_invariants.py wave 5 (generalises Discovery #103).

    Upper bound (trivial): b(K) <= c(K). From any diagram.

    🔴 CORRECTION 2026-08-10 (finding 15 from the Galway laptop, cause found here):
      This used to read 'det == ±1  ⇒  unknot ⇒ interval [1, 1]'. **That is wrong.**
      det = |Δ(-1)|, and there are plenty of nontrivial knots with Δ ≡ 1 or
      Δ(-1) = ±1 — in our OWN corpus: `conway_cw` (Conway, g=3), `conway_kt`
      (Kinoshita-Terasaka, g=2), `mi_T3_5` = T(3,5) (g=4), `wz_sqrt19`/`wz_sqrt22`
      (g=6) and **`irrational` = T(3,7), our anchor 3 (g=6)**. For T(p,q) with
      p,q both odd, det is always 1 — that hits a whole family, not
      outliers. The helper would have declared every one of these knots an unknot;
      in fact it hit 6 records (cluster S_18_d1_D3_chiral, g=4), because
      only there was `bridge_number` enriched at all. It showed up
      as **interval [3, 1]** — lower bound greater than upper, i.e. a
      number that cannot exist: the Wirtinger sharpening (wave 5) had
      set the 3 correctly, while the 1 came from the false unknot assumption.

      The correct criterion is the **Seifert genus**: g(K) = 0 ⟺ K is the
      unknot. That is a theorem, not a heuristic, and `genus` is in the
      signature anyway — no new wiring needed.

      det = ±1 is useless for that, but it carries the REVERSE statement, and that one
      is sharp: a nontrivial 2-bridge knot K(p, q) has det = p >= 3.
      So det = ±1 **excludes** two bridges ⇒ b(K) >= 3. The error thus
      turns into a sharpening — and it agrees independently with what
      Wirtinger had found for exactly these 6 records (also 3).

    Returns a dict with:
      lower_bound, upper_bound, interval, is_2_bridge_excluded,
      lower_bound_source, upper_bound_source, plus an optional Murasugi note.
    """
    g = int(genus) if genus is not None else None
    c = int(crossing_number) if crossing_number is not None else None

    # Unknot — and ONLY the genus decides that (g(K) = 0 ⟺ unknot).
    if g is not None and g == 0:
        return {
            'lower_bound': 1,
            'upper_bound': 1,
            'interval':    [1, 1],
            'is_2_bridge_excluded': True,  # the unknot is 1-bridge, not 2-bridge
            'lower_bound_source': 'unknot (Seifert genus 0 ⟺ unknot)',
            'upper_bound_source': 'unknot (Seifert genus 0 ⟺ unknot)',
        }

    det1 = det is not None and abs(int(det)) == 1
    if det1:
        lower = 3
        lower_src = ('det = ±1 excludes 2 bridges ⇒ b(K) >= 3 '
                     '(a nontrivial 2-bridge knot K(p,q) has det = p >= 3). '
                     '⚠️ det = ±1 does NOT mean unknot — counterexamples in our own '
                     'corpus: conway_cw, conway_kt, mi_T3_5, irrational = T(3,7).')
    else:
        lower = 2
        lower_src = ('nontrivial knot (Genus > 0) ⇒ b(K) >= 2 '
                     '(Discovery #109: old "Schubert b>=g+1" rule removed)')
    upper = c if c is not None else None
    out = {
        'lower_bound': lower,
        'upper_bound': upper,
        'interval':    [lower, upper] if upper is not None else None,
        # Honestly: without det=±1 we cannot exclude 2-bridge.
        # Refinements via enrich_optional_invariants wave 5.
        'is_2_bridge_excluded': bool(det1),
        'lower_bound_source': lower_src,
        'upper_bound_source': f'trivial: b(K) <= c(K) = {c}' if c is not None else None,
    }
    # Murasugi/Sakuma: for ALTERNATING knots additionally
    # b(K) <= ceil((c+1)/2). This inequality is often quoted in knot theory,
    # but strongly heuristic — we carry it as a
    # PROBABLE upper bound with the (alt_*_) prefix, NOT as a proven one.
    if alternating and c is not None:
        alt_upper = (c + 1) // 2 + ((c + 1) % 2)  # ceil((c+1)/2)
        out['alt_upper_bound_heuristic']        = alt_upper
        out['alt_upper_bound_heuristic_source'] = (
            f'Murasugi/Sakuma heuristic for alternating: b(K) <= ceil((c+1)/2) = {alt_upper}; '
            f'NOT proven in this pipeline, only a plausibility hint')
    return out


def bridge_count_from_signed_gauss(sg):
    """Counts the cyclically connected 'O' blocks = bridges
    in THIS diagram. That is an upper bound for the
    bridge_number of the knot (= min over all diagrams).

    Bridge = a connected strand arc that runs over other strands
    the whole time. In the signed Gauss code: a sequence of 'O…' tokens
    without an intervening 'U…'."""
    types = [t[0] for t in re.findall(r'([OU])\d+[+-]', sg)]
    n = len(types)
    if n == 0:
        return 0
    if 'O' not in types or 'U' not in types:
        return 1
    # Cyclic: count U → O transitions
    bridges = sum(1 for i in range(n)
                  if types[i] == 'O' and types[(i - 1) % n] == 'U')
    return bridges


# ──────────────────────────────────────────────
# ──────────────────────────────────────────────
# STAGE-6 GUARD (approved by Werner 2026-08-23: "build the guard in")
# ──────────────────────────────────────────────
# OCCASION, measured the same day: `rosette_C8_b155` (TARGET c=308) ran TWENTY
# MINUTES against the memory wall — nine earlyoom kills at ~15.8 GB each, each
# time retried, while the box stood at 1 GB free and everything else starved.
# The existing brake in pipeline_run.sh lowers the NUMBER of workers (c>36 ⟹ 2),
# not their SIZE — 2 x 15 GB blow 32 Gi just the same.
#
# TWO LAYERS, both MEASURED rather than set (count over 474 records with
# unknotting_number on 23.08.):
#   · unknotting_number was NEVER obtained above c=64 (largest five: 32,33,44,60,64)
#     ⟹ UNKNOT_C_MAX = 80 lies comfortably above; nothing that historically
#       succeeded is lost, but c=308 is rejected in milliseconds.
#   · largest successful braid length: 201 (median 15, 90% point 53)
#     ⟹ UNKNOT_L_MAX = 260.
#   · STAGE6_MEM_GB caps the ADDRESS SPACE of the spawned process; the
#     pool workers INHERIT the limit, so it applies per worker.
#
# ⚠️ HONEST ABOUT THE RESIDUAL RISK (Footgun 13): against C/C++ cores RLIMIT_AS acts as a
# CRASH (std::bad_alloc ⟹ SIGABRT), not as a MemoryError. A pool worker killed by a
# signal can make `pool.map` hang. Exactly for that reason the
# THRESHOLD is the first layer and the limit only the second: anything above the threshold
# is not even started.
UNKNOT_C_MAX  = int(os.environ.get('UNKNOT_C_MAX', 80))
UNKNOT_L_MAX  = int(os.environ.get('UNKNOT_L_MAX', 260))
STAGE6_MEM_GB = float(os.environ.get('STAGE6_MEM_GB', 6))

MEM_PRELUDE = (
    "import resource as _res, os as _os\n"
    "_gb = float(_os.environ.get('STAGE6_MEM_GB', %r))\n"
    "try:\n"
    "    _b = int(_gb * 1024**3)\n"
    "    _res.setrlimit(_res.RLIMIT_AS, (_b, _b))\n"
    "except Exception:\n"
    "    pass\n"
) % STAGE6_MEM_GB


def _unknot_zu_gross(mass, wert, grenze, was):
    """Labelled rejection instead of a silent gap — an empty cell looks
    like a finding and is none (Disc #674/#680)."""
    return {'unknotting_number': None,
            'unknotting_uebersprungen': True,
            'grund': (f'{was} {wert} > {grenze} ({mass}). The brute force is '
                      f'exponential; measured over 474 records, '
                      f'unknotting_number was never obtained above c=64 or braid length 201. '
                      f'Above the threshold it is not computed '
                      f'but rejected — otherwise the run hits the '
                      f'memory wall (23.08.: 20 min, 9 earlyoom kills).'),
            'schwelle': grenze, 'gemessen': wert,
            'nachholbar_mit': f'{mass}=<larger> as environment variable'}


# 4. UNKNOTTING NUMBER — Sage brute-force braid (Discovery #60)
# ──────────────────────────────────────────────

def unknotting_number_brute_force(braid_word, n_strands, max_subset=4,
                                  processes=12):
    """Brute force over all subset sizes 0, 1, …, max_subset:
    per subset: sign flip at the subset positions in the braid word,
    construct the Knot, check whether trivial (det == 1 ∧ Jones polynomial == 1).
    Distributed over multiprocessing.Pool(processes).

    Scales as 2^n with n = braid_length, but we go through subset sizes
    in ascending order → early abort."""
    if braid_word is not None and len(braid_word) > UNKNOT_L_MAX:
        return _unknot_zu_gross('UNKNOT_L_MAX', len(braid_word), UNKNOT_L_MAX,
                                'braid length')
    code = SEAL_PRELUDE + MEM_PRELUDE + textwrap.dedent(f"""
        import json, sys
        from sage.all import BraidGroup, Knot
        from itertools import combinations
        from multiprocessing import Pool

        braid_word = {braid_word!r}
        n_strands  = {n_strands}
        max_subset = {max_subset}
        N_proc     = {processes}

        B = BraidGroup(n_strands)

        def is_trivial_after_flip(subset):
            modified = list(braid_word)
            for k in subset:
                modified[k] = -modified[k]
            try:
                b = B(modified)
                K = Knot(b)
                # Test 1: Det == 1
                if abs(int(K.determinant())) != 1:
                    return False
                # Test 2: Jones polynomial = 1
                jp = K.jones_polynomial()
                # trivial Jones polynomial (unknot) is 1
                if str(jp).strip() in ('1', 't^0'):
                    return True
                return False
            except Exception:
                return False

        L = len(braid_word)
        print(json.dumps({{'msg': f'Brute-force over braid length {{L}}, subsets k=1..{{max_subset}} (k=0 is decided by the caller from det/Jones in the pflicht JSON)'}}), flush=True)

        u_found = None
        for k in range(1, max_subset + 1):
            subsets = list(combinations(range(L), k))
            if not subsets:
                continue
            with Pool(N_proc) as pool:
                results = pool.map(is_trivial_after_flip, subsets, chunksize=max(1, len(subsets)//(N_proc*8)))
            n_unk = sum(results)
            print(json.dumps({{'k': k, 'n_subsets': len(subsets), 'unknot_count': int(n_unk)}}), flush=True)
            if n_unk > 0:
                u_found = k
                break

        out = {{'unknotting_number': u_found,
                'searched_up_to_subset_size': max_subset,
                'braid_length': L,
                'n_strands': n_strands}}
        if u_found is None:
            out['note'] = f'u > {{max_subset}} or not found — raise max_subset'
        print(json.dumps(out), flush=True)
    """)
    cmd = ['conda', 'run', '-n', 'sage_stable', '--no-capture-output',
           'sage', '-c', code]
    res = _run_with_pgkill(cmd, timeout=3600)
    if res.returncode != 0:
        return {'error': f'rc={res.returncode}', 'stderr': res.stderr[-1500:]}
    out = {'progress': []}
    for line in res.stdout.splitlines():
        line = line.strip()
        if line.startswith('{') and line.endswith('}'):
            try:
                d = json.loads(line)
                if 'unknotting_number' in d:
                    out.update(d)
                else:
                    out['progress'].append(d)
            except Exception:
                continue
    return out


# ──────────────────────────────────────────────
# 4b. UNKNOTTING NUMBER via signed Gauss crossing flip (Regina)
# ──────────────────────────────────────────────

def unknotting_via_signed_gauss_flip(signed_gauss, max_subset=4, processes=12):
    """Brute-force upper bound for u(K) on the MIN-CROSSING diagram via signed Gauss.

    Topologically clean crossing change in signed Gauss format:
      at crossing i swap O <-> U AND + <-> -.
    That is the standard crossing change (Discovery #59 concerns PD label
    rotation, not the signed Gauss flip — the latter preserves classicality,
    see smoke tests of 2026-05-09).

    Triviality test: Regina simplify to size==0 ⇒ unknot. Conservative
    (false negatives possible), but sufficient in most cases for subset flips
    of an existing min-crossing diagram.

    Advantage over the braid brute force (Discovery #93):
      Saturn diagram 12 crossings vs. Sage braid 9 strands × 40 letters.
      2^12 = 4096 subsets in total, with early abort typically <100.
    """
    crossings = re.findall(r'([OU])(\d+)([+-])', signed_gauss)
    if not crossings:
        return {'error': 'parse signed_gauss', 'sg': signed_gauss[:200]}
    indices = sorted({int(idx) for _, idx, _ in crossings})
    n_crossings = len(indices)
    if n_crossings > UNKNOT_C_MAX:
        return _unknot_zu_gross('UNKNOT_C_MAX', n_crossings, UNKNOT_C_MAX,
                                'diagram size')
    code = SEAL_PRELUDE + MEM_PRELUDE + textwrap.dedent(f"""
        import json, re, sys
        from itertools import combinations
        from multiprocessing import Pool
        import regina

        signed_gauss = {signed_gauss!r}
        indices      = {indices!r}
        max_subset   = {max_subset}
        N_proc       = {processes}

        def crossing_change_sg(sg, flip_set):
            out = []
            for m in re.finditer(r'([OU])(\\d+)([+-])', sg):
                typ, idx, sgn = m.group(1), m.group(2), m.group(3)
                if int(idx) in flip_set:
                    typ = 'O' if typ == 'U' else 'U'
                    sgn = '+' if sgn == '-' else '-'
                out.append(typ + idx + sgn)
            return ''.join(out)

        def trivial_test(subset):
            '''Returns 0 = NOT-trivial (det != 1), 1 = simplify→0 (certainly unknot),
                       2 = det==1 ∧ Jones==1 but Regina simplify > 0 (unknot candidate).
            '''
            try:
                sg = crossing_change_sg(signed_gauss, set(subset))
                L = regina.Link.fromSignedGauss(sg)
                if not L.isClassical():
                    return 0
                # det check before simplify (det computation is faster)
                jones = str(L.jones(regina.Algorithm.Treewidth))
                # Regina Jones at x=i = ±det(K), without a Sage round trip:
                # det(K) = |jones(-1)| in the usual convention, but Regina uses x.
                # Instead: simplify and the larger tests there
                L.simplify()
                if L.size() == 0:
                    return 1
                # diagram still has crossings; test via Jones==1
                jp = str(L.jones(regina.Algorithm.Treewidth)).strip()
                if jp in ('1', 't^0'):
                    return 2
                return 0
            except Exception:
                return 0

        L0 = regina.Link.fromSignedGauss(signed_gauss)
        L0.simplify()
        n0 = L0.size()
        print(json.dumps({{'msg': f'Brute-force signed-gauss-flip on {{n0}}-crossing diagram, k=1..{{max_subset}}'}}), flush=True)

        u_found = None
        u_witness = None
        for k in range(1, max_subset + 1):
            subsets = list(combinations(indices, k))
            with Pool(N_proc) as pool:
                results = pool.map(trivial_test, subsets, chunksize=max(1, len(subsets)//(N_proc*8)))
            n_simp     = sum(1 for r in results if r == 1)
            n_jones1   = sum(1 for r in results if r == 2)
            print(json.dumps({{'k': k, 'n_subsets': len(subsets),
                               'unknot_simplify_zero': n_simp,
                               'unknot_jones_only':    n_jones1}}), flush=True)
            if n_simp > 0 or n_jones1 > 0:
                u_found = k
                # Witness: bevorzugt simplify→0, sonst jones==1
                ex_idx = next(i for i, r in enumerate(results) if r == 1) if n_simp > 0 \
                         else next(i for i, r in enumerate(results) if r == 2)
                u_witness = list(subsets[ex_idx])
                print(json.dumps({{'witness_subset': u_witness,
                                   'witness_type': 'simplify_zero' if n_simp > 0 else 'jones_only'}}), flush=True)
                break

        out = {{'unknotting_number_upper_bound': u_found,
               'searched_up_to_subset_size': max_subset,
               'method': 'signed_gauss_crossing_flip',
               'n_crossings_diagram': n0,
               'witness_subset': u_witness}}
        if u_found is None:
            out['note'] = f'no subset up to k={{max_subset}} leads to the unknot via Regina simplify'
        print(json.dumps(out), flush=True)
    """)
    cmd = ['conda', 'run', '-n', 'base', '--no-capture-output',
           'python', '-c', code]
    # Fix 2026-06-04 (Disc #255 lesson A, approved by Werner): catch the timeout instead of letting
    # the stage die — items (0)-(4) were finished in the kf286 run and were lost
    # with the uncaught exception (rc=3, failure move).
    try:
        res = _run_with_pgkill(cmd, timeout=1800)
    except subprocess.TimeoutExpired:
        return {'error': 'timeout', 'timeout_seconds': 1800,
                'note': 'unknotting search aborted after 1800s; result = timeout, no finding'}
    if res.returncode != 0:
        return {'error': f'rc={res.returncode}', 'stderr': res.stderr[-1500:]}
    out = {'progress': []}
    for line in res.stdout.splitlines():
        line = line.strip()
        if line.startswith('{') and line.endswith('}'):
            try:
                d = json.loads(line)
                if 'unknotting_number_upper_bound' in d:
                    out.update(d)
                elif 'witness_subset' in d:
                    out['witness_subset'] = d['witness_subset']
                else:
                    out['progress'].append(d)
            except Exception:
                continue
    return out


# ──────────────────────────────────────────────
# 5. KAUFFMAN POLYNOMIAL via CHA Markov trace (Sage)
# ──────────────────────────────────────────────

def kauffman_via_cha(braid_word, n_strands):
    """Attempt via the CHA Markov trace (Discovery #55). Marin data only for n ≤ 4."""
    if n_strands > 4:
        return {'error': f'n_strands={n_strands} > 4 — Marin data only go up to n=4 '
                          f'(Discovery #55, KeyError {n_strands})',
                'unavailable_in_pipeline': True}
    code = SEAL_PRELUDE + textwrap.dedent(f"""
        import json
        from sage.all import CubicHeckeAlgebra, BraidGroup
        try:
            CHA = CubicHeckeAlgebra({n_strands})
            B = BraidGroup({n_strands})
            b = B({braid_word!r})
            chab = CHA(b)
            mt = chab.formal_markov_trace()
            kauff = mt.specialize_kauffman()
            print(json.dumps({{'kauffman_polynomial': str(kauff)}}))
        except Exception as e:
            print(json.dumps({{'error': str(e)}}))
    """)
    cmd = ['conda', 'run', '-n', 'sage_stable', '--no-capture-output',
           'sage', '-c', code]
    res = _run_with_pgkill(cmd, timeout=600)
    if res.returncode != 0:
        return {'error': f'rc={res.returncode}', 'stderr': res.stderr[-1000:]}
    for line in reversed(res.stdout.strip().splitlines()):
        line = line.strip()
        if line.startswith('{') and line.endswith('}'):
            try: return json.loads(line)
            except: continue
    return {'error': 'no json'}


def kauffman_via_knotpy(pd_code_str):
    """Fallback: Kauffman 2-variable polynomial via knotpy state sum
    (Discovery #106, Discovery #114a). Applies when n_strands > 4 or
    on a CHA error. PD code format as from Regina:
    'PD[X[a,b,c,d], X[...], ...]'.
    Convention note: knotpy's state sum can deliver the mirror polynomial
    (a -> a^-1) relative to KnotInfo. For amphichiral
    knots this does not matter.
    """
    code = textwrap.dedent(f"""
        import json, re
        try:
            import knotpy
            pd_str = {pd_code_str!r}
            inner = pd_str.strip()
            if inner.startswith('PD['):
                inner = inner[3:]
                if inner.endswith(']'):
                    inner = inner[:-1]
            xings = re.findall(r'X\\[([\\d,\\s]+)\\]', inner)
            if not xings:
                print(json.dumps({{'error': 'no crossings parsed', 'pd_head': pd_str[:200]}}))
            else:
                crossings = [tuple(int(x.strip()) for x in c.split(',')) for c in xings]
                K = knotpy.from_pd_notation(crossings)
                F = knotpy.kauffman(K)
                s = str(F)
                amphi_check = None
                try:
                    import sympy
                    a, z = sympy.symbols('a z')
                    F_mirror = F.subs(a, 1/a)
                    diff = sympy.simplify(F - F_mirror)
                    amphi_check = (diff == 0)
                except Exception:
                    amphi_check = None
                print(json.dumps({{
                    'kauffman_polynomial': s,
                    'kauffman_polynomial_n_chars': len(s),
                    'kauffman_polynomial_n_crossings_input': len(crossings),
                    'kauffman_polynomial_method': 'knotpy.kauffman state-sum (Discovery #106 / #114a Fallback)',
                    'kauffman_polynomial_amphi_symmetry_check': amphi_check,
                    'kauffman_polynomial_convention_note': 'knotpy state-sum; for chiral knots possibly the mirror relative to KnotInfo (Discovery #106)',
                }}))
        except Exception as e:
            import traceback
            print(json.dumps({{'error': str(e), 'traceback': traceback.format_exc()[-1500:]}}))
    """)
    # Disc #222/#286 fix (Werner 2026-06-15): NOT 'conda run' (starts the inner Python in its own
    # process group -> pgkill/timeout does not apply -> 2^c orphan hangs for 37 min). Direct base Python =
    # the child, pgkill kills it, the 600s timeout WORKS -> Kauffman always attempted, cleanly killed
    # where too slow, NO hang. (base = miniforge3 root, has knotpy.)
    cmd = [site_config.PYTHON_BASE, '-c', code]
    try:
        res = _run_with_pgkill(cmd, timeout=600)
    except subprocess.TimeoutExpired:
        return {'error': 'knotpy_timeout_600s'}
    if res.returncode != 0:
        return {'error': f'rc={res.returncode}', 'stderr': res.stderr[-1500:]}
    for line in reversed(res.stdout.strip().splitlines()):
        line = line.strip()
        if line.startswith('{') and line.endswith('}'):
            try: return json.loads(line)
            except: continue
    return {'error': 'no json', 'stdout': res.stdout[-1000:]}


# ──────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────

def main():
    # Hole B (OOM 2026-07-01): pass an external kill through to the conda/sage pool grandchildren,
    # otherwise the multiprocessing.Pool workers survive as RAM-eating orphans.
    install_signal_killer()
    parser = argparse.ArgumentParser()
    parser.add_argument('input_json', help='pflicht_invarianten.json to extend')
    parser.add_argument('--output', default=None,
                        help='Output path (default: input with _complete suffix)')
    parser.add_argument('--max-unknotting-subset', type=int, default=4)
    parser.add_argument('--processes', type=int, default=12)
    args = parser.parse_args()

    # OOM cap (env override, backward compatible: unset -> the --processes arg stays).
    # The generated regina pool (line ~571: fromSignedGauss/homfly/jones/pd per t_step)
    # grew to 5-9GB/worker for densified high-c knots -> 12 workers = OOM.
    # PFLICHT_EXTRA_PROCESSES throttles ONLY this multiprocessing stage (hard-rule-6
    # exception = documented OOM path). Set by r3_reresolve.py via R3_PFLICHT_EXTRA.
    _env_pe = os.environ.get("PFLICHT_EXTRA_PROCESSES")
    if _env_pe:
        args.processes = max(1, int(_env_pe))

    with open(args.input_json) as f:
        report = json.load(f)

    output = args.output or args.input_json.replace('.json', '_complete.json')

    # Input values
    a    = report['pflichtliste']['ALGEBRA']
    d    = report['pflichtliste']['DIAGRAMM']
    h    = report['pflichtliste']['HYPERBOLISCH']
    sg   = report['signed_gauss']
    pd   = d.get('pd_code')
    if not pd:
        print("ERROR: pd_code missing in the input.", file=sys.stderr)
        sys.exit(2)
    braid_word = d.get('braid_word')
    if isinstance(braid_word, str):
        braid_word = json.loads(braid_word)
    n_strands  = d.get('braid_n_strands') or d.get('braid_index')

    print(f"[EXTRA] Input: {args.input_json}")
    print(f"[EXTRA]   crossing_number_min = {d.get('crossing_number_min')}")
    print(f"[EXTRA]   alternating         = {d.get('alternating')}")
    print(f"[EXTRA]   amphicheiral        = {h.get('is_amphicheiral')}")
    print(f"[EXTRA]   braid: n={n_strands}, len={d.get('braid_length')}")

    # 0. MFW lower bound for braid_index from HOMFLY (Discovery #99)
    print(f"\n[EXTRA] (0) braid_index — MFW lower bound from HOMFLY …")
    homfly = a.get('homfly_polynomial')
    if homfly and n_strands:
        mfw = mfw_braid_index_bound(homfly)
        if mfw:
            d['braid_index_lower_bound']    = mfw['lower_bound']
            d['braid_index_upper_bound']    = n_strands
            d['braid_index_interval']       = [mfw['lower_bound'], n_strands]
            d['braid_index_mfw_v_spread']   = mfw['spread']
            d['braid_index_mfw_note']       = (
                f"MFW theorem: b(K) >= {mfw['lower_bound']} from the HOMFLY v-spread "
                f"{mfw['spread']} (v-exp ∈ [{mfw['min_v_exp']}, {mfw['max_v_exp']}]); "
                f"upper bound = current Sage braid strands = {n_strands} "
                f"(NOT guaranteed minimal — Markov reduction is NP-hard, "
                f"Discovery #56). If lower_bound == upper_bound, b(K) "
                f"is determined exactly.")
            print(f"  v-exp-Spread = {mfw['spread']}  (variable '{mfw['variable_used']}')")
            print(f"  Interval: {mfw['lower_bound']} <= b(K) <= {n_strands}")
            if mfw['lower_bound'] == n_strands:
                print(f"  → b(K) = {n_strands} EXACT (MFW bound attained sharply)")
        else:
            print(f"  MFW parsing failed on homfly_polynomial")
    else:
        print(f"  Skipped (homfly or n_strands missing)")

    # 1. Bridge number — trivial bounds (Discovery #109: old g+1 rule removed)
    print(f"\n[EXTRA] (1) bridge_number — trivial bounds (sharpening in enrich wave 5) …")
    g_seif = a.get('genus') or a.get('knot_floer_homology', {}).get('seifert_genus')
    c_min  = d.get('crossing_number_min')
    det    = a.get('determinant_sage') or a.get('determinant')
    altern = d.get('alternating', False)
    bnb = bridge_number_bounds(genus=g_seif, crossing_number=c_min,
                                det=det, alternating=altern)
    d['bridge_number_lower_bound']  = bnb['lower_bound']
    d['bridge_number_upper_bound']  = bnb['upper_bound']
    d['bridge_number_interval']     = bnb['interval']
    d['bridge_number_is_2_excluded'] = bnb['is_2_bridge_excluded']
    d['bridge_number_lower_source'] = bnb['lower_bound_source']
    d['bridge_number_upper_source'] = bnb['upper_bound_source']
    if 'alt_upper_bound_heuristic' in bnb:
        d['bridge_number_alt_heuristic']        = bnb['alt_upper_bound_heuristic']
        d['bridge_number_alt_heuristic_source'] = bnb['alt_upper_bound_heuristic_source']
    # Keep the diagram heuristic for illustration (cyclic O blocks);
    # it is trivially = c for alternating knots and thus not meaningful
    # (Discovery #93). Only an illustrative value.
    d['bridge_count_diagram_heuristic'] = bridge_count_from_signed_gauss(sg)
    print(f"  Lower:           b(K) >= {bnb['lower_bound']} ({bnb['lower_bound_source']})")
    print(f"  Trivial upper:   b(K) <= c   = {c_min}")
    print(f"  Interval:        [{bnb['lower_bound']}, {bnb['upper_bound']}]")
    print(f"  2-bridge excluded: {bnb['is_2_bridge_excluded']}")
    print(f"  → sharpening via Wirtinger meridional_rank in enrich_optional_invariants wave 5")
    if 'alt_upper_bound_heuristic' in bnb:
        print(f"  Alt heuristic (alternating): b(K) <= ceil((c+1)/2) = "
              f"{bnb['alt_upper_bound_heuristic']} (NOT proven)")

    # 2. s invariant: amphicheiral ⇒ 0 (Disc #75); otherwise KnotJob -sgr (Disc #111)
    print(f"\n[EXTRA] (2) s_invariant …")
    if h.get('is_amphicheiral'):
        a['s_invariant'] = 0
        a['s_invariant_note'] = 'trivially 0 because of amphichirality (Discovery #75)'
        a['s_invariant_method'] = 'amphicheiral_symmetry'
        print(f"  s_invariant = 0 (amphicheiral, Discovery #75)")
    else:
        t0 = time.time()
        sres = s_invariant_via_knotjob_sgr(pd, name=report.get('name', 'knot'))
        report['raw_subprocess_outputs']['knotjob_sgr'] = sres
        if sres.get('s_invariant') is not None:
            a['s_invariant'] = sres['s_invariant']
            a['s_invariant_note'] = (
                f"KnotJob -sgr ({JAVA_BIN}): '"
                f"{sres['raw_stdout'].strip().splitlines()[-1] if sres.get('raw_stdout') else '?'}'"
            )
            a['s_invariant_method'] = 'knotjob_sgr'
            print(f"  → {time.time()-t0:.1f}s")
            print(f"  s_invariant = {sres['s_invariant']} (KnotJob -sgr)")
        else:
            a['s_invariant'] = None
            a['s_invariant_note'] = (
                f"KnotJob -sgr failed: {sres.get('parse_error', 'unknown')}"
            )
            a['s_invariant_method'] = 'knotjob_sgr_failed'
            print(f"  → {time.time()-t0:.1f}s")
            print(f"  s_invariant = None (KnotJob -sgr error: {sres.get('parse_error','?')})")

    # 3. Khovanov via KnotJob
    print(f"\n[EXTRA] (3) khovanov_homology via KnotJob …")
    t0 = time.time()
    khov = khovanov_via_knotjob(pd, name=report.get('name', 'knot'))
    a['khovanov_homology'] = {
        'polynomial':                 khov.get('polynomial'),
        'terms_t_q':                  khov.get('terms_t_q'),
        'delta_diagonals':            khov.get('delta_diagonals'),
        'n_distinct_diagonals':       khov.get('n_distinct_diagonals'),
        'width':                      khov.get('width'),
        'thin':                       khov.get('thin'),
        'torsion':                    khov.get('torsion'),
        'torsion_orders':             khov.get('torsion_orders'),
        'has_torsion':                khov.get('has_torsion'),
        'z2_thin_only':               khov.get('z2_thin_only'),
        'reduced_polynomial':         khov.get('reduced_polynomial'),
        'reduced_terms_t_q':          khov.get('reduced_terms_t_q'),
        'reduced_delta_diagonals':    khov.get('reduced_delta_diagonals'),
        'reduced_n_distinct_diagonals': khov.get('reduced_n_distinct_diagonals'),
        'reduced_width':              khov.get('reduced_width'),
        'reduced_thin':               khov.get('reduced_thin'),
    }
    print(f"  → {time.time()-t0:.1f}s")
    if khov.get('polynomial'):
        print(f"  Khovanov: {khov['polynomial'][:200]}{'...' if len(khov.get('polynomial',''))>200 else ''}")
        print(f"  width={khov.get('width')}, thin={khov.get('thin')}")
    else:
        print(f"  KnotJob output: {khov.get('raw_stdout','')[-300:]}")
    report['raw_subprocess_outputs']['knotjob'] = khov

    # 4. length_spectrum_short via length_spectrum_alt_gen (Discovery #103)
    print(f"\n[EXTRA] (4) length_spectrum_short via length_spectrum_alt_gen …")
    t0 = time.time()
    isosig = h.get('canonical_retriangulation') or report.get('raw_subprocess_outputs', {}).get('snappy', {}).get('isosig')
    # Fix 2026-06-04 (same bug class as unknotting, Disc #255 lesson A,
    # approved by Werner): catch the timeout instead of letting the stage die (the tr15 run broke here).
    try:
        lsp = length_spectrum_via_canonize(pd, n_geos=10, timeout=600, isosig=isosig)
    except subprocess.TimeoutExpired:
        lsp = {'error': 'timeout', 'timeout_seconds': 600,
               'method': 'timeout', 'duration_s': 600,
               'note': 'length_spectrum aborted after 600s; result = timeout'}
    h['length_spectrum_short']     = lsp.get('length_spectrum_short')
    h['length_spectrum_method']    = lsp.get('method')
    h['length_spectrum_duration_s'] = lsp.get('duration_s')
    print(f"  → {time.time()-t0:.1f}s")
    if lsp.get('length_spectrum_short'):
        print(f"  method={lsp.get('method')}, "
              f"{lsp.get('n_geodesics_raw', len(lsp['length_spectrum_short']))} geodesics")
        for g in lsp['length_spectrum_short'][:5]:
            if 'length_real' in g:
                print(f"    L = {g['length_real']:.6f} + {g['length_imag']:+.6f}i  "
                      f"(mult {g['multiplicity']})")
            else:
                print(f"    L = {g.get('length', '?')}  (mult {g.get('mult', '?')})")
    else:
        errs = {k: v for k, v in lsp.items() if k.startswith('err')}
        print(f"  Length spectrum not obtainable; errors: {errs}")
    report['raw_subprocess_outputs']['length_spectrum'] = lsp

    # 5. unknotting_number — brute force via signed Gauss crossing flip
    # (Discovery #100; replaces the braid-based brute force of Discovery #60,
    # because a 12-crossing diagram has 2^12 subsets instead of 2^40 and the
    # min-crossing diagram almost always realises u(K).)
    print(f"\n[EXTRA] (5) unknotting_number signed Gauss crossing flip "
          f"(Pool({args.processes}), max_subset={args.max_unknotting_subset}) …")
    # Pre-check: det != 1 ⇒ NOT trivial (lower bound u >= 1);
    # det == 1 ∧ Jones == 1 ⇒ trivial (u = 0).
    det = a.get('determinant_sage') or a.get('determinant')
    jones = (a.get('jones_polynomial') or '').strip()
    if det is not None and abs(int(det)) == 1 and jones in ('1', 't^0'):
        print(f"  Pre-check: |det|=1 ∧ Jones=1 ⇒ knot trivial, u = 0")
        unk = {'unknotting_number_upper_bound': 0,
               'unknotting_number_lower_bound': 0,
               'shortcut': 'pre_check_det_jones'}
        a['unknotting_number'] = 0
        a['unknotting_number_upper_bound'] = 0
        a['unknotting_number_lower_bound'] = 0
        report['raw_subprocess_outputs']['unknotting'] = unk
    elif sg and (a.get('signature') is not None
                 and abs(int(a['signature'])) / 2 > args.max_unknotting_subset):
        # Fix 2026-06-04 (Disc #255 lesson B, approved by Werner): futility gate.
        # u(K) >= |sig|/2 (own measurement, no lookup) — if the lower bound
        # already lies above the search depth, the negative result is predetermined:
        # skip the search instead of burning a timeout (kf286: sig=18 => u>=9 >> 4).
        sig_lb = abs(int(a['signature'])) // 2 + (abs(int(a['signature'])) % 2 > 0)
        print(f"  Pre-check: |sig|/2 = {abs(int(a['signature']))/2} > max_subset={args.max_unknotting_subset}"
              f" ⇒ search provably fruitless, skip (u >= {sig_lb})")
        a['unknotting_number_upper_bound'] = None
        a['unknotting_number_lower_bound'] = sig_lb
        a['unknotting_method'] = 'skipped_signature_bound'
        a['unknotting_number_searched_to'] = 0
        report['raw_subprocess_outputs']['unknotting'] = {
            'shortcut': 'signature_bound_exceeds_search_depth',
            'unknotting_number_lower_bound': sig_lb,
            'searched_up_to_subset_size': 0}
    elif sg:
        u_lower = 0 if (det is None) else 1  # det != 1 ⇒ u >= 1
        t0 = time.time()
        unk = unknotting_via_signed_gauss_flip(sg,
                                                max_subset=args.max_unknotting_subset,
                                                processes=args.processes)
        u_upper = unk.get('unknotting_number_upper_bound')
        a['unknotting_number_upper_bound'] = u_upper
        a['unknotting_number_lower_bound'] = u_lower
        a['unknotting_witness_subset']     = unk.get('witness_subset')
        a['unknotting_method']             = unk.get('method')
        a['unknotting_number_searched_to'] = unk.get('searched_up_to_subset_size')
        # Discovery #114c: unknotting_number = upper_bound given a stage-1 witness
        # ----------------------------------------------------------------------
        # If stage 1 (Regina simplify→0) delivers a witness subset AND all
        # k < u_upper were searched exhaustively without a hit, then
        # diagram_u(K) = u_upper is provable on the current diagram.
        # The global u(K) could be < u_upper on a DIFFERENT diagram
        # (Regina simplify false-negative rate + non-minimal diagrams),
        # hence document the strict_lower marker.
        # Convention parallel to bridge_index_berge_indicated (Discovery #109).
        if u_upper is None:
            a['unknotting_number'] = None
            a['unknotting_number_note'] = (
                f'u > {args.max_unknotting_subset} on this diagram — '
                f'raise max_subset or check another minimal diagram')
            a['unknotting_number_strict_lower'] = u_lower
            a['unknotting_number_strict_upper'] = None
        else:
            # Check whether k=u_lower..u_upper-1 are all without a hit
            progs = unk.get('progress', []) or []
            ks_searched = sorted({p['k'] for p in progs if 'k' in p})
            ks_required = list(range(u_lower, u_upper))
            ks_with_hits = {p['k'] for p in progs
                            if (p.get('unknot_simplify_zero', 0) > 0
                                or p.get('unknot_jones_only', 0) > 0)}
            exhaustive_below_upper = (
                set(ks_required).issubset(set(ks_searched))
                and not (set(ks_required) & ks_with_hits)
            )

            if u_lower == u_upper:
                # Unique value
                a['unknotting_number']              = u_upper
                a['unknotting_number_strict_lower'] = u_lower
                a['unknotting_number_strict_upper'] = u_upper
                a['unknotting_number_assignment_basis'] = (
                    f'lower == upper == {u_upper}: eindeutig.')
            elif exhaustive_below_upper and unk.get('witness_subset'):
                # Discovery #114c: assign with strict marker
                a['unknotting_number']              = u_upper
                a['unknotting_number_strict_lower'] = u_lower
                a['unknotting_number_strict_upper'] = u_upper
                a['unknotting_number_indicated']    = True
                a['unknotting_number_assignment_basis'] = (
                    f'Discovery #114c: Stage-1-Witness {unk.get("witness_subset")} '
                    f'(Regina simplify→0, topologically reliable — no false '
                    f'positives) at k={u_upper}; k={ks_required} all exhaustive '
                    f'without hit (stage 1 + stage 2). diagram_u(K) = {u_upper} '
                    f'provable on the c-minimal diagram. global u(K) could '
                    f'be < {u_upper} on another diagram, hence strict_lower '
                    f'= {u_lower} stays formal, as not strictly excluded.')
            else:
                # Lower != upper, but not tested exhaustively
                a['unknotting_number']              = None
                a['unknotting_number_strict_lower'] = u_lower
                a['unknotting_number_strict_upper'] = u_upper
                a['unknotting_number_assignment_basis'] = (
                    f'k={ks_required} not completely tested '
                    f'(searched k={ks_searched}, with hits={sorted(ks_with_hits)}); '
                    f'no Discovery #114c auto-assign.')

            a['unknotting_number_note'] = (
                f'{u_lower} <= u(K) <= {u_upper}; Witness-Subset '
                f'{unk.get("witness_subset")} gives the unknot via crossing flip; '
                f'k={list(range(u_lower, u_upper))} without hit.')
        print(f"  → {time.time()-t0:.1f}s")
        for prog in unk.get('progress', []):
            if 'k' in prog:
                print(f"    k={prog['k']}: {prog['n_subsets']} subsets, "
                      f"simplify→0={prog.get('unknot_simplify_zero', 0)}, "
                      f"jones=1={prog.get('unknot_jones_only', 0)}")
        if u_upper is not None:
            print(f"  u(K) ∈ [{u_lower}, {u_upper}] "
                  f"(Witness {unk.get('witness_subset')})")
        else:
            print(f"  u(K) > {args.max_unknotting_subset} on this diagram")
        report['raw_subprocess_outputs']['unknotting'] = unk
    else:
        print(f"\n[EXTRA] (5) unknotting_number — no signed_gauss, skipped")

    # 6. Kauffman polynomial — CHA Markov trace (n<=4) primary, knotpy as
    #    fallback (Discovery #114a). For orbit (n=7) the CHA route ran into a
    #    KeyError, knotpy delivered in 2.4s.
    print(f"\n[EXTRA] (6) kauffman_polynomial via CHA Markov trace …")
    kauff_cha = None
    kauff_knotpy = None
    kauff_value = None
    kauff_method = None
    if braid_word and n_strands:
        t0 = time.time()
        # Fix 2026-06-04 (same bug class, Disc #255 lesson A, approved by Werner):
        # catch the timeout -> the knotpy fallback applies, the stage does not die.
        try:
            kauff_cha = kauffman_via_cha(braid_word, n_strands)
        except subprocess.TimeoutExpired:
            kauff_cha = {'error': 'timeout (600s)'}
        print(f"  CHA → {time.time()-t0:.1f}s")
        if kauff_cha.get('kauffman_polynomial'):
            kauff_value  = kauff_cha['kauffman_polynomial']
            kauff_method = 'CHA-Markov-Trace (Sage CubicHeckeAlgebra, n<=4, Discovery #55)'
            print(f"  CHA success: {kauff_value[:200]}")
        else:
            print(f"  CHA NOT available: {kauff_cha.get('error')}")

    if kauff_value is None and pd:
        print(f"[EXTRA] (6b) kauffman_polynomial — Fallback knotpy (Discovery #114a)")
        t0 = time.time()
        kauff_knotpy = kauffman_via_knotpy(pd)
        print(f"  knotpy → {time.time()-t0:.1f}s")
        if kauff_knotpy.get('kauffman_polynomial'):
            kauff_value  = kauff_knotpy['kauffman_polynomial']
            kauff_method = kauff_knotpy.get('kauffman_polynomial_method',
                                             'knotpy.kauffman state-sum')
            print(f"  knotpy success: {kauff_knotpy['kauffman_polynomial_n_chars']} chars, "
                  f"amphi_symm={kauff_knotpy.get('kauffman_polynomial_amphi_symmetry_check')}")
        else:
            print(f"  knotpy NOT available: {kauff_knotpy.get('error')}")

    if kauff_value is not None:
        a['kauffman_polynomial']                  = kauff_value
        a['kauffman_polynomial_pipeline_status']  = 'computed'
        a['kauffman_polynomial_method']           = kauff_method
        if kauff_knotpy and kauff_method and 'knotpy' in kauff_method:
            a['kauffman_polynomial_n_chars']      = kauff_knotpy.get('kauffman_polynomial_n_chars')
            a['kauffman_polynomial_n_crossings_input'] = kauff_knotpy.get('kauffman_polynomial_n_crossings_input')
            a['kauffman_polynomial_amphi_symmetry_check'] = kauff_knotpy.get('kauffman_polynomial_amphi_symmetry_check')
            a['kauffman_polynomial_convention_note']      = kauff_knotpy.get('kauffman_polynomial_convention_note')
            a['kauffman_polynomial_cha_error']    = (kauff_cha or {}).get('error')
    else:
        a['kauffman_polynomial']                 = None
        a['kauffman_polynomial_pipeline_status'] = 'not_computable'
        a['kauffman_polynomial_note']            = (
            f"CHA: {(kauff_cha or {}).get('error', 'no input')}; "
            f"knotpy: {(kauff_knotpy or {}).get('error', 'no input')}"
        )
        # Disc #218 fix A (2026-05-19 ~16:35): persist the knotpy error explicitly
        # so that enrich_optional_invariants wave 4 does not repeat
        # the same kauffman_via_knotpy call sequentially.
        # Before this fix: pflicht_extra timeout + enrich timeout =
        # 2 * 600s = 20+ min hang on c=24+ chiral knots
        # (W32 run 6 kf376, 2026-05-18).
        if kauff_knotpy and kauff_knotpy.get('error'):
            a['kauffman_polynomial_knotpy_error'] = kauff_knotpy.get('error')
        if kauff_cha and kauff_cha.get('error'):
            a['kauffman_polynomial_cha_error']    = kauff_cha.get('error')

    # Summary
    print(f"\n{'='*60}")
    print(f"ADDED MANDATORY INVARIANTS")
    print(f"{'='*60}")
    print(f"  braid_index_interval       : "
          f"[{d.get('braid_index_lower_bound')}, {d.get('braid_index_upper_bound')}] "
          f"(MFW-Spread {d.get('braid_index_mfw_v_spread')})")
    print(f"  bridge_number_interval     : "
          f"[{d.get('bridge_number_lower_bound')}, {d.get('bridge_number_upper_bound')}] "
          f"(Schubert lower; 2-bridge excluded: {d.get('bridge_number_is_2_excluded')})")
    print(f"  s_invariant                : {a.get('s_invariant')}")
    print(f"  khovanov_homology width    : {a.get('khovanov_homology', {}).get('width')}")
    print(f"  khovanov_homology thin     : {a.get('khovanov_homology', {}).get('thin')}")
    print(f"  length_spectrum_short      : "
          f"{len(h.get('length_spectrum_short') or [])} geodesics")
    print(f"  unknotting_number          : {a.get('unknotting_number')}")
    print(f"  kauffman_polynomial        : "
          f"{'available' if a.get('kauffman_polynomial') else a.get('kauffman_polynomial_note')}")

    with open(output, 'w') as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\n[EXTRA] Writing {output}")


if __name__ == '__main__':
    main()
