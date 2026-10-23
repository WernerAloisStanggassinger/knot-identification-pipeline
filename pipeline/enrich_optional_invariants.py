#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
enrich_optional_invariants.py
=============================
Enriches a pflicht_invarianten_*.json with optional invariants from
optionale_invarianten.txt. Five waves:

  Wave 1 (--phase trivial): derivable from existing mandatory fields,
    no subprocess. E.g. amphicheiral, invertible, prime, components,
    volume_real_part, injectivity_radius, ...

  Wave 2 (--phase sage): needs a Sage subprocess for Vassiliev v2/v3,
    multivariable_alexander, seifert_matrix, knot_group_generators.

  Wave 3 (--phase snappy): needs a snappy subprocess for cusp_area,
    cusp_volume, meridian/longitude_translation.

  Wave 4 (--phase kauffman): kauffman_polynomial via knotpy
    (Discovery #106 — replaces the CHA Markov-trace n<=4 limit from pflicht_extra).

  Wave 5 (--phase bridge): meridional_rank_upper_bound from
    Wirtinger knot_group_generators_simplified, plus a Berge-indicated
    tightening of bridge_number (generalises Discovery #103, Saturn).

Default: all waves; can be reduced with --skip-{sage,snappy,kauffman,bridge}.

Writes into the NEW top-level block 'optional_invariants' in the JSON
(NOT into pflichtliste, so that the mandatory/optional separation stays clear).

Usage:
  python enrich_optional_invariants.py <pflicht_*.json> [--in-place]
                                                       [--skip-sage]
                                                       [--skip-snappy]
                                                       [--skip-kauffman]
                                                       [--skip-bridge]
"""
import argparse
import site_config
import json
import re
import subprocess
import sys
import textwrap
import time
from pathlib import Path

# Disc #222: subprocess.run + conda run -> orphaned inner python on timeout.
# pipeline_subprocess.run_capture_with_group_kill kills the whole process group.
from pipeline_subprocess import (
    run_capture_with_group_kill as _run_capture_with_group_kill,
    install_signal_killer,
)

SEAL_PRELUDE = (
    "import sys\n"
    "sys.path.insert(0, '@FREEZE_DIR@')\n"
    "from pipeline_seal import seal_all\n"
    "seal_all(verbose=False)\n".replace('@FREEZE_DIR@', site_config.FREEZE_DIR)
)


# ──────────────────────────────────────────────
# Wave 1: trivial extracts (no subprocess)
# ──────────────────────────────────────────────

def derive_trivial_invariants(report):
    a = report['pflichtliste']['ALGEBRA']
    d = report['pflichtliste']['DIAGRAMM']
    h = report['pflichtliste']['HYPERBOLISCH']
    raw = report.get('raw_subprocess_outputs', {})
    sn  = raw.get('snappy', {})
    out = {}

    # === DIAGRAMM/NOTATION ===
    out['braid_length']           = d.get('braid_length')
    out['braid_notation']         = str(d.get('braid_word'))
    out['gauss_notation']         = d.get('gauss_code_reduced')
    out['pd_notation_vector']     = d.get('pd_code')
    out['signed_gauss']           = d.get('signed_gauss_reduced') or report.get('signed_gauss')
    out['oriented_gauss']         = d.get('oriented_gauss')
    out['dt_code']                = d.get('dt_code')
    out['knot_sig']               = d.get('knot_sig')
    out['writhe']                 = d.get('writhe')
    out['frequency']              = d.get('frequency')
    out['frequency_ratio']        = d.get('frequency_ratio')
    out['adequate']               = d.get('alternating')  # alternating ⇒ adequate (Lickorish-Thistlethwaite)
    _hv = h.get('hyperbolic_volume')
    out['hyperbolic']             = (isinstance(_hv, (int, float)) and _hv > 1e-6)
    out['bridge_index']           = d.get('bridge_number_interval')   # interval, not a unique value
    out['braid_index_interval']   = d.get('braid_index_interval')

    # === TOPOLOGY/CLASSIFICATION ===
    out['amphicheiral']           = h.get('is_amphicheiral')
    out['negative_amphicheiral']  = None  # handedness-specific, Saturn is both
    out['positive_amphicheiral']  = None
    out['reversible']             = h.get('is_invertible_knot')
    out['invertible']             = h.get('is_invertible_knot')
    out['symmetry_type']          = ('amphicheiral_invertible'
                                     if h.get('is_amphicheiral') and h.get('is_invertible_knot')
                                     else 'asymmetric_or_other')
    # quasialternating: alternating ⇒ QA (Manolescu-Ozsvath thin homology)
    # Saturn alternating + Khovanov-thin ⇒ QA.
    out['quasialternating']       = bool(d.get('alternating') and a.get('khovanov_homology', {}).get('thin'))
    # torus_knot T(p,q): genus = (p-1)(q-1)/2; all torus knots are non-amphichiral (except the unknot).
    # Saturn is amphichiral ⇒ NOT a torus knot.
    out['torus_knot']             = False if h.get('is_amphicheiral') else None
    out['torus_knot_parameters']  = None
    # satellite: hyperbolic knots are NOT satellites (by the Thurston classification).
    out['satellite']              = False if out['hyperbolic'] else None
    out['satellite_parameters']   = None
    out['cable']                  = False if out['hyperbolic'] else None
    # Two-bridge: only 'False' if explicitly excluded, otherwise 'None'
    # (the pipeline has no positive 2-bridge proof path; conservatively "unknown").
    out['two_bridge']             = False if d.get('bridge_number_is_2_excluded') else None
    # Prime: hyperbolic ⇒ prime (Thurston trichotomy: connected sums have essential
    # tori in the complement, so they are not hyperbolic. Contrapositive: hyperbolic ⇒ prime).
    # Non-hyperbolic knots (torus/satellite) need an extra argument; conservatively None.
    alex = (a.get('alexander_polynomial') or '').strip()
    out['prime']                  = True  if out.get('hyperbolic') else None
    out['composite']              = False if out.get('hyperbolic') else None
    out['mutation']               = None  # mutation partner unknown without a table lookup
    out['mutant_knot']            = None

    # Quasi-positive: positive braid (all generators positive) ⇒ quasi-positive.
    # The Saturn Sage braid has negative generators, so it is not obviously quasi-positive.
    bw = d.get('braid_word') or []
    if isinstance(bw, str):
        try: bw = json.loads(bw)
        except: bw = []
    out['quasi_positive_braid']   = bool(bw) and all(g > 0 for g in bw)
    out['positive']               = out['quasi_positive_braid']  # positive = positive braid

    # === ALGEBRA/POLYNOMIALS ===
    kh = a.get('khovanov_homology', {})
    out['khovanov_polynomial']                       = kh.get('polynomial')
    out['khovanov_unreduced_integral_polynomial']    = kh.get('polynomial')
    out['khovanov_reduced_polynomial']               = kh.get('reduced_polynomial')
    out['khovanov_torsion']                          = kh.get('torsion')
    out['khovanov_torsion_orders']                   = kh.get('torsion_orders')
    out['khovanov_width']                            = kh.get('width')
    out['khovanov_thin']                             = kh.get('thin')
    out['rasmussen_s']                               = a.get('s_invariant')
    out['ozsvath_szabo_tau']                         = a.get('tau_invariant')
    out['nu_invariant']                              = a.get('knot_floer_homology', {}).get('nu')
    out['knot_floer_tau']                            = a.get('knot_floer_homology', {}).get('tau')
    out['knot_floer_nu']                             = a.get('knot_floer_homology', {}).get('nu')
    out['knot_floer_epsilon']                        = a.get('knot_floer_homology', {}).get('epsilon')
    out['hat_hfk_rank']                              = a.get('knot_floer_homology', {}).get('total_rank')
    # HFK thickness: 1 for an L-space knot, 1 for thin, > 1 for non-thin
    out['knot_floer_thickness']                      = (1 if a.get('knot_floer_homology', {}).get('total_rank')
                                                       and kh.get('thin') else None)
    out['l_space_knot']                              = a.get('knot_floer_homology', {}).get('L_space_knot')
    out['seifert_genus']                             = a.get('knot_floer_homology', {}).get('seifert_genus')
    out['_3_genus']                                  = a.get('genus')
    out['_4_genus']                                  = a.get('smooth_4_genus')
    # g_4^top <= g_4^smooth (every smooth surface is also topologically locally flat) and
    # g_4^top >= |signature|/2. Formerly equality stood here -- a claim, not a theorem
    # (2026-07-23, Disc #556). A value only if the bounds coincide, otherwise None + bounds.
    out['topological_4_genus']                       = a.get('topological_4_genus')
    out['topological_four_genus']                    = a.get('topological_4_genus')
    out['topological_4_genus_lower_bound']           = a.get('topological_4_genus_lower_bound')
    out['topological_4_genus_upper_bound']           = a.get('topological_4_genus_upper_bound')
    out['smooth_4_genus_lower_bound']                = a.get('smooth_4_genus_lower_bound')
    out['smooth_4_genus_upper_bound']                = a.get('smooth_4_genus_upper_bound')
    out['smooth_4_genus_interval']                   = a.get('smooth_4_genus_interval')
    # smooth slice ⇒ topologically slice. Since 2026-07-23 a['slice'] is three-valued
    # (True only for the unknot / False proven not slice / None undetermined) — None must NOT
    # become False, otherwise "we don't know" turns into a claim again.
    out['topologically_slice']                       = a.get('slice')
    out['unknotting_number_lower_bound']             = a.get('unknotting_number_lower_bound')
    out['unknotting_number_upper_bound']             = a.get('unknotting_number_upper_bound')
    out['unknotting_witness_subset']                 = a.get('unknotting_witness_subset')
    # Concordance order: slice ⇒ order 1 in the concordance group (= trivial).
    # `is True` is deliberate: with slice=None the order is unknown, not 1 (2026-07-23).
    out['concordance_order']                         = (1 if a.get('slice') is True else None)
    # Concordance genus = g_4. None means undetermined, not 0 (2026-07-23).
    out['concordance_genus']                         = a.get('smooth_4_genus')

    # === VASSILIEV / FINITE-TYPE ===
    # v_2 = c_2 (Casson invariant, the Conway polynomial coefficient of z²)
    # Saturn Conway: 1 - 6z² + 9z⁴ ⇒ coefficient of z² is -6 ⇒ v_2 = -6.
    conway_str = a.get('conway_polynomial', '')
    m_z2 = re.search(r'([+-]?\s*\d+)\s*z\^?2(?!\d)', conway_str.replace(' ', ''))
    if m_z2:
        try:
            v2_val = int(m_z2.group(1).replace(' ', '') or '1')
            out['vassiliev_invariant_2'] = v2_val
            out['v2'] = v2_val
        except Exception:
            pass
    # v_3: from the Jones polynomial derivative at t=1, J''(1)/6 or convention.
    # Saturn is amphichiral ⇒ v_3 = 0 (v_3 is of odd order
    # and vanishes on amphichiral knots).
    if h.get('is_amphicheiral'):
        out['vassiliev_invariant_3'] = 0
        out['v3'] = 0

    # === NUMERICAL INVARIANTS ===
    # Dealternating number: 0 for alternating knots.
    out['dealternating_number']   = 0 if d.get('alternating') else None
    # Turaev genus: 0 for alternating knots (Lowrance-Spyropoulos).
    out['turaev_genus']           = 0 if d.get('alternating') else None

    # === HYPERBOLIC/GEOMETRY ===
    out['volume_real_part']       = h.get('hyperbolic_volume')
    out['volume_imaginary_part']  = 0.0  # knot complements have real volume
    out['num_tetrahedra']         = sn.get('num_tetrahedra')
    out['solution_type']          = sn.get('solution_type')
    out['orientability']          = True  # knot complements are orientable
    # injectivity_radius = half of Re(shortest geodesic)
    ls_short = h.get('length_spectrum_short')
    if ls_short:
        first = ls_short[0]
        if 'length_real' in first:
            out['injectivity_radius']     = abs(first['length_real']) / 2.0
            out['shortest_geodesic_real'] = first['length_real']
            out['shortest_geodesic_imag'] = first.get('length_imag')

    # === LINK-SPECIFIC (trivial for a knot) ===
    out['components']             = 1
    out['nullity']                = 0
    out['splitting_number']       = None  # n/a for a knot

    return out


# ──────────────────────────────────────────────
# Wave 2: Sage subprocess for Vassiliev (verification),
#          multivariable_alexander, seifert_matrix
# ──────────────────────────────────────────────

def sage_extras(braid_word, n_strands):
    code = SEAL_PRELUDE + textwrap.dedent(f"""
        import json
        from sage.all import BraidGroup, Knot, var
        try:
            B = BraidGroup({n_strands})
            b = B({braid_word!r})
            K = Knot(b)
            out = {{}}
            # Multivariable Alexander - for a knot identical to the standard Alexander
            try:
                out['multivariable_alexander'] = str(K.alexander_polynomial())
            except Exception as e:
                out['multivariable_alexander_err'] = str(e)
            # Seifert-Matrix
            try:
                S = K.seifert_matrix()
                out['seifert_matrix']         = str(S)
                out['seifert_matrix_rank']    = int(S.rank())
                out['seifert_matrix_nrows']   = int(S.nrows())
            except Exception as e:
                out['seifert_matrix_err'] = str(e)
            # knot group generators: number of minimal Wirtinger generators
            # (lower bound for bridge_number)
            try:
                G = K.fundamental_group()
                out['knot_group_generators_initial'] = int(G.ngens())
                G_simp = G.simplification_isomorphism().codomain()
                out['knot_group_generators_simplified'] = int(G_simp.ngens())
            except Exception as e:
                out['knot_group_err'] = str(e)
            print(json.dumps(out, default=str))
        except Exception as e_top:
            print(json.dumps({{'error_top': str(e_top)}}, default=str))
    """)
    cmd = ['conda', 'run', '-n', 'sage_stable', '--no-capture-output',
           'sage', '-c', code]
    try:
        res = _run_capture_with_group_kill(cmd, timeout=300)
    except subprocess.TimeoutExpired:
        return {'error': 'sage_timeout_300s'}
    if res.returncode != 0:
        return {'error': f'rc={res.returncode}', 'stderr': res.stderr[-1500:]}
    for line in reversed(res.stdout.strip().splitlines()):
        line = line.strip()
        if line.startswith('{') and line.endswith('}'):
            try: return json.loads(line)
            except: continue
    return {'error': 'no json', 'stdout': res.stdout[-1000:]}


# ──────────────────────────────────────────────
# Wave 3: snappy extras (cusp area, volume etc.)
# ──────────────────────────────────────────────

def snappy_extras(isosig):
    code = SEAL_PRELUDE + textwrap.dedent(f"""
        import json, snappy
        try:
            M = snappy.Manifold({isosig!r})
            out = {{}}
            out['volume'] = float(M.volume())
            try:
                ca = M.cusp_areas()
                out['cusp_areas']      = [float(x) for x in ca]
                out['cusp_area_total'] = float(sum(ca))
            except Exception as e:
                out['cusp_areas_err'] = str(e)
            try:
                # cusp_neighborhood() starts at displacement 0 -> volume(0) is the same constant
                # 0.324759526419 for EVERY one-cusp knot and does not measure the knot
                # (Disc #745, 648 records with one single value). The maximum is only
                # there after set_displacement(stopping_displacement()).
                cn = M.cusp_neighborhood()
                _v0 = float(cn.volume(0))
                cn.set_displacement(cn.stopping_displacement())
                _vmax = float(cn.volume(0))
                _vol = float(out['volume'])
                # canary: the maximum cannot be smaller than the displacement-0 value
                # and not larger than the total volume.
                if not (_vmax == _vmax) or _vmax <= 0:
                    out['cusp_volume_err'] = 'impossible maximum %r' % _vmax
                elif _vmax < _v0 - 1e-9:
                    out['cusp_volume_err'] = 'maximum %r smaller than displacement-0 %r' % (_vmax, _v0)
                elif _vmax > _vol + 1e-6:
                    out['cusp_volume_err'] = 'cusp volume %r larger than total volume %r' % (_vmax, _vol)
                else:
                    out['cusp_volume_max'] = _vmax
                    out['cusp_volume_at_displacement_0'] = _v0
                    out['cusp_volume_max_via'] = 'set_displacement(stopping_displacement())'
            except Exception as e:
                out['cusp_volume_err'] = str(e)
            try:
                ti = M.cusp_translations()
                out['cusp_translations'] = [
                    {{'meridian': str(t[0]), 'longitude': str(t[1])}}
                    for t in ti
                ]
            except Exception as e:
                out['cusp_translations_err'] = str(e)
            try:
                out['homology_str'] = str(M.homology())
            except Exception as e:
                out['homology_err'] = str(e)
            try:
                out['num_cusps']   = int(M.num_cusps())
            except Exception as e:
                out['num_cusps_err'] = str(e)
            print(json.dumps(out, default=str))
        except Exception as e_top:
            print(json.dumps({{'error_top': str(e_top)}}, default=str))
    """)
    cmd = ['conda', 'run', '-n', 'snappy_env', '--no-capture-output',
           'python', '-c', code]
    try:
        res = _run_capture_with_group_kill(cmd, timeout=120)
    except subprocess.TimeoutExpired:
        return {'error': 'snappy_timeout_120s'}
    if res.returncode != 0:
        return {'error': f'rc={res.returncode}', 'stderr': res.stderr[-1500:]}
    for line in reversed(res.stdout.strip().splitlines()):
        line = line.strip()
        if line.startswith('{') and line.endswith('}'):
            try: return json.loads(line)
            except: continue
    return {'error': 'no json', 'stdout': res.stdout[-1000:]}


# ──────────────────────────────────────────────
# Wave 4: Kauffman polynomial via knotpy (Discovery #106)
# ──────────────────────────────────────────────

def kauffman_via_knotpy(pd_code_str):
    """Kauffman two-variable polynomial via the knotpy state sum.
    PD code format as from Regina: 'PD[X[a,b,c,d], X[...], ...]'.
    Convention note: knotpy's state sum can deliver a mirror polynomial
    (a -> a^-1) relative to KnotInfo (see the Discovery #106
    trefoil example). For amphichiral knots this does not arise.
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
                # amphichirality symmetry check: F(a) == F(1/a) ?
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
                    'kauffman_polynomial_method': 'knotpy.kauffman state-sum (Discovery #106)',
                    'kauffman_polynomial_amphi_symmetry_check': amphi_check,
                    'kauffman_polynomial_convention_note': 'knotpy state-sum; for chiral knots possibly the mirror relative to KnotInfo (Discovery #106)',
                }}))
        except Exception as e:
            import traceback
            print(json.dumps({{'error': str(e), 'traceback': traceback.format_exc()[-1500:]}}))
    """)
    cmd = ['conda', 'run', '-n', 'base', '--no-capture-output',
           'python', '-c', code]
    try:
        res = _run_capture_with_group_kill(cmd, timeout=600)
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
# Wave 5: bridge-number tightening via the Wirtinger meridional_rank
# (generalises Discovery #103, Saturn)
# ──────────────────────────────────────────────

def derive_meridional_rank_and_bridge(report, optional):
    """From Wirtinger-simplified generators (knot_group_generators_simplified
    from wave 2) derive the upper bound for meridional_rank and
    tighten bridge_number via the Berge conjecture.

    Theorems:
      Boileau-Zieschang: bridge_number(K) >= meridional_rank(K)
      Berge (proven c<=12, conjectured c>=13):
        bridge_number(K) = meridional_rank(K)

    Wirtinger heuristic (as in Discovery #103, Saturn):
      Sage's simplification_isomorphism().codomain().ngens() is
      an upper bound for rank(pi_1(S^3 \\ K)). Empirically, for
      small alternating knots = meridional_rank (the Saturn case).
      We treat the result as meridional_rank_upper_bound.
    """
    out = {}
    mr_upper = optional.get('knot_group_generators_simplified')
    if not isinstance(mr_upper, int):
        return {'meridional_rank_note': 'no knot_group_generators_simplified available (wave 2 Sage skipped or failed)'}

    d = report['pflichtliste'].get('DIAGRAMM', {})
    a = report['pflichtliste'].get('ALGEBRA', {})
    c_min = d.get('crossing_number_min', 999) or 999
    det_K = a.get('determinant', 1) or 1
    # 🔴 CORRECTION 2026-08-10 (finding 15, sweep of the same cause):
    # non-triviality is decided by the SEIFERT GENUS (g > 0 ⟺ not the
    # unknot), not by det > 1. det = |Δ(-1)| is unfit for this: Conway,
    # Kinoshita-Terasaka, every T(p,q) with p,q odd (among them our anchor
    # `irrational` = T(3,7)) have det = 1 and are highly non-trivial.
    g_K = a.get('genus')
    if g_K is None:
        g_K = (a.get('knot_floer_homology') or {}).get('seifert_genus')
    nontrivial = (g_K is not None and int(g_K) > 0) or (det_K and abs(int(det_K)) > 1)

    out['meridional_rank_upper_bound']         = mr_upper
    # Trivial lower bound: nontrivial knot ⇒ mr >= 2 (genus > 0 OR det > 1)
    out['meridional_rank_lower_bound']         = 2 if nontrivial else 1
    out['meridional_rank_method']              = 'Sage simplification_isomorphism().codomain().ngens() (Wirtinger heuristic as in Discovery #103)'

    # Berge status
    berge_proven = (isinstance(c_min, int) and c_min <= 12)
    out['berge_conjecture_proven_for_this_K']  = berge_proven
    out['berge_conjecture_status']             = (
        f'PROVEN (c={c_min}<=12, Boileau-Boyer-Zieschang)'
        if berge_proven
        else f'CONJECTURED only (c={c_min}>12)'
    )

    # ====== FRESH CALCULATION of the bridge bounds from first principles ======
    # (do NOT inherit the possibly wrong pflicht_extra bound — see Discovery #109:
    #  pflicht_extra sets 'Schubert b>=g+1', which is NOT a general theorem; counterexample
    #  the (5,2) torus knot with g=2, b=2.)
    fresh_lo = 2 if nontrivial else 1     # (2026-08-10: genus instead of det, see above)
    if optional.get('two_bridge') is False:
        fresh_lo = max(fresh_lo, 3)
    # Boileau-Zieschang: bridge >= meridional_rank >= mr_lower
    fresh_lo = max(fresh_lo, out['meridional_rank_lower_bound'])
    fresh_hi = int(c_min)

    # Detect an inconsistency with the old bound
    cur = optional.get('bridge_index')
    inherited_lo, inherited_hi = None, None
    if isinstance(cur, list) and len(cur) == 2 and all(isinstance(x, (int, float)) for x in cur):
        inherited_lo, inherited_hi = int(cur[0]), int(cur[1])

    inconsistency_note = None
    if inherited_lo is not None and inherited_lo > mr_upper:
        inconsistency_note = (
            f'Previous lower bound {inherited_lo} (probably from pflicht_extra '
            f'"Schubert b>=g+1") > meridional_rank_upper {mr_upper}. '
            f'This inequality is NOT a general theorem (counterexample: '
            f'(5,2) torus knot g=2, b=2). Overwriting with fresh_lo={fresh_lo} '
            f'from first principles. See Discovery #109.'
        )
        cur_lo = fresh_lo
    else:
        cur_lo = max(fresh_lo, inherited_lo if inherited_lo is not None else 0)
    cur_hi = inherited_hi if inherited_hi is not None else fresh_hi

    # Upper bound via Berge: bridge = mr <= mr_upper
    if berge_proven:
        new_hi = min(cur_hi, mr_upper)
        out['bridge_index']                    = [cur_lo, new_hi]
        out['bridge_index_upper_source']       = (
            f'Berge theorem (c={c_min}<=12, PROVEN): b = mr <= {mr_upper}'
        )
    else:
        # Caution: not proven, hence do NOT narrow the interval,
        # but document the Berge-indicated value.
        out['bridge_index']                    = [cur_lo, cur_hi]
        out['bridge_index_berge_indicated']    = [cur_lo, mr_upper]
        out['bridge_index_upper_source']       = (
            f'Trivial b<=c={c_min} (Berge conjecture for c={c_min}>12 not yet proven; '
            f'Berge-indicated upper bound: {mr_upper})'
        )
    out['bridge_index_lower_source'] = (
        f'fresh_lo={fresh_lo} = max(det>1=>2, 2-bridge-ausschluss=>3, '
        f'Boileau-Zieschang b>=mr>={out["meridional_rank_lower_bound"]})'
    )
    if inconsistency_note:
        out['bridge_index_inconsistency_with_pflicht_extra'] = inconsistency_note

    return out


# ──────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────

def main():
    # Hole B (OOM 2026-07-01): pass an external kill on to the conda/sage/knotpy grandchildren.
    install_signal_killer()
    ap = argparse.ArgumentParser()
    ap.add_argument('json_path')
    ap.add_argument('--in-place', action='store_true')
    ap.add_argument('--skip-sage', action='store_true')
    ap.add_argument('--skip-snappy', action='store_true')
    ap.add_argument('--skip-kauffman', action='store_true',
                    help='Skip wave 4 (Kauffman via knotpy, Discovery #106)')
    ap.add_argument('--skip-bridge', action='store_true',
                    help='Skip wave 5 (bridge_number via Wirtinger mr, Discovery #103)')
    args = ap.parse_args()

    p = Path(args.json_path)
    with p.open() as f:
        report = json.load(f)

    optional = report.get('optional_invariants', {})

    # === Wave 1: trivial ===
    print('=== Wave 1: trivial extracts ===')
    t0 = time.time()
    triv = derive_trivial_invariants(report)
    optional.update(triv)
    print(f'  {len(triv)} fields in {time.time()-t0:.2f}s')
    populated = sum(1 for v in triv.values() if v is not None and v != [])
    print(f'  populated: {populated} / {len(triv)}')

    # === Wave 2: Sage ===
    if not args.skip_sage:
        print('\n=== Wave 2: Sage extras ===')
        t0 = time.time()
        bw = report['pflichtliste']['DIAGRAMM'].get('braid_word')
        if isinstance(bw, str):
            try: bw = json.loads(bw)
            except: bw = None
        n_str = report['pflichtliste']['DIAGRAMM'].get('braid_n_strands')
        if bw and n_str:
            sage_out = sage_extras(bw, n_str)
            print(f'  {time.time()-t0:.1f}s, keys: {sorted(sage_out.keys())}')
            optional.update({k: v for k, v in sage_out.items() if not k.startswith('error')})
            if 'error' in sage_out or 'error_top' in sage_out:
                print(f'  Sage error: {sage_out.get("error") or sage_out.get("error_top")}')
        else:
            print('  skipped (no braid_word).')

    # === Wave 3: snappy ===
    if not args.skip_snappy:
        print('\n=== Wave 3: snappy extras ===')
        t0 = time.time()
        iso = report['pflichtliste']['HYPERBOLISCH'].get('canonical_retriangulation')
        if iso:
            sn_out = snappy_extras(iso)
            print(f'  {time.time()-t0:.1f}s, keys: {sorted(sn_out.keys())}')
            optional.update({k: v for k, v in sn_out.items() if not k.startswith('error')})
            if 'error' in sn_out or 'error_top' in sn_out:
                print(f'  snappy error: {sn_out.get("error") or sn_out.get("error_top")}')
        else:
            print('  skipped (no isosig).')

    # === Wave 4: Kauffman via knotpy (Discovery #106) ===
    # Discovery #114b: additionally fill the mandatory slot as well (idempotent),
    # if it is still empty there (or unreliable from a CHA failure).
    if not args.skip_kauffman:
        print('\n=== Wave 4: Kauffman via knotpy (Discovery #106) ===')
        # Disc #218 fix B (2026-05-19 ~16:35): idempotency pre-check before
        # the expensive kauffman_via_knotpy. Wave 4 was meant to be "idempotent
        # for older JSONs" (Disc #114b), but the implementation
        # was not — the compute was unconditional, only the write was
        # protected. For c=24+ chiral knots without a CHA hit,
        # first pflicht_invarianten_extra called knotpy.kauffman → 600s timeout,
        # then this stage AGAIN → another 600s = 20+ min hang
        # (W32 run 6 kf376, 2026-05-18). Disc #218 documents the fix.
        alg = report['pflichtliste'].get('ALGEBRA', {})
        kauff_already   = alg.get('kauffman_polynomial')
        knotpy_err_prev = str(alg.get('kauffman_polynomial_knotpy_error') or '')
        if kauff_already:
            n_chars = alg.get('kauffman_polynomial_n_chars') or len(str(kauff_already))
            print(f'  skipped [Disc #218 fix]: pflicht_extra has already '
                  f'delivered a Kauffman polynomial ({n_chars} chars, '
                  f'method={alg.get("kauffman_polynomial_method")}).')
        elif 'knotpy_timeout' in knotpy_err_prev:
            print(f'  skipped [Disc #218 fix]: pflicht_extra has already '
                  f'documented a knotpy error ({knotpy_err_prev}); '
                  f'calling it again would trigger the same hang sequentially.')
            optional['kauffman_polynomial_skip_reason'] = (
                f'pflicht_extra knotpy: {knotpy_err_prev}')
        else:
            t0 = time.time()
            pd_code = report['pflichtliste']['DIAGRAMM'].get('pd_code')
            if pd_code:
                kf_out = kauffman_via_knotpy(pd_code)
                print(f'  {time.time()-t0:.1f}s, keys: {sorted(kf_out.keys())}')
                if 'error' in kf_out:
                    print(f'  knotpy error: {kf_out.get("error")}')
                    # Store the error info as a note anyway
                    optional['kauffman_polynomial_error'] = kf_out.get('error')
                else:
                    optional.update({k: v for k, v in kf_out.items()
                                     if not k.startswith('error')})
                    kp_str = kf_out.get('kauffman_polynomial', '')
                    print(f'  kauffman_polynomial: {len(kp_str)} chars, '
                          f'amphi_symm_check={kf_out.get("kauffman_polynomial_amphi_symmetry_check")}')

                    # Discovery #114b: also fill the mandatory slot if it is still
                    # empty. pflicht_invarianten_extra.py has had its own knotpy
                    # fallback since Discovery #114a, but wave 4 is an
                    # additional safeguard (idempotent for older JSONs).
                    if not alg.get('kauffman_polynomial'):
                        alg['kauffman_polynomial']                  = kf_out.get('kauffman_polynomial')
                        alg['kauffman_polynomial_pipeline_status']  = 'computed'
                        alg['kauffman_polynomial_method']           = kf_out.get('kauffman_polynomial_method')
                        alg['kauffman_polynomial_n_chars']          = kf_out.get('kauffman_polynomial_n_chars')
                        alg['kauffman_polynomial_n_crossings_input'] = kf_out.get('kauffman_polynomial_n_crossings_input')
                        alg['kauffman_polynomial_amphi_symmetry_check'] = kf_out.get('kauffman_polynomial_amphi_symmetry_check')
                        alg['kauffman_polynomial_convention_note']      = kf_out.get('kauffman_polynomial_convention_note')
                        alg['kauffman_polynomial_source_section']       = (
                            'enrich_optional_invariants wave 4 (Discovery #114b): '
                            'mandatory slot was empty, filled with the knotpy value.'
                        )
                        print(f'  Discovery #114b: mandatory slot pflichtliste.ALGEBRA.kauffman_polynomial filled as well.')
            else:
                print('  skipped (no pd_code).')

    # === Wave 5: bridge number via Wirtinger mr (Discovery #103) ===
    if not args.skip_bridge:
        print('\n=== Wave 5: bridge number via Wirtinger mr (Discovery #103) ===')
        t0 = time.time()
        br_out = derive_meridional_rank_and_bridge(report, optional)
        print(f'  {time.time()-t0:.2f}s, keys: {sorted(br_out.keys())}')
        optional.update(br_out)
        if 'meridional_rank_upper_bound' in br_out:
            print(f'  meridional_rank ∈ [{br_out["meridional_rank_lower_bound"]}, '
                  f'{br_out["meridional_rank_upper_bound"]}]')
            print(f'  bridge_index: {br_out.get("bridge_index")} '
                  f'({br_out.get("berge_conjecture_status")})')
            if 'bridge_index_berge_indicated' in br_out:
                print(f'    Berge-indicated: {br_out["bridge_index_berge_indicated"]}')

    report['optional_invariants'] = optional

    out_path = p if args.in_place else p.with_suffix('.optional.json')
    with out_path.open('w') as f:
        json.dump(report, f, indent=2, default=str)

    # Summary
    print(f'\n{"="*60}')
    print(f'OPTIONAL_INVARIANTS — {len(optional)} fields in total')
    print(f'{"="*60}')
    populated = sum(1 for v in optional.values()
                    if v is not None and v != [] and v != '')
    print(f'  populated: {populated} / {len(optional)}')
    print(f'  null/empty: {len(optional) - populated}')
    print(f'\nWritten: {out_path}')


if __name__ == '__main__':
    main()
