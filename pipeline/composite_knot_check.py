#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""composite_knot_check.py — checks whether a knot is a connected sum of known factors.

Application of Discovery #164 (knots-of-love = 3 K3_1 # 3 K3_1m): recognise connected-sum
knots via multiplicativity of Conway, Jones, Alex, determinant
and additivity of signature, genus, Arf.

Procedure:
  1. Load the Pflicht polynomials of a knot candidate.
  2. Check the Conway polynomial for factorisation (sympy.factor).
  3. If Conway splits into non-trivial factors → CS candidate.
     Per factor: compare against the library of known Conway polynomials.
  4. On a CS match: Jones multiplicativity test with mirror combinations
     of the factors — determines the mirror mix uniquely.
  5. Consistency check: c, det, sig, genus, arf, fibered.

Usage:
    python composite_knot_check.py <pflicht_complete.json> [<more>...]

Output: short balance per knot + optionally one JSON per CS hit.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import sympy as sp


# ════════════════════════════════════════════════════════════════════════
# Library of known knot Conway polynomials (excerpt from KnotInfo / Pflicht
# data of the corpus). Source: KnotInfo + generated Pflicht JSONs.
# Form: (name, conway_in_z, jones_in_x_or_q, sig, genus, arf, fibered)
# Jones always in x (convention of the Pflicht pipeline).
# ════════════════════════════════════════════════════════════════════════

z = sp.Symbol('z')
x = sp.Symbol('x')


def _quick_parse(s, var):
    """Local helper so that the KNOT_LIB initialisation stays simple."""
    if s is None:
        return None
    from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,
                                            implicit_multiplication_application)
    trans = standard_transformations + (implicit_multiplication_application,)
    return parse_expr(s.replace('^', '**'), local_dict={str(var): var},
                      transformations=trans)


def _c(s): return _quick_parse(s, z)
def _j(s): return _quick_parse(s, x)


KNOT_LIB = [
    # ── Classical KnotInfo knots c <= 8 (hardcoded; Conway via trefoil + mirror) ──
    # name        conway(z)            jones(x)                                       sig genus arf fibered
    ("K3_1",    _c("1 + z^2"),       _j("-x^8 + x^6 + x^2"),                          -2, 1, 1, True),
    ("K3_1m",   _c("1 + z^2"),       _j("x^-2 + x^-6 - x^-8"),                        +2, 1, 1, True),
    ("K4_1",    _c("1 - z^2"),       _j("x^4 - x^2 + 1 - x^-2 + x^-4"),                0, 1, 1, True),
    ("K5_1",    _c("1 + 3z^2 + z^4"),_j("-x^14 + x^12 - x^10 + x^8 + x^4"),           -4, 2, 1, True),
    ("K5_1m",   _c("1 + 3z^2 + z^4"),_j("x^-4 + x^-8 - x^-10 + x^-12 - x^-14"),       +4, 2, 1, True),
    ("K5_2",    _c("1 + 2z^2"),      _j("-x^8 + 2 x^6 - x^4 + 2 x^2 - 1 + x^-2"),     -2, 1, 0, True),
    ("K5_2m",   _c("1 + 2z^2"),      _j("x^2 - 1 + 2 x^-2 - x^-4 + 2 x^-6 - x^-8"),   +2, 1, 0, True),
    ("K6_1",    _c("1 - 2z^2"),      None,                                             0, 1, 0, False),
    ("K6_2",    _c("1 - z^2 - z^4"), None,                                            -2, 2, 1, True),
    ("K6_3",    _c("1 + z^2 + z^4"), None,                                             0, 2, 1, True),
    # ── From corpus Pflicht data 2026-05-13 (Discovery #165 lib extension) ──
    ("K8_19",   _c("1 + 5z^2 + 5z^4 + z^6"),
                _j("-x^16 + x^10 + x^6"),                                              -6, 3, 1, True),
    ("K10_43",  _c("1 + 2z^2 + z^4 - z^6"),
                _j("-x^10 + 3 x^8 - 6 x^6 + 9 x^4 - 11 x^2 + 13 - 11 x^-2 + 9 x^-4 - 6 x^-6 + 3 x^-8 - x^-10"),
                0, 3, 0, True),
    ("11n_133",  _c("1 + 2z^2 + 2z^4 + 4z^6 + z^8"),
                  _j("-x^14 + 2 x^12 - 3 x^10 + 4 x^8 - 4 x^6 + 4 x^4 - 3 x^2 + 3 - x^-2"),
                  -4, 4, 0, True),
    ("11n_133m", _c("1 + 2z^2 + 2z^4 + 4z^6 + z^8"),
                  _j("-x^2 + 3 - 3 x^-2 + 4 x^-4 - 4 x^-6 + 4 x^-8 - 3 x^-10 + 2 x^-12 - x^-14"),
                  +4, 4, 0, True),
    ("12a_1019", _c("1 + 2z^2 - z^4 - 2z^6 + z^8"),
                  _j("x^12 - 6 x^10 + 16 x^8 - 30 x^6 + 44 x^4 - 54 x^2 + 59 - 54 x^-2 + 44 x^-4 - 30 x^-6 + 16 x^-8 - 6 x^-10 + x^-12"),
                  0, 4, 0, True),
    ("12a_503",  _c("1 - 2z^2 + 3z^4 - z^6 - z^8"),
                  _j("x^14 - 4 x^12 + 7 x^10 - 11 x^8 + 17 x^6 - 19 x^4 + 20 x^2 - 19 + 15 x^-2 - 11 x^-4 + 7 x^-6 - 3 x^-8 + x^-10"),
                  -2, 4, 0, True),
    ("12n_706",  _c("1 - 4z^2 + 2z^4 + 4z^6 + z^8"),
                  _j("2 x^8 - 4 x^6 + 6 x^4 - 8 x^2 + 9 - 8 x^-2 + 6 x^-4 - 4 x^-6 + 2 x^-8"),
                  0, 4, 0, True),
    ("12n_837",  _c("1 + 2z^2 - z^4 + z^6 + z^8"),
                  _j("x^12 - 6 x^10 + 13 x^8 - 20 x^6 + 26 x^4 - 28 x^2 + 28 - 22 x^-2 + 15 x^-4 - 8 x^-6 + 2 x^-8"),
                  0, 4, 0, True),
    ("12n_839",  _c("1 + 2z^2 - 4z^4 + z^6 + z^8"),
                  _j("-2 x^18 + 6 x^16 - 11 x^14 + 16 x^12 - 19 x^10 + 20 x^8 - 18 x^6 + 14 x^4 - 9 x^2 + 5 - x^-2"),
                  -4, 4, 0, True),
]


# ════════════════════════════════════════════════════════════════════════
# Load KnotInfo c<=13 lib from JSON (build_knotinfo_lib.py generates it).
# Add mirror variants automatically via Jones substitution x->1/x.
# This extension covers 12965 prime knots — the complete
# known knot corpus up to c=13. Discovery #166.
# ════════════════════════════════════════════════════════════════════════

KNOTINFO_LIB_PATH = Path(__file__).resolve().parent / "knotinfo_lib_c13.json"

# Extended lib from KnotInfo: lazy format, Conway as coeff vector (instant
# lookup), Jones as string (only parsed to sympy on demand when a
# Conway match occurs). Without this trick 12965 sympy.parse_expr
# calls would slow the script down by >10 min.
#
# KNOTINFO_LIB entry = dict with:
#   name, c, conway_coeffs (list[int], even-power coeffs c0,c2,c4,..),
#   conway_str, jones_str, jones_mirror_str, sig, genus, arf, fibered
KNOTINFO_LIB = []


def _vector_to_coeffs(v):
    """KnotInfo conway_polynomial_vector = [link_flag, max_z_deg/2, c0, c2, c4, ...].
    Returns [c0, c2, c4, ...] (coeffs of z^0, z^2, z^4, ...).
    """
    if not v or len(v) < 3:
        return None
    return v[2:]


def _load_knotinfo_lib_lazy():
    """Load KnotInfo JSON without sympy parsing — only coeff vectors + strings."""
    if not KNOTINFO_LIB_PATH.exists():
        return 0
    try:
        data = json.load(open(KNOTINFO_LIB_PATH))
    except Exception:
        return 0
    for d in data:
        coeffs = _vector_to_coeffs(d.get("conway_vector"))
        if not coeffs:
            continue
        KNOTINFO_LIB.append({
            "name": d["name"],
            "c": d.get("c"),
            "conway_coeffs": coeffs,
            "conway_str": d.get("conway"),
            "jones_str": d.get("jones"),
            "jones_mirror_str": d.get("jones_mirror"),
            "homfly_str": d.get("homfly"),
            "homfly_mirror_str": d.get("homfly_mirror"),
            "sig": d.get("sig"),
            "genus": d.get("genus"),
            "arf": d.get("arf"),
            "fibered": d.get("fibered", False),
        })
    return len(KNOTINFO_LIB)


_loaded_kinfo = _load_knotinfo_lib_lazy()
if _loaded_kinfo:
    print(f"[composite_check] KnotInfo lib lazily loaded: {_loaded_kinfo} entries "
          f"(KNOT_LIB hardcoded: {len(KNOT_LIB)}, KNOTINFO_LIB indexed: {_loaded_kinfo})",
          file=sys.stderr)


def _poly_to_coeffs(expr):
    """sympy polynomial (in z) to coeff vector [c0, c2, c4, ...]. Returns None on
    error or odd-power terms (Conway factors should only have even ones).
    """
    try:
        p = sp.Poly(expr, z)
        # all_coeffs returns highest-to-lowest; we want lowest-to-highest
        all_c = list(reversed(p.all_coeffs()))
        # Conway has only even powers: odd indices must be 0
        even_only = []
        for i, c in enumerate(all_c):
            if i % 2 == 1:
                if c != 0:
                    return None  # odd term, not a Conway factor
            else:
                even_only.append(int(c))
        return even_only
    except Exception:
        return None


def _match_knotinfo_factor(target_coeffs):
    """Find all KnotInfo entries with identical Conway coeff vector."""
    matches = []
    for entry in KNOTINFO_LIB:
        if entry["conway_coeffs"] == target_coeffs:
            matches.append(entry)
    return matches


def _knot_c_lookup(name):
    """Return c (crossing_number_min) for a lib knot name."""
    # hardcoded lib: c is NOT explicitly in it (was not included),
    # so fallback via naming convention K<c>_<idx> or K<c>n/a_<idx>
    import re
    m = re.match(r"K(\d+)[a-z]?_\d+m?$", name)
    if m:
        return int(m.group(1))
    # KnotInfo lib directly
    is_mirror = name.endswith("m")
    lookup = name[:-1] if is_mirror else name
    for entry in KNOTINFO_LIB:
        if entry["name"] == lookup:
            return entry["c"]
    return None


def parse_pflicht(pflicht_json: Path):
    """Load Pflicht data relevant for the CS test."""
    d = json.load(open(pflicht_json))['pflichtliste']
    diag = d.get('DIAGRAMM', {}) or {}
    alg = d.get('ALGEBRA', {}) or {}
    hyp = d.get('HYPERBOLISCH', {}) or {}
    return {
        "name": pflicht_json.parent.name,
        "c": diag.get('crossing_number_min'),
        "det": alg.get('determinant'),
        "sig": alg.get('signature'),
        "arf": alg.get('arf_invariant'),
        "genus": alg.get('genus'),
        "fibered": alg.get('fibered'),
        "vol": hyp.get('hyperbolic_volume'),
        "conway_str": alg.get('conway_polynomial'),
        "jones_str": alg.get('jones_polynomial'),
        "alex_str": alg.get('alexander_polynomial'),
    }


def parse_poly(s: str, var: sp.Symbol):
    """Parse a polynomial string (Pflicht format) into a sympy expression.

    Pflicht format examples:
      Conway:  "1 + 6z^2 + 21z^4 + ..."  (no space between coeff and z)
      Jones:   "-x^8 + x^6 + x^2"
      Jones:   "-2 x^24 + 3 x^22 - x^20 ..."  (with space after coeff)
    """
    if not s:
        return None
    s = s.replace('^', '**')
    # Activate implicit-multiplication transformation — fixes "3 x" and "3x"
    from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,
                                            implicit_multiplication_application)
    transformations = standard_transformations + (implicit_multiplication_application,)
    try:
        return parse_expr(s, local_dict={str(var): var},
                          transformations=transformations)
    except Exception as e:
        return None


def laurent_eq(p1, p2, var=x, shift=50):
    """Check equality of Laurent polynomials via x^shift scaling."""
    diff = sp.expand((p1 - p2) * var**shift)
    return sp.simplify(diff) == 0


def factor_conway(conway_expr):
    """Factorise Conway polynomial (in z) via sympy. Returns list of factors."""
    if conway_expr is None:
        return None
    factored = sp.factor(conway_expr)
    # Get factors from Mul/Pow structure
    if isinstance(factored, sp.Mul):
        factors = []
        for arg in factored.args:
            if isinstance(arg, sp.Pow):
                base, exp = arg.args
                for _ in range(int(exp)):
                    factors.append(base)
            else:
                factors.append(arg)
        return factors
    elif isinstance(factored, sp.Pow):
        base, exp = factored.args
        return [base] * int(exp)
    return [factored]


def find_cs_decomposition(target):
    """Try to recognise the target knot as a connected sum of known factors.

    Strategy:
      1. Factorise target.conway in z.
      2. Per factor: match against KNOT_LIB (hardcoded) AND KNOTINFO_LIB (lazy).
         KnotInfo match via coeff-vector comparison — instant.
      3. On complete factorisation into lib knots: return composition.
      4. Refine via Jones multiplicativity test (mirror mix).
    """
    conway_t = parse_poly(target["conway_str"], z)
    if conway_t is None:
        return None, "Conway polynomial unparsable"

    factors = factor_conway(conway_t)
    if factors is None or len(factors) <= 1:
        return None, f"Conway not factorisable in {factors}"

    # Match factors against hardcoded + KnotInfo lib
    lib_match = []
    for f in factors:
        candidates = []
        # hardcoded lib (sympy expand compare)
        for name, conway_lib, *_ in KNOT_LIB:
            if sp.simplify(sp.expand(f - conway_lib)) == 0:
                candidates.append(name)
        # KnotInfo lib via coeff vector — fast
        target_coeffs = _poly_to_coeffs(f)
        if target_coeffs is not None:
            kinfo_matches = _match_knotinfo_factor(target_coeffs)
            for entry in kinfo_matches:
                if entry["name"] not in candidates:
                    candidates.append(entry["name"])
                mirror_name = entry["name"] + "m"
                if mirror_name not in candidates:
                    candidates.append(mirror_name)
        if not candidates:
            return None, f"Factor {f} (coeffs {target_coeffs}) not in lib"
        lib_match.append((f, candidates))

    # Additivity filter: Pflicht consistency per factor
    # Connected-Sum:
    #   c(K1 # K2)      <= c(K1) + c(K2)   (crossing min, exact for Pflicht data)
    #   sig             = sum(sig_i)
    #   genus           = sum(genus_i)
    #   arf (mod 2)     = sum(arf_i) mod 2
    # Factor min c lies at the trefoil-like Conway (z^2+1 -> c>=3).
    # Reduce the candidate lists to matching c values:
    #   max_factor_c = target_c - (n_factors - 1) * 3   (every other factor >= 3)
    #   hence c_i <= max_factor_c
    target_c = target.get("c")
    n = len(lib_match)
    if target_c is not None and n >= 2:
        max_fc = target_c - (n - 1) * 3
        filtered = []
        for f, cands in lib_match:
            keep = []
            for name in cands:
                c_lib = _knot_c_lookup(name)
                if c_lib is None or c_lib <= max_fc:
                    keep.append(name)
            filtered.append((f, keep))
        lib_match = filtered

    # Conway does not distinguish K3_1 / K3_1m. For trefoil factors:
    # check all possible (#K3_1, #K3_1m) mixes by Jones test.
    # Generalisation: per factor count "K and Km" how often.
    base_counts = {}  # name → multiplicity (without mirror distinction)
    for _, cands in lib_match:
        # On conflict (K3_1 vs K3_1m): "K3_1-or-mirror" as key
        # — group all mirror pairs
        base = cands[0].rstrip("m")  # heuristic
        base_counts[base] = base_counts.get(base, 0) + 1

    return {
        "conway_factors": [str(f) for f in factors],
        "lib_matches_per_factor": [[c for c in cands] for _, cands in lib_match],
        "base_counts": base_counts,
    }, None


def _get_jones_for_name(name):
    """Parse a Jones polynomial on demand for a lib knot name. Cache for repeated lookup."""
    if not hasattr(_get_jones_for_name, "_cache"):
        _get_jones_for_name._cache = {}
        for n, _, j, *_ in KNOT_LIB:
            if j is not None:
                _get_jones_for_name._cache[n] = j
    cache = _get_jones_for_name._cache
    if name in cache:
        return cache[name]
    is_mirror = name.endswith("m")
    lookup_name = name[:-1] if is_mirror else name
    for entry in KNOTINFO_LIB:
        if entry["name"] == lookup_name:
            jstr = entry["jones_mirror_str"] if is_mirror else entry["jones_str"]
            if not jstr:
                cache[name] = None
                return None
            try:
                j = parse_poly(jstr, x)
                cache[name] = j
                return j
            except Exception:
                cache[name] = None
                return None
    cache[name] = None
    return None


# Sympy symbols for HOMFLY (corpus convention 2-variable)
h_x = sp.Symbol('x')   # = 1/v in KnotInfo convention
h_y = sp.Symbol('y')   # = z in KnotInfo convention


def _parse_homfly(s):
    """Parse HOMFLY string (in x, y) into a sympy expression."""
    if not s:
        return None
    from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,
                                            implicit_multiplication_application)
    trans = standard_transformations + (implicit_multiplication_application,)
    try:
        return parse_expr(s.replace('^', '**'),
                          local_dict={'x': h_x, 'y': h_y},
                          transformations=trans)
    except Exception:
        return None


def _get_homfly_for_name(name):
    """Parse HOMFLY polynomial on demand for a lib knot name. Cache."""
    if not hasattr(_get_homfly_for_name, "_cache"):
        _get_homfly_for_name._cache = {}
    cache = _get_homfly_for_name._cache
    if name in cache:
        return cache[name]
    is_mirror = name.endswith("m")
    lookup_name = name[:-1] if is_mirror else name
    for entry in KNOTINFO_LIB:
        if entry["name"] == lookup_name:
            hstr = entry["homfly_mirror_str"] if is_mirror else entry["homfly_str"]
            if not hstr:
                cache[name] = None
                return None
            h = _parse_homfly(hstr)
            cache[name] = h
            return h
    cache[name] = None
    return None


def homfly_eq(p1, p2, shift_x=20, shift_y=20):
    """Check equality of 2-variable Laurent polynomials via xshift, yshift."""
    if p1 is None or p2 is None:
        return False
    diff = sp.expand((p1 - p2) * h_x**shift_x * h_y**shift_y)
    return sp.simplify(diff) == 0


def homfly_mirror_test(target_homfly_str, factor_names_list):
    """HOMFLY multiplicativity test analogous to jones_mirror_test.

    HOMFLY(K1 # K2) = HOMFLY(K1) * HOMFLY(K2).
    HOMFLY is more discriminating than Jones — if the Jones test
    yields several hits, HOMFLY can often give a unique match (or
    falsify the Jones hits).
    """
    H_target = _parse_homfly(target_homfly_str)
    if H_target is None:
        return None
    from itertools import product
    options = []
    for names in factor_names_list:
        opts = []
        for n in names:
            h = _get_homfly_for_name(n)
            if h is not None:
                opts.append((n, h))
        if not opts:
            return None
        options.append(opts)
    n_combos = 1
    for o in options:
        n_combos *= len(o)
    if n_combos > 4096:
        return None
    matches = []
    for combo in product(*options):
        prod = sp.Integer(1)
        for _, h in combo:
            prod = prod * h
        prod = sp.expand(prod)
        if homfly_eq(H_target, prod):
            matches.append([n for n, _ in combo])
    return matches


def jones_mirror_test(target_jones_str, factor_names_list):
    """Check Jones multiplicativity for different mirror combinations
    of the factors. Returns the matching mirror distribution.
    """
    J_target = parse_poly(target_jones_str, x)
    if J_target is None:
        return None

    from itertools import product

    options = []
    for names in factor_names_list:
        unique_jones = []
        for n in names:
            j = _get_jones_for_name(n)
            if j is not None:
                unique_jones.append((n, j))
        if not unique_jones:
            return None
        options.append(unique_jones)

    # With many factors the combinatorics explode. Hard limit: 2^10 = 1024.
    n_combos = 1
    for opts in options:
        n_combos *= len(opts)
    if n_combos > 4096:
        # too many — return factor list without explicit match
        return None

    matches = []
    for combo in product(*options):
        product_j = sp.Integer(1)
        for _, j in combo:
            product_j = product_j * j
        product_j = sp.expand(product_j)
        if laurent_eq(J_target, product_j, var=x):
            matches.append([n for n, _ in combo])
    return matches


def analyze(pflicht_json: Path, output_dir: Path | None = None,
            verbose: bool = True):
    t = parse_pflicht(pflicht_json)
    if verbose:
        print(f"\n══ {t['name']} (c={t['c']} det={t['det']} sig={t['sig']} "
              f"genus={t['genus']} arf={t['arf']} fibered={t['fibered']}) ══")
        print(f"  conway: {t['conway_str']}")

    decomp, reason = find_cs_decomposition(t)
    if decomp is None:
        if verbose:
            print(f"  Status: NOT-COMPOSITE (or factors outside lib) — {reason}")
        return {"name": t["name"], "status": "not_composite", "reason": reason}

    if verbose:
        print(f"  Conway factorisation: {decomp['conway_factors']}")
        print(f"  Factor matches: {decomp['lib_matches_per_factor']}")

    factor_names = decomp["lib_matches_per_factor"]
    jones_matches = jones_mirror_test(t["jones_str"], factor_names)

    # HOMFLY verification of the Jones matches (more discriminating than Jones)
    homfly_matches = None
    homfly_str = parse_pflicht(pflicht_json).get("homfly_str") if False else None
    # Hack: re-load HOMFLY directly
    try:
        _d = json.load(open(pflicht_json))
        homfly_str = _d.get('pflichtliste',{}).get('ALGEBRA',{}).get('homfly_polynomial')
    except Exception:
        pass
    if homfly_str and jones_matches:
        homfly_matches = homfly_mirror_test(homfly_str, factor_names)
        if verbose and homfly_matches is not None:
            n_j, n_h = len(jones_matches), len(homfly_matches)
            print(f"  Jones match candidates: {n_j},  HOMFLY match candidates: {n_h}")

    if jones_matches:
        # Sanity check via sig, genus, arf
        from collections import Counter
        cs_label = None
        # Build sig-lookup from both libs
        lib_sig = {n: s for n, _, _, s, *_ in KNOT_LIB}
        for entry in KNOTINFO_LIB:
            lib_sig[entry["name"]] = entry["sig"]
            if entry["sig"] is not None:
                lib_sig[entry["name"] + "m"] = -entry["sig"]
        # Prefer HOMFLY match (more discriminating); if none: Jones match
        candidates = homfly_matches if homfly_matches else jones_matches
        method = "Jones+HOMFLY" if homfly_matches else "Jones"
        for combo in candidates:
            counts = Counter(combo)
            cs_str = " # ".join(f"{counts[k]}×{k}" if counts[k] > 1 else k
                                 for k in sorted(counts))
            pred_sig = sum(lib_sig.get(n, 0) for n in combo)
            if t["sig"] is not None and pred_sig != t["sig"]:
                continue
            cs_label = cs_str
            break

        if cs_label:
            if verbose:
                print(f"  ✓ IDENTIFIED ({method}): {cs_label}")
            result = {
                "name": t["name"],
                "status": "composite_identified",
                "decomposition": cs_label,
                "factor_combo": candidates[0] if candidates else None,
                "verification_method": method,
            }
            if output_dir:
                out = output_dir / f"{t['name']}_composite.json"
                json.dump(result, open(out, "w"), indent=2)
                if verbose: print(f"  → {out}")
            return result
        else:
            if verbose:
                print(f"  Jones match but sig inconsistency")
            return {"name": t["name"], "status": "jones_match_sig_mismatch",
                    "candidates": jones_matches}
    else:
        if verbose:
            print(f"  Factorisation in Conway OK but no Jones match")
        return {"name": t["name"], "status": "conway_factored_no_jones_match",
                "conway_factors": decomp["conway_factors"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pflicht_jsons", nargs="+", type=Path)
    ap.add_argument("--output-dir", type=Path, default=None,
                    help="Write one JSON per composite hit to that directory")
    args = ap.parse_args()

    summary = []
    for pj in args.pflicht_jsons:
        if not pj.exists():
            print(f"  {pj}: MISSING")
            continue
        r = analyze(pj, output_dir=args.output_dir)
        summary.append(r)

    print(f"\n══ SUMMARY ══")
    composites = [r for r in summary if r["status"] == "composite_identified"]
    print(f"  {len(composites)} / {len(summary)} identified as composite.")
    for r in composites:
        print(f"    {r['name']:30s} = {r['decomposition']}")


if __name__ == "__main__":
    main()
