# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
pipeline_crosscheck.py — database cross-check for one mandatory-invariants (Pflicht) JSON.

EXPLICITLY NOT blind-test compliant: this script is allowed to look up knot tables
(KnotInfo via database_knotinfo, Fremlin cache, snappy census). The Pflicht JSONs
stay untouched; the output is a SEPARATE <knot>_crosscheck.json.

IMPORTANT: NO seal_all() import here — crosscheck is allowed table lookups. snappy
is called via a subprocess in `snappy_env`, where it offers unsealed Manifold
construction and .identify().

Usage:
  python pipeline_crosscheck.py <pflicht_complete.json> [--output <crosscheck.json>]

Match paths:
  (a) KnotInfo (complete for c<=13): (c, det, sig) index, then polynomial filter
      with a variable-agnostic normaliser (Discovery #134: pipeline x vs
      KnotInfo t).
  (b) Fremlin DHF list (37 knots 3_1..8_21 + 12a_1202 + 15331): via
      Rolfsen name match (from KnotInfo) + Jones match (15331 explicitly).
  (c) snappy census (phase 2b): OrientableCuspedCensus.identify(M) via
      the triangulation isosig from raw_subprocess_outputs.snappy.isosig. Covers
      LinkExteriors / CensusKnots / HTLinkExteriors (c=14..16). Match criterion:
      volume tolerance 1e-9 (Discovery #129: isosig is unstable under D6/D3
      symmetry → volume primary, isosig secondary as tie-breaker).

Output: dict with input_invariants, knotinfo_match, fremlin_match,
        snappy_census_match, summary, novelty_status.
"""
from __future__ import annotations
import argparse
import site_config
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

KI_INDEX = Path(site_config.KNOTINFO_INDEX)
KI_RECORDS = Path(site_config.KNOTINFO_RECORDS)
FREMLIN_CACHE = Path(site_config.FREMLIN_CACHE)


def normalize_poly(s: str) -> str:
    """Canonicalise a polynomial string for a bit match: remove all whitespace,
    unify multiplications, drop a leading plus sign."""
    if not s: return ""
    t = str(s)
    t = re.sub(r"\s+", "", t)
    t = t.replace("*", "")
    t = t.lstrip("+")
    return t


def jones_substitute_t_to_x_squared(jones_t: str) -> str:
    """Legacy coefficient sequence (very blunt!). DEPRECATED since Discovery
    #134 — only kept as a fallback. Prefer poly_canonical().
    """
    if not jones_t: return ""
    coeffs = re.findall(r"([+-]?\d+)", jones_t)
    return "|".join(coeffs)


# ----------------------------------------------------------------------------
# ROBUST POLYNOMIAL NORMALISER (Discovery #134)
# ----------------------------------------------------------------------------
# Problem: pipeline Regina gives Jones in variable x (half-variable convention),
# KnotInfo in variable t. For Jones, classically x^2 == t. The old
# coefficient-only sequence was too blunt (Discovery #134 / celtic-for-playing
# + double-6: 14 KnotInfo candidates for c=12,det=49,sig=0).
#
# Solution: parse a Laurent polynomial in ONE variable into a sorted
# (exp_relative, coef) list. exp_relative := exp - min_exp normalises
# Laurent monomial shifts. For the Jones match, scaling
# (exponent_scale=2) is tried in addition (x^N <-> t^(N/2)).
# ----------------------------------------------------------------------------

# Regex: recognises terms of the form [+/-] [coef] [* var [^exp | **exp]] and also
# fractions `coef/var^exp`. We parse token by token.
_TERM_RE = re.compile(
    r"""(?P<sign>[+-]?)\s*
        (?:(?P<coef>\d+)\s*\*?\s*)?         # optional coefficient
        (?:                                  # either a normal variable term
            (?P<var>[a-zA-Z])
            (?:\s*\*\*\s*|\s*\^\s*)?         # ** or ^ or directly
            (?P<expnum>-?\d+)?
        |
            (?P<bare>)                        # or a bare constant
        )""",
    re.VERBOSE,
)


def _split_terms(s: str):
    """Splits a polynomial string into terms. Handles both
    'a*x^2 - b*x + c' and 'a*x**2 - b/x^2 + c'.
    Also the KnotInfo convention 't^(-8)' (parentheses around negative exponents)."""
    if not s:
        return []
    # Remove whitespace
    t = re.sub(r"\s+", "", str(s))
    # Unwrap parenthesised exponents: x^(-N) → x^-N, x^(+N) → x^N, x^(N) → x^N
    t = re.sub(r"(\^|\*\*)\((-?\+?\d+)\)", r"\1\2", t)
    # Fractions 1/x^N or C/x^N → C*x^-N; also C*z**N/a**M → C*z**N*a**-M
    # Step by step: every '/var(^|**)N' becomes '*var^-N'
    def _fraction(m):
        var, op, exp = m.group(1), m.group(2), m.group(3)
        return f"*{var}^-{exp}"
    t = re.sub(r"/([a-zA-Z])(\*\*|\^)(\d+)", _fraction, t)
    # Fraction without exponent: /x → *x^-1
    t = re.sub(r"/([a-zA-Z])(?![\^\*\d])", r"*\1^-1", t)
    # Split terms at '+' or '-' (not at the start, not directly after ^/** )
    # We mark term boundaries, then split.
    t = re.sub(r"(?<!\^)(?<!\*\*)(?<=[^\^])([+-])", r"|\1", t)
    # If the first character is not a sign, prepend '+'
    t = t.lstrip("|")
    if not t.startswith(("+", "-")):
        t = "+" + t
    # Split
    parts = [p for p in t.split("|") if p]
    return parts


def _parse_term(term: str):
    """Parse a single term such as '+2*x^4' or '-x^-2' or '+9' or
    '+a^4*z^2' or '2x^8' (no *). Returns (coef, ((var, exp),...)).
    """
    t = term.strip()
    sign = 1
    if t.startswith("+"):
        t = t[1:]
    elif t.startswith("-"):
        sign = -1
        t = t[1:]
    coef = 1
    var_exp = {}  # var -> exp
    # Token scan: recognise either a pure number, or a var(^|**)(-?\d+)?
    # Patterns are tried in order.
    token_re = re.compile(
        r"(?P<num>\d+)"
        r"|(?P<var>[a-zA-Z])(?:\*\*|\^)?(?P<exp>-?\d+)?"
        r"|(?P<sep>\*)"
    )
    pos = 0
    while pos < len(t):
        m = token_re.match(t, pos)
        if not m:
            return None
        if m.group("num") is not None:
            coef *= int(m.group("num"))
        elif m.group("var") is not None:
            v = m.group("var")
            e = int(m.group("exp")) if m.group("exp") is not None else 1
            var_exp[v] = var_exp.get(v, 0) + e
        # sep is ignored
        pos = m.end()
    return sign * coef, tuple(sorted(var_exp.items()))


def _gcd(a, b):
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a


def poly_canonical(s: str, normalize_shift: bool = True, gcd_rescale: bool = True):
    """Parse a polynomial string into a canonical (exp, coef) list in ONE
    variable. Multivariable polynomials (HOMFLY, Kauffman) -> None.

    normalize_shift=True subtracts min_exp -> invariant under
    Laurent monomial shift (Sage often gives t-polynomials without normalisation;
    the pipeline normalises to the balanced Jones).

    gcd_rescale=True divides all (shifted) exponents by their GCD ->
    invariant under the variable substitution x -> t^(1/k). This makes
    pipeline Jones in x with exponents {0,2,4,...,16} (Regina convention) match
    KnotInfo Jones in t with exponents {0,1,2,...,8} (Discovery #134).

    Returns: tuple of (exp_int, coef_int) sorted by exp, or None on parse fail
    or multivariable.
    """
    if not s:
        return None
    parts = _split_terms(s)
    coef_by_exp = {}
    for term in parts:
        parsed = _parse_term(term)
        if parsed is None:
            return None
        coef, var_exp = parsed
        if coef == 0:
            continue
        if len(var_exp) == 0:
            exp = 0
        elif len(var_exp) == 1:
            _, exp = var_exp[0]
        else:
            # Multivariable → not covered by this function
            return None
        coef_by_exp[exp] = coef_by_exp.get(exp, 0) + coef
    # Drop zero coefficients
    coef_by_exp = {e: c for e, c in coef_by_exp.items() if c != 0}
    if not coef_by_exp:
        return tuple()
    if normalize_shift:
        min_e = min(coef_by_exp.keys())
        coef_by_exp = {(e - min_e): c for e, c in coef_by_exp.items()}
    if gcd_rescale and len(coef_by_exp) > 1:
        exps_nonzero = [e for e in coef_by_exp.keys() if e != 0]
        if exps_nonzero:
            g = exps_nonzero[0]
            for e in exps_nonzero[1:]:
                g = _gcd(g, e)
            if g > 1:
                coef_by_exp = {(e // g): c for e, c in coef_by_exp.items()}
    return tuple(sorted(coef_by_exp.items()))


def jones_canonical_with_mirror(s: str):
    """Returns (canonical_form, mirror_form) for Jones. Mirror corresponds to
    t -> t^-1 (or x -> x^-1): negate all exponents and shift again."""
    base = poly_canonical(s, normalize_shift=True)
    if base is None:
        return None, None
    if not base:
        return base, base
    # Mirror: negate exponents, then shift again
    mirrored_raw = [(-e, c) for (e, c) in base]
    min_e = min(e for e, _ in mirrored_raw)
    mirror = tuple(sorted((e - min_e, c) for (e, c) in mirrored_raw))
    return base, mirror


def jones_equal(pipeline_jones: str, knotinfo_jones: str) -> bool:
    """Robust Jones comparison. Tries (a) direct canonicalisation,
    (b) mirror comparison. Works variable-agnostically (x vs t).

    Note: the pipeline convention is x with exponents twice as large as KnotInfo's t
    (classically x^2 = t). Since poly_canonical normalises the exponent shift AND
    keeps the coefficient order, the form for {x^8,x^6,...,x^-8} is
    identical to {t^4,t^3,...,t^-4} *as an ordered-by-exponent sequence*: 9 terms with the
    same coefficients in the same order.
    """
    p = poly_canonical(pipeline_jones)
    k_can, k_mir = jones_canonical_with_mirror(knotinfo_jones)
    if p is None or k_can is None:
        return False
    return p == k_can or p == k_mir


def alex_equal(pipeline_alex: str, knotinfo_alex: str) -> bool:
    """Alexander polynomial comparison. Alexander is only determined up to
    multiplication by +/- t^k → we compare the canonicalised form (shift-normalised)
    AND the mirror (exponent negation).
    """
    p = poly_canonical(pipeline_alex)
    k_can, k_mir = jones_canonical_with_mirror(knotinfo_alex)
    if p is None or k_can is None:
        return False
    if p == k_can or p == k_mir:
        return True
    # Additionally: Alexander may be globally multiplied by -1
    neg = tuple((e, -c) for (e, c) in p)
    neg_sorted = tuple(sorted(neg))
    if neg_sorted == k_can or neg_sorted == k_mir:
        return True
    return False


def load_mandatory(path: Path) -> dict:
    d = json.loads(path.read_text())
    pl = d["pflichtliste"]
    D, A, H = pl["DIAGRAMM"], pl["ALGEBRA"], pl["HYPERBOLISCH"]
    snappy_raw = d.get("raw_subprocess_outputs", {}).get("snappy", {}) or {}
    # genus: prefer 3_genus, else seifert_genus, else genus
    genus = A.get("3_genus") or A.get("seifert_genus") or A.get("genus")
    try:
        genus = int(genus) if genus is not None else None
    except (TypeError, ValueError):
        genus = None
    return {
        "name": d.get("name"),
        "c": D.get("crossing_number_min"),
        "alt": D.get("alternating"),
        "det": A.get("determinant"),
        "sig": A.get("signature"),
        "arf": A.get("arf_invariant"),
        "fib": A.get("fibered"),
        "genus": genus,
        "vol": H.get("hyperbolic_volume"),
        "sym": H.get("symmetry_group"),
        "amphi": H.get("is_amphicheiral"),
        "jones_str": str(A.get("jones_polynomial", "")),
        "alex_str":  str(A.get("alexander_polynomial", "")),
        "homfly_str": str(A.get("homfly_polynomial", "")),
        "isosig": snappy_raw.get("isosig") or H.get("canonical_retriangulation"),
        "pd_list_str": d.get("pd_list_str", ""),
        "_source": str(path),
    }


def _norm_sym(s):
    """Normalise symmetry-group strings for comparison: remove all whitespace,
    the plus notation '+' stays. 'Z/2+ Z/2' == 'Z/2 + Z/2' == 'Z/2+Z/2'."""
    if s is None: return None
    return "".join(str(s).split())


def _filter_by_aux_invariants(inv: dict, names: list, records: dict, mirror: bool = False) -> tuple[list, dict]:
    """Additionally filter candidates by fibered/genus/arf/alt/sym_group.
    Records are only checked against inv when BOTH values are known
    (None lets everything through — silent DB entries must not
    be filtered out).

    For a mirror match, invariant properties (fibered, genus, arf, alt,
    sym_group) are compared as equal — they are all mirror-invariant.

    sym_group (Discovery #174) disambiguates mutant pairs such as 13a_1231/1237
    whose algebraic invariants (Alex, Jones, HOMFLY, Khovanov,
    theta) are all identical but whose full_symmetry_group differs.

    Returns (kept, info dict with filter statistics).
    """
    fib = inv.get("fib")
    genus = inv.get("genus")
    arf = inv.get("arf")
    alt = inv.get("alt")
    sym = _norm_sym(inv.get("sym"))
    kept = []
    fib_drops = 0
    g_drops = 0
    arf_drops = 0
    alt_drops = 0
    sym_drops = 0
    for n in names:
        r = records.get(n, {})
        if fib is not None and r.get("fibered") is not None and bool(fib) != bool(r["fibered"]):
            fib_drops += 1
            continue
        if genus is not None and r.get("genus") is not None and int(genus) != int(r["genus"]):
            g_drops += 1
            continue
        if arf is not None and r.get("arf") is not None and int(arf) != int(r["arf"]):
            arf_drops += 1
            continue
        if alt is not None and r.get("alt") is not None and bool(alt) != bool(r["alt"]):
            alt_drops += 1
            continue
        if sym is not None and r.get("full_sym_group"):
            r_sym = _norm_sym(r["full_sym_group"])
            if r_sym is not None and sym != r_sym:
                sym_drops += 1
                continue
        kept.append(n)
    info = {
        "kept": len(kept),
        "dropped_fibered": fib_drops,
        "dropped_genus": g_drops,
        "dropped_arf": arf_drops,
        "dropped_alternating": alt_drops,
        "dropped_sym_group": sym_drops,
    }
    return kept, info


def _polynomial_filter(cand: list, inv: dict, records: dict) -> tuple[list, list]:
    """Chains the Jones and Alex filters. Returns (after_jones, after_jones_and_alex).
    If jones_str/alex_str is empty or no candidate has the polynomial, the list
    passes through (all candidates survive the respective step).
    """
    j_filtered = list(cand)
    if inv.get("jones_str"):
        j = [n for n in cand if jones_equal(inv["jones_str"], records.get(n, {}).get("jones", ""))]
        if j:
            j_filtered = j
    ja_filtered = list(j_filtered)
    if len(j_filtered) > 1 and inv.get("alex_str"):
        a = [n for n in j_filtered if alex_equal(inv["alex_str"], records.get(n, {}).get("alex", ""))]
        if a:
            ja_filtered = a
    return j_filtered, ja_filtered


def _try_match_branch(cand: list, inv: dict, records: dict, mirror: bool) -> dict:
    """Tries the filter pipeline (poly + aux) on a candidate list.
    Returns a dict with intermediate results + optional best_match + confidence + status.
    If not resolvable: empty best_match, status='ambiguous after ...'.
    """
    info = {}
    if not cand:
        return {"best_match": None, "status": "no candidate"}

    aux_filtered, aux_info = _filter_by_aux_invariants(inv, cand, records, mirror=mirror)
    info["aux_filter"] = aux_info
    # If the aux filter throws everything away: do not reduce, keep the original.
    work = aux_filtered if aux_filtered else cand
    if not aux_filtered:
        info["aux_filter_empty"] = True

    j_filtered, ja_filtered = _polynomial_filter(work, inv, records)
    if len(work) != len(cand) or len(j_filtered) != len(work):
        info["after_aux"] = list(work)
        info["after_jones"] = list(j_filtered)
    if len(ja_filtered) != len(j_filtered):
        info["after_alex"] = list(ja_filtered)

    final = ja_filtered
    tag = "mirror " if mirror else ""

    if len(final) == 1:
        levels = []
        if any(aux_info.get(k) for k in ("dropped_fibered", "dropped_genus", "dropped_arf", "dropped_alternating")):
            levels.append("aux")
        if inv.get("jones_str"):
            levels.append("jones")
        if len(j_filtered) > 1 and inv.get("alex_str"):
            levels.append("alex")
        levels_s = "+".join(["c", "det", "sig"] + levels) if levels else "c+det+sig"
        info["best_match"] = final[0]
        info["confidence"] = f"high ({tag}c+det+sig+{'+'.join(levels) if levels else 'unique'})"
        info["status"] = "matched (mirror)" if mirror else "matched"
        if mirror:
            info["mirror_of"] = final[0]
        return info

    if not final:
        info["status"] = f"{tag}no candidate survives filters"
        return info

    # Several survivors — if the original length was 1 and it was not filtered out
    if len(cand) == 1:
        info["best_match"] = cand[0]
        info["confidence"] = f"medium ({tag}c+det+sig only)"
        info["status"] = "matched (mirror)" if mirror else "matched"
        if mirror:
            info["mirror_of"] = cand[0]
        return info

    info["status"] = f"{tag}ambiguous ({len(final)} candidates after filters)"
    info["surviving_candidates"] = list(final)
    return info


def knotinfo_match(inv: dict, idx: dict, records: dict) -> dict:
    c, det, sig = inv["c"], inv["det"], inv["sig"]
    out = {"queried_c_det_sig": [c, det, sig]}

    if c is None or det is None:
        out["status"] = "skipped (missing c or det)"
        return out
    # Unknot special case: c=0, det=1 → K0_1 regardless of sig
    if c == 0 and det == 1:
        out["candidates_by_c_det_sig"] = ["K0_1"]
        out["best_match"] = "K0_1"
        out["confidence"] = "high (Unknot c=0,det=1)"
        out["status"] = "matched (Unknot)"
        return out
    if sig is None:
        out["status"] = "skipped (signature missing — c={}, det={})".format(c, det)
        return out

    # Direct index lookup
    key = f"{c}|{det}|{sig}"
    cand = idx.get(key, [])
    out["candidates_by_c_det_sig"] = cand

    # If c > 13: KnotInfo does not cover this completely
    if c > 13:
        out["status"] = f"out of KnotInfo range (c={c} > 13, KnotInfo only fully covers c<=13)"
        out["coverage_warning"] = True
        out["partial_candidates"] = cand or []
        return out

    # Main branch: non-mirror candidates
    direct = _try_match_branch(cand, inv, records, mirror=False) if cand else {"status": "no candidate"}
    out["direct"] = direct
    if direct.get("best_match"):
        out["best_match"] = direct["best_match"]
        out["confidence"] = direct["confidence"]
        out["status"] = direct["status"]
        return out

    # Mirror branch: -sig (the cross-check is also tried when the direct branch
    # is ambiguous — direct candidates with the same sig=-x can e.g.
    # all be fibered=False, while the actual knot is fibered=True
    # and only lies in the mirror class +sig; Discovery #wave6 kf73 bug).
    key_m = f"{c}|{det}|{-sig}"
    cand_m = idx.get(key_m, [])
    if cand_m:
        out["mirror_candidates_by_c_det_neg_sig"] = cand_m
        mirror_res = _try_match_branch(cand_m, inv, records, mirror=True)
        out["mirror"] = mirror_res
        if mirror_res.get("best_match"):
            out["best_match"] = mirror_res["best_match"]
            out["confidence"] = mirror_res["confidence"]
            out["status"] = mirror_res["status"]
            out["mirror_of"] = mirror_res["best_match"]
            return out

    # Neither direct nor mirror unique
    direct_count = len(direct.get("surviving_candidates", cand)) if cand else 0
    mirror_count = len(out.get("mirror", {}).get("surviving_candidates", cand_m)) if cand_m else 0
    if not cand and not cand_m:
        out["status"] = "no KnotInfo candidate"
    else:
        out["status"] = (
            f"ambiguous (direct={direct_count}, mirror={mirror_count})"
        )
    return out


def fremlin_match(inv: dict, ki_result: dict, fremlin: list) -> dict:
    out = {"status": "no match"}
    best = ki_result.get("best_match")
    if best:
        # Strip the K prefix: K3_1 → 3_1
        bare = best.lstrip("K")
        for r in fremlin:
            if r["knot_id"] == bare or r.get("title_id") == bare.replace("_", ""):
                out["status"] = "matched via knotinfo-name"
                out["knot_id"] = r["knot_id"]
                out["url"] = r["source_url"]
                out["fremlin_symmetries"] = r.get("symmetries_raw")
                return out
    # Direct Jones match (15331) via the robust normaliser
    if inv["jones_str"]:
        for r in fremlin:
            if "jones_t_raw" not in r:
                continue
            if jones_equal(inv["jones_str"], r["jones_t_raw"]):
                out["status"] = "matched via jones-polynomial"
                out["knot_id"] = r["knot_id"]
                out["url"] = r["source_url"]
                out["fremlin_symmetries"] = r.get("symmetries_raw")
                return out
    return out


# ----------------------------------------------------------------------------
# snappy census match (phase 2b, Discovery #133/#134)
# ----------------------------------------------------------------------------
# Called via a subprocess in `snappy_env`: snappy is not installed in the current
# pipeline env, and the current env has pipeline_seal active (not in the
# crosscheck script, admittedly — but snappy always lives in snappy_env). So we call
# conda run -n snappy_env python -c '...' and parse JSON from stdout.
#
# Match criteria (prioritised):
#   1. volume tolerance 1e-9 between pflicht.vol and the identify() manifold's vol
#   2. isosig match as an additional tie-breaker (Discovery #129 warns: under
#      D6/D3 symmetry the isosig can be unstable, hence not a hard match
#      criterion)
# ----------------------------------------------------------------------------

SNAPPY_SUBPROCESS_TEMPLATE = r"""
import json, re, sys
try:
    import snappy
except Exception as e:
    print(json.dumps({"ok": False, "error": "snappy import failed: " + str(e)}))
    sys.exit(0)

isosig = __ISOSIG__
expected_vol = __VOL__
expected_c = __C__
tol = __TOL__
# Extract the crossing count from snappy names: "12n_706", "12n_487", "11a_367"
_NAME_C_RE = re.compile(r"^[KL]?(\d+)[an]")

def _name_to_c(name):
    m = _NAME_C_RE.match(str(name))
    return int(m.group(1)) if m else None

result = {"ok": True, "isosig_used": isosig, "expected_vol": expected_vol}
try:
    M = snappy.Manifold(isosig)
    # Disc #327 (census gate): randomise into a geometric shape BEFORE identify(). Otherwise
    # identify() misses a known census manifold -> that is how wz_sqrt43 = K10_934 (a known
    # ~c10 in a c=36 diagram) slipped through wrongly as novel. census_sweep did this, the crosscheck did not.
    for _ in range(30):
        try:
            if M.solution_type() == 'all tetrahedra positively oriented':
                break
            M.randomize()
        except Exception:
            break
    result["solution_type"] = str(M.solution_type())
    result["constructed_vol"] = float(M.volume())
    result["constructed_isosig"] = M.triangulation_isosig()
    if expected_vol is not None:
        result["vol_diff"] = abs(result["constructed_vol"] - expected_vol)
        result["vol_match"] = result["vol_diff"] < tol
    else:
        result["vol_diff"] = None
        result["vol_match"] = None
    # identify() checks against all census DBs internally (OrientableCuspedCensus,
    # LinkExteriors, CensusKnots, HTLinkExteriors).
    try:
        ids = M.identify()
        result["identify"] = [str(x) for x in ids]
    except Exception as e:
        result["identify"] = []
        result["identify_error"] = str(e)
    # HTLinkExteriors volume scan for c=14..16 if identify is empty.
    # ONLY for hyperbolic knots: at vol≈0 (torus knots, composites with
    # non-hyperbolic factors) every non-hyperbolic census
    # entry matches, the scan is worthless. Discovery #169.
    if (not result["identify"]
            and expected_vol is not None
            and expected_vol > 1e-6):
        ht_matches = []
        ht_dropped_by_c = 0
        try:
            for L in snappy.HTLinkExteriors:
                try:
                    if L.num_cusps() != 1:
                        continue
                    v = float(L.volume())
                except Exception:
                    continue
                if abs(v - expected_vol) < tol:
                    name = L.name()
                    # c filter: only candidates with the same crossing count
                    if expected_c is not None:
                        cand_c = _name_to_c(name)
                        if cand_c is not None and cand_c != expected_c:
                            ht_dropped_by_c += 1
                            continue
                    ht_matches.append({"name": name, "vol": v})
                    if len(ht_matches) >= 8:
                        break
            result["ht_volume_scan"] = ht_matches
            if ht_dropped_by_c:
                result["ht_dropped_by_c"] = ht_dropped_by_c
        except Exception as e:
            result["ht_volume_scan_error"] = str(e)
    elif expected_vol is not None and expected_vol <= 1e-6:
        result["ht_volume_scan_skipped"] = "vol≈0 (non-hyperbolic; volume match meaningless)"
except Exception as e:
    result["ok"] = False
    result["error"] = "Manifold construction or identify failed: " + str(e)

print(json.dumps(result))
"""


def snappy_census_match(inv: dict, env_name: str = "snappy_env", tol: float = 1e-9, timeout_s: int = 300) -> dict:
    """For pipeline output with an isosig: call snappy in a snappy_env subprocess,
    fetch the identify() result and the volume match.

    Returns: dict with status, snappy outputs, best_match_name (if unique).
    """
    out = {"status": "skipped (no isosig in pflicht-JSON)"}
    isosig = inv.get("isosig")
    if not isosig:
        return out

    # Locate the conda binary
    conda = shutil.which("conda") or site_config.CONDA_BIN
    if not Path(conda).exists():
        out["status"] = f"skipped (conda not found at {conda})"
        return out

    expected_vol = inv.get("vol")
    expected_c = inv.get("c")
    vol_repr = repr(expected_vol) if isinstance(expected_vol, (int, float)) else "None"
    c_repr = repr(expected_c) if isinstance(expected_c, int) else "None"
    code = (SNAPPY_SUBPROCESS_TEMPLATE
            .replace("__ISOSIG__", repr(isosig))
            .replace("__VOL__", vol_repr)
            .replace("__C__", c_repr)
            .replace("__TOL__", repr(tol)))

    try:
        proc = subprocess.run(
            [conda, "run", "-n", env_name, "python", "-c", code],
            capture_output=True, text=True, timeout=timeout_s,
        )
    except subprocess.TimeoutExpired:
        out["status"] = f"timeout after {timeout_s}s"
        return out
    except Exception as e:
        out["status"] = f"subprocess error: {e}"
        return out

    # Stdout can contain lines like 'pipeline_seal: ...' — we look for the
    # last JSON line.
    stdout = (proc.stdout or "").strip()
    if not stdout:
        out["status"] = "no stdout"
        out["stderr_tail"] = (proc.stderr or "")[-500:]
        return out
    payload = None
    for line in reversed(stdout.splitlines()):
        line = line.strip()
        if line.startswith("{") and line.endswith("}"):
            try:
                payload = json.loads(line)
                break
            except json.JSONDecodeError:
                continue
    if payload is None:
        out["status"] = "no JSON in subprocess stdout"
        out["stdout_tail"] = stdout[-500:]
        out["stderr_tail"] = (proc.stderr or "")[-500:]
        return out

    if not payload.get("ok"):
        out["status"] = "snappy subprocess error"
        out["error"] = payload.get("error")
        return out

    out["isosig_used"] = payload.get("isosig_used")
    out["constructed_vol"] = payload.get("constructed_vol")
    out["expected_vol"] = payload.get("expected_vol")
    out["vol_diff"] = payload.get("vol_diff")
    out["vol_match"] = payload.get("vol_match")
    out["constructed_isosig"] = payload.get("constructed_isosig")
    out["isosig_match"] = (
        payload.get("constructed_isosig") == payload.get("isosig_used")
        if payload.get("constructed_isosig") and payload.get("isosig_used") else None
    )
    identify = payload.get("identify") or []
    out["identify_raw"] = identify
    if "identify_error" in payload:
        out["identify_error"] = payload["identify_error"]
    if "ht_volume_scan" in payload:
        out["ht_volume_scan"] = payload["ht_volume_scan"]
    if "ht_volume_scan_error" in payload:
        out["ht_volume_scan_error"] = payload["ht_volume_scan_error"]
    if "ht_volume_scan_skipped" in payload:
        out["ht_volume_scan_skipped"] = payload["ht_volume_scan_skipped"]
    if "ht_dropped_by_c" in payload:
        out["ht_dropped_by_c"] = payload["ht_dropped_by_c"]

    # Identify names look like "12n_706(0,0)" or "K12n706(0,0)" — we
    # extract the knot name + normalise to the classical form (Discovery #175):
    #   c<=10: keep the K prefix (K3_1, K8_19)
    #   c>=11: drop the K prefix + underscore between letter and number (12n_706)
    knot_names = []
    _CLASSICAL_RE = re.compile(r"^K?(\d+)([an])_?(\d+)(m?)$")
    for raw in identify:
        m = re.match(r"^([A-Za-z]?\d+[an]?_?\d+|[KL]\d+[an]\d+|[a-z]\d+_\d+|[mstv]\d+)(?:\(.*\))?$", raw)
        bare = m.group(1) if m else raw.split("(")[0]
        # Normalise to the classical form
        km = _CLASSICAL_RE.match(bare)
        if km:
            c = int(km.group(1))
            if c <= 10:
                bare_classical = bare if bare.startswith("K") else "K" + bare
            else:
                bare_classical = f"{km.group(1)}{km.group(2)}_{km.group(3)}{km.group(4)}"
            knot_names.append(bare_classical)
        else:
            knot_names.append(bare)
    out["identified_names"] = knot_names

    # Confidence evaluation
    if out.get("vol_match") and knot_names:
        if len(knot_names) == 1:
            out["best_match"] = knot_names[0]
            out["confidence"] = "high (vol<{} + snappy.identify)".format(tol)
        else:
            out["best_match"] = knot_names[0]
            out["confidence"] = f"medium (vol-match + {len(knot_names)} identify candidates)"
            out["alt_matches"] = knot_names[1:]
        out["status"] = "matched"
    elif knot_names:
        out["best_match"] = knot_names[0]
        out["confidence"] = "low (identify-only, vol-diff={})".format(out.get("vol_diff"))
        out["status"] = "matched (vol-mismatch)"
    elif out.get("ht_volume_scan"):
        out["status"] = "ht_volume_only"
        out["ht_volume_candidates"] = out["ht_volume_scan"]
    else:
        out["status"] = "no snappy-census match"
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pflicht_json", type=Path)
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--snappy-env", default="snappy_env",
                    help="conda env name for the snappy subprocess (default: snappy_env)")
    ap.add_argument("--skip-snappy", action="store_true",
                    help="skips the snappy census subprocess")
    args = ap.parse_args()

    if not KI_INDEX.exists():
        raise SystemExit(f"KnotInfo index missing: {KI_INDEX}")
    idx = json.loads(KI_INDEX.read_text())
    records = json.loads(KI_RECORDS.read_text())
    fremlin = json.loads(FREMLIN_CACHE.read_text()) if FREMLIN_CACHE.exists() else []

    inv = load_mandatory(args.pflicht_json)
    ki = knotinfo_match(inv, idx, records)
    fr = fremlin_match(inv, ki, fremlin)
    if args.skip_snappy:
        sc = {"status": "skipped (--skip-snappy)"}
    else:
        sc = snappy_census_match(inv, env_name=args.snappy_env)

    # Consolidated status (Discovery #134: snappy census for c>13 as an
    # additional match path)
    sc_best = sc.get("best_match") if sc.get("status", "").startswith("matched") else None

    if ki.get("best_match"):
        mirror_tag = " (mirror)" if "mirror" in ki.get("status", "") else ""
        labels = [f"{ki['best_match']}{mirror_tag} (KnotInfo)"]
        if fr.get("status", "").startswith("matched"):
            labels.append(f"Fremlin {fr['knot_id']}")
        if sc_best:
            labels.append(f"snappy {sc_best}")
        summary = "Identified as " + " = ".join(labels)
        if fr.get("status", "").startswith("matched") or sc_best:
            nov = "known (multi-DB)"
        else:
            nov = "known (KnotInfo)"
    elif ki.get("status", "").startswith("out of KnotInfo"):
        if sc_best:
            extras = []
            if fr.get("status", "").startswith("matched"):
                extras.append(f"Fremlin {fr['knot_id']}")
            extras_s = (" = " + " = ".join(extras)) if extras else ""
            summary = f"Out of KnotInfo range (c={inv['c']}); identified via snappy as {sc_best}{extras_s}"
            nov = "known (snappy-Census)"
        elif fr.get("status", "").startswith("matched"):
            summary = f"Out of KnotInfo range (c={inv['c']}); Fremlin {fr['knot_id']} match"
            nov = "possibly known (Fremlin only)"
        else:
            # Hyperbolicity guard (Werner 2026-06-12): a NON-hyperbolic knot
            # (vol≈0) is a torus / satellite / connected sum — a KNOWN family, NEVER
            # a novelty (doctrine: novelty = hyperbolic prime + no match). The census
            # (snappy/Fremlin) is hyperbolic-focused, which is why classical
            # torus knots (T(3,7), T(6,7)) and connected sums miss it — they are NOT new.
            _v = inv.get("vol")
            if not (isinstance(_v, (int, float)) and _v > 0.5):
                summary = (f"Out of KnotInfo range (c={inv['c']}); non-hyperbolic (vol≈0) — "
                           f"torus / satellite / connected sum: a KNOWN family, NOT a novelty")
                nov = "known (non-hyperbolic family: torus / satellite / connected sum)"
            else:
                summary = f"Out of KnotInfo range (c={inv['c']}); no Fremlin/snappy match — NOVELTY CANDIDATE"
                nov = "novelty candidate (c > 13, no Fremlin, no snappy-Census, hyperbolic)"
    elif ki.get("status", "").startswith("ambiguous"):
        if sc_best:
            summary = f"KnotInfo ambiguous ({len(ki.get('candidates_by_c_det_sig', []))} cands); snappy identifies {sc_best}"
            nov = "known (snappy disambiguation)"
        else:
            summary = f"Ambiguous: {len(ki.get('candidates_by_c_det_sig', []))} KnotInfo candidates"
            nov = "known (ambiguous)"
    else:
        if sc_best:
            extras = []
            if fr.get("status", "").startswith("matched"):
                extras.append(f"Fremlin {fr['knot_id']}")
            extras_s = (" = " + " = ".join(extras)) if extras else ""
            summary = f"snappy identifies {sc_best}{extras_s}"
            nov = "known (snappy-Census)"
        elif fr.get("status", "").startswith("matched"):
            summary = f"No KnotInfo match; Fremlin {fr['knot_id']} match"
            nov = "possibly known (Fremlin only)"
        else:
            summary = "No DB match"
            nov = "novelty candidate"

    result = {
        "name": inv["name"],
        "source_pflicht_json": inv["_source"],
        "input_invariants": {
            "c": inv["c"], "det": inv["det"], "sig": inv["sig"],
            "fib": inv["fib"], "alt": inv["alt"], "amphi": inv["amphi"],
            "vol": inv["vol"], "sym": inv["sym"],
            "isosig": inv.get("isosig"),
        },
        "knotinfo_match": ki,
        "fremlin_match": fr,
        "snappy_census_match": sc,
        "summary": summary,
        "novelty_status": nov,
    }

    out_path = args.output or args.pflicht_json.with_name(
        args.pflicht_json.stem.replace("_pflicht_invarianten_complete", "") + "_crosscheck.json"
    )
    out_path.write_text(json.dumps(result, indent=2, ensure_ascii=False))
    print(f"{inv['name']:<22} | {summary}  →  {out_path.name}")


if __name__ == "__main__":
    main()
