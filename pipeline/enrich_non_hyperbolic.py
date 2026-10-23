#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
enrich_non_hyperbolic.py
========================
Standalone enrichment: for non-hyperbolic knots (vol ≈ 0)
test two algorithmic hypotheses and write the result as a new
block 'non_hyperbolic_classification' into the pflicht_*_complete.json.

Background (Discoveries #142 + #143, 2026-05-12):
  - knots-of-love (c=18, det=729=3^6, vol≈0) was identified as a 6-fold connected
    sum of trefoils via Alex = (t^2 - t + 1)^6, bit-exact
    against the pipeline Alex.
  - irrational (c=14, det=1, vol=0, 3-strand braid) was identified as
    mirror(T(3,7)) via the standard torus Alex formula.
  - Both tests are algebraically very cheap + very informative.
    They get past the snappy limitation (vol≈0 → snappy gives no
    further information) without using table lookups.

Doctrine preserved:
  - The pipeline seal (#64) stays untouched — all tests are
    purely algorithmic, via Alex polynomial comparison against standard
    formulas or hard-coded Alex polynomials of small prime knots.
  - Writes a NEW block 'non_hyperbolic_classification' into
    pflicht_complete.json, leaves all mandatory fields
    untouched.

Logic:
  Step 0: vol check. vol > vol_threshold → hyperbolic → prime
          (Thurston trichotomy, Discovery #128) → skip with a marker.
  Step 1 (torus): iterate over coprime (p,q) with (p-1)*(q-1) =
                  deg(alex). Compare bit-exactly with the standard T(p,q)
                  Alex polynomial. Mirror via target == T(p,q)(1/t)*t^deg.
  Step 2 (connected sum): check whether target == (alex(K))^k for a
                          known small prime knot K and k >= 2.
                          Table: 3_1, 4_1, 5_1, 5_2, 6_1, 6_2, 6_3,
                          7_1, 8_19 (T(3,4)).
                          Plus a full sympy factorisation as an
                          indicator.

Usage:
  python enrich_non_hyperbolic.py <pflicht_*.json> [--in-place]
"""
import argparse
import json
import re
import sys
from math import gcd
from pathlib import Path


# ---------------------------------------------------------------------------
# Alex polynomials of known prime knots (small c, prime, blind-test compliant).
#
# Doctrine audit (Discovery #144 / #64): these polynomials are STANDARD TEXTBOOK
# CONSTANTS of knot theory (Rolfsen table / Bar-Natan). No
# KnotInfo database pull. The values can be derived from classical sources (Murasugi,
# Adams, Lickorish) and verified against standard formulas (T(p,q)).
# The pipeline seal is preserved — we compare the
# computed mandatory Alex against these CONSTANTS, not against a
# dynamic DB query.
# ---------------------------------------------------------------------------
def _known_prim_alex(t):
    """Returns dict {name: sympy polynomial in t}."""
    return {
        "3_1":  t**2 - t + 1,
        "4_1":  t**2 - 3*t + 1,
        "5_1":  t**4 - t**3 + t**2 - t + 1,        # = T(2,5)
        "5_2":  2*t**2 - 3*t + 2,
        "6_1":  2*t**2 - 5*t + 2,
        "6_2":  -t**4 + 3*t**3 - 3*t**2 + 3*t - 1,
        "6_3":  t**4 - 3*t**3 + 5*t**2 - 3*t + 1,
        "7_1":  t**6 - t**5 + t**4 - t**3 + t**2 - t + 1,  # = T(2,7)
        "7_2":  3*t**2 - 5*t + 3,
        "8_19": t**6 - t**5 + t**3 - t + 1,        # = T(3,4)
    }


# ---------------------------------------------------------------------------
# Alex string → sympy polynomial in t
# ---------------------------------------------------------------------------
def parse_alex_string(s, t):
    """Parses the Alex string from the pipeline pflicht (x as variable, often without *)."""
    if not s:
        return None
    s = str(s).strip()
    # 'x' → 't', replace whitespace between number and variable by '*'
    s = s.replace("x", "t")
    # Pipeline format: "1 t^2", "2 t^11" — insert * where needed
    s = re.sub(r"(\d)\s+([a-zA-Z])", r"\1*\2", s)
    s = re.sub(r"(\d)\s*\^", r"\1**", s)
    s = s.replace("^", "**")
    try:
        import sympy as sp
        return sp.expand(sp.sympify(s, locals={"t": t}))
    except Exception as e:
        return None


# ---------------------------------------------------------------------------
# T(p,q) torus-knot Alex polynomial
# ---------------------------------------------------------------------------
def torus_alex(p, q, t):
    """Δ(T(p,q)) = (t^pq - 1)(t - 1) / ((t^p - 1)(t^q - 1)).

    Returns the polynomial as a sympy expression. None if the division
    is not exact.
    """
    import sympy as sp
    if gcd(p, q) != 1:
        return None
    num = (t**(p*q) - 1) * (t - 1)
    den = (t**p - 1) * (t**q - 1)
    q_poly, r_poly = sp.div(num, den, t)
    if r_poly != 0:
        return None
    return sp.expand(q_poly)


def torus_min_crossings(p, q):
    """Minimal crossing number of T(p,q) (p, q coprime)."""
    return min(p * (q - 1), q * (p - 1))


def torus_genus(p, q):
    return (p - 1) * (q - 1) // 2


# ---------------------------------------------------------------------------
# Torus test: iterate over all coprime (p,q) with (p-1)(q-1) = deg(alex)
# ---------------------------------------------------------------------------
def try_torus_match(target_alex, t, alex_deg, c_pipeline, s_pipeline):
    """Check all coprime (p,q) with (p-1)(q-1)=alex_deg for a bit match."""
    import sympy as sp
    if not isinstance(alex_deg, int) or alex_deg < 2:
        return None
    target_mirror = sp.expand(target_alex.subs(t, 1/t) * t**alex_deg)

    candidates = []
    for p in range(2, alex_deg + 2):
        q_target = (alex_deg // (p - 1)) + 1
        if (p - 1) * (q_target - 1) != alex_deg:
            continue
        if gcd(p, q_target) != 1:
            continue
        if p > q_target:
            continue  # only p ≤ q
        a = torus_alex(p, q_target, t)
        if a is None:
            continue
        direct_match = (a == target_alex)
        mirror_match = (a == target_mirror)
        if direct_match or mirror_match:
            # Mirror detection via Alex alone is unreliable (the torus Alex
            # is palindromic). Sign discrimination via the s_invariant:
            # standard T(p,q) has s = +(p-1)(q-1). s_pipeline < 0 ⇒ mirror.
            s_torus = (p - 1) * (q_target - 1)
            mirror_via_alex_only = mirror_match and not direct_match
            mirror_via_signature = (s_pipeline is not None and
                                    isinstance(s_pipeline, (int, float)) and
                                    s_pipeline < 0)
            candidates.append({
                "p":       p,
                "q":       q_target,
                "alex_direct_match": direct_match,
                "alex_mirror_match": mirror_match,
                "mirror_via_signature": mirror_via_signature,
                "is_mirror_variant":   mirror_via_alex_only or mirror_via_signature,
                "c_torus": torus_min_crossings(p, q_target),
                "genus_torus": torus_genus(p, q_target),
                "s_torus": s_torus,
                "s_pipeline": s_pipeline,
            })
    return candidates


# ---------------------------------------------------------------------------
# Connected-sum test: k-th powers of known prime knots
# ---------------------------------------------------------------------------
def try_connected_sum_power(target_alex, t, det_pipeline):
    """Check whether target = (Alex(K))^k for a small prime knot K."""
    import sympy as sp
    candidates = []
    for name, alex_K in _known_prim_alex(t).items():
        alex_K_expand = sp.expand(alex_K)
        det_K = abs(int(alex_K_expand.subs(t, -1)))
        if det_K == 0:
            continue
        # Which k fit the det balance?
        if det_pipeline is None or det_pipeline <= 1:
            ks = range(2, 8)
        else:
            ks = []
            for k in range(2, 8):
                if det_K**k == det_pipeline:
                    ks.append(k)
            if not ks:
                continue
        for k in ks:
            candidate = sp.expand(alex_K_expand**k)
            if candidate == target_alex:
                candidates.append({
                    "factor_name":  name,
                    "factor_alex":  str(alex_K_expand),
                    "k":            k,
                    "det_factor":   det_K,
                    "det_check":    det_K**k,
                })
    return candidates


def try_full_factorization(target_alex, t):
    """Sympy.factor + heuristic: each factor against known prime knots.

    If ALL irreducible factors match known prime-knot Alex polynomials
    (heterogeneous connected-sum decomposition such as 3_1 # 4_1 # ...), this gives
    a genuine connected-sum identification, not only a partial one.
    """
    import sympy as sp
    factored = sp.factor(target_alex)
    factors_str = str(factored)
    factor_list = sp.factor_list(target_alex)
    # factor_list = (coeff, [(poly, multiplicity), ...])
    matched = []
    unmatched = []
    for poly, mult in factor_list[1]:
        poly_exp = sp.expand(poly)
        matched_name = None
        for name, alex_K in _known_prim_alex(t).items():
            alex_K_exp = sp.expand(alex_K)
            if poly_exp == alex_K_exp:
                matched_name = name
                break
            # Mirror (t -> 1/t * t^deg)
            deg = sp.Poly(alex_K_exp, t).degree()
            mirror = sp.expand(alex_K_exp.subs(t, 1/t) * t**deg)
            if poly_exp == mirror:
                matched_name = name + " (mirror)"
                break
        if matched_name:
            matched.append({"factor": str(poly_exp), "mult": mult,
                            "matched_as": matched_name})
        else:
            unmatched.append({"factor": str(poly_exp), "mult": mult})

    # Fully matched (heterogeneous connected sum)
    is_full_cs = (not unmatched) and (len(matched) >= 1)
    cs_decomp = None
    if is_full_cs:
        parts = []
        for m in matched:
            if m["mult"] == 1:
                parts.append(m["matched_as"])
            else:
                parts.append(f"{m['mult']}×{m['matched_as']}")
        cs_decomp = " # ".join(parts)

    return {
        "factored_str":         factors_str,
        "matched":              matched,
        "unmatched":            unmatched,
        "leading_coeff":        str(factor_list[0]),
        "is_full_connected_sum": is_full_cs,
        "connected_sum_decomposition": cs_decomp,
    }


# ---------------------------------------------------------------------------
# Disc #323 — connected-sum gate before "vol>threshold ⇒ prime"
# ---------------------------------------------------------------------------
def _alex_reducible_into_knot_factors(target, t):
    """True if the Alexander polynomial splits into >=2 non-trivial factors
    (deg>=1, not a bare x shift). NECESSARY condition for a connected sum
    (Δ_{K1#K2}=Δ_{K1}·Δ_{K2}). Cheap pre-filter: only if it fires is the
    expensive geometric re-check worthwhile."""
    if target is None:
        return False
    try:
        import sympy as sp
        const, facs = sp.factor_list(sp.expand(target))
    except Exception:
        return False
    nontriv = 0
    for base, mult in facs:
        bp = sp.Poly(base, t)
        if bp.degree() == 0:
            continue
        if bp == sp.Poly(t, t):   # bare t factor (Laurent normalisation)
            continue
        nontriv += mult
    return nontriv >= 2


def _robust_hyperbolic(pd_code, n_tries=40):
    """Build the knot from pd_code -> snappy.exterior -> randomize loop.
    Returns True if 'all tetrahedra positively oriented' is ever reached
    (= GENUINE hyperbolic structure, knot prime), otherwise False.
    Returns None if snappy/spherogram/pd is not available (verification impossible).
    Disc #323: a connected sum is NEVER hyperbolic — snappy does report a
    volume (sum of the JSJ pieces), but finds no positively oriented solution."""
    if not pd_code or not isinstance(pd_code, str):
        return None
    try:
        import re as _re
        import snappy
        import spherogram
    except Exception:
        return None
    try:
        xs = _re.findall(r'X\[([0-9,\s]+)\]', pd_code)
        pdl = [tuple(int(tt) for tt in x.split(',')) for x in xs]
        if not pdl:
            return None
        M = spherogram.Link(pdl).exterior()
        for _ in range(n_tries):
            try:
                M.randomize()
            except Exception:
                pass
            if M.solution_type() == 'all tetrahedra positively oriented':
                return True
        return False
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pflicht_json", type=Path)
    ap.add_argument("--in-place", action="store_true",
                    help="Patch directly into the input file.")
    ap.add_argument("--vol-threshold", type=float, default=1e-3,
                    help="vol < threshold → non-hyperbolic (default 1e-3).")
    args = ap.parse_args()

    if not args.pflicht_json.exists():
        sys.exit(f"ERROR: {args.pflicht_json} missing")

    d = json.loads(args.pflicht_json.read_text())
    pl = d.get("pflichtliste", {})
    alg = pl.get("ALGEBRA", {})
    diag = pl.get("DIAGRAMM", {})
    hyp = pl.get("HYPERBOLISCH", {})

    alex_str = alg.get("alexander_polynomial")
    vol      = hyp.get("hyperbolic_volume")
    det      = alg.get("determinant")
    c        = diag.get("crossing_number_min")
    s_inv    = alg.get("s_invariant")
    genus    = alg.get("genus")
    fibered  = alg.get("fibered")
    name     = d.get("name", args.pflicht_json.stem)

    print(f"[enrich_non_hyp] {name}: vol={vol}, det={det}, c={c}")

    result = {
        "vol_threshold":  args.vol_threshold,
        "vol":            vol,
        "is_hyperbolic":  None,
        "torus_match":    None,
        "connected_sum_match": None,
        "alex_factorization": None,
        "decision":       None,
    }

    # Step 0: hyperbolicity  (Disc #323 connected-sum gate)
    run_nonhyp_tests = False
    if vol is None:
        # No volume does NOT mean "not hyperbolic" — it means UNMEASURED.
        # Until 2026-08-23 the run ended here with decision="no_vol_field" and
        # run_nonhyp_tests=False: NOTHING AT ALL was checked. For lack of an alternative
        # the labeler turned this into "known_nonhyp" — an unmeasured knot counted
        # as a known one. Found on pipeline_failures/gitter_d3_k3_s1p5_g0p25/myknot
        # (c=48, det=19^8), whose stage 6 died on 'contains degenerate tetrahedra';
        # it really is an 8-fold connected sum, but nobody knew.
        # The algebraic tests (torus / connected sum) need NO volume, only
        # the Alexander polynomial — and _robust_hyperbolic decides independently of it.
        genuine = _robust_hyperbolic(diag.get("pd_code"))
        result["robust_hyperbolicity"] = genuine
        if genuine is True:
            result["is_hyperbolic"] = True
            result["decision"] = "hyperbolic_verified_without_vol_field"
            print("[enrich_non_hyp]   vol missing, BUT snappy finds a positively oriented structure → hyperbolic (recompute vol).")
        elif genuine is False:
            result["is_hyperbolic"] = False
            run_nonhyp_tests = True
            result["decision"] = "no_vol_field_snappy_nonhyperbolic"
            print("[enrich_non_hyp]   vol missing, snappy finds NO positively oriented structure → non-hyperbolic. Torus + connected-sum tests follow.")
        else:
            result["decision"] = "no_vol_field_unverifiable"
            print("[enrich_non_hyp]   vol missing and snappy check not possible → OPEN (no finding, no classification).")
    elif vol > args.vol_threshold:
        # Disc #323: vol>threshold does NOT prove hyperbolicity — snappy sums the
        # JSJ pieces of a connected sum into one volume (that is how 12 su_* K#K slipped
        # through wrongly as novel_hyperbolic). Suspicious only if the Alexander
        # polynomial splits into >=2 knot factors (necessary for #). Then re-check the
        # hyperbolic structure robustly (snappy + randomize), so that an unfavourably
        # stored triangulation does not wrongly demote a GENUINE hyperbolic knot.
        try:
            import sympy as sp
        except ImportError:
            sys.exit("ERROR: sympy not available — enrich_non_hyperbolic needs sympy")
        t = sp.symbols("t")
        target = parse_alex_string(alex_str, t)
        if target is not None and _alex_reducible_into_knot_factors(target, t):
            genuine = _robust_hyperbolic(diag.get("pd_code"))
            result["robust_hyperbolicity"] = genuine
            if genuine is False:
                result["is_hyperbolic"] = False
                run_nonhyp_tests = True
                print(f"[enrich_non_hyp]   vol={vol:.4f}>{args.vol_threshold}, BUT reducible Alex + NO positively oriented snappy solution (Disc #323) → non-hyperbolic. Torus+CS tests follow.")
            elif genuine is True:
                result["is_hyperbolic"] = True
                result["decision"] = "hyperbolic_prim_verified_reducible_alex"
                print(f"[enrich_non_hyp]   vol={vol:.4f}>{args.vol_threshold}, reducible Alex, but snappy confirms a positively oriented structure → prime.")
            else:  # None: not verifiable → conservatively the old assumption (Disc #128)
                result["is_hyperbolic"] = True
                result["decision"] = "hyperbolic_prim_unverified_reducible_alex"
                print(f"[enrich_non_hyp]   vol={vol:.4f}>{args.vol_threshold}, reducible Alex, snappy check not possible → conservatively prime (Disc #128).")
        else:
            result["is_hyperbolic"] = True
            result["decision"] = "hyperbolic_prim_irreducible_alex"
            print(f"[enrich_non_hyp]   vol={vol:.4f} > {args.vol_threshold}, Alex irreducible → hyperbolic ⇒ prime (Discovery #128).")
    else:
        result["is_hyperbolic"] = False
        run_nonhyp_tests = True
        print(f"[enrich_non_hyp]   vol={vol} ≤ {args.vol_threshold} → non-hyperbolic. Torus + connected-sum tests follow.")

    if run_nonhyp_tests:
        try:
            import sympy as sp
        except ImportError:
            sys.exit("ERROR: sympy not available — enrich_non_hyperbolic needs sympy")

        t = sp.symbols("t")
        target = parse_alex_string(alex_str, t)
        if target is None:
            result["decision"] = "alex_parse_failed"
            print(f"[enrich_non_hyp]   ERROR: Alex parse failed: {alex_str[:80]}")
        else:
            alex_deg = sp.Poly(target, t).degree()
            print(f"[enrich_non_hyp]   Alex parsed, deg = {alex_deg}")

            # 1) Torus test
            torus_cands = try_torus_match(target, t, alex_deg, c, s_inv)
            if torus_cands:
                result["torus_match"] = torus_cands
                result["decision"] = "torus_knot"
                for tc in torus_cands:
                    label = f"T({tc['p']},{tc['q']})"
                    if tc["is_mirror_variant"]:
                        label = f"mirror({label}) = T({tc['p']},-{tc['q']})"
                    print(f"[enrich_non_hyp]   ✓ TORUS: {label}  (Alex bit-exact)")
                    print(f"[enrich_non_hyp]     c_torus={tc['c_torus']}, genus={tc['genus_torus']}, |s|={tc['s_torus']}")
                    print(f"[enrich_non_hyp]     vs. pipeline: c={c}, genus={genus}, s={s_inv}")
                    if tc["is_mirror_variant"]:
                        print(f"[enrich_non_hyp]     Mirror indicator: s_pipeline={s_inv} < 0 (standard T convention s>0)")
            else:
                # 2) Connected-sum test
                cs_cands = try_connected_sum_power(target, t, det)
                if cs_cands:
                    result["connected_sum_match"] = cs_cands
                    result["decision"] = "connected_sum_power"
                    for cs in cs_cands:
                        print(f"[enrich_non_hyp]   ✓ CONNECTED SUM: {cs['k']} × {cs['factor_name']}  (Alex bit-exact)")
                        print(f"[enrich_non_hyp]     det = {cs['det_factor']}^{cs['k']} = {cs['det_check']}")
                else:
                    # 3) Full factorisation as an indicator
                    fact = try_full_factorization(target, t)
                    result["alex_factorization"] = fact
                    if fact.get("is_full_connected_sum"):
                        # Heterogeneous connected sum — all factors matched
                        result["decision"] = "connected_sum_heterogeneous"
                        result["connected_sum_decomposition"] = fact["connected_sum_decomposition"]
                        print(f"[enrich_non_hyp]   ✓ CONNECTED SUM (heterogeneous): {fact['connected_sum_decomposition']}")
                        print(f"[enrich_non_hyp]     Factor match: {[m['matched_as'] for m in fact['matched']]}")
                    elif fact["matched"]:
                        result["decision"] = "partial_factorization"
                        print(f"[enrich_non_hyp]   ! Partial factor match: {[m['matched_as'] for m in fact['matched']]}")
                        if fact["unmatched"]:
                            print(f"[enrich_non_hyp]     unmatched factors: {[m['factor'][:40] for m in fact['unmatched']]}")
                    else:
                        # Disc #327 satellite gate: non-hyperbolic + not torus + no
                        # connected sum. Trichotomy (prime knot = hyperbolic XOR torus XOR satellite):
                        # a non-hyperbolic prime knot with IRREDUCIBLE Alexander is a
                        # SATELLITE (known family, NEVER novelty). Reducible -> composite (factors
                        # merely not named). Both are known_nonhyp anyway via is_hyperbolic=False.
                        if _alex_reducible_into_knot_factors(target, t):
                            result["decision"] = "composite_unnamed_factors"
                            result["satellite_or_composite"] = "composite"
                        else:
                            result["decision"] = "satellite_candidate"
                            result["satellite_or_composite"] = "satellite"
                        print(f"[enrich_non_hyp]   - non-hyp, not torus, no CS -> {result['decision']} "
                              f"(known family, NO novelty; trichotomy). Factorisation: {fact['factored_str'][:80]}")

    # In-place patch
    d["non_hyperbolic_classification"] = result
    if args.in_place:
        args.pflicht_json.write_text(json.dumps(d, indent=2, ensure_ascii=False, default=str))
        print(f"[enrich_non_hyp]   in-place patch -> {args.pflicht_json}")
    else:
        print(json.dumps(result, indent=2, default=str))


if __name__ == "__main__":
    main()
