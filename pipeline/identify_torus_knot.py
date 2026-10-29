# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
identify_torus_knot.py — algebraic identification of torus knots T(p,q).

A knot T(p,q) with p,q coprime and p,q>=2 has characteristic values:
  c   = min(p*(q-1), q*(p-1))
  det = 1
  vol = 0 (Seifert fibre manifold)
  fibered = True
  genus = (p-1)*(q-1)/2
  Alex(t) = (t^(p*q) - 1) * (t - 1) / ((t^p - 1) * (t^q - 1))

Mirror T(p,q) has the same values except sig (sign flipped).

Usage:
  conda run -n gnottheorie_env python identify_torus_knot.py <pflicht_complete.json> \
      [--output <torus_check.json>] [--max-pq 30]

Output: torus_check.json with:
  is_torus_knot (bool), pq (tuple if match), mirror (bool),
  algebra_check (compares all fields)
  hypothesis_chain (list of the (p,q) tried)

Discovery #171.
"""
from __future__ import annotations
import argparse
import json
import math
import re
import sys
from pathlib import Path


def _torus_alex_coeffs(p: int, q: int) -> list[int]:
    """Alexander polynomial coefficients of T(p,q) as a list, lowest degree first.

    Alex(t) = (t^(p*q) - 1)(t - 1) / ((t^p - 1)(t^q - 1))

    We compute this via polynomial division over Z (all factors are integer
    polynomials). The resulting polynomial has degree (p-1)*(q-1).
    """
    deg = (p - 1) * (q - 1)
    # Numerator: (t^(p*q) - 1) * (t - 1)
    num = [0] * (p * q + 2)
    num[0] = 1   # constant
    num[1] = -1  # -t
    num[p * q] = -1  # -t^pq
    num[p * q + 1] = 1  # t^(pq+1)
    # Sanity: (t-1)(t^pq - 1) = t^(pq+1) - t^pq - t + 1
    # Check: num = [1, -1, 0, ..., 0, -1, 1]

    # Denominator: (t^p - 1)(t^q - 1) = t^(p+q) - t^p - t^q + 1
    den = [0] * (p + q + 1)
    den[0] = 1
    den[p] = -1
    den[q] = -1
    den[p + q] = 1
    # if p == q: den[p] = -2 would be needed, but we assume gcd(p,q)=1
    # hence p != q (except for (1,1), which we exclude).

    # Polynom-Division (Z[t])
    quot = _poly_divide(num, den)
    if quot is None:
        return []
    # trim to deg
    while len(quot) > deg + 1 and quot[-1] == 0:
        quot.pop()
    return quot


def _poly_divide(num: list[int], den: list[int]) -> list[int] | None:
    """Z[t] polynomial division. Returns quotient or None if not divisible."""
    num = list(num)
    deg_d = len(den) - 1
    while deg_d > 0 and den[deg_d] == 0:
        deg_d -= 1
        den = den[:deg_d + 1]
    if deg_d == 0:
        # constant denominator
        if den[0] == 0:
            return None
        return [v // den[0] for v in num]
    quot = [0] * max(1, len(num) - deg_d)
    while True:
        # remove leading zeros at the end
        while len(num) > 0 and num[-1] == 0:
            num.pop()
        if len(num) - 1 < deg_d:
            break
        deg_n = len(num) - 1
        lead_n = num[-1]
        lead_d = den[-1]
        if lead_n % lead_d != 0:
            return None
        shift = deg_n - deg_d
        coef = lead_n // lead_d
        quot[shift] += coef
        # subtract coef * den * t^shift
        for i, dv in enumerate(den):
            num[shift + i] -= coef * dv
    if any(v != 0 for v in num):
        return None
    return quot


def _parse_alex_to_coeffs(s: str) -> list[int] | None:
    """Parse an Alexander polynomial string like 't^12 - t^11 + ... + 1' or
    'x^12 - x^11 + ...' to a coefficient list, lowest degree first.
    Negative exponents are shifted to min=0."""
    if not s:
        return None
    s = re.sub(r"\s+", "", s)
    # binary '-' after variable/exponent → ' + -' for easier split
    s = re.sub(r"(?<=[^\^\*])-", "+-", s)
    # split
    if s.startswith("+"):
        s = s[1:]
    terms = [t for t in s.split("+") if t]
    coef_by_exp: dict[int, int] = {}
    var_re = re.compile(r"^([+-]?\d*)\*?([a-zA-Z])?(?:\^(-?\d+))?$")
    for term in terms:
        m = var_re.match(term)
        if not m:
            # More complex terms, e.g. with explicit '*'
            # Try generic: (sign)(coef)(*var^exp)?
            sign = 1
            t = term
            if t.startswith("-"): sign = -1; t = t[1:]
            elif t.startswith("+"): t = t[1:]
            mm = re.match(r"^(\d*)\*?([a-zA-Z])\^?(-?\d+)?$", t)
            if not mm:
                if re.match(r"^\d+$", t):
                    coef_by_exp[0] = coef_by_exp.get(0, 0) + sign * int(t)
                    continue
                return None
            coef = int(mm.group(1)) if mm.group(1) else 1
            exp = int(mm.group(3)) if mm.group(3) else 1
            coef_by_exp[exp] = coef_by_exp.get(exp, 0) + sign * coef
            continue
        coef_s, var, exp_s = m.group(1), m.group(2), m.group(3)
        sign = 1
        if coef_s in ("-",):
            sign = -1; coef_s = ""
        elif coef_s in ("+",):
            coef_s = ""
        coef = int(coef_s) if coef_s else 1
        if var is None:
            exp = 0
        else:
            exp = int(exp_s) if exp_s else 1
        coef_by_exp[exp] = coef_by_exp.get(exp, 0) + sign * coef
    if not coef_by_exp:
        return None
    min_e = min(coef_by_exp.keys())
    shifted = {e - min_e: c for e, c in coef_by_exp.items() if c != 0}
    if not shifted:
        return None
    max_e = max(shifted.keys())
    return [shifted.get(i, 0) for i in range(max_e + 1)]


def _norm_alex(coeffs: list[int]) -> tuple[int, ...]:
    """Canonicalisation: mirror = reverse. We take the lexicographically smaller one."""
    a = tuple(coeffs)
    b = tuple(reversed(coeffs))
    # global sign flip as well: -coeffs
    candidates = [a, b, tuple(-x for x in a), tuple(-x for x in b)]
    return min(candidates)


def _coprime(p: int, q: int) -> bool:
    return math.gcd(p, q) == 1


def find_torus_knot(c: int, det: int, vol: float, fibered: bool,
                    alex_coeffs: list[int], max_pq: int = 30) -> dict:
    """Tries all (p,q) with p<q, coprime, c=min(p(q-1),q(p-1))."""
    result = {
        "is_torus_knot": False,
        "pq": None,
        "alex_match": False,
        "preconditions": {"vol_zero": vol < 1e-6, "fibered": bool(fibered),
                          "det_info": det},
        "tried": [],
    }
    # All T(p,q) are non-hyperbolic (Seifert fibre manifold) and fibered.
    # Det varies: det=1 only for p,q both odd; for p=2,q=2k+1: det=2k+1.
    # Hence det is not a pre-filter — the Alexander polynomial comparison discriminates.
    if vol >= 1e-6 or not fibered:
        result["status"] = "preconditions not met (need vol≈0 + fibered=True)"
        return result
    target_norm = _norm_alex(alex_coeffs)
    for p in range(2, max_pq + 1):
        for q in range(p + 1, max_pq + 1):
            if not _coprime(p, q):
                continue
            c_pq = min(p * (q - 1), q * (p - 1))
            if c_pq != c:
                continue
            cand_coeffs = _torus_alex_coeffs(p, q)
            cand_norm = _norm_alex(cand_coeffs)
            match = cand_norm == target_norm
            result["tried"].append({"p": p, "q": q, "c": c_pq,
                                    "alex_match": match,
                                    "expected_genus": (p - 1) * (q - 1) // 2})
            if match:
                result["is_torus_knot"] = True
                result["pq"] = [p, q]
                result["alex_match"] = True
                result["status"] = f"matched T({p},{q})"
                return result
    if result["tried"]:
        result["status"] = f"c={c} matched {len(result['tried'])} (p,q) candidates, but no Alex match"
    else:
        result["status"] = f"no (p,q) with c={c} in the search range p,q<={max_pq}"
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pflicht_json", type=Path)
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--max-pq", type=int, default=40)
    args = ap.parse_args()

    d = json.loads(args.pflicht_json.read_text())
    pl = d["pflichtliste"]
    A, H = pl.get("ALGEBRA", {}), pl.get("HYPERBOLISCH", {})
    D = pl.get("DIAGRAMM", {})
    c = D.get("crossing_number_min")
    det = A.get("determinant")
    vol = H.get("hyperbolic_volume")
    fibered = A.get("fibered")
    alex_str = str(A.get("alexander_polynomial", ""))
    alex_coeffs = _parse_alex_to_coeffs(alex_str)

    out = {
        "name": d.get("name"),
        "source_pflicht_json": str(args.pflicht_json),
        "inputs": {
            "c": c, "det": det, "vol": vol, "fibered": fibered,
            "alex_str": alex_str, "alex_coeffs": alex_coeffs,
        },
    }
    if alex_coeffs is None or c is None or det is None:
        out["status"] = "skipped (missing alex_coeffs / c / det)"
    else:
        res = find_torus_knot(c, det, vol or 0, bool(fibered), alex_coeffs, args.max_pq)
        out.update(res)

    out_path = args.output or args.pflicht_json.with_name(
        args.pflicht_json.stem.replace("_pflicht_invarianten_complete", "") + "_torus_check.json"
    )
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False))
    name = out.get("name") or args.pflicht_json.stem
    if out.get("is_torus_knot"):
        print(f"{name:<22} | T{tuple(out['pq'])} (matched alex)  →  {out_path.name}")
    else:
        print(f"{name:<22} | {out.get('status','?')}  →  {out_path.name}")


if __name__ == "__main__":
    main()
