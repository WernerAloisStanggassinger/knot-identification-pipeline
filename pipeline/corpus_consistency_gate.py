#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Knot-Identification Pipeline — Werner Alois Stanggassinger & Claude Code (Anthropic).
"""corpus_consistency_gate.py — theorem-based data gate (machine 183, 2026-06-12).

Checks every corpus knot against UNCONDITIONAL knot theorems. Every violation = a data
error (wrong sig/genus/det/Alexander). Found on 2026-06-12: S_24_d25_Z3_chiral 4_genus
underestimated (g4=1, forced g4=2; Disc/error log).

DELIBERATELY ONLY UNCONDITIONAL THEOREMS (no „amphi⟹det=□" — that holds only for NEGATIVE-amphi,
Hartley-Kawauchi; the strong form produces false positives, see machine 178/179).

Usage:
    python3 corpus_consistency_gate.py [--db PATH] [--strict]
    # library: from corpus_consistency_gate import check_corpus -> (n_viol, violations)
Exit code 1 on violation (for pre-commit gate), 0 otherwise.
"""
from __future__ import annotations
import json
import site_config
import sys
import argparse

try:
    import sympy as sp
    from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,
                                            implicit_multiplication_application)
    _HAVE_SYMPY = True
    _T = sp.symbols("t")
    _TR = standard_transformations + (implicit_multiplication_application,)
except Exception:
    _HAVE_SYMPY = False

DEFAULT_DB = "@PIPELINE_DIR@/corpus_db.json".replace('@PIPELINE_DIR@', site_config.PIPELINE_DIR)


def _num(x):
    return x if isinstance(x, (int, float)) else None


def _alex(r):
    if not _HAVE_SYMPY:
        return None
    raw = str(r.get("alexander") or r.get("alex") or "")
    if not raw or raw == "None":
        return None
    raw = raw.replace("^", "**").replace("x", "t")
    try:
        return parse_expr(raw, transformations=_TR, local_dict={"t": _T})
    except Exception:
        return None


def check_record(r):
    """Returns list of (check, detail) violations for ONE record."""
    out = []
    kid = r.get("distinct_knot_id") or r.get("stem")
    det, sig = _num(r.get("det")), _num(r.get("sig"))
    g3, g4 = _num(r.get("3_genus")), _num(r.get("4_genus"))

    if det is not None and int(abs(det)) % 2 == 0:
        out.append(("det_ungerade", f"det={det} gerade"))
    if sig is not None and int(sig) % 2 != 0:
        out.append(("sig_gerade", f"sig={sig} ungerade"))
    if sig is not None and g4 is not None and abs(sig) / 2 > g4 + 1e-9:
        # additionally: if g3 is known, g4 is uniquely forced
        fix = f" (g4 erzwungen = {int(abs(sig)//2)}..{int(g3)})" if g3 is not None else ""
        out.append(("slice_bound", f"|sig|/2={abs(sig)/2} > g4={g4}{fix}"))
    if g3 is not None and g4 is not None and g4 > g3 + 1e-9:
        out.append(("g4_le_g3", f"g4={g4} > g3={g3}"))
    if sig is not None and g3 is not None and abs(sig) > 2 * g3 + 1e-9:
        out.append(("sig_genus_bound", f"|sig|={abs(sig)} > 2*g3={2*g3}"))
    # Fox-Milnor: slice (g4=0) ⟹ det = Quadrat ∧ sig = 0 (notwendig). Verstoss ⟹ g4=0 falsch.
    if g4 is not None and g4 == 0:
        from math import isqrt
        d = abs(int(det)) if det is not None else None
        if d is not None and isqrt(d) ** 2 != d:
            out.append(("foxmilnor_slice_det", f"g4=0 (slice) but det={d} not a square → g4≠0"))
        if sig is not None and sig != 0:
            out.append(("foxmilnor_slice_sig", f"g4=0 (slice) aber sig={sig}≠0 → g4≠0"))

    # Arf-Invariante ↔ Determinante mod 8: Arf=0 ⟺ det≡±1 mod 8, Arf=1 ⟺ det≡±3 mod 8.
    arf = r.get("arf")
    if arf in (0, 1) and det is not None:
        dm = abs(int(det)) % 8
        exp = 0 if dm in (1, 7) else (1 if dm in (3, 5) else None)
        if exp is not None and arf != exp:
            out.append(("arf_det_mod8", f"Arf={arf} aber det≡{dm} mod 8 (erwartet Arf={exp})"))

    a = _alex(r)
    if a is not None:
        d1 = a.subs(_T, 1)
        try:
            if abs(int(d1)) != 1:
                out.append(("alex_at_1", f"Δ(1)={d1} ≠ ±1"))
        except Exception:
            pass
        if det is not None:
            try:
                dm1 = abs(int(a.subs(_T, -1)))
                if dm1 != abs(int(det)):
                    out.append(("det_eq_alex_m1", f"|Δ(−1)|={dm1} ≠ det={abs(int(det))}"))
            except Exception:
                pass
    return [(kid, c, d) for c, d in out]


def check_corpus(db_path=DEFAULT_DB, only_distinct=True):
    db = json.load(open(db_path))
    knots = db["knots"]
    if only_distinct:
        knots = [k for k in knots if k.get("is_distinct_representative")]
    violations = []
    for r in knots:
        violations.extend(check_record(r))
    return len(violations), violations


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default=DEFAULT_DB)
    ap.add_argument("--all", action="store_true", help="all records, not only distinct ones")
    args = ap.parse_args()
    n, viol = check_corpus(args.db, only_distinct=not args.all)
    if not viol:
        print("✓ Consistency gate GREEN — no theorem violation.")
        return 0
    print(f"⚠ {n} VIOLATION(S) (= data errors):")
    for kid, c, d in viol:
        print(f"  {c:16} {kid}: {d}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
