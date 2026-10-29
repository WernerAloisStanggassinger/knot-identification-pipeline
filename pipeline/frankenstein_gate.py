#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""frankenstein_gate.py — same-knot guard as a pre-commit gate (Disc #251/#252/#253).

Occasion (2026-06-04, Werner's go "mach das gate" (EN: "make the gate")): the Phase-C complex showed
that a record can be stitched together from a diagram layer (dt/braid/PD of knot A) and an algebra layer
(det/genus/alex of object B) ("Frankenstein").
The full-corpus audit (Disc #253) measured all existing records — this
module prevents mechanically that FUTURE Frankensteins ever get into corpus_db
or through the referee gate again.

Three layers:
  TIER 0 (free, pure stdlib — runs ALWAYS, in build_corpus_db as a HARD gate
          and in crosscheck_consistency as a standard check). Bit-tested
          impossibility flags from r3_phantom_diag (caught 7/7 phantoms):
            - cmin > dt length                  (minimal c above an existing diagram)
            - genus > (dt-1)/2                  (Seifert bound)
            - KFH total_rank < det              (rank >= det always)
            - KFH seifert_genus != ALGEBRA.genus
          plus two new free checks (algebra-internal consistency):
            - |Delta(-1)| != det                (Alexander vs determinant)
            - deg(Delta) > 2*genus              (Seifert bound)
            - determinant != determinant_sage   (two-witness divergence)
            - cmin > raw PD crossing number     (pd_list_str)
  TIER 1 (--pd, needs regina ~ miniforge base): det(stored PD after simplify)
          == ALGEBRA.det — the detector that caught kf86/kf373 (Disc #252).
  TIER 2 (--guard-new, GNOTTHEORIE_ENV: snappy+spherogram+numpy): full same-knot
          guard, REWRITE 2026-07-02 — dense reconstruction @N=2048 (vol converges
          only there, resolution lesson Disc #253) -> tone (snappy vol) + point eval
          (det/sig via spherogram seifert_matrix @ t=-1, sage-free, NO regina).
          vol match certifies the projection; Frankenstein = vol matches, algebra
          does not. Result cached in the registry (one run per fseries hash+det).

Registry: frankenstein_guard_registry.json (next to corpus_db.json).
  stem -> {fseries_sha256, det, vol, deg_alex, verdict, source, checked_at}
  Scope = machine/own records outside bekannt/singletons_unknots (audit scope
  Disc #253; bekannt = identified via table cross-check, unknots =
  documented degenerations). crosscheck_consistency FAILs if a
  scope record has no fresh registry entry (fseries hash + det match)
  -> new wave records MUST go through --guard-new before the gate turns green.

CLI:
  python3 frankenstein_gate.py --scan              # Tier 0 over the whole corpus (read-only)
  ~/miniforge3/bin/python frankenstein_gate.py --scan --pd    # + Tier 1
  ~/miniforge3/envs/gnottheorie_env/bin/python frankenstein_gate.py --guard-new [--jobs N]
                                                   # Tier 2 (tone+point eval): new/changed
                                                   # records + registry write. GNOTTHEORIE_ENV!
  python3 frankenstein_gate.py --seed-registry     # once: seal the existing stock
                                                   # (provenance = full audit Disc #253)
  python3 frankenstein_gate.py --registry-status   # coverage report
"""
from __future__ import annotations
import argparse
import site_config
import hashlib
import json
import re
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(site_config.ROOT)
HERE = Path(__file__).resolve().parent
CORPUS_DB = HERE / "corpus_db.json"
REGISTRY = HERE / "frankenstein_guard_registry.json"
TOOLS = ROOT / "claude_setup" / "tools"

# Same-knot guard (rewrite 2026-07-02, successor of Disc #391): tone (snappy vol) +
# point eval (det via spherogram seifert_matrix @ t=-1, sage-free) instead of PD-ladder
# census. Resolution lesson (Nyquist/Disc #253): vol converges to the true value only at N=2048
# (N<=1536 = wrong/simpler knots) -> fixed dense reconstruction.
# Seed 0 primary; vol match certifies that the projection hit the right knot;
# otherwise fallback seeds. Runs in GNOTTHEORIE_ENV (snappy+spherogram+numpy),
# NOT base — det needs NO regina/symbolic Alexander any more (that was the 142 s killer).
GUARD_N = 2048
GUARD_SEEDS = [0, 9001, 9002, 9003]
VOL_TOL = 1e-4
# 🔴 2026-08-24: the earlier one-liner "DET_PRIME = 1000003  # >> every corpus det"
# was an ASSUMPTION with an expiry date — and it has expired. gitter_d3_k3_s1p5_g0p25
# has det = 19^8 = 16,983,563,041, SEVENTEEN-THOUSAND TIMES the modulus. The guard
# then reported "det=487911 (should be 16983563041)" and would have condemned the record as
# Frankenstein — yet 19^8 mod 1000003 = 512092, signed remainder
# -487911: the reported number WAS the right determinant, only folded.
# Cure: NOT raising the modulus (the next determinant is larger again),
# but comparing modularly on BOTH sides — and over several primes,
# so that a random collision (1e-6 with one) is practically excluded (1e-18).
DET_PRIME = 1000003         # backward compatibility; comparison runs via:
DET_PRIMES = (1000003, 999983, 1000033)

# Registry scope: bekannt = table-identified (own evidence chain),
# singletons_unknots = documented sweep degenerations (Werner 2026-06-04).
SCOPE_EXCLUDED_CLUSTERS = {"bekannt", "singletons_unknots", "pipeline_failures",
                           "singletons"}


# ───────────────────────── Helpers (pure stdlib) ─────────────────────────────

def _int(v):
    """Pipeline fields come as int OR str depending on the age of the stage."""
    if v is None or isinstance(v, bool):
        return None
    try:
        return int(str(v).strip())
    except (ValueError, TypeError):
        return None


def _seq_len(v):
    """Length of dt_code / braid_word — str ('30 36 -34 …') or list."""
    if v is None:
        return None
    if isinstance(v, (list, tuple)):
        return len(v) or None
    n = len(re.findall(r"-?\d+", str(v)))
    return n or None


def alex_eval_minus1(s):
    """|Delta(-1)| from the stored Alexander string — regex parser,
    NO sympify (error-log doctrine 2026-06-04: implicit multiplication).
    Format: 'x^8 - 16 x^7 + 90 x^6 - … + 1' (also t variable, also '1')."""
    if s is None:
        return None
    s = str(s).replace("**", "^").replace("*", "").replace(" ", "")
    s = s.replace("−", "-")  # Unicode minus
    var = "x" if "x" in s else ("t" if "t" in s else None)
    if var is None:
        return abs(_int(s)) if _int(s) is not None else None
    if var == "t":
        s = s.replace("t", "x")
    total = 0
    pos = 0
    for m in re.finditer(r"([+-]?)(\d*)(x(?:\^(-?\d+))?)?", s):
        if m.start() != pos or not m.group(0):
            if m.start() >= len(s):
                break
            return None  # unparseable piece -> check abstains
        pos = m.end()
        sign, coeff, xpart, exp = m.groups()
        if not coeff and not xpart:
            continue
        c = int(coeff) if coeff else 1
        if sign == "-":
            c = -c
        k = (int(exp) if exp else 1) if xpart else 0
        total += c * (-1 if k % 2 else 1)
        if pos >= len(s):
            break
    return abs(total)


def alex_degree_span(s):
    """Degree span (max-min exponent) of the Alexander string — for Laurent-safe
    deg(Delta) <= 2*genus."""
    if s is None:
        return None
    s = str(s).replace("**", "^")
    exps = [int(e) for e in re.findall(r"x\^(-?\d+)|t\^(-?\d+)", s)
            for e in e if e] if ("x" in s or "t" in s) else []
    # bare x/t (exponent 1) and constant (exponent 0):
    s_nospace = s.replace(" ", "")
    if re.search(r"[xt](?!\^)", s_nospace):
        exps.append(1)
    if re.search(r"(^|[+-])\d+($|[+-])", s_nospace):
        exps.append(0)
    if not exps:
        return 0 if _int(s) is not None else None
    return max(exps) - min(exps)


def pd_str_count(s):
    """Crossing number of the raw pd_list_str (same regex as r3_phantom_diag)."""
    if not s:
        return None
    n = len(re.findall(
        r"[\[(]\s*\d+\s*,\s*\d+\s*,\s*\d+\s*,\s*\d+\s*[\])]", str(s)))
    return n or None


def _kfh(alg):
    k = alg.get("knot_floer_homology")
    if isinstance(k, dict):
        return k
    if isinstance(k, str):
        try:
            import ast
            v = ast.literal_eval(k)
            return v if isinstance(v, dict) else {}
        except (ValueError, SyntaxError):
            return {}
    return {}


# ───────────────────────── TIER 0 — free gates ─────────────────────────────

def tier0_flags(complete: dict) -> list[str]:
    """Impossibility/consistency flags from a loaded _complete.json.
    Empty = clean. Every flag = Frankenstein evidence (Disc #251/#252 criteria
    + Alexander self-consistency). Missing fields: the check abstains."""
    flags = []
    pf = complete.get("pflichtliste") or {}
    diag = pf.get("DIAGRAMM") or {}
    alg = pf.get("ALGEBRA") or {}

    det = _int(alg.get("determinant"))
    det_sage = _int(alg.get("determinant_sage"))
    genus = _int(alg.get("genus"))
    cmin = _int(diag.get("crossing_number_min"))
    dt_len = _seq_len(diag.get("dt_code"))
    kfh = _kfh(alg)
    kfh_rank = _int(kfh.get("total_rank"))
    kfh_genus = _int(kfh.get("seifert_genus"))
    alex = alg.get("alexander_polynomial")

    # (1)+(2) diagram vs algebra layer (r3_phantom_diag, caught 7/7)
    if dt_len and cmin is not None and cmin > dt_len:
        flags.append(f"cmin {cmin} > dt {dt_len} (impossible)")
    if dt_len and genus is not None and genus > (dt_len - 1) / 2:
        flags.append(f"genus {genus} > (dt-1)/2 = {(dt_len - 1) / 2} (impossible)")
    # (3)+(4) KFH vs Algebra
    if kfh_rank and det and kfh_rank < det:
        flags.append(f"KFH total_rank {kfh_rank} < det {det} "
                     f"(impossible: rank >= det) -> KFH of ANOTHER knot")
    if kfh_genus is not None and genus is not None and kfh_genus != genus:
        flags.append(f"KFH genus {kfh_genus} != ALGEBRA genus {genus}")
    # (5) Alexander self-consistency
    if alex is not None and det is not None:
        d_alex = alex_eval_minus1(alex)
        if d_alex is not None and d_alex != det:
            flags.append(f"|Delta(-1)| = {d_alex} != det {det} "
                         f"(Alexander belongs to ANOTHER knot)")
    # (6) Seifert bound deg(Delta) <= 2*genus
    if alex is not None and genus is not None:
        span = alex_degree_span(alex)
        if span is not None and span > 2 * genus:
            flags.append(f"deg(Delta) {span} > 2*genus {2 * genus} (impossible)")
    # (7) two-witness divergence
    if det is not None and det_sage is not None and det != det_sage:
        flags.append(f"determinant {det} != determinant_sage {det_sage}")
    # (8) cmin vs raw PD
    pd_n = pd_str_count(complete.get("pd_list_str"))
    if pd_n and cmin is not None and cmin > pd_n:
        flags.append(f"cmin {cmin} > raw-PD {pd_n} (impossible)")
    return flags


def tier0_check_file(path) -> list[str]:
    """Wrapper for build_corpus_db: path -> flags ([] = clean, also on
    read error — structure duty is the schema validator's business)."""
    try:
        return tier0_flags(json.load(open(path)))
    except (OSError, json.JSONDecodeError, ValueError):
        return []


# ───────────────────────── TIER 1 — stored PD vs algebra (regina) ────────────

def tier1_flag(complete: dict):
    """det(stored PD after simplify) == ALGEBRA.det (Disc #252-Detektor).
    Needs regina; None = clean or not checkable."""
    # Regina runs in a process of its own (Werner 2026-08-29 „stell um" (EN: "switch it over")): Regina is
    # GPL v2+, a call is undisputedly free, an import in the same process is not.
    # Tier 0 thereby stays stdlib-pure — now even stricter than before.
    import os as _os, sys as _sys
    _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
    from regina_bridge import single as _regina_einer
    pd_s = complete.get("pd_list_str")
    pf = complete.get("pflichtliste") or {}
    det = _int((pf.get("ALGEBRA") or {}).get("determinant"))
    if not pd_s or det is None:
        return None
    pd = [list(map(int, t)) for t in re.findall(
        r"[\[(]\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*[\])]", str(pd_s))]
    if not pd:
        return None
    _r = _regina_einer({"pd": pd}, """
L = regina.Link.fromPD([list(map(int, t)) for t in auf['pd']])
L.simplify()
ax = L.alexander()
deg = ax.degree()
return {'det': abs(sum(int(str(ax[i])) * ((-1) ** i) for i in range(deg + 1)))}
""")
    if 'det' not in _r:
        return None                      # bridge failed -> as "not checkable"
    pd_det = _r['det']
    if pd_det != det:
        return (f"det(stored PD) = {pd_det} != ALGEBRA.det {det} "
                f"(Frankenstein: diagram and algebra = two objects)")
    return None


# ───────────────────────── TIER 2 — same-knot guard (fseries ladder) ─────────

def pick_fseries(folder: Path, stem: str):
    """Dense stem fseries; NEVER the rounded _sculpture_ variant
    (PD lesson 06-03, same rule as r3_corpus_objectaudit)."""
    cands = [p for p in Path(folder).glob("*.fseries")
             if not p.name.startswith("_sculpture")]
    if not cands:
        return None
    exact = [p for p in cands if p.stem == stem]
    return (exact or cands)[0]


def fseries_sha(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _det_mod_prime(M, p=None):
    """det of an integer matrix mod p, as a signed remainder. Point eval
    (Disc #248) instead of bignum.
    ⚠️ The remainder is the exact determinant ONLY if |det| < p/2. Whoever compares it
    with a raw determinant without this condition measures themselves
    — see det_stimmt() below and the block at DET_PRIMES."""
    p = p or DET_PRIME
    A = [[int(x) % p for x in row] for row in M]
    n = len(A)
    det = 1
    for i in range(n):
        piv = next((r for r in range(i, n) if A[r][i] % p), None)
        if piv is None:
            return 0
        if piv != i:
            A[i], A[piv] = A[piv], A[i]
            det = -det
        inv = pow(A[i][i], p - 2, p)
        det = det * A[i][i] % p
        for r in range(i + 1, n):
            f = A[r][i] * inv % p
            if f:
                A[r] = [(A[r][c] - f * A[i][c]) % p for c in range(n)]
    det %= p
    return det - p if det > p // 2 else det


def _signed(r, p):
    return r - p if r > p // 2 else r


def det_matches(S, rec_det):
    """Does the determinant of the reconstructed projection hit that of the record?

    The comparison is MODULAR ON BOTH SIDES, over several primes. Thus
    the check is independent of how large the determinant is — exactly
    the property the old single modulus lacked.
    Returns: (hits, displayed_det) — displayed_det is the exact number
    if it fits under the smallest bound, otherwise the text 'congruent'.
    """
    if rec_det is None:
        return False, None
    remainders = []
    for p in DET_PRIMES:
        r = _det_mod_prime(S, p) % p
        remainders.append((p, r))
        if r != (rec_det % p) and r != ((-rec_det) % p):
            # Sign with the SAME prime the remainder comes from — otherwise
            # the error message shows a number that exists nowhere.
            p0, r0 = remainders[0]
            return False, _signed(r0, p0)
    p0 = DET_PRIMES[0]
    exact_val = abs(_signed(remainders[0][1], p0)) if abs(rec_det) < p0 // 2 else None
    return True, (exact_val if exact_val is not None else f"congruent mod {len(DET_PRIMES)} primes")


VOL_RANDOMIZE = 30          # like the calibrated ladder in certify_hyperbolicity.py


def _volume_escalated(E, log=None):
    """Volume of a manifold — with ESCALATION instead of a single attempt.

    🔴 Why (Werner 2026-08-27: „den guard reparieren" (EN: "repair the guard")): up to here it read
        vol_gueltig = str(E.solution_type()).startswith("all tetrahedra positively oriented")
        if not vol_gueltig: vol = None
    — **one** attempt, and final on failure. For large knots the
    first triangulation is almost never positively oriented; the guard then reported
    NO_DATA and thereby measured the mood of the first triangulation, not the knot.
    Caught on 27.08. at `UEBER_k3_s3` (247 crossings): the guard said
    `vol=None`, while `certify_hyperbolicity.py` **rigorously proved**
    the same manifold hyperbolic with `ManifoldHP + randomize#4`
    (interval arithmetic). det and |sig| already hit in the guard.

    The same ladder as there — Manifold → ManifoldHP → randomize → canonical —,
    measured only by `solution_type`, because no Sage is available here
    (`verify_hyperbolicity` needs it). So this is **weaker than a
    proof** and only means: the triangulation is not degenerate.

    ⚠️ The ladder must NOT turn a degenerate manifold into a valid one.
    Exactly that is what the counter-case in `--self-test-vol` checks: a
    connected sum stays without volume even after full escalation.

    Returns: (vol_or_None, path).
    """
    def gut(M):
        try:
            if not str(M.solution_type()).startswith("all tetrahedra positively oriented"):
                return None
            return float(M.volume())
        except Exception:
            return None

    v = gut(E)
    if v is not None:
        return v, "Manifold"
    MH = None
    try:
        MH = E.high_precision()
        v = gut(MH)
        if v is not None:
            return v, "ManifoldHP"
    except Exception:
        MH = None
    if MH is not None:
        for i in range(1, VOL_RANDOMIZE + 1):
            try:
                MH.randomize()
            except Exception:
                break
            v = gut(MH)
            if v is not None:
                return v, f"ManifoldHP+randomize#{i}"
        try:
            v = gut(MH.canonical_retriangulation())
            if v is not None:
                return v, "canonical_retriangulation"
        except Exception:
            pass
    return None, "fully escalated, no valid triangulation"


def _compare(actual, expected_, applies):
    """Three-valued instead of two-valued: 'trifft' | 'weicht_ab' | 'kein_sollwert'
    | 'nicht_gerechnet'.

    🔴 Why (Werner 2026-08-27: „zieh das nach" (EN: "bring that up to date")): up to here
    `det_stimmt(S, None)` simply returned `False` — **"I have nothing to compare"
    could not be told from "it contradicts itself".** On 27.08. this reported
    `UEBER_k3_s4` as FRANKENSTEIN although volume and signature
    hit: regina's symbolic Alexander polynomial had run into the
    120 s latch in stage 5, hence the field `determinant` was missing — the VALUE
    was receipted three times over, only under other keys.
    An alarm that screams even when the datum is missing is soon overlooked.
    (The signature was always treated correctly: `rec_sig is None or …`.
    The asymmetry was the error.)
    """
    if expected_ is None:
        return "kein_sollwert"
    if actual is None:
        return "nicht_gerechnet"
    return "trifft" if applies else "weicht_ab"


def _verdict(vol_ok, rec_vol, vol, det_st, sig_st):
    """The decision, factored out of the expensive machinery — so that it can
    be tested without a single knot (`--self-test-verdict`)."""
    contradiction = "weicht_ab" in (det_st, sig_st)
    receipted = sum(1 for x in (det_st, sig_st) if x == "trifft")
    if vol_ok:
        if contradiction:
            return "FRANKENSTEIN"          # geometry matches, algebra contradicts
        if receipted >= 1:
            return "OK"
        # Volume hits, but NO algebraic quantity could be checked.
        # That is weaker than OK and therefore gets a name of its own —
        # a green over an empty check set would be blindness (footgun 32).
        return "OK_NUR_GEOMETRIE"
    if rec_vol is None:
        if vol is not None:
            return "WIDERSPRUCH_HYPERBOLISCH"
        if contradiction:
            return "NO_DATA"
        if receipted >= 1:
            return "OK_ALGEBRAISCH"
    return "NO_DATA"


def self_test_verdict():
    """Counter-cases for the verdict logic. Purely logical, no computation —
    hence it runs in seconds and is actually run."""
    cases = [
        ("A everything hits",                 (True,  9.0, 9.0, "trifft",       "trifft"),        "OK"),
        ("B det contradicts",             (True,  9.0, 9.0, "weicht_ab",    "trifft"),        "FRANKENSTEIN"),
        ("C |sig| contradicts",           (True,  9.0, 9.0, "trifft",       "weicht_ab"),     "FRANKENSTEIN"),
        ("D det MISSING, sig hits",        (True,  9.0, 9.0, "kein_sollwert","trifft"),        "OK"),
        ("E both target values missing",       (True,  9.0, 9.0, "kein_sollwert","kein_sollwert"), "OK_NUR_GEOMETRIE"),
        ("F record without vol, curve with",   (False, None, 9.0, "trifft",      "trifft"),        "WIDERSPRUCH_HYPERBOLISCH"),
        ("G both without vol, algebra ok",   (False, None, None,"trifft",      "trifft"),        "OK_ALGEBRAISCH"),
        ("H both without vol, algebra empty", (False, None, None,"kein_sollwert","kein_sollwert"),"NO_DATA"),
        ("I both without vol, contradiction",  (False, None, None,"weicht_ab",   "trifft"),        "NO_DATA"),
        ("J vol missed",                 (False, 9.0, 3.0, "trifft",       "trifft"),        "NO_DATA"),
    ]
    bad = 0
    for name, arg, expected_ in cases:
        actual = _verdict(*arg)
        gut = actual == expected_
        bad += not gut
        print(f"  {'✅' if gut else '🔴'} {name:32s} -> {actual:26s} expected {expected_}")
    print(f"\n  {'ALL GREEN' if not bad else f'{bad} FAILED'}")
    return 0 if not bad else 1


def self_test_vol():
    """Counter-cases for the volume escalation. WITHOUT them the ladder is decoration —
    one that eventually squeezes a volume out of ANY manifold would be
    worse than the single attempt from before.

      A hyperbolic (4_1)                    -> volume found
      B COUNTER-CASE connected sum 3_1#3_1  -> even fully escalated NO volume
      C COUNTER-CASE torus knot T(2,5)      -> even fully escalated NO volume
    If B or C fails, the ladder turns degeneracy into validity.
    """
    import spherogram, snappy  # noqa: F401
    # Braid closures instead of PD lists: spherogram builds them directly, and regina
    # is not in this environment (gnottheorie_env).
    cases = [
        ("A 4_1      hyperbolic",        [1, -2, 1, -2],          True),
        ("B 3_1#3_1  connected sum",    [1, 1, 1, 2, 2, 2],      False),
        ("C T(2,5)   torus knot",         [1, 1, 1, 1, 1],         False),
        ("D 5_2      hyperbolic",        [1, 1, 1, 2, -1, 2],     True),
    ]
    bad = 0
    for name, word, expected_vol in cases:
        L = spherogram.Link(braid_closure=word)
        L.simplify(mode="global")
        vol, weg = _volume_escalated(L.exterior())
        gut = (vol is not None) == expected_vol
        bad += not gut
        print(f"  {'✅' if gut else '🔴'} {name:32s} vol={('%.6f' % vol) if vol else 'None':>12}  "
              f"[{weg}]   expected {'volume' if expected_vol else 'NO volume'}")
    print(f"\n  {'ALL GREEN' if not bad else f'{bad} FAILED'}")
    return 0 if not bad else 1


def _tier2_work(task):
    """One record: dense reconstruction @GUARD_N -> tone (snappy vol) + point eval
    (det/sig via spherogram seifert_matrix @ t=-1). vol match certifies the
    projection (invariant, projection-independent); Frankenstein = vol matches,
    algebra does not. Runs in the worker (gnottheorie_env: snappy+spherogram+numpy)."""
    stem, fs_path, rec_det, rec_deg, rec_vol, rec_sig = task
    import types
    import numpy as np
    for _n, _b in (("complex", complex), ("float", float),
                   ("int", int), ("bool", bool)):
        if not hasattr(np, _n):
            setattr(np, _n, _b)
    # r3_polygon_regina_det imports regina at top level (only numpy functions used) -> stub it
    sys.modules.setdefault("regina", types.ModuleType("regina"))
    sys.path.insert(0, str(TOOLS))
    sys.path.insert(0, str(HERE))
    from build_curve_json import parse_fseries, reconstruct_curve
    from r3_polygon_regina_det import rot, pd_from_polygon
    import spherogram
    import snappy  # noqa: F401  — spherogram.Link.exterior() needs snappy in the namespace

    offset, a, b = parse_fseries(fs_path)
    pts = reconstruct_curve(offset, a, b, GUARD_N)
    best = {"verdict": "NO_DATA", "record_det": rec_det, "record_vol": rec_vol,
            "note": "no projection reproduces stored vol"}
    for seed in GUARD_SEEDS:
        try:
            pd, C = pd_from_polygon(pts @ rot(seed).T)
            if not pd:
                continue
            L = spherogram.Link(pd)
            L.simplify(mode="global")
            E = L.exterior()
            vol, vol_route = _volume_escalated(E)
            # snappy PRINTS a number even for non-hyperbolic knots and marks it
            # invalid only via solution_type. Whoever does not ask that
            # reads a degenerate number as a volume. This is exactly how this
            # guard displayed gitter_d3_k3_s1p5_g0p25 (8-fold connected sum,
            # det 19^8) on 24.08. with "vol=92.606934" — the number belongs to
            # 'contains degenerate tetrahedra' and means: NO volume.
            # Cf. Disc #687, the same error class one level deeper.
            # (the check now sits in `_volumen_eskaliert`, together with the ladder)
        except Exception:
            continue
        vol_ok = (rec_vol is not None and vol is not None
                  and abs(vol - rec_vol) < VOL_TOL)
        det = sig = None
        det_ok = False
        try:
            V = np.array(L.seifert_matrix(), dtype=object)
            S = (V + V.T).tolist()
            det_ok, det = det_matches(S, rec_det)
            ev = np.linalg.eigvalsh(np.array(S, dtype=float))
            sig = int((ev > 1e-8).sum() - (ev < -1e-8).sum())
        except Exception:
            pass
        det_st = _compare(det, rec_det, det_ok)
        sig_st = _compare(None if sig is None else abs(sig),
                           None if rec_sig is None else abs(rec_sig),
                           sig is not None and rec_sig is not None
                           and abs(sig) == abs(rec_sig))
        sig_ok = sig_st in ("trifft", "kein_sollwert")
        info = {"seed": seed, "raw_crossings": C, "vol_weg": vol_route, "simpl_crossings": len(L.crossings),
                "vol": None if vol is None else round(vol, 6),
                "vol_gueltig": vol is not None,
                "det": det, "abs_sig": None if sig is None else abs(sig),
                "record_det": rec_det, "record_vol": rec_vol, "record_sig": rec_sig,
                "vol_ok": vol_ok, "det_ok": det_ok, "sig_ok": sig_ok,
                "det_status": det_st, "sig_status": sig_st}
        verdict = _verdict(vol_ok, rec_vol, vol, det_st, sig_st)
        if vol_ok:
            info["verdict"] = verdict
            return stem, info
        if rec_vol is None:
            # The record has NO volume — for torus knots and
            # connected sums that is the truth, not a gap. A
            # vol match can never exist here; whoever checks only that leaves
            # the whole known-non-hyperbolic class (41 records) uncertified
            # forever. So the algebraic route — but NARROW:
            if vol is not None:
                # Record says "no volume", the curve yields a valid one.
                # That is a contradiction and is reported, not waved through.
                info["verdict"] = "WIDERSPRUCH_HYPERBOLISCH"
                return stem, info
            if verdict == "OK_ALGEBRAISCH":
                # Both sides say "not hyperbolic", and det AND |sig| of the
                # independently reconstructed projection hit the record.
                # Weaker than a vol certificate and therefore a verdict OF ITS OWN —
                # it must never be confused with "OK".
                info["verdict"] = verdict
                return stem, info
        best = {**info, "verdict": "NO_DATA"}   # vol miss: next seed
    return stem, best


# ───────────────────────── Corpus / Registry ─────────────────────────────────

def load_corpus():
    return json.load(open(CORPUS_DB, encoding="utf-8"))


def in_scope(rec) -> bool:
    if rec.get("category") not in ("machine", "own"):
        return False
    if (rec.get("cluster_folder") or "") in SCOPE_EXCLUDED_CLUSTERS:
        return False
    return True


def rec_deg_from_alex(rec):
    degs = [int(m) for m in re.findall(r"x\^(\d+)", str(rec.get("alex") or ""))]
    return max(degs) if degs else None


def load_registry() -> dict:
    if REGISTRY.exists():
        return json.loads(REGISTRY.read_text())
    return {"_meta": {}, "records": {}}


def save_registry(reg):
    reg["_meta"]["updated_at"] = datetime.now().isoformat()
    REGISTRY.write_text(json.dumps(reg, indent=1, sort_keys=True))


def registry_coverage(db=None, reg=None):
    """Per scope record: covered (fseries hash + det match the registry) or not.
    -> (covered:list, uncovered:list[(stem, reason)]). Used by the crosscheck gate
    — pure stdlib, no regina."""
    db = db or load_corpus()
    reg = reg or load_registry()
    recs = reg.get("records", {})
    covered, uncovered = [], []
    for r in db["knots"]:
        if not in_scope(r):
            continue
        stem = r.get("stem")
        folder = Path(r["pflicht_path"]).parent
        fs = pick_fseries(folder, stem)
        e = recs.get(stem)
        if e is None:
            uncovered.append((stem, "no registry entry (new record? -> --guard-new)"))
            continue
        if fs is None:
            uncovered.append((stem, "no .fseries in the folder"))
            continue
        if e.get("fseries_sha256") != fseries_sha(fs):
            uncovered.append((stem, "fseries changed since the guard run -> --guard-new"))
            continue
        if e.get("det") != r.get("det"):
            uncovered.append((stem, f"det drift: Registry {e.get('det')} vs corpus "
                                    f"{r.get('det')} -> --guard-new"))
            continue
        ev, rv = e.get("vol"), r.get("vol")
        if ev is not None and rv is not None and abs(ev - rv) > VOL_TOL:
            uncovered.append((stem, f"vol drift: Registry {ev} vs corpus {rv} -> --guard-new"))
            continue
        # OK_NUR_GEOMETRIE counts as checked: the volume is a
        # projection-INDEPENDENT invariant and it hit. It is, however,
        # weaker than OK (no algebraic quantity was checkable) and
        # therefore carries a name of its own — never merge.
        if e.get("verdict") not in ("OK", "SEEDED", "OK_ALGEBRAISCH",
                                    "OK_NUR_GEOMETRIE"):
            uncovered.append((stem, f"guard verdict {e.get('verdict')}"))
            continue
        covered.append(stem)
    return covered, uncovered


# ───────────────────────── CLI modes ──────────────────────────────────────────

def cmd_scan(with_pd=False):
    db = load_corpus()
    n_checked = n_flagged = 0
    for r in db["knots"] + db.get("validation_knots", []):
        p = Path(r["pflicht_path"])
        if not p.exists():
            print(f"  ! {r.get('stem')}: pflicht_path missing: {p}")
            continue
        d = json.load(open(p))
        flags = tier0_flags(d)
        if with_pd:
            try:
                f1 = tier1_flag(d)
                if f1:
                    flags.append(f1)
            except Exception as e:
                flags.append(f"Tier-1 not checkable: {type(e).__name__}: {e}")
        n_checked += 1
        if flags:
            n_flagged += 1
            print(f"FLAG  {r.get('stem')}  [{r.get('cluster_folder')}]")
            for f in flags:
                print(f"      🚩 {f}")
    print(f"\nTier-0{'+ Tier-1' if with_pd else ''}: {n_checked} records checked, "
          f"{n_flagged} flagged")
    return 1 if n_flagged else 0


def cmd_seed():
    """Seal the existing stock: provenance = full-corpus audit 2026-06-04 (Disc #253,
    'all 268 measured or documented-excluded'). Records WHICH
    RESULTS file carries the stem, where findable."""
    dur = ROOT / "claude_setup" / "tmp_durable"
    sources = {}
    for fn in ("r3_corpus_objectaudit_RESULTS.json", "r3_s18_cmin_RESULTS.json",
               "r3_phantom_diag_RESULTS.json", "r3_phaseC_object_RESULTS.json"):
        p = dur / fn
        if not p.exists():
            continue
        d = json.loads(p.read_text())
        items = d if isinstance(d, list) else (
            [{"stem": k} for k in d] if isinstance(d, dict) else [])
        for it in items:
            s = it.get("stem") if isinstance(it, dict) else None
            if s:
                sources.setdefault(s, fn)
    db = load_corpus()
    reg = load_registry()
    reg["_meta"].setdefault(
        "seeded",
        "2026-06-04 full-corpus audit Disc #253 — all existing records measured "
        "or documented-excluded; hashes = fseries state at seed time")
    n = 0
    for r in db["knots"]:
        if not in_scope(r):
            continue
        stem = r.get("stem")
        if stem in reg["records"]:
            continue
        fs = pick_fseries(Path(r["pflicht_path"]).parent, stem)
        if fs is None:
            print(f"  ! {stem}: no .fseries — not seeded")
            continue
        # Stem match in the RESULTS files is partly prefix-free (kf86 vs kf86_sh-…)
        src = sources.get(stem) or next(
            (v for k, v in sources.items() if k and (stem.startswith(k) or k.startswith(stem))),
            None)
        reg["records"][stem] = {
            "fseries_sha256": fseries_sha(fs),
            "det": r.get("det"),
            "deg_alex": rec_deg_from_alex(r),
            "verdict": "SEEDED",
            "source": src or "korpus-voll-audit-2026-06-04 (Disc #253)",
            "checked_at": datetime.now().isoformat(),
        }
        n += 1
    save_registry(reg)
    print(f"Registry seeded: {n} new entries, total {len(reg['records'])}")


def cmd_guard_new(jobs):
    """Tier 2 for all uncovered scope records, in parallel (full-throttle default)."""
    import multiprocessing as mp
    db = load_corpus()
    reg = load_registry()
    _, uncovered = registry_coverage(db, reg)
    todo = []
    by_stem = {r.get("stem"): r for r in db["knots"]}
    for stem, why in uncovered:
        r = by_stem.get(stem)
        if r is None:
            continue
        fs = pick_fseries(Path(r["pflicht_path"]).parent, stem)
        if fs is None:
            print(f"  ! {stem}: no .fseries ({why}) — guard not possible")
            continue
        todo.append((stem, str(fs), r.get("det"), rec_deg_from_alex(r),
                     r.get("vol"), r.get("sig")))
    if not todo:
        print("Nothing to guard — all scope records covered.")
        return 0
    print(f"same-knot guard (tone+point eval) for {len(todo)} records "
          f"(N={GUARD_N}, Seeds {GUARD_SEEDS}, {jobs} workers) …")
    rc = 0
    with mp.Pool(jobs) as pool:
        for stem, res in pool.imap_unordered(_tier2_work, todo):
            r = by_stem[stem]
            fs = pick_fseries(Path(r["pflicht_path"]).parent, stem)
            reg["records"][stem] = {
                "fseries_sha256": fseries_sha(fs),
                "det": r.get("det"),
                "vol": r.get("vol"),
                "deg_alex": rec_deg_from_alex(r),
                "verdict": res["verdict"],
                "detail": {k: v for k, v in res.items() if k != "verdict"},
                "source": "frankenstein_gate --guard-new (Ton+Punkt-Eval)",
                "checked_at": datetime.now().isoformat(),
            }
            mark = ("✓" if res["verdict"] in ("OK", "OK_ALGEBRAISCH")
                    else "◐" if res["verdict"] == "OK_NUR_GEOMETRIE" else "🚩")
            print(f"  {mark} {stem}: {res['verdict']}  "
                  f"vol={res.get('vol','?')}(soll {res.get('record_vol','?')}) "
                  f"det={res.get('det','?')}(soll {res.get('record_det','?')}) "
                  f"|sig|={res.get('abs_sig','?')}")
            if res["verdict"] not in ("OK", "OK_ALGEBRAISCH", "OK_NUR_GEOMETRIE"):
                rc = 1
            save_registry(reg)  # after every record — a crash loses nothing
    return rc


def cmd_status():
    covered, uncovered = registry_coverage()
    print(f"Registry coverage: {len(covered)} covered, {len(uncovered)} open")
    for stem, why in uncovered:
        print(f"  - {stem}: {why}")
    return 1 if uncovered else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--scan", action="store_true", help="Tier 0 over the whole corpus")
    ap.add_argument("--pd", action="store_true", help="with --scan also Tier 1 (regina)")
    ap.add_argument("--seed-registry", action="store_true")
    ap.add_argument("--guard-new", action="store_true")
    ap.add_argument("--self-test-verdict", action="store_true",
                    help="Counter-cases for the verdict logic — purely logical, seconds")
    ap.add_argument("--self-test-vol", action="store_true",
                    help="Counter-cases for the volume escalation (Disc: guard repair 27.08.)")
    ap.add_argument("--registry-status", action="store_true")
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--guard-n", type=int, default=None,
                    help="Points of the curve reconstruction (default GUARD_N=2048). "
                         "Set higher if the reconstruction changes the KNOT: "
                         "on 24.08. gitter_d3_k3_s1p5_g0p25 (8-fold "
                         "connected sum, 8 closely adjacent lobes) gave at 2048 "
                         "det=487911 instead of 19^8 — the pipeline itself had scanned this "
                         "curve with 131072 frames and seen 237..701 crossings "
                         "per frame. 2048 is an ASSUMPTION, not a measurement.")
    args = ap.parse_args()
    if args.guard_n:
        globals()["GUARD_N"] = args.guard_n
        print(f"[guard] reconstruction with N={args.guard_n} instead of the default 2048",
              flush=True)
    if args.self_test_verdict:
        sys.exit(self_test_verdict())
    if args.self_test_vol:
        sys.exit(self_test_vol())
    if args.scan:
        sys.exit(cmd_scan(with_pd=args.pd))
    if args.seed_registry:
        cmd_seed()
        sys.exit(0)
    if args.guard_new:
        sys.exit(cmd_guard_new(args.jobs))
    if args.registry_status:
        sys.exit(cmd_status())
    ap.print_help()


if __name__ == "__main__":
    main()
