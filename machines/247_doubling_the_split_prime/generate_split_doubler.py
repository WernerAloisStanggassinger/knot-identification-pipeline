# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 247 — doubling_the_split_prime (cyclotomic module machine / split-prime spectrum).

Generates for primes (p, ℓ) the periodic pretzel knot P(ℓ^p) = P(ℓ,…,ℓ) [p bands].
⚠️ CORRECTION 2026-06-13: the GROUP H₁(Σ₂) is ALWAYS (Z/ℓ)^{p−1}⊕Z/pℓ (symmetry-trivial, does not see
ℓ mod p — 3173/3173 confirmed). The Naik phenomenon sits in the Z/p MODULE: the (Z/ℓ)^{p−1} stage
= F_ℓ[x]/Φ_p(x) splits like ℓ in ℤ[ζ_p] (ℓ≡1 split → p−1 characters; inert → 1 block; partial
in between). ℓ mod p STEERS the spectrum. „Verdoppler" (doubler) = historical name.

Three layers:
  h1_guarantee(p, ℓ)      — GROUP (Smith, always p−1) + MODULE pointer.
  module_decomposition(p,ℓ) — THE Naik spectrum: Φ_p mod ℓ factorization (split/inert/partial).
  pretzel_embedding(p,ℓ)  — C_p-symmetric 3D realization (p pillars of ℓ half-twists each + closure).
  Corpus anchor: P(7,7,7) form [7,21]/det 147 = corpus S_18_d147_Z3 (same FORM, different knot c≠).

Math derivation + bets: 247_split_prim_verdoppler.txt.
This file = construction; pipeline identification stays blind (Hard-Rule 1).
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
YARD = HERE.parent
sys.path.insert(0, str(YARD))
from machine_common import (  # noqa: E402
    densify, embedding_report, polyline_to_fseries,
    count_projection_crossings)

try:
    import sympy as sp
    from sympy.matrices.normalforms import smith_normal_form
    _HAVE_SYMPY = True
except Exception:  # pragma: no cover
    _HAVE_SYMPY = False


# ───────────────────────────── math core (proven) ──────────────────────────
def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def plumbing_matrix(p: int, ell: int):
    """Star plumbing of the pretzel double cover: central knot (framing 0)
    + p arms (framing ℓ). coker = H₁(Σ₂(P(ℓ^p)))."""
    if not _HAVE_SYMPY:
        raise RuntimeError("sympy required for h1_guarantee (conda env base/gnottheorie).")
    n = p + 1
    M = sp.zeros(n, n)
    for i in range(1, n):
        M[0, i] = 1
        M[i, 0] = 1
        M[i, i] = ell
    return M


def _ord_mod(ell: int, p: int):
    """Multiplicative order of ℓ mod p (None if p | ℓ = ramified)."""
    r = ell % p
    if r == 0:
        return None
    o, cur = 1, r
    while cur != 1:
        cur = (cur * r) % p
        o += 1
    return o


def module_decomposition(p: int, ell: int) -> dict:
    """The ACTUAL Naik content: the Z/p module structure of the (Z/ℓ)^{p−1} stage
    = F_ℓ[x]/Φ_p(x), splits like ℓ in ℤ[ζ_p] (= Φ_p mod ℓ factorizes).
    Delivers type + block count/degree. ℓ mod p steers the spectrum, NOT the group rank."""
    f = _ord_mod(ell, p)
    if f is None:
        return {"split_type": "ramified", "ord_p_ell": None,
                "n_blocks": None, "block_degree": None, "n_characters": 0}
    n_blocks = (p - 1) // f
    typ = ("split" if f == 1 else "inert" if f == p - 1 else "partial")
    # explicit factorization of Φ_p mod ℓ as a cross-check
    _X = sp.symbols('x')
    facs = sp.Poly(sp.cyclotomic_poly(p, _X), _X, modulus=ell).factor_list()[1]
    degs = sorted(g.degree() for g, m in facs for _ in range(m))
    return {"split_type": typ, "ord_p_ell": f,
            "n_blocks": n_blocks, "block_degree": f,
            "n_characters": n_blocks if f == 1 else 0,
            "phi_p_mod_ell_degrees": degs,
            "consistent": (degs == [f] * n_blocks)}


def h1_guarantee(p: int, ell: int) -> dict:
    """H₁(Σ₂(P(ℓ^p))) — GROUP (Smith) + MODULE (Naik spectrum).
    ⚠️ Correction 2026-06-13: the group ℓ-rank is ALWAYS p−1 (symmetry-trivial,
    does NOT see the decomposition of ℓ). The Naik phenomenon sits in the Z/p module (see
    module_decomposition). 3173/3173 exhaustively confirmed (p≤71, ℓ≤997)."""
    M = plumbing_matrix(p, ell)
    D = smith_normal_form(M)
    facs = sorted(abs(int(D[i, i])) for i in range(D.rows) if D[i, i] != 0)
    nontrivial = [f for f in facs if f > 1]
    det = abs(int(M.det()))
    pred_det = p * ell ** (p - 1)
    ell_rank = sum(1 for f in facs if f % ell == 0)
    return {
        "p": p, "ell": ell,
        "det": det,
        "det_predicted_p_l_pow_pm1": pred_det,
        "det_ok": det == pred_det,
        "invariant_factors": facs,
        "H1_nontrivial": nontrivial,
        # GROUP — always p−1, symmetry-trivial (NOT Naik):
        "group_ell_rank": ell_rank,
        "group_ell_rank_always_pm1": ell_rank == p - 1,
        # MODULE — THE Naik spectrum (steerable via ℓ mod p):
        "module": module_decomposition(p, ell),
    }


# ───────────────────────── geometric realization ─────────────────────────
def pretzel_embedding(p: int, ell: int, *, R: float = 1.0,
                      rho_frac: float = 0.45, arc_bump: float = 0.35,
                      pts_per_band: int = 240, pts_per_arc: int = 24
                      ) -> np.ndarray:
    """C_p-symmetric embedding of P(ℓ^p) as a closed polyline.

    p twist pillars (azimuth 2πk/p, radius R), two antipodal strands each,
    rotating over z∈[−1,1] by ℓ·π (ℓ half-twists). Tangential cap ports;
    cyclic top/bottom arcs connect band k with k+1 (mod p).
    Free of self-intersection as long as rho < R·sin(π/p).
    """
    if p % 2 == 0 or ell % 2 == 0:
        raise ValueError(f"P({ell}^{p}): p={p}, ℓ={ell} must both be odd (otherwise no knot).")
    rho = rho_frac * R * np.sin(np.pi / p)   # guaranteed < R·sin(π/p)

    def frame(k):
        phi = 2 * np.pi * k / p
        rhat = np.array([np.cos(phi), np.sin(phi), 0.0])
        that = np.array([-np.sin(phi), np.cos(phi), 0.0])
        return rhat, that

    def strand(k, s):
        """Strand s∈{0,1} of band k, from z=−1 (bottom) to z=+1 (top)."""
        rhat, that = frame(k)
        t = np.linspace(0.0, 1.0, pts_per_band)
        z = -1.0 + 2.0 * t
        psi = ell * np.pi * t + s * np.pi
        off = (np.cos(psi)[:, None] * that[None, :]
               + np.sin(psi)[:, None] * rhat[None, :]) * rho
        return R * rhat[None, :] + off + z[:, None] * np.array([0, 0, 1.0])[None, :]

    def arc(pA, pB, bump):
        """Smooth arc pA→pB, slightly arched in z (bump>0 = upward)."""
        mid = 0.5 * (pA + pB)
        mid = mid.copy()
        mid[2] += bump
        t = np.linspace(0.0, 1.0, pts_per_arc)[:, None]
        # quadratic Bézier through (pA, ctrl, pB) with ctrl such that apex ~mid
        ctrl = 2 * mid - 0.5 * (pA + pB)
        return ((1 - t) ** 2 * pA[None, :]
                + 2 * (1 - t) * t * ctrl[None, :]
                + t ** 2 * pB[None, :])

    # strand polylines per band: strand0 B(+)→T(−), strand1 B(−)→T(+)
    S = {(k, s): strand(k, s) for k in range(p) for s in (0, 1)}
    # ports (coordinates = strand end points, exact):
    #   strand0: bot=(k,'B','+')  top=(k,'T','-')
    #   strand1: bot=(k,'B','-')  top=(k,'T','+')
    port = {}
    for k in range(p):
        port[(k, 'B', '+')] = S[(k, 0)][0]
        port[(k, 'T', '-')] = S[(k, 0)][-1]
        port[(k, 'B', '-')] = S[(k, 1)][0]
        port[(k, 'T', '+')] = S[(k, 1)][-1]

    # edges as (portA, portB, points[portA..portB]):
    edges = []
    for k in range(p):
        edges.append(((k, 'B', '+'), (k, 'T', '-'), S[(k, 0)]))   # strand0
        edges.append(((k, 'B', '-'), (k, 'T', '+'), S[(k, 1)]))   # strand1
    for k in range(p):
        kn = (k + 1) % p
        a, b = (k, 'T', '+'), (kn, 'T', '-')
        edges.append((a, b, arc(port[a], port[b], +arc_bump)))   # top arc
        a, b = (k, 'B', '+'), (kn, 'B', '-')
        edges.append((a, b, arc(port[a], port[b], -arc_bump)))   # bottom arc

    # incidence: every port lies on exactly 2 edges (1 strand + 1 arc)
    inc = {}
    for idx, (a, b, _) in enumerate(edges):
        inc.setdefault(a, []).append(idx)
        inc.setdefault(b, []).append(idx)
    assert all(len(v) == 2 for v in inc.values()), "Port-Inzidenz != 2 — Closure kaputt"

    # walk: traverse one closed component
    used = set()
    start_port = (0, 'B', '+')
    cur_port = start_port
    out = []
    for _ in range(len(edges)):
        ei = next((e for e in inc[cur_port] if e not in used), None)
        if ei is None:
            break
        used.add(ei)
        a, b, pts = edges[ei]
        if cur_port == a:
            out.append(pts)
            cur_port = b
        else:
            out.append(pts[::-1])
            cur_port = a
    if len(used) != len(edges):
        raise RuntimeError(f"P({ell}^{p}): walk closed only {len(used)}/{len(edges)} "
                           f"edges ⟹ several components (not a knot).")
    if cur_port != start_port:
        raise RuntimeError(f"P({ell}^{p}): walk not closed.")
    # join, drop duplicate seam points
    P = np.vstack([seg if i == 0 else seg[1:] for i, seg in enumerate(out)])
    return P


# ───────────────────────────────── generation ────────────────────────────────
def make_knot(p: int, ell: int, out_root: Path) -> tuple[dict, Path | None]:
    name = f"sd_p{p}_l{ell}"
    entry = {"name": name, "machine": "247_doubling_the_split_prime",
             "p": p, "ell": ell, "pretzel": f"P({ell}^{p})"}
    if not (_is_prime(p) and _is_prime(ell)):
        entry["rejected"] = "p,ℓ müssen prim sein"
        return entry, None
    if p % 2 == 0 or ell % 2 == 0 or ell == p:
        entry["rejected"] = f"p={p},ℓ={ell}: beide ungerade prim + ℓ≠p nötig (sonst kein Knoten)"
        return entry, None
    # ℓ mod p free: the machine realizes the WHOLE spectrum (split/inert/partial);
    # ℓ≡1 = maximal Naik case, but NOT forced (correction 2026-06-13).
    entry["h1"] = h1_guarantee(p, ell)
    P = pretzel_embedding(p, ell)
    # side projection (x/y) = diagram view ≈ p·ℓ twist crossings + a few
    # arc crossings. The z axis looks at the twist helices END-ON
    # (each turn projects onto the same circle → apparent crossings) and is
    # NOT a meaningful diagram view — hence min over x/y.
    n_side = min(count_projection_crossings(P, axis=0),
                 count_projection_crossings(P, axis=1))
    entry["projection_crossings_side"] = int(n_side)
    entry["crossings_expected_p_times_l"] = p * ell
    Pd = densify(P)
    emb = embedding_report(Pd)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = out_root / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = None
    for nt in (256, 384, 512, 768, 1024):
        try:
            rep = polyline_to_fseries(
                P, folder / f"{name}.fseries", knot_name=name,
                source="maschinen_ideen/247", n_terms=nt, nyquist_max_rel=0.02)
            break
        except RuntimeError:
            continue
    if rep is None:
        entry["rejected"] = "Nyquist rot bis n_terms=1024"
        return entry, folder
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


# ──────────────────────────────── commands ─────────────────────────────────
ORACLE = {"p": 3, "ell": 7, "det": 147, "facs": [1, 1, 7, 21], "group_ell_rank": 2,
          "corpus": "S_18_d147_Z3 (H₁=[7,21]) — gleiche Verschlingungsform, ANDERER Knoten (c≠)"}
# Diverse spectrum cases: split (max), inert, partial — shows: group constant, module varies.
SMOKE_CASES = [(3, 7, "split"), (3, 5, "inert"), (5, 11, "split"),
               (5, 19, "partial"), (7, 29, "split"), (7, 3, "inert"), (7, 2, "partial")]


def cmd_demo(_args):
    print("[DEMO] Smoke test machine 247 (cyclotomic module machine / split-prime spectrum)")
    ok = 0
    total = 0

    # S1: group H₁ against corpus linking-form oracle (NOT knot identity)
    total += 1
    g = h1_guarantee(ORACLE["p"], ORACLE["ell"])
    s1 = (g["det"] == ORACLE["det"]
          and g["invariant_factors"] == ORACLE["facs"]
          and g["group_ell_rank"] == ORACLE["group_ell_rank"]
          and g["det_ok"])
    print(f"  (S1) P(7,7,7) group vs corpus form {ORACLE['corpus']}: "
          f"det {g['det']} · factors {g['invariant_factors']} · "
          f"7-rank {g['group_ell_rank']} {'✓ SAME FORM' if s1 else '✗'}")
    ok += s1

    # S2: group ALWAYS p−1 (symmetry-trivial); MODULE carries the Naik spectrum (split/inert/partial)
    total += 1
    s2 = True
    for p, ell, exp_type in SMOKE_CASES:
        g = h1_guarantee(p, ell)
        m = g["module"]
        good = (g["det_ok"] and g["group_ell_rank_always_pm1"]
                and m["split_type"] == exp_type and m.get("consistent", False))
        s2 &= good
        print(f"  (S2) P({ell}^{p}): det {g['det']} · GROUP ℓ-rank {g['group_ell_rank']}(=p−1 always) "
              f"· MODULE {m['split_type']}: {m['n_blocks']}×degree{m['block_degree']} "
              f"(ord={m['ord_p_ell']}) {'✓' if good else '✗'}")
    ok += s2

    # S3: geometric realization P(7,7,7) — one component + embedding
    total += 1
    try:
        P = pretzel_embedding(3, 7)
        n_side = min(count_projection_crossings(P, axis=0),
                     count_projection_crossings(P, axis=1))
        emb = embedding_report(densify(P))
        # side view ≈ p·ℓ (+ few arc crossings); tolerance generous,
        # the topology lies in trace+H₁, not in the projection count.
        sane = (3 * 7) <= n_side <= 2 * (3 * 7)
        s3 = emb["embedded"] and sane
        print(f"  (S3) P(7,7,7) embedding: {len(P)} points · side-projection crossings "
              f"{n_side} (≈ p·ℓ={3*7} + arcs) · d_min_rel {emb['d_min_rel']} "
              f"· embedded {'✓' if emb['embedded'] else '✗'} {'✓' if s3 else '✗'}")
    except Exception as e:
        s3 = False
        print(f"  (S3) embedding ERROR: {e}")
    ok += s3

    # S4: knot condition — even ℓ OR even p ⟹ ValueError (no knot)
    total += 1
    s4 = True
    for p, ell in [(2, 7), (3, 4)]:
        try:
            pretzel_embedding(p, ell)
            s4 = False
        except ValueError:
            pass
    print(f"  (S4) knot condition (p,ℓ odd enforced): {'✓' if s4 else '✗'}")
    ok += s4

    print(f"[DEMO] {ok}/{total} smoke checks passed.")
    return 0 if ok == total else 1


def cmd_h1(args):
    g = h1_guarantee(args.p, args.ell)
    print(json.dumps(g, indent=2, ensure_ascii=False))
    return 0


def cmd_build(args):
    out_root = Path(args.out) if args.out else (HERE / "new_knots")
    entry, folder = make_knot(args.p, args.ell, out_root)
    print(json.dumps(entry, indent=2, ensure_ascii=False))
    if folder:
        print(f"→ {folder}")
    return 0


def main():
    ap = argparse.ArgumentParser(description="Machine 247 — split-prime doubler")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo").set_defaults(func=cmd_demo)
    h = sub.add_parser("h1"); h.add_argument("p", type=int); h.add_argument("ell", type=int)
    h.set_defaults(func=cmd_h1)
    b = sub.add_parser("build"); b.add_argument("p", type=int); b.add_argument("ell", type=int)
    b.add_argument("--out", default=None); b.set_defaults(func=cmd_build)
    args = ap.parse_args()
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()
