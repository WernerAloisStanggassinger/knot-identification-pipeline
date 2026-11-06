# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 01 — knotting_oscillations (Lissajous/Fourier generator).

Generates knots directly from frequency triples. Core insight: the .fseries IS
this world — cos(n·t+φ) = cosφ·cos(nt) − sinφ·sin(nt) is EXACTLY a
coefficient row. No FFT, no resampling error: the generator writes
the coefficients analytically.

Modes:
  --plan                 parameter space: pairwise coprime triples up to
                         --max-freq, parity classes + bet predictions.
  --make NX NY NZ        one triple, phase grid (--phases k), writes
                         ready-to-fire knot folders + MANIFEST.json.
  --demo                 smoke test (roundtrip exactness, embedding gate,
                         singular-phase defence).

Bets (bit-fixed BEFORE every pipeline run, lit check 2026-06-04 HOLDS):
  all frequencies ODD       ⇒ strongly-plus-amphichiral ⇒ Δ=F², det = odd
                              square ⇒ det ≡ 1 mod 8 [Hartley-Kawauchi 1979]
  one frequency EVEN        ⇒ 2-periodic [Lamm]
  EVERY Lissajous           ⇒ Arf = 0 ⇒ det ≡ 1 mod 8 [BHJS/Murasugi]

Doctrine: do NOT fire before Monday — this script writes only .fseries folders
to out_fseries/, NO pipeline runs, NO corpus ingestion. When firing:
fixed sign convention (Disc #257) + POLYGON_PD=1 witness +
frankenstein_gate --guard-new per record.
"""
from __future__ import annotations
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "the_knotfamily_machine"))
from fseries_io import write_fseries, reconstruct_from_fourier  # noqa: E402

MIN_TERMS = 20          # format convention (saturn): at least 20 rows
EMBED_SAMPLES = 4096    # polygon resolution for the embedding gate
EMBED_MIN_REL = 0.003   # min. non-neighbour distance / diameter (gate)


# ---------------------------------------------------------------- core: exact

def lissajous_coeffs(freqs, phases, amps=(1.0, 1.0, 1.0), n_terms=None):
    """Exact Fourier coefficients of a Lissajous triple.

    x_d(t) = A_d·cos(n_d·t + φ_d)  ⇒  a[n_d-1,d] = A_d·cosφ_d,
                                       b[n_d-1,d] = −A_d·sinφ_d.
    """
    n_terms = max(MIN_TERMS, max(freqs)) if n_terms is None else n_terms
    a = np.zeros((n_terms, 3))
    b = np.zeros((n_terms, 3))
    for d, (n, phi, A) in enumerate(zip(freqs, phases, amps)):
        a[n - 1, d] = A * math.cos(phi)
        b[n - 1, d] = -A * math.sin(phi)
    return a, b


def harmonic_coeffs(spec, n_terms=None):
    """General harmonic mode (the generalisation of the machine).

    spec: list of terms {"axis":0|1|2, "n":int, "amp":float, "phase":float}
    — finitely many harmonics per direction (Fourier (p,q,r) world).
    """
    max_n = max(t["n"] for t in spec)
    n_terms = max(MIN_TERMS, max_n) if n_terms is None else n_terms
    a = np.zeros((n_terms, 3))
    b = np.zeros((n_terms, 3))
    for t in spec:
        d, n, A, phi = t["axis"], t["n"], t.get("amp", 1.0), t.get("phase", 0.0)
        a[n - 1, d] += A * math.cos(phi)
        b[n - 1, d] += -A * math.sin(phi)
    return a, b


# ------------------------------------------------------------ Embedding-Gate

def min_nonadjacent_dist(P, exclude=8):
    """Minimal distance of non-adjacent polygon vertices (cyclic), chunked.

    Returns (d_min, i, j). exclude = arc neighbourhood that does not count.
    """
    n = len(P)
    d_min, arg = np.inf, (0, 0)
    chunk = 256
    idx = np.arange(n)
    for s in range(0, n, chunk):
        block = P[s:s + chunk]                       # (c,3)
        D = np.linalg.norm(block[:, None, :] - P[None, :, :], axis=2)  # (c,n)
        ring = np.abs(idx[None, :] - idx[s:s + chunk, None])
        ring = np.minimum(ring, n - ring)
        D[ring <= exclude] = np.inf
        k = np.unravel_index(np.argmin(D), D.shape)
        if D[k] < d_min:
            d_min, arg = float(D[k]), (s + int(k[0]), int(k[1]))
    return d_min, arg


def embedding_report(a, b, n_samples=EMBED_SAMPLES):
    """Gate: is the curve embedded? Returns a dict with d_min_rel + verdict."""
    P = reconstruct_from_fourier(a, b, n_samples=n_samples)
    diam = float(np.linalg.norm(P.max(axis=0) - P.min(axis=0)))
    d_min, (i, j) = min_nonadjacent_dist(P)
    rel = d_min / diam if diam > 0 else 0.0
    return {
        "n_samples": n_samples,
        "diameter": round(diam, 6),
        "d_min": round(d_min, 6),
        "d_min_rel": round(rel, 6),
        "closest_pair": [int(i), int(j)],
        "embedded": bool(rel >= EMBED_MIN_REL),
    }


# ------------------------------------------------------- parameter space / plan

def coprime_triples(max_freq):
    """Pairwise coprime triples 2 ≤ nx < ny < nz ≤ max_freq."""
    out = []
    for nx in range(2, max_freq + 1):
        for ny in range(nx + 1, max_freq + 1):
            if math.gcd(nx, ny) != 1:
                continue
            for nz in range(ny + 1, max_freq + 1):
                if math.gcd(nx, nz) == 1 and math.gcd(ny, nz) == 1:
                    out.append((nx, ny, nz))
    return out


def parity_class(freqs):
    n_even = sum(1 for f in freqs if f % 2 == 0)
    return "all_odd" if n_even == 0 else f"{n_even}_even"


def predictions(freqs):
    """The bets — bit-fixed before every measurement (lit check 06-04)."""
    p = {
        "arf": 0,
        "det_mod_8": 1,
        "source": "Lit-Check 2026-06-04 (BHJS/Lamm/Hartley-Kawauchi/Murasugi)",
    }
    if parity_class(freqs) == "all_odd":
        p.update({
            "strongly_plus_amphichiral": True,
            "det_is_odd_square": True,
            "alexander_is_perfect_square": True,
        })
    else:
        p.update({"period_2": True})
    return p


def plan(max_freq):
    triples = coprime_triples(max_freq)
    rows = []
    for t in triples:
        rows.append({
            "freqs": list(t),
            "parity": parity_class(t),
            "predictions": predictions(t),
        })
    n_odd = sum(1 for r in rows if r["parity"] == "all_odd")
    return {
        "max_freq": max_freq,
        "n_triples": len(rows),
        "n_all_odd": n_odd,
        "n_one_even": len(rows) - n_odd,
        "triples": rows,
    }


# ----------------------------------------------------------------- generation

def make_knot(freqs, phases, out_root, amps=(1.0, 1.0, 1.0)):
    """One Lissajous knot: ready-to-fire folder + manifest entry.

    Returns (manifest_entry, fseries_path|None). Rejected (not
    embedded) phases produce NO folder, only the manifest entry.
    """
    nx, ny, nz = freqs
    name = (f"lj_{nx}_{ny}_{nz}"
            f"_ph{phases[0]:.3f}_{phases[1]:.3f}_{phases[2]:.3f}"
            ).replace(".", "p")
    a, b = lissajous_coeffs(freqs, phases, amps)
    emb = embedding_report(a, b)
    entry = {
        "name": name,
        "machine": "01_knotting_oscillations",
        "mode": "lissajous",
        "freqs": list(freqs),
        "phases": [round(p, 6) for p in phases],
        "amps": list(amps),
        "parity": parity_class(freqs),
        "predictions": predictions(freqs),
        "embedding": emb,
        "coefficients_exact": True,
    }
    if not emb["embedded"]:
        entry["rejected"] = "singular/near-singular phase (embedding gate)"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    path = write_fseries(folder / f"{name}.fseries", a, b,
                         knot_name=name,
                         source="maschinen_ideen/01_knotting_oscillations")
    # Roundtrip truth: read back, compare coefficients bit-wise (6 dec.)
    back = np.loadtxt(path, comments="%")
    a2 = back[:, [0, 2, 4]]
    b2 = back[:, [1, 3, 5]]
    rt = float(max(np.abs(a2 - a).max(), np.abs(b2 - b).max()))
    entry["roundtrip_max_err"] = rt
    if rt > 5e-7:
        raise RuntimeError(f"{name}: roundtrip error {rt} > format resolution")
    return entry, path


def phase_raster(k):
    """k² phase pairs (φx, φy) in the open grid, φz = 0 (t-shift convention)."""
    vals = [(i + 1) * math.pi / (k + 1) for i in range(k)]
    return [(px, py, 0.0) for px in vals for py in vals]


# ----------------------------------------------------------------------- CLI

def cmd_plan(args):
    P = plan(args.max_freq)
    out = HERE / f"parameter_plan_maxfreq{args.max_freq}.json"
    out.write_text(json.dumps(P, indent=1))
    print(f"[PLAN] {P['n_triples']} pairwise coprime triples ≤ {args.max_freq}: "
          f"{P['n_all_odd']} all_odd (amphi bet) + {P['n_one_even']} with an even "
          f"frequency (2-periodic bet) → {out.name}")
    for r in P["triples"][:12]:
        print("   ", r["freqs"], r["parity"])
    if P["n_triples"] > 12:
        print(f"    … +{P['n_triples']-12} more in the JSON")


def cmd_make(args):
    out_root = HERE / "out_fseries"
    freqs = (args.freqs[0], args.freqs[1], args.freqs[2])
    entries, made = [], 0
    for ph in phase_raster(args.phases):
        e, p = make_knot(freqs, ph, out_root)
        entries.append(e)
        made += p is not None
    man_path = HERE / "MANIFEST.json"
    manifest = json.loads(man_path.read_text()) if man_path.exists() else {"knots": []}
    known = {k["name"] for k in manifest["knots"]}
    manifest["knots"] += [e for e in entries if e["name"] not in known]
    man_path.write_text(json.dumps(manifest, indent=1))
    rej = len(entries) - made
    print(f"[MAKE] {freqs} {parity_class(freqs)}: {made} knots written, "
          f"{rej} phases rejected by the embedding gate → out_fseries/ + MANIFEST.json")


def cmd_demo(args):
    print("[DEMO] Smoke test machine 01")
    # 1) Exactness: coefficients analytic vs. FFT of the reconstructed curve
    a, b = lissajous_coeffs((3, 4, 7), (0.7, 0.2, 0.0))
    P = reconstruct_from_fourier(a, b, n_samples=2048)
    spec = np.fft.rfft(P, axis=0)
    a_fft = (2.0 / 2048) * np.real(spec[1:MIN_TERMS + 1])
    b_fft = (-2.0 / 2048) * np.imag(spec[1:MIN_TERMS + 1])
    err = max(np.abs(a_fft - a).max(), np.abs(b_fft - b).max())
    print(f"  (1) analytic == FFT re-measurement: max|Δ| = {err:.2e}  "
          f"{'✓' if err < 1e-12 else '✗ FAIL'}")
    assert err < 1e-12
    # 2) Embedding gate: generic phase embedded
    emb = embedding_report(a, b)
    print(f"  (2) (3,4,7) φ=(0.7,0.2,0): d_min_rel={emb['d_min_rel']}  "
          f"embedded={emb['embedded']}  {'✓' if emb['embedded'] else '✗ FAIL'}")
    assert emb["embedded"]
    # 3) Singularity guard: φ=(0,0,0) is point-symmetry singular (curve doubling)
    a0, b0 = lissajous_coeffs((3, 5, 7), (0.0, 0.0, 0.0))
    emb0 = embedding_report(a0, b0)
    print(f"  (3) (3,5,7) φ=0 (singular): d_min_rel={emb0['d_min_rel']}  "
          f"embedded={emb0['embedded']}  {'✓' if not emb0['embedded'] else '✗ FAIL'}")
    assert not emb0["embedded"]
    # 4) Write + roundtrip (in a tmp subfolder, no manifest entry)
    e, p = make_knot((3, 4, 7), (0.7, 0.2, 0.0), HERE / "demo_smoke")
    print(f"  (4) write+readback: roundtrip_max_err={e['roundtrip_max_err']:.1e} "
          f"→ {p.relative_to(HERE)}  ✓")
    print("[DEMO] 4/4 green — machine ready to fire (firing only after Monday).")


def main():
    ap = argparse.ArgumentParser(description="Machine 01 — knotting_oscillations")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("plan");  p1.add_argument("--max-freq", type=int, default=12)
    p2 = sub.add_parser("make")
    p2.add_argument("freqs", type=int, nargs=3, metavar=("NX", "NY", "NZ"))
    p2.add_argument("--phases", type=int, default=3, help="k → k² phase grid")
    sub.add_parser("demo")
    args = ap.parse_args()
    {"plan": cmd_plan, "make": cmd_make, "demo": cmd_demo}[args.cmd](args)


if __name__ == "__main__":
    main()
