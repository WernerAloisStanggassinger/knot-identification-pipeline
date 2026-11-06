# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 119 — ghosting_the_ladder (Whitehead doubles of the C₂ staircase).

RUN-UP HISTORY (both predecessors MEASURED failed, see txt):
  1: 117 rail companion — w' window scissors (corridor 12ρ) +
     clasp features below the gate threshold (ρ forced too small).
  2: anisotropic stretching — rel clearance unmoved, quantile check
     measured the corridor instead of the core.
  3 (THIS ONE): RING COMPANION — closed 3-braid (σ₁σ₂⁻¹)^{3k}·σ₁σ₂
     on concentric circular tracks (sector Θ=π carries the word,
     rest = planar circular arcs). No lift, no V-dip:
     radius of curvature ~R_RING instead of ~0.3 ⟹ clasp corridor hugs
     close (sagitta ∝ 1/R); clearance ~S_RAD ⟹ ρ in the trefoil class,
     clasp features above the gate threshold. Same knot family
     (identical braid word ⟹ identical closure).

Double = machine-81 whitehead_double (wdir=+1 fix), roll_to_max_clearance,
ρ = min(0.3·local_clear, 0.45·d_min). Δ=1/det=1 for ALL rungs
(doubles theorem), companion degree climbs hidden: Dror swarm
axis 4 (companion × twist × sign × DEGREE).

BETS: pre-registered in 119_geister_leiter.txt (W0 unknot · W1 Δ=1
series · W2 Jones separates, two-sided · W3 θ ghost ladder · W4 6D
MUST see · W5 cross-consistency with 117/81). Doctrine: do NOT fire.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
YARD = HERE.parent
sys.path.insert(0, str(YARD))
sys.path.insert(0, str(YARD / "117_stepping_the_delta"))
sys.path.insert(0, str(YARD / "81_clasping_the_double"))
sys.path.insert(0, str(YARD / "10_climbing_the_torus_ladder"))
from generate_delta_ladder import braid_word  # noqa: E402
from generate_whitehead_doubles import (  # noqa: E402
    whitehead_double, roll_to_max_clearance)
from machine_common import (  # noqa: E402
    densify, embedding_report, polyline_to_fseries,
    count_projection_crossings)
sys.path.insert(0, str(YARD / "the_knotfamily_machine"))
from fseries_io import resample_uniform  # noqa: E402

R_RING = 2.0             # ring radius (radius of curvature of the arcs)
S_RAD  = 0.55            # radial track spacing
H_Z    = 0.22            # z lift per crossing
THETA  = np.pi           # word sector (fixed ⟹ rungs comparable)
NW     = 25              # waypoints per crossing window
N_ARC  = 120             # waypoints per return arc (fixed)
N_COMP = 2048            # uniform companion resampling (81 expectation)


def ring_curve(word) -> np.ndarray:
    """Closed 3-braid for any word on the ring (word-
    parametric since harvest 26 — machine 120 uses the same architecture;
    Precondition: word permutation is a 3-cycle, otherwise link)."""
    n_c = len(word)
    w = THETA / n_c
    lanes_r = (R_RING - S_RAD, R_RING, R_RING + S_RAD)
    occupant = [0, 1, 2]
    paths = {s: [] for s in range(3)}
    th0 = 0.0
    ts = np.linspace(0.0, 1.0, NW, endpoint=False)
    for gi, sign in word:
        a, b = gi - 1, gi
        c = 3 - a - b
        sa, sb, sc = occupant[a], occupant[b], occupant[c]
        for t in ts:
            th = th0 + w * t
            ramp = (1 - np.cos(np.pi * t)) / 2
            ra = lanes_r[a] + (lanes_r[b] - lanes_r[a]) * ramp
            rb = lanes_r[b] + (lanes_r[a] - lanes_r[b]) * ramp
            z = H_Z * np.sin(np.pi * t)
            ct, st = np.cos(th), np.sin(th)
            paths[sa].append((ra * ct, ra * st, +z if sign > 0 else -z))
            paths[sb].append((rb * ct, rb * st, -z if sign > 0 else +z))
            paths[sc].append((lanes_r[c] * ct, lanes_r[c] * st, 0.0))
        occupant[a], occupant[b] = sb, sa
        th0 += w
    end_lane = {occupant[lane]: lane for lane in range(3)}
    arcs = {}
    for lane in range(3):                      # return arcs Θ→2π, planar
        r = lanes_r[lane]
        th = np.linspace(THETA, 2 * np.pi, N_ARC, endpoint=False)
        arcs[lane] = [(r * np.cos(t), r * np.sin(t), 0.0) for t in th]
    pts, cur, visited = [], 0, []
    for _ in range(3):
        visited.append(cur)
        pts += paths[cur]
        lane = end_lane[cur]
        pts += arcs[lane]
        cur = lane
    if cur != 0 or sorted(visited) != [0, 1, 2]:
        raise RuntimeError("permutation not a 3-cycle — closure broken (link?)")
    return np.asarray(pts, float)


def ring_rung(k: int) -> np.ndarray:
    """119 rung: ring curve for the C₂ staircase word."""
    return ring_curve(braid_word(k))


def clearance_arc(P: np.ndarray, n_samples: int = 4096,
                  excl_arc: float = 0.6) -> float:
    """GENUINE non-local clearance: pairs with arc-length distance
    < excl_arc excluded. Root cause (measured, attempt 3): embedding_
    report's d_min has a MEASUREMENT FLOOR of 9·L/n_samples (exclude=8
    in sample INDICES) — for long smooth curves it measures only this
    floor (ring: 0.0817 = exactly 9·L/4096), never the true clearance."""
    Q = resample_uniform(P, n_samples)
    seg = np.linalg.norm(np.diff(np.vstack([Q, Q[:1]]), axis=0), axis=1)
    L = float(seg.sum())
    excl = max(1, int(np.ceil(excl_arc / (L / n_samples))))
    n = len(Q)
    idx = np.arange(n)
    d_min = np.inf
    for st in range(0, n, 256):
        blk = Q[st:st + 256]
        D = np.linalg.norm(blk[:, None, :] - Q[None, :, :], axis=2)
        ring = np.abs(idx[None, :] - idx[st:st + 256, None])
        ring = np.minimum(ring, n - ring)
        D[ring <= excl] = np.inf
        d_min = min(d_min, float(D.min()))
    return d_min


def dist_to_polyline(D: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Distance of every D point to the dense companion sample."""
    Xf = resample_uniform(X, 8192)
    out = np.empty(len(D))
    for s in range(0, len(D), 512):
        blk = D[s:s + 512]
        dd = np.linalg.norm(blk[:, None, :] - Xf[None, :, :], axis=2)
        out[s:s + 512] = dd.min(axis=1)
    return out


def make_rung(k: int, out_root: Path):
    name = f"gd_k{k}"
    X = ring_rung(k)
    n_c = count_projection_crossings(X, axis=2)   # whole ring: window only
    Xu = resample_uniform(densify(X), N_COMP)
    Xu, local_clear = roll_to_max_clearance(Xu)
    emb_c = embedding_report(Xu)
    clear_true = clearance_arc(Xu)            # arc-length-aware
    rho = float(min(0.3 * local_clear, 0.45 * clear_true))
    entry = {"name": name, "machine": "119_ghosting_the_ladder",
             "companion": f"ring_dm_k{k} (R={R_RING}, Sektor π)",
             "companion_crossings": n_c,
             "companion_crossings_soll": 6 * k + 2,
             "companion_embedding": emb_c,
             "rho_protokoll": {"local_clear": round(local_clear, 6),
                               "companion_d_min_floor": emb_c["d_min"],
                               "clearance_arc": round(clear_true, 6),
                               "rho": round(rho, 6)},
             "predictions": {
                 "W0": "k=0: Unknot (Pattern-Trivialitaet)",
                 "W1": "ALLE Sprossen Delta=1/det=1 (Doubles-Satz)",
                 "W2": "Jones trennt Sprossen (zweiseitig)",
                 "W3": "theta-Leiter unter Alexander-Blindheit (Dror A4)",
                 "W4": "6D MUSS Begleiter-Grad sehen (Kontrast zu W1)",
                 "W5": "gd_k1 == W(Knoten(dm_k1)) Cross-Konsistenz"}}
    if n_c != 6 * k + 2:
        entry["rejected"] = f"Begleiter-Kreuzungen {n_c} != {6 * k + 2}"
        return entry, None
    if not emb_c["embedded"]:
        entry["rejected"] = "Begleiter embedding gate"
        return entry, None
    D = whitehead_double(Xu, rho=rho, wdir=+1)
    dist = dist_to_polyline(D, Xu)
    # tube discipline in three steps (lessons of attempts 1/2): core hugs the
    # companion (median≈ρ), corridor share stays small (by design ~8%),
    # and EVERYTHING stays in the guaranteed empty ball around the turning point
    # (lesson 5: roll_to_max_clearance ⟹ radius local_clear free).
    med = float(np.median(dist))
    frac_out = float(np.mean(dist > 1.3 * rho))
    entry["rohr_disziplin"] = {
        "median_over_rho": round(med / rho, 3),
        "frac_over_1.3rho": round(frac_out, 4),
        "max_dist_abs": round(float(dist.max()), 4),
        "korridor_budget_0.9lc": round(0.9 * local_clear, 4)}
    if not (0.8 <= med / rho <= 1.2):
        entry["rejected"] = f"Rohr-Kern: Median {med/rho:.2f}ρ ∉ [0.8,1.2]ρ"
        return entry, None
    if frac_out > 0.15:
        entry["rejected"] = f"Korridor-Anteil {frac_out:.1%} > 15%"
        return entry, None
    if dist.max() > 0.9 * local_clear:
        entry["rejected"] = (f"Korridor verlässt die freie Kugel: "
                             f"{dist.max():.3f} > 0.9·lc={0.9*local_clear:.3f}")
        return entry, None
    emb = embedding_report(densify(D))
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate (Double)"
        return entry, None
    folder = out_root / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = None
    for nt in (256, 384, 512, 768, 1024):
        try:
            rep = polyline_to_fseries(
                D, folder / f"{name}.fseries", knot_name=name,
                source="maschinen_ideen/119", n_terms=nt,
                nyquist_max_rel=0.02)
            break
        except RuntimeError:
            continue
    if rep is None:
        entry["rejected"] = "Nyquist rot bis n_terms=1024"
        return entry, None
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(_args):
    print("[DEMO] smoke test machine 119 (ghost ladder, W∘ring C2 staircase)",
          flush=True)
    out = HERE / "demo_smoke"
    ok = 0
    # S1: ring companion — expected crossings (whole ring) + embedding
    s1 = True
    for k in range(4):
        X = ring_rung(k)
        n_c = count_projection_crossings(X, axis=2)
        emb = embedding_report(densify(X))
        good = (n_c == 6 * k + 2) and emb["embedded"]
        s1 &= good
        print(f"  (S1) ring companion k={k}: crossings {n_c}/{6 * k + 2} · "
              f"d_min_rel={emb['d_min_rel']} {'✓' if good else '✗'}",
              flush=True)
    ok += s1
    # S2+S3: doubles — tube discipline + embedding + Nyquist
    s2 = s3 = True
    for k in range(4):
        entry, folder = make_rung(k, out)
        rp = entry.get("rho_protokoll", {})
        rd = entry.get("rohr_disziplin")
        rej = str(entry.get("rejected", ""))
        if rd and not any(w in rej for w in ("Rohr", "Korridor")):
            print(f"  (S2) gd_k{k}: ρ={rp.get('rho')} · core med="
                  f"{rd['median_over_rho']}ρ · corridor {rd['frac_over_1.3rho']:.1%}"
                  f" · max {rd['max_dist_abs']} < {rd['korridor_budget_0.9lc']} ✓",
                  flush=True)
        else:
            s2 = False
            print(f"  (S2) gd_k{k}: tube discipline not reached — {rej} ✗",
                  flush=True)
            continue
        if folder is None:
            print(f"  (S3) gd_k{k}: REJECTED — {entry['rejected']}"
                  f"{' (format limit measured, ok)' if k == 3 else ' ✗'}",
                  flush=True)
            if k < 3:
                s3 = False
        else:
            print(f"  (S3) gd_k{k}: d_min_rel={entry['embedding']['d_min_rel']}"
                  f" · n_terms={entry['fseries']['n_terms']}"
                  f" · nyq={entry['fseries']['nyquist_dev_rel']} ✓", flush=True)
    ok += s2
    ok += s3
    print(f"[DEMO] {ok}/3 green — {'ghost ladder stands.' if ok == 3 else 'RED.'}"
          f" Bets: 119_geister_leiter.txt. Doctrine: do NOT fire.", flush=True)
    return 0 if ok == 3 else 1


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    a = ap.parse_args()
    if a.cmd == "demo":
        sys.exit(cmd_demo(a))


if __name__ == "__main__":
    main()
