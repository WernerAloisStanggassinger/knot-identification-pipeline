# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 49 — bending_space_conformally (Möbius/inversion zoo).

Sphere inversion P → c + r²(P−c)/|P−c|² on stock knots: knot type
stays (inversion = homeomorphism of R³∪∞, curve stays finite
as long as c outside), geometry becomes WILDLY different — the perfect
embedding zoo for the 6D forgetting question (27 sister, conformal
instead of pumped) + sculpturally independent forms.
BET EXACT: center outside ⇒ same knot (up to mirror
convention) ⇒ SONO normal form MUST collapse onto the same point.
FLIP-OVER CLAUSE: center near the curve inflates outliers → center
distance rule + embedding/Nyquist gates. Doctrine: do NOT ignite.
"""
from __future__ import annotations
import site_config
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402
from the_knotfamily_machine.fseries_io import reconstruct_from_fourier  # noqa: E402

ROOT = Path(site_config.ROOT)
CENTER_GAP_MIN = 0.25   # min |P−c| / diameter — center distance rule
N_SAMPLES = 8192


def load_fseries(path):
    """Robust: skip %-comments, 6 columns → a(:,3), b(:,3)."""
    C = np.loadtxt(path, comments="%", ndmin=2)
    if C.shape[1] != 6:
        raise ValueError(f"{path}: {C.shape[1]} columns (expected 6)")
    return C[:, [0, 2, 4]], C[:, [1, 3, 5]]


def invert(P, center, radius=None):
    """Sphere inversion; radius default = RMS distance (keeps scale civil)."""
    d = P - np.asarray(center, float)
    n2 = (d ** 2).sum(axis=1)
    if radius is None:
        radius = float(np.sqrt(np.sqrt(n2).mean() ** 2))
    return np.asarray(center, float) + (radius ** 2) * d / n2[:, None], radius


def center_gap(P, center):
    diam = float(np.linalg.norm(P.max(0) - P.min(0)))
    return float(np.linalg.norm(P - np.asarray(center, float),
                                axis=1).min() / diam)


def pick_source():
    """Stock source: first readable bekannt machine whose RAW curve
    itself passes the embedding gate (build lesson: many machine curves
    lie just BELOW the gate margin — then the inversion is not the
    culprit, but the source). Examples vary — deliberately NOT always
    the same anchor knot."""
    for pat in ["bekannt/11n_133__kf100_sh+0.5_tw+2/*.fseries",
                "bekannt/11n_133__kf58_helix_twist_t3.0/*.fseries",
                "bekannt/*/*.fseries"]:
        hits = sorted(Path(site_config.GENERATED_KNOTS_DIR).glob(pat))
        for h in hits:
            try:
                a, b = load_fseries(h)
                P = reconstruct_from_fourier(a, b, n_samples=2048)
                if embedding_report(P, n_samples=2048)["embedded"]:
                    return h
            except Exception:
                continue
    raise FileNotFoundError("no gate-worthy bekannt .fseries found")


def make_knot(src_path, center_rel, out_root, tag):
    a, b = load_fseries(src_path)
    P = reconstruct_from_fourier(a, b, n_samples=N_SAMPLES)
    diam = float(np.linalg.norm(P.max(0) - P.min(0)))
    center = P.mean(axis=0) + np.asarray(center_rel, float) * diam
    gap = center_gap(P, center)
    name = f"inv_{Path(src_path).parent.name[:24]}_{tag}"
    entry = {"name": name, "machine": "49_bending_space_conformally",
             "source": str(src_path), "center_rel": list(map(float, center_rel)),
             "center_gap_rel": round(gap, 4),
             "predictions": {
                 "wette_exakt": "Zentrum außerhalb ⇒ gleicher Knoten — "
                                "SONO-Normalform muss kollabieren (Zündung)",
                 "zweck": "Einbettungs-Zoo für 6D-Vergessens-Frage"}}
    if gap < CENTER_GAP_MIN:
        entry["rejected"] = f"Zentrum-Gate (gap={gap:.3f} < {CENTER_GAP_MIN})"
        return entry, None
    Q, radius = invert(P, center)
    entry["inv_radius"] = round(radius, 6)
    emb = embedding_report(Q)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    # inversion generates high harmonics → terms generous
    rep = polyline_to_fseries(Q, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/49",
                              n_terms=max(160, 2 * a.shape[0]),
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 49 (inversion zoo)")
    rng = np.random.default_rng(5)
    # (1) involution exact: invert∘invert = id (same c, r)
    P = rng.normal(size=(2048, 3)) + np.array([4.0, 0, 0])
    c = np.zeros(3)
    Q, r = invert(P, c, radius=2.0)
    P2, _ = invert(Q, c, radius=2.0)
    err = np.abs(P2 - P).max()
    ok1 = err < 1e-12
    print(f"  (1) Involution: max|Δ|={err:.2e}  {'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) conformal anchor: circle → exact circle (radius residual)
    t = np.linspace(0, 2 * np.pi, 4096, endpoint=False)
    circ = np.column_stack([3 + np.cos(t), np.sin(t), 0 * t])
    Qc, _ = invert(circ, np.zeros(3), radius=1.5)
    ctr = Qc.mean(axis=0)
    # circle centre via least squares (mean is not enough — inversion
    # distributes points unevenly): linearize |q−m|²=R²
    A = np.column_stack([2 * Qc, np.ones(len(Qc))])
    sol, *_ = np.linalg.lstsq(A, (Qc ** 2).sum(axis=1), rcond=None)
    m = sol[:3]
    R = np.linalg.norm(Qc - m, axis=1)
    resid = float(R.max() - R.min()) / float(R.mean())
    ok2 = resid < 1e-12
    print(f"  (2) circle→circle: radius residual={resid:.2e}  "
          f"{'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) stock knots inverted + gates green
    src = pick_source()
    e3, f3 = make_knot(src, (0.9, 0.35, 0.55), HERE / "demo_smoke", "far")
    ok3 = f3 is not None
    print(f"  (3) zoo specimen from {Path(src).parent.name}: "
          f"gap={e3['center_gap_rel']}, "
          f"nyq={e3.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok3 else '✗ FAIL: ' + str(e3.get('rejected'))}")
    assert ok3
    # (4) center gate fires for a near center
    e4, f4 = make_knot(src, (0.05, 0.02, 0.0), HERE / "demo_smoke", "near")
    ok4 = f4 is None and "Zentrum-Gate" in str(e4.get("rejected", ""))
    print(f"  (4) center rule: {e4.get('rejected')}  "
          f"{'✓' if ok4 else '✗ FAIL'}")
    assert ok4
    print("[DEMO] 4/4 green — inversion machine ready to fire "
          "(embedding zoo for the 6D forgetting question).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
