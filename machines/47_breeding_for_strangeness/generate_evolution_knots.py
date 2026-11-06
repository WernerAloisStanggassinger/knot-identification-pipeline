# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 47 — breeding_for_strangeness (evolution machine).

Genetic algorithm on fseries coefficients; fitness = DISTANCE
from the stock (novelty search!) via cheap proxies: writhe (27 Gauss),
ACN estimator (same double-integral kernel, |·|), embedding margin.
The machine actively searches for what the corpus does NOT have. Mutations =
coefficient noise + the comb effect of machine 32; crossover = the hybrid of machine 31
(xy from A, z from B). Stock = corpus samples + a growing archive.
⚠ RISK NO. 1 (honest, demanded in the founding txt): proxy foreignness
can breed GEOMETRY foreignness instead of TOPOLOGY foreignness — whether „foreign in the
proxy" = „foreign in the knot" is decided only by ignition (pipeline).
Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
import site_config
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "27_pumping_the_twist_account"))
sys.path.insert(0, str(HERE.parent / "31_crossing_the_fan"))
sys.path.insert(0, str(HERE.parent / "32_running_through_effects_pedals"))
sys.path.insert(0, str(HERE.parent))
from generate_writhe_pumped import writhe  # noqa: E402
from generate_hybrid_knots import hybrid  # noqa: E402
from generate_effect_knots import comb  # noqa: E402
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402
from the_knotfamily_machine.fseries_io import reconstruct_from_fourier  # noqa: E402

ROOT = Path(site_config.ROOT)
N_TERMS = 16
N_PROXY = 384          # cheap sampling for the GA loop
PROXY_SCALE = np.array([6.0, 30.0, 0.02])   # grobe Spannen Wr/ACN/Marge


def acn_proxy(P, n=N_PROXY):
    """Average-crossing-number estimator: |Gauss kernel| instead of signed
    (gleicher Doppelintegral-Kern wie 27er-Writhe)."""
    from the_knotfamily_machine.fseries_io import resample_uniform
    Q = resample_uniform(P, n)
    m = (Q + np.roll(Q, -1, axis=0)) / 2.0
    t = np.roll(Q, -1, axis=0) - Q
    d = m[:, None, :] - m[None, :, :]
    r3 = np.linalg.norm(d, axis=2) ** 3
    idx = np.arange(n)
    ring = np.abs(idx[None, :] - idx[:, None])
    ring = np.minimum(ring, n - ring)
    mask = ring > 2
    cross = np.cross(t[:, None, :], t[None, :, :])
    num = np.abs(np.einsum("ijk,ijk->ij", d, cross))
    with np.errstate(divide="ignore", invalid="ignore"):
        val = np.where(mask, num / r3, 0.0)
    return float(val.sum()) / (4 * np.pi)


def proxies(a, b):
    """(Wr, ACN, margin) or None if not embedded (gate first)."""
    P = reconstruct_from_fourier(a, b, n_samples=N_PROXY)
    rep = embedding_report(P, n_samples=N_PROXY)
    if not rep["embedded"]:
        return None
    return np.array([writhe(P, n=N_PROXY), acn_proxy(P),
                     rep["d_min_rel"]])


def novelty(v, bestand):
    """Distance to the nearest stock point, proxy-normalized."""
    if not bestand:
        return np.inf
    D = [np.linalg.norm((v - w) / PROXY_SCALE) for w in bestand]
    return float(min(D))


def load_bestand(max_n=6):
    """Corpus samples: the first embeddable known-table machine curves (cheap, small n)."""
    out, meta = [], []
    for f in sorted((Path(site_config.GENERATED_KNOTS_DIR) / "bekannt").glob(
            "*/*.fseries"))[:60]:
        try:
            C = np.loadtxt(f, comments="%", ndmin=2)
            if C.shape[1] != 6:
                continue
            v = proxies(C[:, [0, 2, 4]], C[:, [1, 3, 5]])
            if v is not None:
                out.append(v)
                meta.append(f.parent.name)
            if len(out) >= max_n:
                break
        except Exception:
            continue
    return out, meta


def random_genome(rng):
    k = np.arange(1, N_TERMS + 1, dtype=float)
    A = k ** -rng.uniform(1.0, 2.0)
    ph = rng.uniform(0, 2 * np.pi, size=(N_TERMS, 3))
    a = A[:, None] * np.cos(ph)
    b = A[:, None] * np.sin(ph)
    a[0], b[0] = (1, 0, 0), (0, 1, 0)     # circle carrier
    return a, b


def mutate(a, b, rng, sigma=0.08):
    if rng.random() < 0.25:               # 32er-Comb-Mutation
        try:
            return comb(a, b, gamma=rng.uniform(0.1, 0.4),
                        tau=rng.uniform(0.3, 2.8))
        except ValueError:
            pass                          # gain gate — back to noise
    amp = np.abs(a) + np.abs(b) + 1e-3
    return (a + sigma * amp * rng.normal(size=a.shape),
            b + sigma * amp * rng.normal(size=b.shape))


def crossover(gA, gB):
    return hybrid(gA[0], gA[1], gB[0], gB[1])


def evolve(seed=47, pop_size=16, gens=6):
    rng = np.random.default_rng(seed)
    bestand, bestand_meta = load_bestand()
    pop = []
    while len(pop) < pop_size:
        g = random_genome(rng)
        v = proxies(*g)
        if v is not None:
            pop.append((g, v, novelty(v, bestand)))
    history = []
    for gen in range(gens):
        pop.sort(key=lambda x: -x[2])
        history.append(round(pop[0][2], 4))
        elite = pop[:max(2, pop_size // 4)]
        nxt = list(elite)
        while len(nxt) < pop_size:
            pa = elite[rng.integers(len(elite))]
            if rng.random() < 0.3 and len(elite) > 1:
                pb = elite[rng.integers(len(elite))]
                g = crossover(pa[0], pb[0])
            else:
                g = mutate(*pa[0], rng)
            v = proxies(*g)
            if v is None:
                continue                  # Embedding-Gate: stirbt still
            nxt.append((g, v, novelty(v, bestand)))
        # archive pressure: the generation's best becomes part of the stock
        bestand.append(nxt[0][1])
        pop = nxt
    pop.sort(key=lambda x: -x[2])
    return pop[0], history, bestand_meta


def make_winner(best, out_root, tag):
    (a, b), v, nov = best
    P = reconstruct_from_fourier(a, b, n_samples=8192)
    name = f"ev_{tag}"
    entry = {"name": name, "machine": "47_breeding_for_strangeness",
             "proxies": {"writhe": round(float(v[0]), 4),
                         "acn_est": round(float(v[1]), 4),
                         "margin": round(float(v[2]), 5)},
             "novelty_final": round(float(nov), 4),
             "predictions": {
                 "risiko_nr_1": "Proxy-Fremdheit kann Geometrie- statt "
                                "Topologie-Fremdheit sein — Zündung "
                                "entscheidet",
                 "identity": "Zünd-Messung"}}
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate (Voll-Auflösung)"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = polyline_to_fseries(P, folder / f"{name}.fseries", knot_name=name,
                              source="maschinen_ideen/47",
                              n_terms=max(64, 2 * N_TERMS),
                              nyquist_max_rel=0.02)
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke-Test Maschine 47 (Evolutions-Maschine)")
    # (1) Proxy-Eichung: Kreis → Wr≈0, ACN≈0; T(2,3) → |Wr| deutlich
    sys.path.insert(0, str(HERE.parent / "10_climbing_the_torus_ladder"))
    from generate_torus_ladder import torus_coeffs
    a_c = np.zeros((4, 3)); b_c = np.zeros((4, 3))
    a_c[0], b_c[0] = (1, 0, 0), (0, 1, 0)
    v_circ = proxies(a_c, b_c)
    a_t, b_t = torus_coeffs(2, 3)
    v_t23 = proxies(a_t, b_t)
    ok1 = (v_circ is not None and abs(v_circ[0]) < 0.02
           and v_circ[1] < 0.05 and v_t23 is not None
           and abs(v_t23[0]) > 1.5)
    print(f"  (1) Proxy-Eichung: Kreis Wr={v_circ[0]:.4f}/ACN={v_circ[1]:.3f}, "
          f"T(2,3) Wr={v_t23[0]:.3f}/ACN={v_t23[1]:.2f}  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) stock loads from the corpus
    bestand, meta = load_bestand()
    ok2 = len(bestand) >= 3
    print(f"  (2) stock: {len(bestand)} corpus proxy points "
          f"({', '.join(m[:18] for m in meta[:3])}…)  "
          f"{'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) Evolution: Novelty-Fortschritt messbar (Werner-Wette)
    best, history, _ = evolve(seed=47, pop_size=12, gens=5)
    ok3 = history[-1] >= history[0]
    print(f"  (3) Novelty trajectory over generations: {history}  "
          f"{'✓ progress' if ok3 else '✗ FAIL (no progress)'}")
    assert ok3
    # (4) winner built + full gates green
    e4, f4 = make_winner(best, HERE / "demo_smoke", "s47_g5")
    ok4 = f4 is not None
    print(f"  (4) breeding winner: Wr={e4['proxies']['writhe']}, "
          f"ACN={e4['proxies']['acn_est']}, novelty={e4['novelty_final']}, "
          f"nyq={e4.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok4 else '✗ FAIL: ' + str(e4.get('rejected'))}")
    assert ok4
    print("[DEMO] 4/4 green — evolution machine ready to fire. "
          "⚠ Risiko Nr. 1 bleibt offen ausgewiesen: fremd-im-Proxy "
          "≠ proven foreign-in-knot (ignition decides).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
