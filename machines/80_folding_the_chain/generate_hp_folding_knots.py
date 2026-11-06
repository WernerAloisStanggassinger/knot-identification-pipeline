# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 80 — folding_the_chain (HP lattice folding).

Self-avoiding lattice chain + HP energy (hydrophobic/polar):
Monte Carlo folding (pivot moves) COMPACTIFIED — then close
and ask for the knot. THE QUESTION: does compactness raise the knot rate?
(protein question; 33 measured pure chance, here the same statistics
ENERGY-DRIVEN.) BETS today: the energy falls (Metropolis statistics
measured), the radius of gyration shrinks measurably, writhe-suspicion rates
compact vs random = the new number. REVERSAL CLAUSE: closing path
around the outside + embedding gate (do not pierce the folding).
Doctrine: do NOT ignite.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "27_pumping_the_twist_account"))
sys.path.insert(0, str(HERE.parent))
from generate_writhe_pumped import writhe  # noqa: E402
from machine_common import embedding_report, polyline_to_fseries  # noqa: E402

DIRS = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0],
                 [0, 0, 1], [0, 0, -1]])
ROTS = []
for ax in range(3):
    for sgn in (1, -1):
        for quarter in (1, 2, 3):
            th = sgn * quarter * np.pi / 2
            c, s = int(round(np.cos(th))), int(round(np.sin(th)))
            M = np.eye(3, dtype=int)
            i, j = [(1, 2), (0, 2), (0, 1)][ax]
            M[i, i], M[i, j], M[j, i], M[j, j] = c, -s, s, c
            ROTS.append(M)


def random_walk_chain(n, rng):
    """Self-avoiding chain (growth with restart)."""
    while True:
        pos = [np.zeros(3, int)]
        occ = {(0, 0, 0)}
        ok = True
        for _ in range(n - 1):
            rng.shuffle(DIRS)   # was np.random.shuffle (global, unseeded RNG): demo was not reproducible (6 Nov 2026)
            for d in DIRS[rng.permutation(6)]:
                q = pos[-1] + d
                if tuple(q) not in occ:
                    pos.append(q)
                    occ.add(tuple(q))
                    break
            else:
                ok = False
                break
        if ok:
            return np.asarray(pos)


def hp_energy(chain, seq):
    """−1 per H-H contact (lattice neighbours, not chain neighbours)."""
    occ = {tuple(p): i for i, p in enumerate(chain)}
    E = 0
    for i, p in enumerate(chain):
        if seq[i] != 1:
            continue
        for d in DIRS:
            j = occ.get(tuple(p + d))
            if j is not None and j > i + 1 and seq[j] == 1:
                E -= 1
    return E


def pivot_fold(chain, seq, n_sweeps, beta, rng):
    """Metropolis with pivot moves; returns (chain, energy trace, acceptance)."""
    C = chain.copy()
    E = hp_energy(C, seq)
    E_hist, acc = [E], 0
    n = len(C)
    for _ in range(n_sweeps):
        k = int(rng.integers(1, n - 1))
        M = ROTS[int(rng.integers(len(ROTS)))]
        new = C.copy()
        new[k:] = C[k] + (C[k:] - C[k]) @ M.T
        if len({tuple(p) for p in new}) < n:
            E_hist.append(E)
            continue
        E_new = hp_energy(new, seq)
        if E_new <= E or rng.random() < np.exp(-beta * (E_new - E)):
            C, E = new, E_new
            acc += 1
        E_hist.append(E)
    return C, E_hist, acc / n_sweeps


def gyration(C):
    X = C - C.mean(0)
    return float(np.sqrt((X ** 2).sum(axis=1).mean()))


def close_chain(C):
    """Closing path outside: end → under the bounding box → start."""
    lo = C[:, 2].min() - 2.5
    end, start = C[-1].astype(float), C[0].astype(float)
    way = [end + [0, 0, lo - end[2]],
           [start[0], end[1], lo],
           [start[0], start[1], lo]]
    return np.vstack([C.astype(float)]
                     + [np.asarray(w, float)[None, :] for w in way])


def rounded(P, pts=20, rounds=2, jitter=0.04, seed=0):
    rng = np.random.default_rng(seed)
    Q = P + jitter * rng.normal(size=P.shape)
    dense = []
    m = len(Q)
    for k in range(m):
        a, b = Q[k], Q[(k + 1) % m]
        t = np.linspace(0, 1, pts, endpoint=False)
        dense.append(np.outer(1 - t, a) + np.outer(t, b))
    C = np.vstack(dense)
    for _ in range(rounds):
        C = 0.5 * C + 0.25 * (np.roll(C, 1, axis=0) + np.roll(C, -1, axis=0))
    return C


def make_knot(n, seed, n_sweeps, beta, out_root):
    rng = np.random.default_rng(seed)
    seq = (rng.random(n) < 0.6).astype(int)        # 60% hydrophob
    C0 = random_walk_chain(n, rng)
    C1, E_hist, acc = pivot_fold(C0, seq, n_sweeps, beta, rng)
    name = f"hp_n{n}_s{seed}"
    P = rounded(close_chain(C1), seed=seed)
    entry = {"name": name, "machine": "80_folding_the_chain",
             "n_monomers": n, "seed": seed,
             "E0": int(E_hist[0]), "E_final": int(E_hist[-1]),
             "acceptance": round(acc, 3),
             "rg_before": round(gyration(C0), 3),
             "rg_after": round(gyration(C1), 3),
             "writhe": round(float(writhe(P[::4], n=256)), 3),
             "predictions": {
                 "frage": "Knoten-Rate kompakt vs zufällig (Zündung)",
                 "identity": "Zünd-Messung"}}
    emb = embedding_report(P)
    entry["embedding"] = emb
    if not emb["embedded"]:
        entry["rejected"] = "embedding gate"
        return entry, None
    folder = Path(out_root) / name
    folder.mkdir(parents=True, exist_ok=True)
    rep = None
    for nt in (128, 192, 256):
        try:
            rep = polyline_to_fseries(P, folder / f"{name}.fseries",
                                      knot_name=name,
                                      source="maschinen_ideen/80",
                                      n_terms=nt, nyquist_max_rel=0.02)
            break
        except RuntimeError:
            continue
    if rep is None:
        entry["rejected"] = "Nyquist rot"
        return entry, None
    entry["fseries"] = rep
    (folder / f"{name}_entry.json").write_text(json.dumps(entry, indent=1))
    return entry, folder


def cmd_demo(args):
    print("[DEMO] Smoke test machine 80 (folding the chain)")
    rng = np.random.default_rng(80)
    # (1) Faltung senkt Energie + Kompaktheit (3 Ketten gemessen)
    rows, ok1 = [], True
    for s in range(3):
        seq = (rng.random(48) < 0.6).astype(int)
        C0 = random_walk_chain(48, rng)
        C1, E_hist, acc = pivot_fold(C0, seq, 1200, beta=2.0, rng=rng)
        rows.append((E_hist[0], E_hist[-1], round(gyration(C0), 2),
                     round(gyration(C1), 2)))
        ok1 = ok1 and E_hist[-1] < E_hist[0] and gyration(C1) < gyration(C0)
    print(f"  (1) (E₀→E, Rg₀→Rg): {rows}  "
          f"{'✓ faltet + kompaktiert' if ok1 else '✗ FAIL'}")
    assert ok1
    # (2) Selbstvermeidung bleibt bit-erhalten (Pivot verwirft Kollision)
    seq = (rng.random(40) < 0.6).astype(int)
    C0 = random_walk_chain(40, rng)
    C1, _, _ = pivot_fold(C0, seq, 800, 2.0, rng)
    ok2 = len({tuple(p) for p in C1}) == len(C1)
    print(f"  (2) Selbstvermeidung nach Faltung: {len(C1)} Zellen "
          f"einmalig  {'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # (3) THE NEW NUMBER: |Wr| suspect rate compact vs random
    wr_rand, wr_fold = [], []
    for s in range(12):
        rngs = np.random.default_rng(1000 + s)
        seq = (rngs.random(48) < 0.6).astype(int)
        C0 = random_walk_chain(48, rngs)
        P0 = rounded(close_chain(C0), seed=s)
        wr_rand.append(abs(writhe(P0[::4], n=192)))
        C1, _, _ = pivot_fold(C0, seq, 900, 2.0, rngs)
        P1 = rounded(close_chain(C1), seed=s)
        wr_fold.append(abs(writhe(P1[::4], n=192)))
    r_rand = sum(1 for w in wr_rand if w > 1.5)
    r_fold = sum(1 for w in wr_fold if w > 1.5)
    print(f"  (3) |Wr|>1.5 suspect: random {r_rand}/12, folded "
          f"{r_fold}/12 (mean {np.mean(wr_rand):.2f}→{np.mean(wr_fold):.2f})"
          f"  ✓ (gemessen — die Protein-Zahl)")
    # (4) folded knot candidate built
    e4 = f4 = None
    for s in range(10):
        e4, f4 = make_knot(48, 2000 + s, 900, 2.0, HERE / "demo_smoke")
        if f4 is not None:
            break
    ok4 = f4 is not None
    print(f"  (4) {e4['name']}: E {e4.get('E0')}→{e4.get('E_final')}, "
          f"Rg {e4.get('rg_before')}→{e4.get('rg_after')}, "
          f"Wr={e4.get('writhe')}, "
          f"nyq={e4.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok4 else '✗ FAIL: ' + str(e4.get('rejected'))}")
    assert ok4
    print("[DEMO] 4/4 green — folding machine ready to fire "
          "(energy-driven knot statistics, protein track).")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo")
    ap.parse_args()
    cmd_demo(None)


if __name__ == "__main__":
    main()
