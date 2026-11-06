# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 131 — pairing_the_theta_blind (13a_1231 / 13a_1237).

THE mutant pair from the Dror/vdV draft (mail_theta_findings_*),
Δ-NOT-trivial (harder than Conway/KT). Loom-4 import from 126.
Braid words from spherogram census (snappy_env, local, documented in
131_theta_blind_paar.txt — construction lookup; identification blind):
  tb_1231 (13a_1231): [-1,2,2,-3,-3,-3,-3,2,1,2,-3,-3,-3,2,2]
  tb_1237 (13a_1237): [1,-2,-2,-2,1,3,1,-2,-2,-2,-2,3,3]

BETS: pre-registered in 131_theta_blind_paar.txt (W0 combinatorics ·
W1 blind ID · W2 Δ-equal-but-≠1 · W3 θ-blind (letter core) ·
W4 6D separation · W5 cross-check against 127). Do NOT fire.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
YARD = HERE.parent
sys.path.insert(0, str(YARD))
sys.path.insert(0, str(YARD / "126_widening_the_loom"))
from generate_loom4 import (  # noqa: E402  (Loom-4 word-agnostic)
    build_station, connector, station_crossings, outside_station,
    word_perm, N_STR)
from machine_common import (  # noqa: E402
    densify, embedding_report, polyline_to_fseries)

WORDS = {
    "tb_1231": [-1, 2, 2, -3, -3, -3, -3, 2, 1, 2, -3, -3, -3, 2, 2],
    "tb_1237": [1, -2, -2, -2, 1, 3, 1, -2, -2, -2, -2, 3, 3],
}
KNOTINFO = {"tb_1231": "13a_1231", "tb_1237": "13a_1237"}
SOURCE_NOTE = "spherogram census K13a1231/K13a1237 (snappy_env, 2026-06-07)"


def braid_word(tag: str):
    return [(abs(x), +1 if x > 0 else -1) for x in WORDS[tag]]


def build_curve(tag: str):
    word = braid_word(tag)
    paths, end_lane = build_station(word)
    pts, cur, visited = [], 0, []
    for _ in range(N_STR):
        visited.append(cur)
        pts += paths[cur]
        lane = end_lane[cur]
        pts += connector(lane)
        cur = lane
    if cur != 0 or sorted(visited) != list(range(N_STR)):
        raise RuntimeError(f"{tag}: permutation not a 4-cycle — no knot")
    return np.asarray(pts, float)


def is_4cycle(p):
    orbit, cur = set(), 0
    for _ in range(N_STR):
        orbit.add(cur)
        cur = p[cur]
    return len(orbit) == N_STR


def make_twin(tag: str, out_root: Path):
    name = tag
    P = build_curve(tag)
    n_c = station_crossings(P)
    soll = len(WORDS[tag])
    entry = {"name": name, "machine": "131_pairing_the_theta_blind",
             "braid_word_raw": WORDS[tag], "source": SOURCE_NOTE,
             "knotinfo_construction_label": KNOTINFO[tag],
             "station_crossings": n_c, "station_crossings_soll": soll,
             "predictions": {
                 "W1": f"Blind-ID: {tag} = {KNOTINFO[tag]} (Crosscheck)",
                 "W2": "Mutanten-Kollision: Delta BEIDE gleich aber !=1 "
                       "(haerter als Conway/KT); Jones/HOMFLY/vol gleich",
                 "W3": "theta(22/7,21/13) trennt NICHT (Brief-Kern, "
                       "mutant_pair_theta_verified) — Geometrie zum Befund",
                 "W4": "6D-Normalformen TRENNEN (6D sieht, wo theta blind)",
                 "W5": "Kreuzcheck 127 (Delta!=1 vs Delta==1, ein theta)"}}
    if n_c != soll:
        entry["rejected"] = f"Kreuzungszahl {n_c} != Soll {soll}"
        return entry, None
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
                source="maschinen_ideen/131", n_terms=nt,
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
    print("[DEMO] Smoke test machine 131 (theta-blind pair 13a_1231/1237)")
    out = HERE / "demo_smoke"
    ok = 0
    # S1: word lengths + 4-cycle
    s1 = (len(WORDS["tb_1231"]) == 15 and len(WORDS["tb_1237"]) == 13)
    perms = {}
    for tag in WORDS:
        p = word_perm(braid_word(tag))
        perms[tag] = p
        four = is_4cycle(p)
        s1 &= four
        print(f"  (S1) {tag} ({KNOTINFO[tag]}): {len(WORDS[tag])} letters"
              f" · permutation {p} · 4-cycle {'✓' if four else '✗'}")
    ok += s1
    # S2: pair geometry
    o1 = outside_station(build_curve("tb_1231"))
    o2 = outside_station(build_curve("tb_1237"))
    same_perm = perms["tb_1231"] == perms["tb_1237"]
    same_out = o1.shape == o2.shape and bool(np.array_equal(o1, o2))
    s2 = same_out if same_perm else True
    print(f"  (S2) permutations {'equal' if same_perm else 'different'}"
          + (f" ⟹ outside bit {'✓' if same_out else '✗'}" if same_perm
             else " ⟹ separated by class (lesson of 124, documented) ✓"))
    ok += s2
    # S3: build both
    s3 = True
    for tag in WORDS:
        entry, folder = make_twin(tag, out)
        if folder is None:
            print(f"  (S3) {tag}: REJECTED — {entry['rejected']} ✗")
            s3 = False
        else:
            print(f"  (S3) {tag}: "
                  f"d_min_rel={entry['embedding']['d_min_rel']}"
                  f" · n_terms={entry['fseries']['n_terms']}"
                  f" · nyq={entry['fseries']['nyquist_dev_rel']} ✓")
    ok += s3
    print(f"[DEMO] {ok}/3 green — "
          f"{'letter pair in place.' if ok == 3 else 'RED.'}"
          f" Bets: 131_theta_blind_paar.txt. Doctrine: do NOT fire.")
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
