# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 127 — pairing_the_conway_twins (Conway/KT mutant pair).

Loom-4 import from 126 (word-agnostic). Braid words from spherogram-
census (snappy_env, 2026-06-07, documented in 127_conway_zwillinge.txt
— construction lookup; pipeline identification stays blind):
  cw (K11n34, Conway):              [-1,2,-1,2,-1,3,-2,-2,-1,3,3]
  kt (K11n42, Kinoshita–Terasaka):  [-1,2,2,-3,-3,2,1,-2,-2,3,-2,3,-2]
DOUBLE VALUE: Dror/vdV θ question on the classic (W3) + redeemed
6D mutant target test (W5, cancelled build stage 3 through the back door).

BETS: preregistered in 127_conway_zwillinge.txt (W0 combinatorics ·
W1 blind ID 11n34/11n42 · W2 mutant collision Δ≡1/Jones/vol ·
W3 θ separation? · W4 genus suspicion · W5 6D separation). Do NOT fire.
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

# Source: spherogram (snappy_env), 2026-06-07 — see founding txt
WORDS = {
    "cw": [-1, 2, -1, 2, -1, 3, -2, -2, -1, 3, 3],
    "kt": [-1, 2, 2, -3, -3, 2, 1, -2, -2, 3, -2, 3, -2],
}
SOURCE_NOTE = "spherogram census K11n34/K11n42 (snappy_env, 2026-06-07)"


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


def make_twin(tag: str, out_root: Path):
    name = f"conway_{tag}"
    P = build_curve(tag)
    n_c = station_crossings(P)
    soll = len(WORDS[tag])
    entry = {"name": name, "machine": "127_pairing_the_conway_twins",
             "braid_word_raw": WORDS[tag], "source": SOURCE_NOTE,
             "station_crossings": n_c, "station_crossings_soll": soll,
             "predictions": {
                 "W1": "Blind-ID: cw=11n34, kt=11n42 (Crosscheck)",
                 "W2": "Mutanten-Kollision: Delta=1 beide, det=1, "
                       "Jones/HOMFLY/vol paarweise identisch",
                 "W3": "theta(22/7,21/13): trennt es das Paar? "
                       "(DIE Dror/vdV-Frage)",
                 "W4": "Genus-VERDACHT (zweiseitig): g(kt) != g(cw)",
                 "W5": "6D-Normalformen der beiden TRENNEN sich "
                       "(eingeloester Mutanten-Zieltest Baustufe 3)"}}
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
                source="maschinen_ideen/127", n_terms=nt,
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
    print("[DEMO] Smoke test machine 127 (Conway twins on Loom-4)")
    out = HERE / "demo_smoke"
    ok = 0
    # S1: word lengths + 4-cycle + source byte fidelity
    s1 = (len(WORDS["cw"]) == 11 and len(WORDS["kt"]) == 13)
    perms = {}
    for tag in ("cw", "kt"):
        p = word_perm(braid_word(tag))
        perms[tag] = p
        cyc = (sorted(p) == list(range(N_STR)) and p != list(range(N_STR)))
        # 4-cycle check: orbit of 0 has length 4
        orbit, cur = set(), 0
        for _ in range(N_STR):
            orbit.add(cur)
            cur = p[cur]
        four = len(orbit) == N_STR
        s1 &= cyc and four
        print(f"  (S1) {tag}: {len(WORDS[tag])} letters · "
              f"permutation {p} · 4-cycle {'✓' if four else '✗'}")
    ok += s1
    # S2: pair geometry — same loom; document the outside comparison
    o_cw = outside_station(build_curve("cw"))
    o_kt = outside_station(build_curve("kt"))
    same_perm = perms["cw"] == perms["kt"]
    same_out = (o_cw.shape == o_kt.shape
                and bool(np.array_equal(o_cw, o_kt)))
    s2 = same_out if same_perm else True
    _perm = 'gleich' if same_perm else 'VERSCHIEDEN'
    _aussen = ('bit-identisch ' + ('✓' if same_out else '✗')) if same_perm \
        else 'klassen-getrennt (dokumentiert, 124er-Lehre) ✓'
    print(f"  (S2) permutations {_perm} ⟹ outside {_aussen}")
    ok += s2
    # S3: embedding + Nyquist BOTH (the pair is the purpose)
    s3 = True
    for tag in ("cw", "kt"):
        entry, folder = make_twin(tag, out)
        if folder is None:
            print(f"  (S3) conway_{tag}: REJECTED — {entry['rejected']} ✗")
            s3 = False
        else:
            print(f"  (S3) conway_{tag}: "
                  f"d_min_rel={entry['embedding']['d_min_rel']}"
                  f" · n_terms={entry['fseries']['n_terms']}"
                  f" · nyq={entry['fseries']['nyquist_dev_rel']} ✓")
    ok += s3
    print(f"[DEMO] {ok}/3 green — {'twins in place.' if ok == 3 else 'RED.'}"
          f" Bets: 127_conway_zwillinge.txt. Doctrine: do NOT fire.")
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
