# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
MACHINE 11 — fusing_mirror_twins (symmetric unions).

Build: 8-strand plat. Left block (slots 0-3) = partial knot K as a
4-plat word (02 continued-fraction machinery), right block (slots 4-7) =
its MIRROR (slot pair i → 6−i, twist sign flipped), seam =
σ4 twists (slots 3,4) ON the mirror axis. Plat closure (07).

BETS (flip clauses in the founding txt):
  - Component count per seam twist number t MEASURED (scan, no parity
    guessed); links rejected.
  - det = det(K)² as CANDIDATE (Lamm, training knowledge) — at firing time over
    several t bit-check. Ribbon ⇒ sig=0, Arf=0, det≡1 mod 8.
  - t=0 anchor: must be K # mirror(K).

Doctrine: do NOT ignite before Monday.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "02_plugging_tangles"))
sys.path.insert(0, str(HERE.parent / "07_closing_differently"))
sys.path.insert(0, str(HERE.parent))
from generate_tangle_knots import (cf_value, cf_to_pq, odd_cf_of,  # noqa: E402
                                   cf_to_braid_word)
from generate_braid_closures import (predict_components,  # noqa: E402
                                     make_knot as braid_make_knot)

N_STRANDS = 8
AXIS_PAIR = 3            # σ4: slots 3,4 — the mirror axis


def mirror_letter(letter):
    """Mirror about the axis between slot 3 and 4: pair i → 6−i, k → −k."""
    i, k = letter
    return (6 - i, -k)


def union_word(cf, t_axis):
    """Symmetric-union word: K block + mirror block + seam twists."""
    x = cf_value(cf)
    p, q = abs(x.numerator), abs(x.denominator)
    if p % 2 == 0:
        raise ValueError(f"Partial knot {cf}: p={p} even = link — "
                         "not a symmetric-union candidate")
    base = cf_to_braid_word(odd_cf_of(p, q))      # slot pairs 0/1 (4-plat)
    left = list(base)
    right = [mirror_letter(l) for l in left]
    word = left + right
    if t_axis:
        word.append((AXIS_PAIR, t_axis))
    return word, (p, q)


def t_scan(cf, t_range=range(-4, 5)):
    """Component count per seam twist (flip clause 3: measure, do not guess)."""
    out = {}
    for t in t_range:
        w, _ = union_word(cf, t)
        out[t] = predict_components(w, N_STRANDS, "plat")
    return out


def make_knot(cf, t_axis, out_root):
    word, (p, q) = union_word(cf, t_axis)
    name = f"su_{'_'.join(map(str, cf))}_t{t_axis}".replace("-", "m")
    entry, folder = braid_make_knot(name, word, N_STRANDS, "plat", 1,
                                    Path(out_root))
    entry["machine"] = "11_fusing_mirror_twins"
    entry["partial_knot"] = {"cf": list(cf), "p": p, "q": q}
    entry["axis_twists"] = t_axis
    entry["predictions"].update({
        "det_candidate": p * p,
        "det_note": f"det = det(K)² = {p*p} — KANDIDAT (Lamm, Lit-Check "
                    "offen; Umdreh-Riss 2), Zünd-Zeit bit-prüfen",
        "ribbon_candidate": True,
        "sig_candidate": 0,
        "det_mod8_candidate": 1,
        "anchor_t0": "t=0 ⇒ K # spiegel(K)" if t_axis == 0 else None,
    })
    return entry, folder


# ----------------------------------------------------------------------- CLI

def cmd_scan(args):
    res = t_scan(args.cf)
    knots = [t for t, c in res.items() if c == 1]
    print(f"[SCAN] Partial {args.cf}: components per seam twist {res}")
    print(f"  knot twists: {knots}")


def cmd_make(args):
    man_path = HERE / "MANIFEST.json"
    man = json.loads(man_path.read_text()) if man_path.exists() else {"knots": []}
    entry, folder = make_knot(args.cf, args.t, HERE / "out_fseries")
    man["knots"] = [k for k in man["knots"] if k["name"] != entry["name"]] + [entry]
    man_path.write_text(json.dumps(man, indent=1))
    print(f"[MAKE] {entry['name']}: det candidate "
          f"{entry['predictions']['det_candidate']} → "
          f"{folder if folder else entry.get('rejected')}")


def cmd_demo(args):
    print("[DEMO] Smoke test machine 11 (symmetric unions)")
    # 1) Mirror map: involutive + sign-flipping
    l = (1, 3)
    ok1 = mirror_letter(mirror_letter(l)) == l and mirror_letter(l) == (5, -3)
    print(f"  (1) mirror involutive, (σ2,+3)→(σ6,−3)  "
          f"{'✓' if ok1 else '✗ FAIL'}")
    assert ok1
    # 2) link partial rejected (p even)
    try:
        union_word([2, 3, 2], 1)        # 16/7
        ok2 = False
    except ValueError:
        ok2 = True
    print(f"  (2) Partial [2,3,2] (p=16) rejected  {'✓' if ok2 else '✗ FAIL'}")
    assert ok2
    # 3) t scan trefoil partial [3]: measure components, find knot twists
    res = t_scan([3])
    knots = [t for t, c in res.items() if c == 1]
    ok3 = len(knots) > 0 and 0 in res
    print(f"  (3) [3]-Scan: {res} → knot twists {knots}  "
          f"{'✓' if ok3 else '✗ FAIL'}")
    assert ok3
    # 4) Build at the first knot twist: det candidate 9 (=3²), all gates
    t0 = sorted(knots, key=abs)[0]
    e, f = make_knot([3], t0, HERE / "demo_smoke")
    ok4 = f is not None and e["predictions"]["det_candidate"] == 9
    print(f"  (4) su_3_t{t0}: det candidate {e['predictions']['det_candidate']} "
          f"(3²) · embedded={e.get('embedding', {}).get('embedded')} · "
          f"Nyquist={e.get('fseries', {}).get('nyquist_dev_rel')}  "
          f"{'✓' if ok4 else '✗ FAIL: ' + str(e.get('rejected'))}")
    assert ok4
    # 5) The 7² point: partial [3,2] (5_2, det 7) → candidate 49
    res2 = t_scan([3, 2])
    knots2 = [t for t, c in res2.items() if c == 1]
    e2, f2 = make_knot([3, 2], sorted(knots2, key=abs)[0], HERE / "demo_smoke")
    ok5 = f2 is not None and e2["predictions"]["det_candidate"] == 49
    print(f"  (5) Partial [3,2]: det candidate "
          f"{e2['predictions']['det_candidate']} (7² — the corpus point!) "
          f"built  {'✓' if ok5 else '✗ FAIL: ' + str(e2.get('rejected'))}")
    assert ok5
    print("[DEMO] 5/5 green — mirror-twin machine ready to fire "
          "(firing after Monday).")


def main():
    ap = argparse.ArgumentParser(description="Machine 11 — symmetric unions")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("scan")
    p1.add_argument("cf", type=int, nargs="+")
    p2 = sub.add_parser("make")
    p2.add_argument("cf", type=int, nargs="+")
    p2.add_argument("--t", type=int, required=True)
    sub.add_parser("demo")
    args = ap.parse_args()
    {"scan": cmd_scan, "make": cmd_make, "demo": cmd_demo}[args.cmd](args)


if __name__ == "__main__":
    main()
