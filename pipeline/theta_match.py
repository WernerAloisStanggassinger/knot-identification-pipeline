# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
theta_match.py — Theta invariant as corpus-internal AND corpus-external
match criterion (Drors-QR-Phase-2, Discovery #172/#173).

Bar-Natan / van der Veen 2509.18456: Theta is a 2-variable polynomial
invariant (in t1, t2), proven strongly distinguishing for knots up to c=15.
Numerical specialisation t1=22/7, t2=21/13 computable for c up to ~600.

Match paths (prioritised):
  1. corpus-internal Theta-num match (corpus siblings)
  2. external DB (knots.sage c<=14): direct match + mirror search
     (theta_A == ±theta_B  →  potentially isotopic or mirror pair)

Call:
  python theta_match.py [--target <stem>] [--corpus-db corpus_db.json]
                        [--external-db knots_sage_theta_db.json]

Output:
  /tmp/theta_match_<target_or_all>.json with match findings.
"""
from __future__ import annotations
import site_config
import argparse
import json
from collections import defaultdict
from pathlib import Path
from fractions import Fraction

CORPUS_DB = Path(__file__).resolve().parent / "corpus_db.json"
EXTERNAL_DB = Path(site_config.THETA_DB)


def _theta_to_fraction(s):
    """Parse theta-num string to Fraction. Returns None on "0" or None or parse error."""
    if s is None: return None
    s = str(s).strip()
    if s in ("0", "", "None"): return None
    if s.startswith("<"): return None  # error markers
    try:
        if "/" in s:
            n, d = s.split("/", 1)
            return Fraction(int(n.strip()), int(d.strip()))
        return Fraction(int(s.strip()))
    except (ValueError, TypeError):
        return None


def _load_external_db(path):
    """Load external DB, index by theta_num (as Fraction). Mirror index in addition."""
    if not path.exists():
        return None, None, None
    raw = json.loads(path.read_text())
    by_theta = {}      # Fraction → list of (name, c)
    by_neg_theta = {}  # for mirror search
    n_skipped = 0
    for name, e in raw.items():
        f = _theta_to_fraction(e.get("theta_at_22_7__21_13"))
        if f is None:
            n_skipped += 1
            continue
        by_theta.setdefault(f, []).append((name, e.get("c")))
        by_neg_theta.setdefault(-f, []).append((name, e.get("c")))
    return raw, by_theta, by_neg_theta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--corpus-db", type=Path, default=CORPUS_DB)
    ap.add_argument("--external-db", type=Path, default=EXTERNAL_DB)
    ap.add_argument("--target", type=str, default=None,
                    help="stem of a corpus knot; default = all singletons.")
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()

    db = json.loads(args.corpus_db.read_text())
    knots = db["knots"]

    # load external DB
    ext_raw, ext_by_theta, ext_by_neg_theta = _load_external_db(args.external_db)
    ext_size = len(ext_raw) if ext_raw else 0

    # Index by theta_full (sym) AND theta_at_22_7__21_13 (num).
    # IMPORTANT: theta_full="0" occurs for amphichiral knots (Theta(K) =
    # Theta(mirror K) = -Theta(K) → Theta = 0). This class is NOT
    # an identification hint but an amphichirality indicator. Filter it out.
    by_theta_sym = defaultdict(list)
    by_theta_num = defaultdict(list)
    n_amphi_theta_zero = 0
    n_theta_error = 0

    def _is_valid_theta_sym(ts):
        if not ts: return False
        s = ts.strip()
        if s == "0": return False
        if s.startswith("<") and s.endswith(">"):
            # error markers like <timeout_theta_main>, <skipped>, etc.
            return False
        return True

    for k in knots:
        ts = k.get("theta_full")
        tn = k.get("theta_at_22_7__21_13")
        if _is_valid_theta_sym(ts):
            by_theta_sym[ts].append(k["stem"])
        elif ts == "0":
            n_amphi_theta_zero += 1
        elif ts and ts.startswith("<"):
            n_theta_error += 1
        if tn is not None and str(tn) != "0":
            by_theta_num[str(tn)].append(k["stem"])

    # balance: are there Theta classes with > 1 member?
    sym_classes = {ts: members for ts, members in by_theta_sym.items() if len(members) > 1}
    num_classes = {tn: members for tn, members in by_theta_num.items() if len(members) > 1}

    # ---- EXTERNAL DB: corpus → knots.sage match ------------------------------
    external_matches = []  # (stem, c, ext_name, kind) — kind ∈ {direct, mirror}
    if ext_by_theta is not None:
        for k in knots:
            f = _theta_to_fraction(k.get("theta_at_22_7__21_13"))
            if f is None:
                continue
            stem = k["stem"]
            c = k.get("c")
            direct = ext_by_theta.get(f, [])
            mirror = ext_by_neg_theta.get(f, [])
            for ext_name, ext_c in direct:
                external_matches.append({"stem": stem, "c": c,
                                          "external": ext_name, "ext_c": ext_c,
                                          "kind": "direct"})
            for ext_name, ext_c in mirror:
                external_matches.append({"stem": stem, "c": c,
                                          "external": ext_name, "ext_c": ext_c,
                                          "kind": "mirror"})

    out = {
        "total_knots_in_db": len(knots),
        "n_with_theta_sym": sum(1 for k in knots if k.get("theta_full")),
        "n_with_theta_num": sum(1 for k in knots if k.get("theta_at_22_7__21_13") is not None),
        "n_amphi_theta_zero": n_amphi_theta_zero,
        "n_theta_error_markers": n_theta_error,
        "n_sym_classes_multi": len(sym_classes),
        "n_num_classes_multi": len(num_classes),
        "external_db_size": ext_size,
        "external_matches": external_matches,
        "n_corpus_with_external_match": len({m["stem"] for m in external_matches}),
    }

    # per Theta-sym class: show members
    out["sym_classes"] = []
    for ts, members in sorted(sym_classes.items(), key=lambda x: -len(x[1])):
        kis = [(stem, next((k.get("knotinfo_id") or k.get("crosscheck_knotinfo")
                            or "—") for k in knots if k["stem"] == stem))
               for stem in members]
        c_set = set()
        for stem in members:
            kk = next(k for k in knots if k["stem"] == stem)
            c_set.add(kk.get("c"))
        out["sym_classes"].append({
            "size": len(members),
            "c_values": sorted(c_set, key=lambda x: (x is None, x)),
            "members": [{"stem": s, "ki": ki} for s, ki in kis],
            "theta_full_head": (ts or "")[:80] + ("..." if ts and len(ts) > 80 else ""),
        })

    # if a target is specified: show it explicitly
    if args.target:
        target = next((k for k in knots if k["stem"] == args.target), None)
        if target is None:
            out["target_error"] = f"stem '{args.target}' not in corpus_db"
        else:
            ts = target.get("theta_full")
            tn = target.get("theta_at_22_7__21_13")
            out["target"] = {
                "stem": target["stem"],
                "c": target.get("c"), "det": target.get("det"),
                "vol": target.get("vol"), "sym": target.get("sym"),
                "theta_at_22_7__21_13": tn,
            }
            if ts:
                same_sym = by_theta_sym.get(ts, [])
                out["target"]["theta_sym_class"] = [s for s in same_sym if s != target["stem"]]
            if tn is not None:
                same_num = by_theta_num.get(str(tn), [])
                out["target"]["theta_num_class"] = [s for s in same_num if s != target["stem"]]

    out_path = args.output or Path(f"/tmp/theta_match_{args.target or 'all'}.json")
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False, default=str))
    print(f"# Theta match balance")
    print(f"  Total corpus knots:         {out['total_knots_in_db']}")
    print(f"  With theta_full (sym):      {out['n_with_theta_sym']}")
    print(f"  With theta_num at 22/7,21/13: {out['n_with_theta_num']}")
    print(f"  Theta-sym classes with >1 member: {out['n_sym_classes_multi']}")
    print(f"  Theta-num classes with >1 member: {out['n_num_classes_multi']}")
    print(f"  External DB (knots.sage):    {ext_size} entries")
    print(f"  Corpus knots with external match: {out['n_corpus_with_external_match']}")
    if external_matches:
        print()
        print(f"## External matches (corpus ↔ knots.sage):")
        # aggregate per stem
        by_stem = defaultdict(list)
        for m in external_matches:
            by_stem[m["stem"]].append(m)
        for stem in sorted(by_stem.keys()):
            ms = by_stem[stem]
            direct = [m for m in ms if m["kind"] == "direct"]
            mirror = [m for m in ms if m["kind"] == "mirror"]
            d_names = ", ".join(m["external"] for m in direct[:5]) + (f" (+{len(direct)-5})" if len(direct) > 5 else "")
            m_names = ", ".join(m["external"] for m in mirror[:5]) + (f" (+{len(mirror)-5})" if len(mirror) > 5 else "")
            parts = []
            if direct: parts.append(f"direct: {d_names}")
            if mirror: parts.append(f"mirror: {m_names}")
            print(f"  {stem:<40}  c={ms[0].get('c')}  → {'; '.join(parts)}")
    print()
    print(f"## Theta-sym classes (top 10):")
    for cls in out["sym_classes"][:10]:
        members_s = ", ".join(m["stem"] + (f"={m['ki']}" if m["ki"] != "—" else "") for m in cls["members"][:5])
        more = f" (+{cls['size']-5})" if cls['size'] > 5 else ""
        print(f"  size={cls['size']:>2}  c={cls['c_values']}  {members_s}{more}")
    if args.target and "target" in out:
        t = out["target"]
        print()
        print(f"## Target: {t['stem']} (c={t['c']}, det={t.get('det')})")
        print(f"  theta_num at (22/7, 21/13) = {t.get('theta_at_22_7__21_13')}")
        if t.get("theta_sym_class"):
            print(f"  THETA-SYM CLASS: {t['theta_sym_class']}")
        if t.get("theta_num_class"):
            print(f"  THETA-NUM CLASS: {t['theta_num_class']}")


if __name__ == "__main__":
    main()
