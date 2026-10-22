# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""
build_knotinfo_index.py — builds data/knotinfo_index/{records.json, index_c_det_sig.json}
from the database_knotinfo package (KnotInfo, C. Livingston and A. H. Moore).
Used ONLY by the cross-check stage 8b, after identification (blind test).

Usage:
  conda run -n sage_stable python data/build_knotinfo_index.py [--out-dir DIR]
  (default DIR = data/knotinfo_index next to this script = the portable KNOTINFO_INDEX folder of site_config)

Fields per record:
  name, c, det, sig, alt (alternating), fibered, genus (three_genus), arf,
  sym_type, full_sym_group, jones, alex.
fibered/alt are normalised to bool ('Y'/'N' → True/False); genus and arf to int
(KnotInfo gives them as strings). None stays None where KnotInfo leaves a value empty.

Reproduction check (22.10.2026, database_knotinfo 2026.10.5): the output is byte-identical to the
index used in the authors' runs (built 14.05.2026): 12966 records, 1792 (c, det, sig) classes.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path


def _yn_to_bool(v):
    if v is None:
        return None
    s = str(v).strip().upper()
    if s in ("Y", "YES", "TRUE", "1"):
        return True
    if s in ("N", "NO", "FALSE", "0"):
        return False
    return None


def _to_int(v):
    if v is None or v == "":
        return None
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "knotinfo_index",
    )
    args = ap.parse_args()

    try:
        from database_knotinfo import link_list
    except ImportError:
        print("database_knotinfo not installed; run inside sage_stable env.", file=sys.stderr)
        sys.exit(2)

    ll = link_list()
    records = {}
    skipped_links = 0
    for r in ll:
        raw_name = r.get("name") or ""
        # KnotInfo records: knots have names like '3_1', '12n_839'.
        # Links have an 'L' prefix (e.g. 'L2a1') — not relevant for the index.
        if raw_name.startswith("L"):
            skipped_links += 1
            continue
        c = _to_int(r.get("crossing_number"))
        det = _to_int(r.get("determinant"))
        sig = _to_int(r.get("signature"))
        if c is None or det is None:
            continue
        # Naming convention of the pipeline (2026-05-14):
        #   c<=10: keep a K prefix (K3_1, K8_19, K10_43)
        #   c>=11: no K prefix (12n_837, 11a_1) — as on the KnotInfo website
        name = ("K" + raw_name) if c <= 10 else raw_name
        records[name] = {
            "name": name,
            "c": c,
            "det": det,
            "sig": sig,
            "alt": _yn_to_bool(r.get("alternating")),
            "fibered": _yn_to_bool(r.get("fibered")),
            "genus": _to_int(r.get("three_genus")),
            "arf": _to_int(r.get("arf_invariant")),
            "sym_type": (r.get("symmetry_type") or "").strip() or None,
            "full_sym_group": (r.get("full_symmetry_group") or "").strip() or None,
            "jones": r.get("jones_polynomial") or "",
            "alex": r.get("alexander_polynomial") or "",
        }

    # Index by (c, det, sig)
    index = {}
    for name, rec in records.items():
        c, det, sig = rec["c"], rec["det"], rec["sig"]
        if sig is None:
            continue
        key = f"{c}|{det}|{sig}"
        index.setdefault(key, []).append(name)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    records_path = args.out_dir / "records.json"
    index_path = args.out_dir / "index_c_det_sig.json"

    records_path.write_text(json.dumps(records, indent=None, ensure_ascii=False))
    index_path.write_text(json.dumps(index, indent=None, ensure_ascii=False))

    n = len(records)
    fib_known = sum(1 for r in records.values() if r["fibered"] is not None)
    g_known = sum(1 for r in records.values() if r["genus"] is not None)
    print(f"Wrote {n} records, {len(index)} (c,det,sig)-classes")
    print(f"  Skipped {skipped_links} Link records")
    print(f"  fibered populated: {fib_known}/{n}")
    print(f"  genus populated:   {g_known}/{n}")
    print(f"  → {records_path}")
    print(f"  → {index_path}")


if __name__ == "__main__":
    main()
