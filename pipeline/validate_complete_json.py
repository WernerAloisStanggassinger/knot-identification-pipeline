#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""Schema validator for pflicht_invarianten_complete.json files.

Two levels:
  - HARD: structural duty — pflichtliste block must exist, with the
          three sections DIAGRAMM / ALGEBRA / HYPERBOLISCH. Error = record
          is NOT ingested by the build.
  - SOFT: type/format warnings on single fields (det as int, alex as str
          or null, etc.). Warnings = ingested, but marked in the log.

Usage:
  - Library: from validate_complete_json import validate_file
              v = validate_file(path)  # -> {ok, errors, warnings, summary}
  - CLI:    python validate_complete_json.py                   # scans all corpus roots
            python validate_complete_json.py --path <file>     # single file
            python validate_complete_json.py --strict          # exit 1 on warnings
            python validate_complete_json.py --json-report     # machine-readable

Design note: deliberately hand-rolled, no jsonschema dependency. Fields are
defined broadly (every knot has different invariants available; composite/
chiral/torus records are filled differently) — the validator checks
existence + type, NOT completeness.
"""
from __future__ import annotations
import argparse
import site_config
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(site_config.ROOT)
CORPUS_ROOTS = [
    "erzeugte_knoten",
    "own_knots",
    "entwicklung_automation",
    "Claude_piplineentwicklung_und_validierung",
]

# type spec: (field_path_tuple, allowed_types, optional)
# allowed_types: tuple of python types; None means "may be null"
FIELD_SPEC: list[tuple[tuple[str, ...], tuple, bool]] = [
    # structural mandatory sections (HARD)
    (("pflichtliste",),                  (dict,),                False),
    (("pflichtliste", "DIAGRAMM"),       (dict,),                False),
    (("pflichtliste", "ALGEBRA"),        (dict,),                False),
    (("pflichtliste", "HYPERBOLISCH"),   (dict,),                False),
    # DIAGRAMM
    (("pflichtliste", "DIAGRAMM", "crossing_number_min"),  (int, type(None)),  True),
    (("pflichtliste", "DIAGRAMM", "pd_code"),              (str, list, type(None)), True),
    (("pflichtliste", "DIAGRAMM", "braid_word"),           (str, list, type(None)), True),
    (("pflichtliste", "DIAGRAMM", "braid_n_strands"),      (int, type(None)),  True),
    (("pflichtliste", "DIAGRAMM", "alternating"),          (bool, type(None)), True),
    # ALGEBRA
    (("pflichtliste", "ALGEBRA", "determinant"),           (int, type(None)),  True),
    (("pflichtliste", "ALGEBRA", "signature"),             (int, type(None)),  True),
    (("pflichtliste", "ALGEBRA", "arf_invariant"),         (int, type(None)),  True),
    (("pflichtliste", "ALGEBRA", "genus"),                 (int, type(None)),  True),
    (("pflichtliste", "ALGEBRA", "fibered"),               (bool, type(None)), True),
    (("pflichtliste", "ALGEBRA", "alexander_polynomial"),  (str, type(None)),  True),
    (("pflichtliste", "ALGEBRA", "jones_polynomial"),      (str, type(None)),  True),
    (("pflichtliste", "ALGEBRA", "homfly_polynomial"),     (str, type(None)),  True),
    (("pflichtliste", "ALGEBRA", "conway_polynomial"),     (str, type(None)),  True),
    (("pflichtliste", "ALGEBRA", "kauffman_polynomial"),   (str, type(None)),  True),
    # HYPERBOLISCH
    (("pflichtliste", "HYPERBOLISCH", "is_amphicheiral"),      (bool, type(None)), True),
    (("pflichtliste", "HYPERBOLISCH", "is_invertible_knot"),   (bool, type(None)), True),
    (("pflichtliste", "HYPERBOLISCH", "hyperbolic_volume"),    (float, int, type(None)), True),
    (("pflichtliste", "HYPERBOLISCH", "symmetry_group"),       (str, type(None)),  True),
    (("pflichtliste", "HYPERBOLISCH", "canonical_retriangulation"), (str, type(None)), True),
]


def _get(d: dict, path: tuple[str, ...]) -> tuple[bool, Any]:
    """Navigiere path, return (exists, value)."""
    cur = d
    for p in path:
        if not isinstance(cur, dict) or p not in cur:
            return (False, None)
        cur = cur[p]
    return (True, cur)


def validate_file(path: Path) -> dict:
    """Validate a single _complete.json file.

    Returns: {
        "ok": bool,                   # False on HARD errors
        "errors": list[str],          # HARD errors (structural duty violated)
        "warnings": list[str],        # SOFT errors (field type, optional)
        "path": str,
    }
    """
    result = {"ok": True, "errors": [], "warnings": [], "path": str(path)}
    try:
        d = json.load(open(path))
    except Exception as e:
        result["ok"] = False
        result["errors"].append(f"JSON-Parse-Fehler: {type(e).__name__}: {e}")
        return result
    if not isinstance(d, dict):
        result["ok"] = False
        result["errors"].append(f"Top-Level ist {type(d).__name__}, erwartet dict")
        return result
    for fpath, allowed, optional in FIELD_SPEC:
        exists, val = _get(d, fpath)
        dotted = ".".join(fpath)
        if not exists:
            if not optional:
                result["ok"] = False
                result["errors"].append(f"PFLICHT-Feld fehlt: {dotted}")
            continue
        if not isinstance(val, allowed):
            msg = f"Type mismatch {dotted}: is {type(val).__name__}, expected one of {[t.__name__ for t in allowed]}"
            if optional:
                result["warnings"].append(msg)
            else:
                result["ok"] = False
                result["errors"].append(msg)
    return result


def scan_corpus() -> list[Path]:
    """Collect all _pflicht_invarianten*_complete.json from all CORPUS_ROOTS."""
    paths = []
    for root in CORPUS_ROOTS:
        rp = ROOT / root
        if not rp.exists():
            continue
        for p in rp.rglob("*_pflicht_invarianten*_complete.json"):
            if any(part.startswith("_sculpture_") for part in p.parts):
                continue  # sculpture symlink
            paths.append(p)
    return paths


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--path", help="Check a single file (instead of the full corpus)")
    ap.add_argument("--strict", action="store_true", help="Exit 1 also on warnings")
    ap.add_argument("--json-report", action="store_true", help="Machine-readable JSON-Output")
    ap.add_argument("--quiet", action="store_true", help="Print only errors/warnings, no pass lines")
    args = ap.parse_args()

    if args.path:
        paths = [Path(args.path)]
    else:
        paths = scan_corpus()

    results = [validate_file(p) for p in paths]
    n_total = len(results)
    n_ok = sum(1 for r in results if r["ok"])
    n_err = n_total - n_ok
    n_warn_files = sum(1 for r in results if r["warnings"])
    n_warn_total = sum(len(r["warnings"]) for r in results)

    if args.json_report:
        # print only fail+warn records (compact)
        bad = [r for r in results if r["errors"] or r["warnings"]]
        print(json.dumps({
            "n_total": n_total,
            "n_ok": n_ok,
            "n_err": n_err,
            "n_warn_files": n_warn_files,
            "n_warn_total": n_warn_total,
            "issues": bad,
        }, indent=2))
    else:
        for r in results:
            if r["errors"]:
                print(f"FAIL  {r['path']}")
                for e in r["errors"]:
                    print(f"    ERROR: {e}")
            elif r["warnings"]:
                if not args.quiet:
                    print(f"WARN  {r['path']}")
                for w in r["warnings"]:
                    print(f"    WARN: {w}")
            else:
                if not args.quiet:
                    print(f"OK    {r['path']}")
        print()
        print(f"== Summary ==")
        print(f"  total:        {n_total}")
        print(f"  ok:           {n_ok}")
        print(f"  errors:       {n_err}")
        print(f"  warn-files:   {n_warn_files}")
        print(f"  warnings:     {n_warn_total}")

    exit_code = 0
    if n_err > 0:
        exit_code = 2
    elif args.strict and n_warn_files > 0:
        exit_code = 1
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
