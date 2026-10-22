# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""
build_fremlin_cache.py — builds data/fremlin_cache/fremlin_cache.json from David Fremlin's pages
"Knots and their symmetries" (https://david.fremlin.de/knots/). The pages are not redistributed here.
Used ONLY by the cross-check stage 8b, after identification (blind test).

Usage:
  python3 data/build_fremlin_cache.py --fetch      # download index.htm + every knot page it links (1 s apart), then parse
  python3 data/build_fremlin_cache.py              # parse pages already in --pages
Options: --pages DIR (default data/fremlin_pages), --out FILE (default data/fremlin_cache/fremlin_cache.json
= the portable FREMLIN_CACHE of site_config).

Extracted per knot page: knot_id, source_url, title, crossing_upper (if stated), fremlin_appears_prime (bool),
symmetries_raw, amphichirality_raw, jones_t_raw (if listed explicitly), date.
The parsing functions are those of the authors' runs except for one widened pattern (Jones polynomial end, marked below);
on the pages of 11.05.2026 the output is byte-identical to the authors' cache (data/README.md).
"""
import argparse
import json
import re
import time
import urllib.request
from pathlib import Path

BASE = "https://david.fremlin.de/knots/"
HERE = Path(__file__).resolve().parent


# The server refuses Python's default user agent (HTTP 403, checked 22.10.2026); we send an honest one of our own.
USER_AGENT = "knot-identification-pipeline/1.0 (data recipe; https://knot-structures.stainlesssteel4u.de/)"


def _get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    return urllib.request.urlopen(req, timeout=60).read()


def fetch(pages: Path):
    pages.mkdir(parents=True, exist_ok=True)
    index = _get(BASE + "index.htm")
    (pages.parent / "fremlin_index.htm").write_bytes(index)
    links = sorted(set(re.findall(r'href="([^"/:]+\.htm)"', index.decode("iso-8859-1"))) - {"table.htm", "index.htm"})
    for name in links:
        time.sleep(1.0)
        (pages / name).write_bytes(_get(BASE + name))
    print(f"fetched {len(links)} pages from {BASE} into {pages}")


def strip_tags(html: str) -> str:
    # the HTML entities we need
    html = html.replace("&pi;", "π").replace("&nbsp;", " ")
    return re.sub(r"<[^>]+>", "", html)


def parse_knot_page(path: Path) -> dict:
    html = path.read_text(encoding="iso-8859-1", errors="replace")
    text = strip_tags(html)
    knot_id = path.stem
    rec = {
        "knot_id": knot_id,
        "source_url": f"https://david.fremlin.de/knots/{knot_id}.htm",
    }

    # title line: "Knot 31" / "Knot 31 Eternal knot" / etc.
    m = re.search(r"Knot\s+(\d[\w_]*)\s*(.*?)\n", text)
    if m:
        rec["title_id"] = m.group(1)
        descr = m.group(2).strip()
        if descr:
            rec["title_descr"] = descr

    # Crossing-upper bound
    m = re.search(r"Crossing number at most\s+(\d+)", text)
    if m:
        rec["crossing_upper"] = int(m.group(1))

    rec["fremlin_appears_prime"] = bool(re.search(r"appears to be prime", text))

    # symmetries line
    m = re.search(r"Symmetr(?:y|ies)\s+([^\n]+)", text)
    if m:
        rec["symmetries_raw"] = re.sub(r"\s+", " ", m.group(1)).strip()

    # Amphichiral
    if "amphichiral" in text:
        m = re.search(r"((?:fully\s+|negative\s+|positive\s+)?amphichiral[^\n,]*)", text)
        if m:
            rec["amphichirality_raw"] = re.sub(r"\s+", " ", m.group(1)).strip()

    # Jones polynomial (Knot 15331 style)
    # Changed from the authors' parser (22.10.2026): the polynomial also ends where a later note begins.
    # Fremlin's note of 24.5.26 on page 15331 follows the polynomial directly, and the old pattern lost it.
    m = re.search(r"Jones polynomial:\s*(.*?)(?:List of knots|table of symmetries|Note added|Back to list)",
                  text, re.DOTALL)
    if m:
        jones = re.sub(r"\s+", " ", m.group(1)).strip()
        # (kept raw: doubled signs and superscripts are not cleaned here)
        rec["jones_t_raw"] = jones

    # date at the end
    m = re.search(r"(\d{1,2}\.\d{1,2}\.\d{2,4})\s*$", text.strip())
    if m:
        rec["date"] = m.group(1)

    return rec


def main():
    ap = argparse.ArgumentParser(description="Build the Fremlin cache used by the cross-check stage.")
    ap.add_argument("--pages", type=Path, default=HERE / "fremlin_pages", help="folder with the downloaded knot pages")
    ap.add_argument("--out", type=Path, default=HERE / "fremlin_cache" / "fremlin_cache.json", help="output JSON")
    ap.add_argument("--fetch", action="store_true", help="download the pages first")
    args = ap.parse_args()
    if args.fetch:
        fetch(args.pages)
    htms = sorted(args.pages.glob("*.htm"))
    records = []
    for p in htms:
        if p.stem in ("questions", "table"):
            continue
        records.append(parse_knot_page(p))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(records, indent=2, ensure_ascii=False))
    print(f"{len(records)} knot records written to {args.out}")


if __name__ == "__main__":
    main()
