#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Knot-Identification Pipeline — Werner Alois Stanggassinger & Claude Code (Anthropic).
"""build_naming_registry.py — the display-name layer (Werner directive 2026-06-16).

PROBLEM (Werner): we use the NAME as IDENTITY, but identity is discovered (computed),
not assigned. For novel knots there is no table authority → we invent param names
(S_39_d7_D7_chiral, m121_ac_k3), which are ugly and drift all the way onto the website.

FIX: separate two layers.
  - CANONICAL KEY  = distinct_knot_id (from corpus_db, derived from the mathematics, join key).
  - DISPLAY NAME   = free label in THIS registry, keyed on distinct_knot_id.
      * known  -> KnotInfo name (auto)
      * novel + sculpture realization -> artist name of the sculpture (auto, anchor like nitro->S_24_d961)
      * novel without a name, symmetry known -> SYSTEMATIC S designation S_<c>_d<det>_<sym>_<chir>
        (Werner 2026-06-21 „bring die s bezeichnungen in ordnung" (EN: put the S designations in order); c taken raw as upper bound;
        artist name remains a later override). Collision (several distinct knots share the
        4-tuple) -> _v<int(vol)> discriminator like the existing _v14/_v27 clusters.
      * novel without a name, symmetry MISSING (sym=None / degenerate) -> status "needs_symmetry"
        (its own work site, flagged in naming_needs_symmetry.txt — cannot be S-named without symmetry).

Corpus/web/cross-references hang on the KEY; the display name is freely changeable without
breaking identity. Renaming for the website never breaks links again.

IDEMPOTENT: an assigned artist name (name_type=artist, status=assigned) is NEVER overwritten —
only provenance/realization count/status are refreshed. Systematic S names are recomputed
deterministically (no drift). Safe to run after every /corpus-rebuild.

S-NAME ENCODING exactly as sort_and_render.fingerprint_from_pflicht (l.139):
    sym = re.sub(r'\\s+','',sym).replace('+','_').replace('/','')   # Z/2->Z2, Z/2 + Z/2->Z2_Z2
    chir = 'amphi' if amphi else ('chiral' if amphi is not None else None)

Usage:
    python3 build_naming_registry.py            # update registry + queue reports
    python3 build_naming_registry.py --check     # check only (exit 1 if needs_name>0 or orphans), write nothing
"""
import json, os, re, sys, datetime
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
CORPUS = os.path.join(HERE, "corpus_db.json")
REGISTRY = os.path.join(HERE, "naming_registry.json")

def load(p, default=None):
    try:
        with open(p) as f: return json.load(f)
    except FileNotFoundError:
        return default

def records(db):
    r = db.get("knots", db)
    return list(r.values()) if isinstance(r, dict) else r

def superscript(s):
    """12n_837 -> 12n₈₃₇ (KnotInfo display)."""
    sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    if "_" in s:
        a, b = s.split("_", 1)
        return a + b.translate(sub)
    return s

def norm_sym(sym_full):
    """Sym token exactly as sort_and_render.fingerprint_from_pflicht: Z/2->Z2, Z/2 + Z/2->Z2_Z2.
    Trivial symmetry group (snappy repr '0' or order 1) -> 'C1' (corpus convention, example
    kf371 = S_24_d175_C1_chiral). Without this mapping '0' fell through the letter check in
    systematic_sname and asymmetric knots (e.g. most braid closures) stayed without an S name
    (Werner 2026-06-24)."""
    if not sym_full or str(sym_full).strip() in ("", "None"):
        return None
    s = re.sub(r"\s+", "", str(sym_full)).replace("+", "_").replace("/", "")
    if s in ("0", "1"):   # snappy: trivial symmetry-group repr -> C1
        return "C1"
    return s or None

def representative(grp):
    """Representative with complete identity fields (c/det/sym)."""
    for r in grp:
        if r.get("c") and r.get("det") and r.get("sym"):
            return r
    return grp[0]

def systematic_sname(rep):
    """S_<c>_d<det>_<sym>_<chir> from corpus fields; c raw as upper bound (Werner 2026-06-21).
    None if not formable (c/det/sym missing or chirality unknown)."""
    c = rep.get("c"); det = rep.get("det")
    sym = norm_sym(rep.get("sym"))
    amphi = rep.get("amphi")
    chir = "amphi" if amphi else ("chiral" if amphi is not None else None)
    if not c or not det or not sym or not chir:
        return None
    # valid group token starts with a letter (Z2, D6, C1, Z2_Z2 ...).
    # 'sym=0' or similar = unclassified/degenerate -> NO S name (-> needs_symmetry).
    if not re.match(r"^[A-Za-z]", sym):
        return None
    return f"S_{c}_d{det}_{sym}_{chir}"

def main():
    check_only = "--check" in sys.argv
    db = load(CORPUS)
    if db is None:
        print("ERROR: corpus_db.json not found", file=sys.stderr); sys.exit(2)
    recs = records(db)
    existing = load(REGISTRY, {"_meta": {}, "names": {}}).get("names", {})

    # group by distinct_knot_id
    groups = {}
    for r in recs:
        did = r.get("distinct_knot_id")
        if not did:
            continue
        groups.setdefault(did, []).append(r)

    names = {}
    needs, auto_known, auto_sculpt, preserved = [], [], [], []
    needs_sym = []                 # novel, but symmetry missing -> cannot be S-named
    cluster_canon = []             # distinct_id already carries a canonical S name -> verbatim
    pending_systematic = {}        # did -> (rep, base_sname, entry)  (raw stems, to be renamed)

    for did, grp in sorted(groups.items()):
        prev = existing.get(did, {})
        prov = sorted({r.get("folder") or r.get("stem") or "?" for r in grp})
        n_real = len(grp)
        ki = next((r.get("knotinfo_id") or r.get("crosscheck_knotinfo") for r in grp
                   if r.get("knotinfo_id") or r.get("crosscheck_knotinfo")), None)
        sculpt = next((r.get("folder") for r in grp if r.get("category") == "own"), None)

        entry = {"display_name": None, "name_type": "provisional", "status": "needs_name",
                 "anchor_sculpture": sculpt, "n_realisierungen": n_real, "provenance": prov}

        # 1) ALWAYS preserve the artist name assigned by Werner (idempotent)
        if prev.get("name_type") == "artist" and prev.get("display_name"):
            entry["display_name"] = prev["display_name"]
            entry["name_type"] = "artist"
            entry["status"] = "assigned"
            entry["anchor_sculpture"] = prev.get("anchor_sculpture", sculpt)
            if prev.get("note"): entry["note"] = prev["note"]
            preserved.append(did)
        # 2) known -> KnotInfo name (auto)
        elif ki:
            entry["display_name"] = superscript(ki)
            entry["name_type"] = "knotinfo"
            entry["status"] = "assigned"
            entry["knotinfo_id"] = ki
            auto_known.append(did)
        # 3) novel + sculpture -> artist name AVAILABLE, but do NOT auto-set as cluster title
        #    (otherwise double title: sculpture page 'nitro' AND cluster page 'nitro'). Werner promotes deliberately.
        elif sculpt:
            entry["display_name"] = None
            entry["name_type"] = "sculpture_anchor"
            entry["status"] = "anchor_available"
            entry["available_name"] = sculpt
            entry["note"] = (f"Artist name of the sculpture '{sculpt}' available — promote as cluster display "
                             f"name? (Werner's decision; otherwise double title with the sculpture page)")
            auto_sculpt.append((did, sculpt))
        # 4) novel without a name
        else:
            # 4a) distinct_id already carries a canonical S name (pipeline-assigned, incl.
            #     _v14/_mut/_sig/_m) -> adopt VERBATIM, NEVER recompute (otherwise reshuffle/drift
            #     that renames the canonical clusters — bug fix 2026-06-21).
            if str(did).startswith("S_"):
                entry["display_name"] = did
                entry["name_type"] = "cluster_canonical"
                entry["status"] = "assigned"
                cluster_canon.append(did)
            # 4b) raw stem -> systematic S designation (if symmetry present), otherwise symmetry work site
            else:
                rep = representative(grp)
                base = systematic_sname(rep)
                if base:
                    # finalised after collision disambiguation further below
                    pending_systematic[did] = (rep, base, entry)
                    names[did] = entry
                    continue
                else:
                    entry["status"] = "needs_symmetry"
                    entry["name_type"] = "provisional"
                    needs_sym.append((did, n_real, prov[0] if prov else "?"))

        names[did] = entry

    # --- Collision disambiguation of the systematic S names ---------------
    # taken = already occupied S names (canonical clusters + KnotInfo not relevant).
    # Raw-stem candidates get the bare name only if the base is NOT occupied by a cluster
    # AND unique; otherwise _v<int(vol)> with tiebreak, always unique against 'taken'.
    systematic = []
    taken = set(cluster_canon)
    by_base = defaultdict(list)
    for did, (rep, base, entry) in pending_systematic.items():
        by_base[base].append(did)
    for base, dids in sorted(by_base.items()):
        if len(dids) == 1 and base not in taken:
            _, _, entry = pending_systematic[dids[0]]
            entry["display_name"] = base
            taken.add(base)
        else:
            # Base occupied (existing cluster = mirror/vol variant) or several candidates
            # -> discriminate by vol (_v<int>), unique against 'taken' (tiebreak letter).
            for did in sorted(dids):
                rep, _, entry = pending_systematic[did]
                vol = rep.get("vol")
                cand = base + (f"_v{int(round(vol))}" if isinstance(vol, (int, float)) else "_v0")
                final = cand; k = 0
                while final in taken:
                    k += 1; final = cand + chr(ord("a") + k - 1)
                entry["display_name"] = final
                taken.add(final)
        for did in dids:
            _, _, entry = pending_systematic[did]
            entry["name_type"] = "systematic"
            entry["status"] = "assigned"
            systematic.append(did)

    # Orphans: registry entries whose distinct_knot_id is no longer in the corpus (drift caught, not silent)
    orphans = [d for d in existing if d not in groups]

    out = {
        "_meta": {
            "purpose": "Display-name layer: distinct_knot_id (canonical) -> display name (free). "
                       "Werner directive 2026-06-16; systematic S names 2026-06-21. "
                       "Identity=key, name=label.",
            "key": "distinct_knot_id",
            "doctrine": "artist name for novel (override); KnotInfo for known; sculpture anchor auto; "
                        "systematic S_<c>_d<det>_<sym>_<chir> when symmetry present (c=upper bound raw). "
                        "Builder NEVER overwrites assigned artist names. Run after /corpus-rebuild.",
            "generated": datetime.date.today().isoformat(),
            "counts": {"total": len(names), "assigned": len(names) - len(needs) - len(needs_sym),
                       "needs_name": len(needs), "needs_symmetry": len(needs_sym),
                       "knotinfo": len(auto_known), "systematic": len(systematic),
                       "cluster_canonical": len(cluster_canon),
                       "artist": sum(1 for e in names.values() if e["name_type"] == "artist"),
                       "sculpture_anchor_available": len(auto_sculpt),
                       "orphans": len(orphans)},
            "orphans": orphans,
        },
        "names": names,
    }

    # Report
    print(f"== Naming registry ==")
    print(f"  distinct knots: {len(names)}")
    print(f"    known (KnotInfo, auto):              {len(auto_known)}")
    print(f"    sculpture anchor (name there, promote?): {len(auto_sculpt)}")
    print(f"    artist name assigned (preserved):    {len(preserved)}")
    print(f"    canonical S name (cluster, verbatim): {len(cluster_canon)}")
    print(f"    systematic S designation (new):      {len(systematic)}")
    print(f"    >>> SYMMETRY MISSING (needs_symmetry): {len(needs_sym)} <<<")
    if needs:
        print(f"    needs_name (unexpected):             {len(needs)}")
    if orphans:
        print(f"  ⚠️ ORPHANS (id no longer in corpus, drift caught): {orphans}")

    if check_only:
        bad = len(needs) > 0 or len(orphans) > 0
        print(f"  --check: {'RED' if bad else 'GREEN'} (needs_name={len(needs)}, orphans={len(orphans)}, "
              f"needs_symmetry={len(needs_sym)} [info, not gate-blocking])")
        sys.exit(1 if bad else 0)

    with open(REGISTRY, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print(f"  written -> {REGISTRY}")

    # Curation queue for Werner (optional artist names — knots now have systematic S names)
    qpath = os.path.join(HERE, "naming_needs_artist_name.txt")
    with open(qpath, "w") as f:
        f.write("# Novel knots — artist name OPTIONAL (co-author decision; override of the systematic S designation).\n")
        f.write("# Eintrag in naming_registry.json: \"<distinct_knot_id>\": {\"display_name\": \"<name>\", \"name_type\": \"artist\", \"status\": \"assigned\"}\n")
        f.write(f"# State {datetime.date.today().isoformat()}: {len(systematic)} with systematic S designation, {len(auto_sculpt)} sculpture anchors.\n\n")
        if auto_sculpt:
            f.write("## SCULPTURE ANCHORS — name available, only promote/confirm (double-title question):\n")
            f.write("#   promoten = in naming_registry.json: name_type:\"artist\", display_name:\"<name>\", status:\"assigned\"\n")
            for did, sc in sorted(auto_sculpt):
                f.write(f"  {did}   <- Skulptur '{sc}'\n")
        f.write("\n## SYSTEMATIC S DESIGNATION assigned — artist name optional (by n_real):\n")
        for did in sorted(systematic, key=lambda d: -names[d]["n_realisierungen"]):
            f.write(f"{names[did]['display_name']}   <- {did}   (n_real={names[did]['n_realisierungen']})\n")
    print(f"  Curation queue -> {qpath}")

    # Symmetry work site: novel without classified symmetry -> no S name possible (Werner 2026-06-21: flag separately)
    spath = os.path.join(HERE, "naming_needs_symmetry.txt")
    with open(spath, "w") as f:
        f.write("# Novel knots WITHOUT classified symmetry (sym=None / degenerate) -> no S designation possible.\n")
        f.write("# Symmetry must first be recomputed (own work site, Werner flag 2026-06-21).\n")
        f.write(f"# {len(needs_sym)} open, state {datetime.date.today().isoformat()}\n\n")
        for did, n, prov in sorted(needs_sym, key=lambda x: -x[1]):
            f.write(f"{did}    (n_real={n}, e.g. {prov})\n")
    print(f"  Symmetry work site -> {spath}")

if __name__ == "__main__":
    main()
