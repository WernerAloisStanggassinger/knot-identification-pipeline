# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""build_corpus_db.py — consolidate all corpus mandatory data into one central DB.

Occasion (2026-05-13 user instruction): doctrine update — the mandatory blind test was the
pipeline validation phase. Now we use every available database
and build a corpus DB of our own that makes all knots generated and
identified so far centrally searchable.

Collects from the four corpus roots:
  erzeugte_knoten/                — machine knots (wave 1-5, polygon sweeps)
  own_knots/                      — user sculptures
  entwicklung_automation/         — pipeline validation knots
  Claude_piplineentwicklung_und_validierung/  — legacy validation knots

per knot (with pflicht_complete.json):
  - identification (path, knot-type prefix if present, category)
  - complete mandatory algebra (jones, alex, conway, det, sig, genus, arf, vol, ...)
  - theta values (sym + num, sym→num fallback)
  - crosscheck status (KnotInfo, snappy, Fremlin)
  - cluster membership (if filed in a c*_d*_* family folder)

Output: pipeline_scripts/corpus_db.json
"""
from __future__ import annotations
import json
import site_config
import re
import sys
from pathlib import Path
from datetime import datetime

ROOT = Path(site_config.ROOT)
CORPUS_ROOTS = [
    "erzeugte_knoten",
    "own_knots",
    "entwicklung_automation",
    "Claude_piplineentwicklung_und_validierung",
]
OUT_PATH = Path(__file__).resolve().parent / "corpus_db.json"


def _category(path: Path) -> str:
    s = str(path)
    if "own_knots" in s: return "own"
    if "erzeugte_knoten" in s: return "machine"
    if "entwicklung_automation" in s: return "validation"
    if "piplineentwicklung_und_validierung" in s: return "validation-legacy"
    return "?"


def _cluster_folder(pfl: Path) -> str | None:
    """Which cluster folder contains this knot (if filed)?

    Stanggassinger classes (S_<c>_d<det>_<sym>_<chir>[_v<vol>][_m]) from Discovery #177.
    Legacy: c<c>_d<det>_<sym>_<chir> (Discovery #152, before the #177 rename).
    """
    parts = pfl.parts
    for p in parts:
        if p.startswith("S_") and "_d" in p:
            return p  # Stanggassinger class
        if p.startswith(("c18_", "c24_", "c30_", "c36_")):
            return p  # legacy cluster (should be empty after #177)
        if p == "bekannt": return "bekannt"
        if p == "pipeline_failures": return "pipeline_failures"
        if p == "singletons": return "singletons"
        if p == "singletons_novelty": return "singletons_novelty"
        if p == "singletons_unknots": return "singletons_unknots"
    return None


def _extract_knotinfo_id(name: str) -> str | None:
    """Extract the KnotInfo type prefix from the folder name (classical convention,
    Discovery #175):
      - c<=10 with K prefix:  K3_1, K3_1m, K8_19, K10_43, K0_1
      - c>=11 without K prefix: 12n_837, 11n_133m, 12n_839m
    """
    # c>=11 (digit + a/n + _ + digits): 12n_837, 11n_133m
    m = re.match(r"(\d+[an]_\d+m?)", name)
    if m:
        return m.group(1)
    # c<=10 (K + digits + _ + digits): K3_1, K10_43, K0_1, K3_1m
    m = re.match(r"(K\d+_\d+m?)", name)
    return m.group(1) if m else None


def _audit_block(d):
    """Newest rigorous hyperbolicity audit block (Disc #538). Prefix search,
    so that later audits (different date) take effect automatically; the youngest wins."""
    keys = sorted(k for k in d if k.startswith("rigorous_hyperbolicity_audit_"))
    return (d.get(keys[-1]) or {}) if keys else {}


def _hyp_flag(d, hyp):
    """is_hyperbolic by priority chain (Disc #538). None = unresolved, never guessed."""
    a = _audit_block(d)
    if a.get("is_hyperbolic") is not None:
        return a["is_hyperbolic"]
    if hyp.get("is_hyperbolic") is not None:
        return hyp["is_hyperbolic"]
    return (d.get("non_hyperbolic_classification") or {}).get("is_hyperbolic")


def _hyp_proof(d, hyp):
    """HOW is the hyperbolicity evidenced? Travels with the flag from #538 on, so that
    'assumed' is never again read as 'found'."""
    a = _audit_block(d)
    if a.get("is_hyperbolic") is not None:
        return a.get("proof")
    if hyp.get("is_hyperbolic") is not None:
        return "pflichtliste_HYPERBOLISCH:" + str(hyp.get("hyperbolicity_reason") or "unspecified")
    dec = (d.get("non_hyperbolic_classification") or {}).get("decision")
    # 'hyperbolic_prim_no_test' carries its own lack of evidence in its name -> report as UNVERIFIED
    if dec in (None, "hyperbolic_prim_no_test"):
        return "UNVERIFIED:" + str(dec)
    return "legacy_gate:" + str(dec)


def _audit_isosig(d):
    """Freshly verified isosig from the audit, if the pipeline did not persist one."""
    ev = (_audit_block(d).get("evidence") or {})
    return ev.get("isosig_fresh") or None


def load_pflicht(pfl: Path):
    """Load mandatory-complete + associated theta/crosscheck files."""
    try:
        d = json.load(open(pfl))
    except Exception as e:
        return None
    pflichtliste = d.get("pflichtliste", {})
    diag = pflichtliste.get("DIAGRAMM", {}) or {}
    alg = pflichtliste.get("ALGEBRA", {}) or {}
    hyp = pflichtliste.get("HYPERBOLISCH", {}) or {}

    folder = pfl.parent
    folder_name = folder.name
    stem = re.sub(r"_pflicht_invarianten.*_complete\.json$", "", pfl.name)

    # Theta data alongside
    theta_sym = folder / f"{stem}_theta.json"
    theta_num = folder / f"{stem}_theta_numeric.json"
    crosscheck = folder / f"{stem}_crosscheck.json"
    torus_check = folder / f"{stem}_torus_check.json"
    # External theta DB (knots.sage c<=14) — load if present
    # Matching is done inline in main(), not per knot

    record = {
        "pflicht_path": str(pfl),
        "folder": folder_name,
        "stem": stem,
        "category": _category(pfl),
        "cluster_folder": _cluster_folder(pfl),
        "knotinfo_id": _extract_knotinfo_id(folder_name),
        # Algebra profile
        "c": diag.get("crossing_number_min"),
        # det fallback (Werner go 2026-06-18, knob 3): Type-a/point-eval knots
        # (c>STAGE5_CMAX) have no symbolic `determinant`, only
        # `determinant_alexander` (Δ(-1)). Otherwise they would come in with det=None
        # → breaks distinctness/naming/Frankenstein gate. det≥1, never 0 → `or` safe.
        "det": alg.get("determinant") or alg.get("determinant_alexander") or alg.get("determinant_sage"),
        "sig": alg.get("signature"),
        "arf": alg.get("arf_invariant"),
        # Disc #816 (20.09.2026): 135 records carried an HFK-exact Seifert genus in
        # ALGEBRA.knot_floer_homology.seifert_genus, but NO ALGEBRA.genus — and thereby fell
        # to genus=None in the corpus. Effect: the genus route of the c bound
        # (c >= 2g + b - 1) and the alternation criterion g == 1/2 bth(Delta) did not take hold
        # there. The fallback is exact, not estimated: HFK DETECTS the
        # genus (Ozsvath-Szabo). [[feedback_a_result_in_a_result_file_is_not_in_the_corpus]]
        "genus": (alg.get("genus")
                  if alg.get("genus") is not None
                  else (alg.get("knot_floer_homology") or {}).get("seifert_genus")),
        "genus_quelle": ("ALGEBRA.genus" if alg.get("genus") is not None
                         else ("HFK.seifert_genus (Ozsvath-Szabo, exakt)"
                               if (alg.get("knot_floer_homology") or {}).get("seifert_genus") is not None
                               else None)),
        "fibered": alg.get("fibered"),
        "alternating": diag.get("alternating"),
        "amphi": hyp.get("is_amphicheiral"),
        "invertibel": hyp.get("is_invertible_knot"),
        "vol": hyp.get("hyperbolic_volume"),
        "sym": hyp.get("symmetry_group"),
        "isosig": hyp.get("canonical_retriangulation") or _audit_isosig(d),
        # Disc #323: authoritative hyperbolicity from the connected-sum gate
        # (enrich_non_hyperbolic). vol alone lies for connected sums (snappy sums
        # the JSJ pieces into one volume) -> label_distinct.is_hyp() reads this field
        # preferentially, fallback vol>0.5 only if it is missing.
        # Disc #538 (2026-07-19): the gate field alone is NOT enough. Its value
        # 'hyperbolic_prim_no_test' literally means "assumed, no test run" and yet sat
        # 335x as is_hyperbolic=true in the corpus (8 of them were connected sums or
        # a satellite). Hence PRIORITY CHAIN + the proof state travels as a field of its own:
        #   1. rigorous_hyperbolicity_audit_* (interval arithmetic / decomposing 2-sphere)
        #   2. pflichtliste.HYPERBOLISCH.is_hyperbolic  (the fresher pipeline block)
        #   3. non_hyperbolic_classification            (legacy gate)
        # RULE: a flag without a carried proof state is not a finding.
        "is_hyperbolic": _hyp_flag(d, hyp),
        "nonhyp_decision": (d.get("non_hyperbolic_classification") or {}).get("decision"),
        "hyperbolicity_proof": _hyp_proof(d, hyp),
        # Polynomials (as strings, for searchable matching)
        "conway": alg.get("conway_polynomial"),
        "jones": alg.get("jones_polynomial"),
        "alex": alg.get("alexander_polynomial"),
        "homfly": alg.get("homfly_polynomial"),
        "kauffman": alg.get("kauffman_polynomial"),
    }

    # Optional-invariants block (phase-2 extension 2026-05-22):
    # HOMFLYPT (a,z) form, genus aliases, bridge, cusp translations.
    # Source: top-level optional_invariants in the _complete.json.
    opt = d.get("optional_invariants", {}) or {}
    if opt.get("homflypt_polynomial"):
        record["homflypt"] = opt.get("homflypt_polynomial")
        record["homflypt_status"] = opt.get("homflypt_polynomial_pipeline_status")
    if opt.get("3_genus") is not None:
        record["3_genus"] = opt.get("3_genus")
    if opt.get("4_genus") is not None:
        record["4_genus"] = opt.get("4_genus")
    if opt.get("bridge_number") is not None:
        record["bridge_number"] = opt.get("bridge_number")
    if opt.get("longitude_translation") is not None:
        record["longitude_translation"] = opt.get("longitude_translation")
    if opt.get("meridian_translation") is not None:
        record["meridian_translation"] = opt.get("meridian_translation")

    # Add theta values
    if theta_sym.exists():
        try:
            ts = json.load(open(theta_sym))
            theta_out = ts.get("theta_output", {}) or {}
            record["theta_full"] = theta_out.get("theta_full")
        except Exception:
            pass
    if theta_num.exists():
        try:
            tn = json.load(open(theta_num))
            record["theta_at_22_7__21_13"] = tn.get("theta")
            record["delta1_at_22_7"] = tn.get("delta1")
        except Exception:
            pass

    # Add torus-knot identification (Discovery #171)
    if torus_check.exists():
        try:
            tc = json.load(open(torus_check))
            if tc.get("is_torus_knot"):
                record["torus_knot"] = True
                record["torus_knot_parameters"] = tc.get("pq")
        except Exception:
            pass

    # Add crosscheck findings
    if crosscheck.exists():
        try:
            cc = json.load(open(crosscheck))
            record["crosscheck_summary"] = cc.get("summary")
            ki = cc.get("knotinfo_match", {})
            ki_status = ki.get("status", "")
            if ki_status.startswith("matched") and ki.get("best_match"):
                record["crosscheck_knotinfo"] = ki.get("best_match")
                if "mirror" in ki_status:
                    record["crosscheck_knotinfo_mirror"] = True
            # Disc #538: a hand-verified `correction` beats EVERYTHING. It arises only
            # when the algorithm demonstrably missed (Jones coincidence without Alexander/
            # genus/primality check -> connected sum read as prime; first mirror
            # 2026-06-11, pulled along corpus-wide 2026-07-19). The folder name then still carries
            # the WRONG table ID -> do NOT derive knotinfo_id from it.
            if ki.get("correction", {}).get("true_identity"):
                record["crosscheck_knotinfo"] = ki["correction"]["true_identity"]
                record["knotinfo_id"] = None
                record["identity_corrected"] = ki["correction"].get("date")
            sn = cc.get("snappy_census_match", {})
            if sn.get("status", "").startswith("matched") and sn.get("best_match"):
                record["crosscheck_snappy"] = sn.get("best_match")
            fr = cc.get("fremlin_match", {})
            if fr.get("status", "").startswith("matched"):
                record["crosscheck_fremlin"] = fr.get("id") or fr.get("name") or fr.get("knot_id")
        except Exception:
            pass

    # Provenance: on-demand-dialled families (novelty_dial, Disc #384) carry a
    # family marker (family_member.json for full records, salvaged_from_failure.json for
    # stage-5+theta salvage). They are genuine corpus knots, but do NOT count into the
    # outreach sweep headline (Werner 2026-07-01: „eigene on-demand-familie" — EN: a family of its own, on demand). Default = sweep.
    record["origin"] = "sweep"
    for marker in ("family_member.json", "salvaged_from_failure.json"):
        mp = folder / marker
        if mp.exists():
            try:
                fam = json.load(open(mp)).get("family")
                if fam:
                    record["origin"] = "on_demand_dial"
                    record["origin_family"] = fam
                    break
            except Exception:
                pass

    return record


def main():
    # Schema validator (pre-commit gate): every _complete.json is checked
    # against the relaxed schema in validate_complete_json.py before it is
    # appended to the DB. HARD errors (structure missing) -> record is NOT ingested,
    # logged to the build output. SOFT warnings -> ingested + logged.
    try:
        from validate_complete_json import validate_file as _validate
    except Exception as _e:
        print(f"  ! Validator not loadable ({_e}), working without pre-commit gate")
        _validate = None
    # Frankenstein gate (Disc #251/#252/#253, wired in 2026-06-04): tier-0
    # impossibility flags (cmin>dt, genus bounds, KFH rank<det,
    # |Delta(-1)|!=det, det!=det_sage, cmin>raw PD). Flag = diagram and
    # algebra layers do not fit together -> record is NOT ingested.
    try:
        from frankenstein_gate import tier0_check_file as _franken
    except Exception as _e:
        print(f"  ! Frankenstein gate not loadable ({_e}), working without")
        _franken = None
    validation_failures: list[tuple[str, list[str]]] = []
    validation_warnings_count = 0
    frankenstein_rejects: list[tuple[str, list[str]]] = []

    db = []
    # Sculpture<->S-cluster cross-link from _sculpture_<stem> symlinks in
    # new_found_knots/<S-cluster>/. Mapping {stem -> S-cluster name}.
    # Background: Werner first draws sculptures by hand with a
    # project name (nitro, ping-pong) as an AutoCAD folder, then extracts
    # the centerline. The sculpture lives mathematically in an S-cluster
    # — own_knots/<stem>/ stays as the original folder, the S-cluster folder
    # gets a _sculpture_<stem> symlink back. In corpus_db the own_knots
    # record is to carry the cluster membership.
    sculpture_cluster_map: dict[str, str] = {}
    nfk_root = ROOT / "erzeugte_knoten" / "new_found_knots"
    if nfk_root.exists():
        for s in nfk_root.glob("S_*/_sculpture_*"):
            if s.is_symlink():
                stem = s.name[len("_sculpture_"):]
                cluster = s.parent.name  # S_24_d961_D6_amphi etc.
                sculpture_cluster_map[stem] = cluster
    # Count exactly ONE _pflicht_invarianten*_complete.json per folder.
    # Background: the canonical-companion pass 2026-05-20 writes a cluster-prefixed
    # variant next to the original in the same folder. Both rglob
    # -> double counting. If both are present: the cluster-prefixed (canonical)
    # wins, the short-name sibling is discarded.
    def _knot_stem(name: str) -> str:
        # "<knot>_pflicht_invarianten_complete.json" -> "<knot>"
        s = name.replace("_pflicht_invarianten_complete.json", "")
        return s.replace("_pflicht_invarianten", "")

    for root in CORPUS_ROOTS:
        root_path = ROOT / root
        if not root_path.exists(): continue
        # Collect per folder
        by_folder: dict[Path, list[Path]] = {}
        for pfl in root_path.rglob("*_pflicht_invarianten*_complete.json"):
            if any(p.startswith("_sculpture_") for p in pfl.parts):
                continue  # sculpture symlink reference, original lies in own_knots/ and is indexed there
            # Every "_"-prefixed DIRECTORY part -> skip: trash/backup/
            # scratch/reexport/pipeline/sweep-work dirs (dedup trashing Disc #233/#234,
            # _n256_pipeline, _reexport, _*_work, _*_backup). Genuine knot folders
            # (own_knots/<name>, S_* clusters, bekannt/, ...) never begin with "_".
            if any(part.startswith("_") for part in pfl.parts[:-1]):
                continue
            # pipeline_failures/ = catch basin of failed runs (stage-9 decision
            # "failure"), NOT corpus stock. A failure folder can nevertheless carry a
            # _complete.json: the core invariants are secured, but stage 6
            # (vol/sym/sig) is missing. Without this filter the build collects it, and the
            # labeler, lacking vol, classes it as "known_nonhyp" — an UNMEASURED
            # knot then counts as a known one. Found 2026-08-23 at
            # pipeline_failures/gitter_d3_k3_s1p5_g0p25/myknot (total 788 -> 789).
            # The way into the corpus leads via salvage (the folder leaves
            # pipeline_failures/ in the process), not via collecting.
            if "pipeline_failures" in pfl.parts[:-1]:
                continue
            by_folder.setdefault(pfl.parent, []).append(pfl)
        # Dedup only for cluster-companion pairs (reorg 2026-05-20): if a
        # file "<cluster>__<short>_pflicht_invarianten_complete.json" and the
        # short variant "<short>_pflicht_invarianten_complete.json" exist in the same
        # folder, discard the cluster-prefixed one — the sister
        # files (_theta_numeric, _crosscheck, _torus_check) are named only under
        # <short>. Other "several _complete in the folder" cases (e.g.
        # own_knots/celtic/ with celtic-for-playing + celtic-trefoil) stay
        # unaffected.
        for folder, files in by_folder.items():
            stems = {f: _knot_stem(f.name) for f in files}
            short_stems = set(stems.values())
            chosen = []
            for f, stem in stems.items():
                if "__" in stem:
                    base = stem.split("__", 1)[1]
                    if base in short_stems and base != stem:
                        continue  # cluster companion, short original exists -> out
                chosen.append(f)
            for c in chosen:
                # Pre-commit gate
                if _validate is not None:
                    v = _validate(c)
                    if not v["ok"]:
                        validation_failures.append((str(c), v["errors"]))
                        continue  # HARD: not ingested
                    if v["warnings"]:
                        validation_warnings_count += len(v["warnings"])
                # Frankenstein gate (tier 0): impossibility flags -> HARD
                if _franken is not None:
                    fl = _franken(c)
                    if fl:
                        frankenstein_rejects.append((str(c), fl))
                        continue  # HARD: not ingested
                rec = load_pflicht(c)
                if rec:
                    if rec.get("category") == "own" and rec.get("cluster_folder") is None:
                        mapped = sculpture_cluster_map.get(rec.get("stem") or "")
                        if mapped:
                            rec["cluster_folder"] = mapped
                            rec["cluster_via_sculpture_symlink"] = True
                    db.append(rec)

    # ---- Duplicate-exclusion gate (Werner 2026-06-19) ----------------------------------
    # Byte-identical .fseries in the SAME cluster = genuine duplicate: the raw m<N> overnight
    # folder stayed lying next to the canonically renamed twin after the promotion, and
    # build_corpus_db counted both (record inflation 625 vs 615, finding 2026-06-19).
    # build_web_export dedups this for serving (_dedup_members, Disc #232) — HERE one
    # stage earlier, so that the count does not inflate in the first place. Keep rule: canonical
    # folder name before raw m<N>_ prefix. GENUINE multiplicity (different parametrisation
    # = different .fseries) stays untouched (multiplicity = evidence). Name conflict
    # (same curve, different crosscheck_knotinfo like cm_k1=12n_533/K10_82) -> warn LOUDLY,
    # do not swallow silently.
    import hashlib
    def _fhash(rec):
        # Read STEM-specifically: a folder can hold SEVERAL knots (e.g.
        # own_knots/celtic/ = celtic-trefoil + celtic-for-playing). glob[0] would hash the
        # same fseries for both -> false duplicate. First <stem>.fseries, otherwise fallback.
        d = Path(rec["pflicht_path"]).parent
        stem = rec.get("stem") or ""
        f = d / f"{stem}.fseries"
        if not f.exists():
            fss = sorted(d.glob("*.fseries"))
            if not fss:
                return None
            f = fss[0]
        return hashlib.md5(f.read_bytes()).hexdigest()
    def _is_raw(rec):
        f = rec.get("folder") or ""
        return bool(re.match(r"^m\d+_", f)) or bool(re.search(r"__m\d+_", f))
    # Key = (cluster, stem, curve-md5): genuine duplicates share ALL three (raw
    # m<N> folder + canonical twin, same stem dm_k3). Different knots (other
    # stem) or different parametrisation (other curve) NEVER collapse.
    _dup_groups: dict[tuple, list] = {}
    _no_hash: list = []
    for rec in db:
        h = _fhash(rec)
        if h is None:
            _no_hash.append(rec); continue
        _dup_groups.setdefault((rec.get("cluster_folder"), rec.get("stem"), h), []).append(rec)
    db = list(_no_hash)
    dedup_dropped: list[tuple] = []
    dedup_name_conflicts: list[tuple] = []
    for (cl, _stem, h), recs in _dup_groups.items():
        if len(recs) == 1:
            db.append(recs[0]); continue
        recs_sorted = sorted(recs, key=lambda r: (_is_raw(r), r.get("folder") or ""))
        keeper = recs_sorted[0]
        db.append(keeper)
        names = {r.get("crosscheck_knotinfo") for r in recs if r.get("crosscheck_knotinfo")}
        if len(names) > 1:
            dedup_name_conflicts.append((cl, keeper.get("folder"), sorted(str(n) for n in names)))
        for r in recs_sorted[1:]:
            dedup_dropped.append((r.get("folder"), keeper.get("folder")))
    if dedup_dropped:
        print(f"  Duplicate exclusion: {len(dedup_dropped)} byte-identical folders collapsed (canonical kept)")
        for dropped, kept in dedup_dropped[:10]:
            print(f"    DROP {dropped}  ==  {kept}")
        if len(dedup_dropped) > 10:
            print(f"    ... +{len(dedup_dropped)-10} more")
    if dedup_name_conflicts:
        print(f"  ⚠️  Duplicate NAME CONFLICT (review needed): {len(dedup_name_conflicts)}")
        for cl, kept, names in dedup_name_conflicts:
            print(f"    {kept}: {names}")

    # External theta DB (knots.sage c<=14, 59,937 knots) — Dror's QR phase 2 (Discovery #173).
    # Match every corpus knot against the external DB via theta_at_22_7__21_13;
    # direct + mirror (-theta) search.
    ext_db_path = ROOT / "external_tools" / "theta_db" / "knots_sage_theta_db.json"
    if ext_db_path.exists():
        from fractions import Fraction
        try:
            ext_raw = json.load(open(ext_db_path))
        except Exception as e:
            ext_raw = {}
            print(f"  ! could not load external theta DB: {e}")
        # Index
        by_t = {}
        for name, e in ext_raw.items():
            s = e.get("theta_at_22_7__21_13")
            if s is None or str(s).strip() in ("0", "", "None"):
                continue
            try:
                if "/" in str(s):
                    n, d_ = str(s).split("/", 1)
                    f = Fraction(int(n.strip()), int(d_.strip()))
                else:
                    f = Fraction(int(str(s).strip()))
            except Exception:
                continue
            by_t.setdefault(f, []).append({"name": name, "c": e.get("c")})
        # Per corpus knot
        for r in db:
            tn = r.get("theta_at_22_7__21_13")
            if tn is None or str(tn) in ("0", "None"):
                continue
            try:
                if "/" in str(tn):
                    n, d_ = str(tn).split("/", 1)
                    f = Fraction(int(n.strip()), int(d_.strip()))
                else:
                    f = Fraction(int(str(tn).strip()))
            except Exception:
                continue
            direct = by_t.get(f, [])
            mirror = by_t.get(-f, [])
            if direct:
                r["external_theta_match_direct"] = [m["name"] for m in direct[:5]]
            if mirror:
                r["external_theta_match_mirror"] = [m["name"] for m in mirror[:5]]

    # Partition: genuine corpus (machine+own) vs validation knots (separate, NOT in the total).
    # Validation copies of a genuine corpus knot (same stem, e.g. double-3=own +
    # validation) are redundant -> drop. Dedup validation among themselves by stem
    # (orbit/saturn lie in validation AND validation-legacy).
    real = [r for r in db if r["category"] in ("machine", "own")]
    real_stems = set(r.get("stem") for r in real)
    validation_all = [r for r in db if r["category"] in ("validation", "validation-legacy")]
    seen_v, validation = set(), []
    for r in validation_all:
        s = r.get("stem")
        if s in real_stems:   # genuine sculpture exists -> validation copy redundant
            continue
        if s in seen_v:       # dedup validation <-> validation-legacy
            continue
        seen_v.add(s); validation.append(r)

    # --- Pull certificates along from the mandatory JSONs (2026-08-13) -------------------
    # WHY HERE: this script writes corpus_db.json FRESH (json.dump, no merge).
    # Fields that stood only by hand in corpus_db were silently gone at the next /corpus-rebuild
    # — exactly that was noticed on 13.08., after c_exact/vol_valid had been
    # entered there. The source of truth is therefore the certificates IN THE RECORD:
    #   _CROSSING_EXACTNESS_*      -> c_exact      (span_t V(K)=cmin  OR reduced alternating)
    #   _HYPERBOLICITY_CERT_*      -> vol_valid + hyperbolicity_proof
    # That way the corpus heals itself on every rebuild instead of rotting.
    def _certificates(path_):
        if not path_ or not Path(path_).exists():
            return (None, None, None)
        try:
            d = json.load(open(path_))
        except Exception:
            return (None, None, None)
        cx = [v for k, v in d.items() if k.startswith("_CROSSING_EXACTNESS") and isinstance(v, dict)]
        c_exact = any(str(v.get("verdikt", "")).upper() == "EXAKT" for v in cx) if cx else None
        hy = [v for k, v in d.items() if k.startswith("_HYPERBOLICITY_CERT") and isinstance(v, dict)]
        if hy:
            proven = any(str(v.get("verdikt", "")).upper().startswith("HYPERBOLISCH BEWIESEN") for v in hy)
            proof = "hyperbolic_rigorous_verify_20260813" if proven else None
            return (c_exact, bool(proven), proof)
        return (c_exact, None, None)

    n_cx = n_vv = 0
    for r in real:
        c_exact, vol_valid, proof = _certificates(r.get("pflicht_path"))
        if c_exact:
            r["c_exact"] = True; n_cx += 1
        if vol_valid is not None:
            r["vol_valid"] = vol_valid; n_vv += 1
            if proof:
                r["hyperbolicity_proof"] = proof
    print("[corpus_db] Zertifikate: c_exact=%d, vol_valid=%d" % (n_cx, n_vv))

    # Stats — meta ONLY from the genuine corpus
    from collections import Counter
    cat_counts = Counter(r["category"] for r in real)
    cluster_counts = Counter(r["cluster_folder"] for r in real if r["cluster_folder"])

    meta = {
        "built_at": datetime.now().isoformat(),
        "total_knots": len(real),
        "n_knots": len(real),
        "distinct_stems": len(set(r.get("stem") for r in real)),
        "categories": dict(cat_counts),
        "cluster_folders": dict(sorted(cluster_counts.items())),
        "validation_knots_count": len(validation),
    }

    out = {"meta": meta, "knots": real, "validation_knots": validation}
    json.dump(out, open(OUT_PATH, "w"), indent=1)
    db = real  # following prints/stats refer to the genuine corpus

    print(f"Wrote {OUT_PATH}")
    print(f"  Total: {len(db)} knots")
    print(f"  Categories: {dict(cat_counts)}")
    print(f"  Cluster: {dict(cluster_counts)}")
    # KnotInfo hits
    n_ki = sum(1 for r in db if r.get("crosscheck_knotinfo"))
    n_sn = sum(1 for r in db if r.get("crosscheck_snappy"))
    n_tk = sum(1 for r in db if r.get("torus_knot"))
    n_ext_theta = sum(1 for r in db
                      if r.get("external_theta_match_direct")
                      or r.get("external_theta_match_mirror"))
    print(f"  KnotInfo hits: {n_ki}, snappy hits: {n_sn}, torus knots: {n_tk}, theta-DB match: {n_ext_theta}")
    # Theta-Coverage
    n_sym = sum(1 for r in db if r.get("theta_full"))
    n_num = sum(1 for r in db if r.get("theta_at_22_7__21_13") is not None)
    print(f"  Theta-Coverage: {n_sym} sym, {n_num} num")
    # Validator report
    print(f"  Pre-Commit-Gate: {len(validation_failures)} HARD-rejected, "
          f"{validation_warnings_count} SOFT-warnings, "
          f"{len(frankenstein_rejects)} Frankenstein-rejected")
    if frankenstein_rejects:
        print(f"  FRANKENSTEIN-rejected files (NOT ingested — diagram/algebra"
              f" = two objects? Clarify, never overwrite):")
        for path, flags in frankenstein_rejects[:10]:
            print(f"    - {path}")
            for f in flags[:4]:
                print(f"        🚩 {f}")
        if len(frankenstein_rejects) > 10:
            print(f"    ... +{len(frankenstein_rejects)-10} more")
    if validation_failures:
        print(f"  HARD-rejected files (NOT ingested):")
        for path, errs in validation_failures[:10]:
            print(f"    - {path}")
            for e in errs[:3]:
                print(f"        {e}")
        if len(validation_failures) > 10:
            print(f"    ... +{len(validation_failures)-10} more")


if __name__ == "__main__":
    main()
