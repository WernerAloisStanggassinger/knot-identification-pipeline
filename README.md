# knot-identification-pipeline

Identifies the knot type of a closed space curve — a sculpture's centreline or a machine-generated curve — directly from its
geometry, given as a Fourier series (`.fseries`). It was built for the stainless-steel knot sculptures of
**Werner Alois Stanggassinger** (sculptor, Wasserburg am Inn, Bavaria) and grew into the tool behind the knot corpus on
https://knot-structures.stainlesssteel4u.de/ (rotatable 3D views, photos and the invariants of every piece).

Werner Alois Stanggassinger in collaboration with Claude Code (Anthropic): the choice of problem, the direction of the work and the
computing hardware are due to the author; the code was written by Claude Code, which orchestrates the computations, working with the
author. The collaboration is stated as a fact, not a disclaimer.
Every invariant is computed by a program, never by a language model: by established programs (among them SnapPy, spherogram,
SageMath, Regina, KnotJob, knot_floer_homology, knotpy, pyknotid and Bar-Natan–van der Veen's Theta) or by code of ours in this
repository. Our own code that computes something itself, rather than passing the work to one of those programs, includes: the
diagram from the curve (stages 1–3, and `polygon_resolve.py`); fallback routes for invariants (`alex_point.py`, the HOMFLY route
`hecke_homfly_fast.py` / `homfly_modular.py` / `homfly_interpolation.py`, checked by `homfly_calibration.py` and
`homfly_crosscheck.py`); the conversion of a PD code into the input of the Theta program (`pd_to_theta_input.py`); the bounds of
stage 7 (braid index from the HOMFLY polynomial by Morton–Franks–Williams, bridge number); the verdict modules `fox_milnor.py`,
`slice_verdict.py`, `four_genus_verdict.py`; the table-free crossing-number bounds `span_lower_bound.py`, `yamada_c_bound.py`,
`turaev_bound.py` (from invariants that the established programs computed); the corpus checks (`frankenstein_gate.py`,
`corpus_consistency_gate.py`, `alexander_normal_form.py`, `theta_mirror.py`, `identify_torus_knot.py`, `composite_knot_check.py`);
the 6D geometry (`sixd_*.py`, `canonical_flow.py`); the clearance functions `clearance_*.py` (used by `build_curve_json.py`, among
the 43 modules); and the tools in `research/` (each folder's README says which
parts are ours, why we built them, and how they were checked). `certify_hyperbolicity.py` is not on this list: the certificate is
SnapPy's `verify_hyperbolicity()`; our part is the order in which it is tried on different triangulations.

## What it does
`pipeline/pipeline_run.sh`, run inside a folder that holds one `.fseries`, goes through these stages:
- **Preflight** — curve sanity, π-symmetry, target crossing count.
1. **Scanner** (GPU) — a camera flight of 131 072 view directions around the curve; every crossing of every projection.
2. **Tracker** — follows each crossing across the views (density clustering in curve parameters).
   `pipeline/worldline_tracker.py` is a second, prototype tracker used by the authors' resolution sweeps: it takes one complete
   projection frame (fewest crossings) with over/under from that frame alone, and uses the other frames only to confirm which
   crossings persist. Not called by `pipeline_run.sh`.
3. **Resolve** — chooses views and builds a signed Gauss code / PD code.
4. **Aggressive simplify** (4a: merge) — Regina and spherogram simplification over many frames; keeps the smallest diagram.
5. **Mandatory invariants** — Alexander, Jones, HOMFLY, Khovanov, signature, determinant, hyperbolic volume and more.
6. **Extra invariants** — Kauffman polynomial and further invariants.
7. **Enrichments** — bridge and braid bounds, unknotting bounds, length spectrum, non-hyperbolic cases (not fatal if one fails).
8. **Summary**; 8a **Theta** (Bar-Natan–van der Veen Θ); 8b **Cross-check** — only *after* identification, comparison with
   KnotInfo, Fremlin and the SnapPy census.
9. **Sort and render** — folder placement and a Blender render (`SKIP_SORT_RENDER=1` skips it).

**6D geometry cross-check** (separate from the stages above). `pipeline/sixd_scan.py` measures what the diagram throws away:
curvature, torsion, writhe and average crossing number of the curve itself, from analytic Fourier derivatives.
`pipeline/canonical_flow.py` first relaxes the curve towards a tight, reproducible shape (the SONO flow, after Pierański's Shrink-On-No-Overlaps algorithm: P. Pierański, in *Ideal Knots*, eds. A. Stasiak,
V. Katritch, L. H. Kauffman, Series on Knots and Everything 19, World Scientific 1998, pp. 20–41), made order-independent
on CPU and GPU by `sixd_robust.py` (`sixd_gpu.py` = optional CuPy fast path); `calibrate_full_corpus_robust.py` runs that over a
whole corpus and survives dying workers (`test_run_jobs_workertod.py` checks it). `sixd_clean_analysis.py` holds the calibrated
validity filter and `sixd_crosscheck.py` compares fingerprints with the algebraic identification. It is a counter-check: it
separates most distinct knots but is blind to mutants, so the invariants decide. The flow takes minutes per knot on a GPU.

**Corpus layer** (what turns many identified knots into the knot corpus). `build_corpus_db.py` collects every record into
`corpus_db.json`, refusing records that fail `validate_complete_json.py` (schema) or tier 0 of `frankenstein_gate.py`
(consistency of the stored diagram data, the algebra and the knot Floer data with each other, for example |Δ(−1)| = det and
deg Δ ≤ 2g; a check whose fields are missing abstains; a record that cannot be read is refused as UNVERIFIED). The gate's tier 1 (the
determinant of the stored PD code, recomputed by Regina, against the stored determinant) and tier 2 (same knot from the curve) run
only when called separately; tier 1 has three outcomes since 3 November 2026 — passed (or nothing to compare: no stored PD code or no determinant), flagged,
or UNVERIFIED when its Regina
subprocess fails (before, such a failure counted as passed) or when the stored PD text is not a complete PD code (every number in
a 4-tuple, labels 1..2n, each exactly twice — all 839 stored PD codes of the authors' corpus are); `frankenstein_gate.py --scan --pd` exits with 1 on any flagged or
UNVERIFIED record (a missing or unreadable record file counts as UNVERIFIED). Regina runs in the Python named by `PYTHON_BASE` (`site_config.py`; under `portable` the current Python), which
must have Regina installed, otherwise every record with a stored PD code and a determinant is UNVERIFIED; `corpus_consistency_gate.py` checks the corpus against unconditional theorems.
`label_distinct_knots.py` decides which records are the same knot and labels each distinct knot (known, novel hyperbolic,
non-hyperbolic, unresolved), using the table-free lower bounds for the crossing number in `span_lower_bound.py`,
`yamada_c_bound.py` and `turaev_bound.py`; `certify_hyperbolicity.py` proves hyperbolicity with interval arithmetic.
The JSON field names these modules write stay German (for example `pflicht_path`, `n_realisierungen`), because they are the field
names of the published corpus files. The records they read come in two forms: the authors' records carry German file and field
names (`<knot>_pflicht_invarianten_complete.json`, `pflichtliste`, `DIAGRAMM`, `HYPERBOLISCH`, …), the pipeline in this release
writes English ones (`<knot>_mandatory_invariants_complete.json`, `mandatory_list`, `DIAGRAM`, `HYPERBOLIC`, …). Since 3 November
2026 `pipeline/record_names.py` lets the corpus layer read both: it finds both file names, and it maps the keys of an English record
back to the German names with `tools/translation_table.json` in reverse (string values stay as written); a German record is passed on
unchanged. (Until then the released
corpus layer could not read the released pipeline's output.) Checked on 3 November 2026: on the authors' 843 records the corpus
layer with this change writes the same `corpus_db.json` as without it (only the build timestamp differs); the same 843 records
translated into the English names with `tools/translation_table.json` give the same `corpus_db.json` again (same records, paths and
file names aside; the order of the records follows the file system); and on the record the released pipeline wrote for
`examples/trefoil.fseries`, `build_corpus_db.py` finds and reads it (c = 3, det = 3, KnotInfo match 3_1), `validate_complete_json.py`
passes it, tiers 0 and 1 of `frankenstein_gate.py` raise no flag (tier 1, with `PYTHON_BASE` set to a Python that has Regina, recomputed
det = 3; re-measured 3 November 2026), `identify_torus_knot.py` gives T(2,3)
and `composite_knot_check.py` runs through.
Three English names also occur untranslated in the authors' records (`proof`, `route`, `exact`) and are not mapped back.
`frankenstein_gate.py --guard-new` (tier 2) needs a module of the authors' that is not in this repository; tiers 0 and 1 run as released.

**Blind test.** During identification no knot table is consulted: `pipeline_seal.py` blocks the table lookups of SnapPy,
spherogram and Sage. Tables are read only in stage 8, to compare.

## Research tools
`research/` holds tools we built for questions beyond the pipeline. Each folder has its own README: what it does, why we built
it instead of using an established program, and how it was checked.
- `research/theta_census_gpu/` — θ of every knot of a census file and of its mirror, exactly modulo a prime on the GPU, as a
  filter that misses no θ-zero (given a correct kernel and diagram conversion, both calibrated); plus an attempt to certify,
  with verified computation (needs SageMath), whether the θ-zeros of hyperbolic knots are chiral or amphichiral.
- `research/ropelength/` — polygonal ropelength checker, a GPU tightener, and an upper bound on the smooth ropelength from a
  polygon (the corner-rounding of Ashton–Cantarella–Piatek–Rawdon).
- `research/volume_census/` — SnapPy volume of every knot of a census file, resumable.
- `research/mutation_search/` — search for a Conway mutation between two given census knots, on their census diagrams and on
  many larger diagrams.
- `research/chirality_witnesses/` — ten chirality witnesses (linking forms of cyclic branched covers, Tristram–Levine, Chern–Simons,
  HOMFLY, knot Floer, Khovanov, writhe and flype tests for alternating knots) and a chirality / amphichirality certificate for
  hyperbolic knots, each with counter-cases.
- `research/jones_gpu/` — an independent second implementation of the Jones polynomial at a point (exact mod p, a tensor network
  on the PD code, CPU reference and GPU batch); slower than Regina's exact Jones polynomial, meant as a cross-check.
- `research/theta_collision/` — which census knots with hyperbolic volumes equal within a tolerance (for example mutants) θ
  tells apart: groups by volume, θ values from `theta_census_gpu/theta_gpu_census.py --values`; a pair counts as separated when
  the values modulo p prove A ≠ B and A ≠ mirror B as oriented knots (given a correct kernel and diagram conversion).
- `research/whitehead_double/` — HOMFLY polynomial of the untwisted Whitehead doubles of a knot (both clasps) and a comparison of
  two knots by it that ignores mirror images; different values prove different knots up to mirror image, equal values prove nothing.
- `research/torsion_polynomial/` — Dunfield–Friedl–Jackson Conjecture 1.7: the hyperbolic torsion polynomial (SnapPy) against genus
  and fiberedness (knot Floer homology) for every hyperbolic knot of a census file, with precision escalation and calibration.
- `research/kauffman_polynomial/` — the two-variable Kauffman polynomial F(a, z) from a PD code, pure Python, exact integers
  (Kauffman's descending-diagram recursion); equal to KnotInfo's F on all 249 prime knots with at most 10 crossings.

## Install
- An NVIDIA GPU with CUDA (scanner and tracker use CuPy).
- Four conda environments with these names, because the stages call them by name: `gnottheorie_env` (scanner, tracker, most
  stages), `base` (Regina, pyknotid, knotpy), `snappy_env` (SnapPy/spherogram subprocesses), `sage_stable` (SageMath).
  `environments/<name>.txt` lists the Python version and every package of each environment as used by the authors (measured
  22.10.2026). These are today's versions: the published results were computed with the versions current at the time,
  since the tools are kept up to date.
- Third-party programs and data that are not included: see `THIRD_PARTY.md`.
- `python3 pipeline/site_config.py` shows where the pipeline will look for everything (`KNOT_PIPELINE_<NAME>` overrides one location).

## Run
```
python3 pipeline/make_camera_trajectory.py --n 131072 --out data/camera_131072     # once
cd my_knot/                                   # a folder with exactly one .fseries
bash /path/to/knot-identification-pipeline/pipeline/pipeline_run.sh
```
`examples/trefoil.fseries` is a small example curve. On the authors' machine (RTX 3070) a typical knot takes about six minutes.
The main result is `<knot>_mandatory_invariants_complete.json` plus `<knot>_summary.md`.

## How this release relates to the code that produced the results
- **Two code trees, on purpose.** "The authors' copy" in this README is a living tree, not only a historical state: the authors
  compute their corpus with their own working copy of the same programs, whose module, file and field names are still German
  (`pflicht_invarianten.py`, `pflichtliste`, `<knot>_pflicht_invarianten_complete.json`). It keeps those names so that the 843
  records already computed — and the data published from them (corpus files, website downloads, Zenodo deposits) — stay readable
  without conversion; those published data are an archive and stay as they are. This repository is the English, public form of
  the same programs. The working copy is not published separately because this code is proven to be the same program up to a
  consistent renaming (22.10. and 27.10.2026, below) plus the changes listed below; every program change since is either made in
  both trees or described below as a change of this release only. The corpus layer here reads records of both kinds (`record_names.py`).
- **Limits found by outside reviews and fixed or stated on 3 November 2026.**
  - *Locations.* Several modules built paths from `ROOT` and the authors' folder layout instead of the names in `site_config.py`;
    with `KNOT_PIPELINE_PROFILE=authors` that coincides, with the default profile `portable` it did not (the corpus database, the
    files next to it that `label_distinct_knots.py` reads, the 6D output folders, and the folder the pipeline writes its knots to).
    Now `sixd_scan.py`, `canonical_flow.py` and `label_distinct_knots.py` use `PIPELINE_DIR` and `SIXD_DIR`, and
    `build_corpus_db.py` and `validate_complete_json.py` also scan `GENERATED_KNOTS_DIR` and count its records as machine records.
    Two new names in `site_config.py`: `BACKUP_DIR` (where `label_distinct_knots.py --apply` copies the database first; `_backups/`
    in the repository, ignored by git, like `output/`) and `THETA_DB` (the external θ database that `build_corpus_db.py` matches
    against; it is derived from Bar-Natan's `knots.sage` (knots up to 14 crossings) and not distributed, so under `portable` the match is skipped
    with a message); `build_corpus_db.py` takes `new_found_knots/` from `GENERATED_KNOTS_DIR`. These path changes are made only in
    this release (the authors' copy keeps its fixed paths, which are the `authors` profile); against the tags
    `proven-equivalent-sixd-20261029` and `proven-equivalent-corpus-20261029` they show as program changes.
    `frankenstein_gate.py` still reads two places of the authors' layout: the module that `--guard-new` needs (above) and the
    result files of the authors' audit of 4 June 2026 for its one-off option `--seed-registry`; neither exists under `portable`.
    Under `authors` every path is the same as before. Checked: with `portable`, the trefoil record written by the released pipeline
    into `output/` goes through `build_corpus_db.py` (c = 3, det = 3, torus knot T(2,3), KnotInfo match 3_1) and
    `label_distinct_knots.py` (1 record, 1 distinct knot). `tools/path_equivalence.py` checks that `ROOT` was replaced cleanly; it does
    not check that both profiles point to the same place.
  - *Certificates in the corpus.* `build_corpus_db.py` takes `c_exact`, `vol_valid` and `hyperbolicity_proof` from blocks
    `_CROSSING_EXACTNESS_*` and `_HYPERBOLICITY_CERT_*` inside a record. The authors write these blocks with a tool of theirs that is
    **not** in this repository. `certify_hyperbolicity.py` writes its results to a list of its own (`--out`), which no module reads.
    On records made only with this repository these three fields stay empty.
  - *No-match verdict and stage reports* (found in a trefoil run of v1.7.0 without SageMath). When the signature was missing,
    `pipeline_crosscheck.py` skipped the KnotInfo lookup and still ended with "novelty candidate", for the trefoil as well; for
    c ≤ 13 there was no hyperbolicity check on that path either. Now a skipped lookup or a missing volume gives "undecided", a
    non-hyperbolic knot (volume ≤ 0.5) "known (non-hyperbolic family)", and a hyperbolic knot without a match only "no DB match
    (hyperbolic)": whether it is a novelty candidate is decided by `label_distinct_knots.py` against the enumeration limits, as
    before. No corpus label comes from this field. `pipeline_run.sh` printed "done" for stages 8a and 8b even after they had failed;
    now it does so only when the script succeeded and its JSON file exists, counts the failed stages and names their number in
    the final line (the exit code stays 0: stages 8a and 8b are not fatal). Both fixes were made in the authors' copy as well.
    A later review (3 November 2026) found two more status errors, fixed in both trees: a SnapPy `identify()` candidate whose
    volume differs from the record's was passed on as a match (`matched (vol-mismatch)` counted as `matched` in
    `pipeline_crosscheck.py`, `build_corpus_db.py` and `sort_and_render.py`); now it is no match (where no other database matched, the verdict is "undecided")
    (none of the authors' 697 crosscheck files contains such a case). And with `--skip-snappy`, or a failed SnapPy subprocess, a
    hyperbolic knot without a match still got "no DB match (hyperbolic)"; now "undecided", as also for census manifolds of equal
    volume without identification. In `frankenstein_gate.py` an unreadable record no longer passes tier 0 and a failed Regina call
    no longer passes tier 1 (both UNVERIFIED, above). Checked: on the authors' 843 records `build_corpus_db.py` with these changes
    writes every field it builds exactly as in the published `corpus_db.json`.
  - *Invariant code, strictness and what is certified* (same review). Fixed in both trees: `alex_point.py` now takes as many
    primes as the CRT lift needs for its coefficient bound (4^n for n crossings: 4 primes up to 61 crossings, more beyond; before,
    always 4 — the largest coefficient in the authors' corpus has 74 bits, the 4-prime lift is unique up to about 123 bits, and on
    the largest diagrams the result is unchanged); `jones_at_i` (determinant from the Jones polynomial) gives no value for an
    unreadable term or a non-square instead of skipping the term or rounding; `fox_milnor.py` refuses an expression whose
    denominator is not ±t^k (it is not a Laurent polynomial); `slice_verdict.py` says when the full Fox–Milnor test could not run;
    `homfly_crosscheck.py` compares the determinant from HOMFLY exactly as an integer and names the two cases that give none (a non-square,
    a negative power of z; before, a floating square root compared with a tolerance). Checked on the authors' corpus: the new `jones_at_i` gives the same value as
    before on all 600 stored Jones polynomials (each equal to the record's determinant), the new Fox–Milnor parser the same verdict
    on all 841 stored Alexander polynomials. **Not certified, by design:** our HOMFLY interpolation (`homfly_modular.py`,
    `homfly_interpolation.py`) stops adding primes when one more prime leaves the CRT result unchanged (`stable`), and takes its
    exponent box from bounds measured on 178 knots against Regina; both are then checked against fresh evaluation points, the
    determinant and the Alexander polynomial. That makes a wrong polynomial very unlikely but is not a proof: `stable` is not a
    coefficient-bound certificate, and the box has no general derivation here.
    Only in this release: stage 8a says plainly that `Theta_numeric.sage` is missing (it is not distributed; `THETA_SYMBOLIC=1`
    uses `Theta.sage`); the summary leaves `chern_simons` empty for a non-hyperbolic knot (the record keeps SnapPy's raw number);
    three remaining German messages are English ("Consistency", "Certificates", and "determined" in
    `unknotting_number_assignment_basis`); `build_corpus_db.py` prints its help for `--help` instead of starting a build (without arguments it builds as before); knotpy, which is not on PyPI,
    is pinned to its GitHub commit in `environments/base.txt`.
  - *Names left German.* The options `--stand` and `--datum` of `unknot_label.py` and `--datum` of `homfly_interpolation.py` were not
    translated, nor `--wand` of `homfly_crosscheck.py`; nor were a number of identifiers, for example `load_pflicht`
    (`build_corpus_db.py`), `parse_pflicht` and `pflicht_json(s)` (`composite_knot_check.py`, `identify_torus_knot.py`), `laden`
    (`span_lower_bound.py`, `yamada_c_bound.py`, `unknot_label.py`), `_regina_einer` (`frankenstein_gate.py`), `erg`
    (`homfly_crosscheck.py`), and in `regina_bridge.py` the argument `koerper` and names such as `auf`, `strikt`, `kehrwert`,
    `tetraeder` inside the Regina code it runs. A few such values reach the output: the cause `homfly_unparsbar`, the verdicts
    `trifft`, `weicht_ab`, `kein_sollwert`, `nicht_gerechnet` of `frankenstein_gate.py` and the key `vol_abweichung` of
    `certify_hyperbolicity.py`. Renaming them would change programs proven equivalent. `tools/map_equivalence.py` checks that a renaming is the same everywhere; it does not check that no German name is left.
  - *Self-tests.* `label_distinct_knots.py --self-check` needs a full corpus and fails without one; `certify_hyperbolicity.py
    --self-check` needs SnapPy inside SageMath; `scan_triage.py` needs pyarrow, `unknot_second_opinion.py` pyknotid.
- The 43 original modules of `pipeline/` (the stages above; the 6D modules, the corpus layer and `worldline_tracker.py` came later,
  see below) were translated from German in two steps. 22.10.2026: comments, docstrings and messages. 23.10.2026:
  also identifiers, module and file names, command-line options (three were missed, see above), environment variables, JSON field names, category labels and
  the remaining texts. `tools/map_equivalence.py` proved, on 23.10. and again on 27.10.2026, for every Python file of the 43 that the
  translated version is the program proven equivalent to the originals (git tag `proven-equivalent-20261027`, see the next points)
  up to one consistent renaming, and wrote that renaming to `tools/translation_table.json` (identifiers, keys, file renames, free
  texts). For `.sh` and `.json` files it only checks that no renamed name survives in them; `pipeline_run.sh` is not compared
  structurally. **On today's tree the gate is red against that tag**, and it is meant to be read with the changes since:
  `python3 -I tools/map_equivalence.py proven-equivalent-20261027 --rename pipeline/pflicht_invarianten.py=pipeline/mandatory_invariants.py …`
  (the seven renames listed under `files` in `tools/translation_table.json`; without them it also reports those seven files as
  missing) reports 49 problems on the tree of this release (3 November 2026): 25 are the German keys that the corpus layer keeps on
  purpose (`pflichtliste`, `DIAGRAMM`, `HYPERBOLISCH`, file-name patterns, and `pflicht_json` / `source_pflicht_json` in the output of
  `identify_torus_knot.py`), two of them in `record_names.py`, which maps English records onto exactly these names; 3 more German
  keys inside the corpus layer (`quelle` in `label_distinct_knots.py`, `fehler` twice in `regina_bridge.py`, which only the corpus
  layer reads); and 6 are program changes after the tag in two files (`pd_to_theta_input.py`, kink guard, 3.11.; `site_config.py`, the
  name `SIXD_DIR`, 29.10., and the names `BACKUP_DIR` and `THETA_DB`, 3.11.); 15 more are the fixes of 3.11. after outside reviews, in eight files (`pipeline_crosscheck.py`, `summary_writer.py`,
  `sort_and_render.py`, `alex_point.py`, `fox_milnor.py`, `homfly_crosscheck.py`, `mandatory_invariants.py`, `slice_verdict.py`),
  all described above; four of the eight messages for `homfly_crosscheck.py` are about identifiers and come from its one changed
  comparison line (`ok = (dh == det_rec)` instead of a tolerance). The gate has been red against this tag since the corpus layer came in (release v1.1.0). The files added
  on 29.10. have their own proofs, described below. The deliberate counter-cases the gates were tested with (25 and 22,
  named in the commit messages) are not in this repository. Since 3 November 2026 the gate refuses an unknown or empty revision instead of reporting it as
  equivalent and writing an empty translation table.
  The originals were proven equivalent on 22.10.2026 (tag `proven-equivalent`). On 27.10.2026 every one of the 43 originals was
  compared again with that state: 42 differ only in texts (their comments and messages were translated in place on 25.10.), one
  changed as a program: `build_curve_json.py` labels each curve with its licence (sculpture curves all rights reserved, machine-knot
  curves CC-BY-NC-4.0). Tag `proven-equivalent-20261027` = `proven-equivalent` plus exactly that change. The gate was
  counter-tested with 25 deliberate changes it must catch or pass; the built-in self-tests of six modules give the same
  results before and after (45 of 45 checks).
- **3 November 2026, one program change in the pipeline:** `pipeline/pd_to_theta_input.py` now refuses a PD code with a
  Reidemeister-I kink (a crossing whose 4-tuple repeats a label) with `ValueError`, because θ came out wrong on such diagrams (99 of 160 enlarged test diagrams, all of them with a kink of the
  kind described below; details in `research/theta_census_gpu/README.md`). Its self-test also checks the refusal. Otherwise the module is unchanged: it is EQUIVALENT
  (`tools/strict_equivalence.py`) to the authors' copy, which carries the same guard. `theta_crosscheck.py` and
  `theta_numeric_crosscheck.py` call it without catching the error, so on a kinked diagram they now stop with that message
  instead of possibly writing a wrong θ; in a pipeline run that is a failed stage 8a, which is reported and counted, and the run
  goes on (stages 8a and 8b are not fatal). The cause was found the same day by an outside reviewer: the θ kernels set the matrix entries
  of a crossing with `=` instead of adding them (`+=`, as `Theta.m` does), which goes wrong exactly at a crossing with i = j + 1, i.e.
  at a kink; the conversion itself is right (details and tests in `research/theta_census_gpu/README.md`, where `theta_modp.py` is
  fixed). Both Sage kernels that stage 8a can call set the entries the same way: `Theta.sage`, the Sage code of
  Bar-Natan and van der Veen (`THETA_SYMBOLIC=1`; tested on enlarged diagrams of 4_1: the 3 with a crossing i = j + 1 wrong, the 3
  without right), and `Theta_numeric.sage` (default; our modification of `Theta.sage`, not in this repository), which took over its
  matrix set-up (30 enlarged diagrams: the 15 with i = j + 1 wrong; after changing it to `+=` on 3.11.2026, none). So the guard stays,
  and it protects both routes. The pipeline's diagrams are simplified first: none of the 834 PD diagrams of the 843 records of
  the authors' corpus has a kink, nor do the 5 diagrams of its validation knots (checked 3 November 2026; the other 9 records
  are unknots without a diagram; 834 + 5 = the 839 PD codes named in `research/theta_census_gpu/README.md`). The research tools
  that use the same conversion on Burton's census diagrams are covered in `research/theta_census_gpu/README.md`: census diagrams are
  minimal, so they carry no kink.
- **Results the authors published before this release (knot corpus, website downloads, Zenodo deposits) use the German
  names.** `tools/translation_table.json` maps every name this code reads or writes, for example `pflichtliste` → `mandatory_list`,
  `pflicht_invarianten_complete.json` → `mandatory_invariants_complete.json`, `beweis` → `proof`, `grund` → `cause`. In
  `pipeline/proof_register.json` additionally `belegt_von` → `receipted_by` and `_doku` → `_doc`.
- Comments refer to the authors' internal notes (for example "Discovery #N", the authors' error log, notes and correspondence
  folders); these notes are not public. Where a comment states a result, the code or a test in this repository is what supports it.
- A few names stay German on purpose because they name things outside the code: the authors' folder names
  (`erzeugte_knoten/`, `bekannt/`), the material `Edelstahl` in the Blender template, and historical file names.
- All absolute paths moved into `pipeline/site_config.py`. With `KNOT_PIPELINE_PROFILE=authors` every location is the authors' own.
- `tools/strict_equivalence.py` and `tools/path_equivalence.py` prove, file by file, that each of the 43 released modules is the same
  program as the authors' copy at the time of the proof (abstract syntax trees compared after the allowed changes; both gates
  counter-tested with deliberate changes that they must detect). The originals themselves are not public; the proof was run by the
  authors on 22.10.2026 for all 43 modules. Later program changes, each made in the authors' copy as well:
  `build_curve_json.py` (27.10.2026, licence label of each curve), `pd_to_theta_input.py` (3.11.2026, kink guard) and
  `pipeline_crosscheck.py` and `sort_and_render.py` (3.11.2026, no-match verdict and SnapPy status), `alex_point.py`, `fox_milnor.py`,
  `homfly_crosscheck.py`, `mandatory_invariants.py` and `slice_verdict.py` (3.11.2026, strictness and CRT bound), all described above; the stage driver `pipeline_run.sh` (23.10.2026; a shell script, outside the Python proof): under `set -e` the stage-2 tracker
  ladder (minimum cell count 2000, 1000, …, 100) stopped at the first failing step instead of trying the next, so from 24.08. to
  23.10.2026 a run whose tracker failed at 2000 ended there; it changes which runs finish, not the values of a run that finishes;
  and again on 3.11.2026, stages 8a/8b report success only when they succeeded.
  Only in this release, not in the authors' copy: the location changes of 3.11.2026 (under *Locations* above); the corpus layer reads English records since 3.11.2026 (`record_names.py`, above),
  `theta_crosscheck.py` writes its note "pathological Cs order" (e5218b0) and, since 3.11.2026, its message for a missing PD code in English; and, from 3.11.2026, `summary_writer.py`
  leaves `chern_simons` empty for a non-hyperbolic knot, `pipeline_run.sh` names a missing `Theta_numeric.sage`, three messages are
  English (`build_corpus_db.py`, `mandatory_invariants_extra.py`, `enrich_unknotting_signed_gauss.py`) and `build_corpus_db.py`
  prints its help for `--help` instead of starting a build. Results the authors computed before such
  a change came from the code before it. The kink guard only adds a refusal: on a diagram it accepts, the conversion is unchanged.
  `pipeline_run.sh`, `mem_watchdog.sh` and `proof_register.json` were translated in the same steps as the 43 modules; they are not
  Python, so the equivalence gates do not cover them.
- The eight 6D modules (`canonical_flow.py`, the five `sixd_*.py`, `calibrate_full_corpus_robust.py`, `test_run_jobs_workertod.py`;
  the clearance modules `clearance_numpy.py` and `clearance_segment.py`, used by `build_curve_json.py`, are among the 43) were added on 29.10.2026 the same way: comments, docstrings and messages translated and paths moved into
  `site_config.py` (`SIXD_DIR`), proven EQUIVALENT to the originals by `strict_equivalence.py` and `path_equivalence.py` (tag
  `proven-equivalent-sixd-20261029`); then identifiers, one command-line option and the remaining texts translated, proven by
  `map_equivalence.py` against that tag, renaming in `tools/translation_table_sixd.json` (for example `eich_vollkorpus_robust.py`
  → `calibrate_full_corpus_robust.py`, `--max-versuche` → `--max-attempts`). The only change to the earlier files is the new
  name `SIXD_DIR` in `site_config.py`. On the example `examples/trefoil.fseries` the released `sixd_scan.py` writes the same
  fingerprint file as the original, byte for byte.
- The corpus layer (15 modules) was added on 29.10.2026 the same way: stage 1 proven EQUIVALENT by `strict_equivalence.py` and
  `path_equivalence.py` (tag `proven-equivalent-corpus-20261029`); stage 2 translated identifiers, seven command-line options and
  seven messages, proven by `map_equivalence.py`, renaming in `tools/translation_table_corpus.json` (for example
  `span_untere_schranke.py` → `span_lower_bound.py`, `--eichen` → `--calibrate-census`). JSON field names were left unchanged on
  purpose. Built from the same records, the released `build_corpus_db.py` writes the same `corpus_db.json` as the original
  (all 843 records; only the build timestamp differs), and a dry run of the released `label_distinct_knots.py` prints exactly
  the same report as the original.
- `worldline_tracker.py` was added on 29.10.2026 too (already English; its one path moved into `site_config.py`): `path_equivalence.py`
  EQUIVALENT to the original; on the crossings of `examples/trefoil.fseries` both write the same result file, byte for byte.
- `tests/test_scanner_vs_numpy.py` checks the GPU scanner against an independent NumPy brute force (22.10.2026: 261 reference
  crossings, 0 differing).

## License
GNU General Public License v3.0 or later (`LICENSE.txt`, SPDX `GPL-3.0-or-later`).

## Citation
See `CITATION.cff`. Each GitHub release is archived on Zenodo with its own DOI.
Author: Werner Alois Stanggassinger, ORCID [0009-0005-8569-8317](https://orcid.org/0009-0005-8569-8317) — all deposits and preprints are listed there.
