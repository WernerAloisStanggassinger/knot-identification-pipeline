# knot-identification-pipeline

Identifies the knot type of a closed space curve — a sculpture's centreline or a machine-generated curve — directly from its
geometry, given as a Fourier series (`.fseries`). It was built for the stainless-steel knot sculptures of
**Werner Alois Stanggassinger** (sculptor, Wasserburg am Inn, Bavaria) and grew into the tool behind the knot corpus on
https://knot-structures.stainlesssteel4u.de/ (rotatable 3D views, photos and the invariants of every piece).

Stated as a fact: the code was written by Claude Code (Anthropic), which orchestrates the computations; the sculptures, the
geometries and the questions are Werner Alois Stanggassinger's, and so is the machine it runs on. The mathematics comes from the
programs it calls (SnapPy, spherogram, SageMath, Regina, KnotJob, Theta), not from a language model.

## What it does
`pipeline/pipeline_run.sh`, run inside a folder that holds one `.fseries`, goes through these stages:
- **Preflight** — curve sanity, π-symmetry, target crossing count.
1. **Scanner** (GPU) — a camera flight of 131 072 view directions around the curve; every crossing of every projection.
2. **Tracker** — follows each crossing across the views (density clustering in curve parameters).
3. **Resolve** — chooses views and builds a signed Gauss code / PD code.
4. **Aggressive simplify** (4a: merge) — Regina and spherogram simplification over many frames; keeps the smallest diagram.
5. **Mandatory invariants** — Alexander, Jones, HOMFLY, Khovanov, signature, determinant, hyperbolic volume and more.
6. **Extra invariants** — Kauffman polynomial and further invariants.
7. **Enrichments** — bridge and braid bounds, unknotting bounds, length spectrum, non-hyperbolic cases (not fatal if one fails).
8. **Summary**; 8a **Theta** (Bar-Natan–van der Veen Θ); 8b **Cross-check** — only *after* identification, comparison with
   KnotInfo, Fremlin and the SnapPy census.
9. **Sort and render** — folder placement and a Blender render (`SKIP_SORT_RENDER=1` skips it).

**Blind test.** During identification no knot table is consulted: `pipeline_seal.py` blocks the table lookups of SnapPy,
spherogram and Sage. Tables are read only in stage 8, to compare.

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
- Every file in `pipeline/` was translated from German in two steps. 22.10.2026: comments, docstrings and messages. 23.10.2026:
  also identifiers, module and file names, command-line options, environment variables, JSON field names, category labels and
  the remaining texts. `tools/map_equivalence.py` proves for every file that the translated version is the program proven
  equivalent to the originals (git tag `proven-equivalent-20261027`, see the next points) up to one consistent renaming,
  and writes that renaming to `tools/translation_table.json` (identifiers, keys, file renames, free texts).
  The originals were proven equivalent on 22.10.2026 (tag `proven-equivalent`). On 27.10.2026 every one of the 43 originals was
  compared again with that state: 42 differ only in texts (their comments and messages were translated in place on 25.10.), one
  changed as a program: `build_curve_json.py` labels each curve with its licence (sculpture curves all rights reserved, machine-knot
  curves CC-BY-NC-4.0). Tag `proven-equivalent-20261027` = `proven-equivalent` plus exactly that change. The gate was
  counter-tested with 25 deliberate changes it must catch or pass; the built-in self-tests of six modules give the same
  results before and after (45 of 45 checks).
- **Results the authors published before this release (knot corpus, website downloads, Zenodo deposits) use the German
  names.** `tools/translation_table.json` maps every name this code reads or writes, for example `pflichtliste` → `mandatory_list`,
  `pflicht_invarianten_complete.json` → `mandatory_invariants_complete.json`, `beweis` → `proof`, `grund` → `cause`. In
  `pipeline/proof_register.json` additionally `belegt_von` → `receipted_by` and `_doku` → `_doc`.
- A few names stay German on purpose because they name things outside the code: the authors' folder names
  (`erzeugte_knoten/`, `bekannt/`), the material `Edelstahl` in the Blender template, and historical file names.
- All absolute paths moved into `pipeline/site_config.py`. With `KNOT_PIPELINE_PROFILE=authors` every location is the authors' own.
- `tools/strict_equivalence.py` and `tools/path_equivalence.py` prove, file by file, that each released module is the same program
  as the one used in the authors' runs (abstract syntax trees compared after the allowed changes; both gates counter-tested with
  deliberate changes that they must detect). The originals themselves are not public; the proof was run by the authors on
  22.10.2026 for all 43 modules.
- `tests/test_scanner_vs_numpy.py` checks the GPU scanner against an independent NumPy brute force (22.10.2026: 261 reference
  crossings, 0 differing).

## License
GNU General Public License v3.0 or later (`LICENSE.txt`, SPDX `GPL-3.0-or-later`).

## Citation
See `CITATION.cff`. Each GitHub release is archived on Zenodo with its own DOI.
Author: Werner Alois Stanggassinger, ORCID [0009-0005-8569-8317](https://orcid.org/0009-0005-8569-8317) — all deposits and preprints are listed there.
