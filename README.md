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
The main result is `<knot>_pflicht_invarianten_complete.json` ("mandatory invariants") plus `<knot>_summary.md`.

## How this release relates to the code that produced the results
- Every file in `pipeline/` was translated from German: comments, docstrings and message texts only. Identifiers, JSON field
  names and data strings stayed as they were, so **many field names and some output texts are German** (for example
  `pflicht_invarianten` = mandatory invariants, `beweis` = proof, `fehler` = error). Translating them would change the output
  format, so it is left for a later version.
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
