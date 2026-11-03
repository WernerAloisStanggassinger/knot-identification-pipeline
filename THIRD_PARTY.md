# Third-party programs and data

This repository contains only our own code, plus two start polygons derived from ridgerunner's example files (see its row). The pipeline calls the programs and reads the data below; they are **not included**.
Get each one from its authors, then point the pipeline at it (`pipeline/site_config.py`: put it at the default location, or set
the environment variable `KNOT_PIPELINE_<NAME>`). `python3 pipeline/site_config.py` lists every location the pipeline will use.

Libraries installed with pip/conda (SnapPy, spherogram, SageMath, Regina, pyknotid, knotpy, NumPy, SciPy, SymPy, pandas, CuPy, pyarrow …)
are listed per conda environment in `environments/`, not here. knotpy (Boštjan Gabrovšek) is not on PyPI: it is installed from GitHub, https://github.com/bgabrovsek/knotpy, commit 8ce3ad2 (pinned in `environments/base.txt`).

| What | Authors | Where to get it | Used by | site_config name |
|---|---|---|---|---|
| **Theta** (`Theta.sage`): the two-variable invariant Θ = (Δ, θ) | Dror Bar-Natan, Roland van der Veen | Sage implementation: https://www.rolandvdv.nl/Theta/ · paper and related files: http://drorbn.net/Theta · paper arXiv:2509.18456 | `theta_crosscheck.py` (stage 8a only with `THETA_SYMBOLIC=1`; the full polynomial, small diagrams) | `THETA_DIR` (folder containing `Theta.sage`) |
| **Theta, numerical** (`Theta_numeric.sage`) | our modification of `Theta.sage`, following the paper's Implementation section on numerical specialisation (T1 = 22/7, T2 = 21/13) | not included yet: it is derived from the authors' code, whose licence we have not seen stated | `theta_numeric_crosscheck.py` (stage 8a by default; without this file stage 8a reports an error and is not run — set `THETA_SYMBOLIC=1` to use `Theta.sage` instead) | `THETA_DIR` |
| **`Rot[]` of `Theta.m`** (rotation numbers and crossing list of a PD code) | Dror Bar-Natan, Roland van der Veen (`Theta.m`, http://drorbn.net/Theta) | — | `pipeline/pd_to_theta_input.py` and `research/theta_census_gpu/pd_to_theta_input.py` are a line-by-line translation of `Rot[]` (`Theta.m` lines 22–40) into Python. We have not seen a licence stated for `Theta.m` (the arXiv listing of 2509.18456 carries arXiv's non-exclusive distribution licence only). The GPL notice in the headers of these two files covers our translation and additions, not the original; we have not yet asked the authors. | — |
| **khoca**: sl(N) and Khovanov homology | Lukas Lewark, Sebastian Oehms (package metadata of khoca 1.5) | pip `khoca` 1.5 (GPLv2+), https://people.math.ethz.ch/~llewark/khoca.php, in an environment of its own (the authors use `khoca_env`; it is not among `environments/*.txt`) | `research/chirality_witnesses/khovanov_witness.py` | — |
| **KnotJob** (`KnotJob.jar`): Khovanov homology, s-invariants | Dirk Schütz | https://www.maths.dur.ac.uk/users/dirk.schuetz/knotjob.html (GNU GPL v3) | `mandatory_invariants.py`, `mandatory_invariants_extra.py` | `KNOTJOB_JAR`, `JAVA_BIN` |
| **KnotInfo** index (`index_c_det_sig.json`, `records.json`) | KnotInfo (C. Livingston, A. H. Moore) via the `database_knotinfo` package | built locally from `database_knotinfo`: `data/build_knotinfo_index.py` (see `data/README.md`) | `pipeline_crosscheck.py` — **after** identification only, never during it (blind test) | `KNOTINFO_INDEX`, `KNOTINFO_RECORDS` |
| **Fremlin** cache (`fremlin_cache.json`) | David Fremlin, "Knots and their symmetries" | https://david.fremlin.de/knots/ — built locally: `data/build_fremlin_cache.py --fetch` (see `data/README.md`) | `pipeline_crosscheck.py` — after identification only | `FREMLIN_CACHE` |
| **Burton's knot census** (`<c>{a,n}-{hyp,torus,satellite}.csv.bz2`, DT codes of all prime knots up to 19 crossings) | Benjamin A. Burton, *The next 350 million knots* (SoCG 2020) | https://people.smp.uq.edu.au/BenjaminBurton/knots/ (the files we downloaded) | `research/theta_census_gpu/`, `research/volume_census/`, `research/mutation_search/`, `research/theta_collision/` (through the first two), `research/whitehead_double/` (census names), `research/su2_abelian_surgery/` (census names), `research/torsion_polynomial/` | — (paths are command-line arguments; `BURTON_CENSUS_DIR` for the mutation search, the Whitehead doubles, `certify_knots.py` and `calibrate.py`) |
| **ridgerunner** 2.2.2 (example polygons `3.1.vect`, `10.124.vect`; the tightener we compare with) | Jason Cantarella, Michael Piatek (its AUTHORS file); method: Ashton, Cantarella, Piatek, Rawdon, arXiv:1002.1723 | its distribution (GNU GPL v3; `share/ridgerunner/`) | `research/ropelength/calibration/` (two start polygons made from its example files) | — |
| **Blender** | Blender Foundation | https://www.blender.org/ | `sort_and_render.py` (renders, optional) | `BLENDER`, `RENDER_TEMPLATE` |

Our own files that are needed but generated rather than shipped:
- `camera_131072.npz` (the camera flight, 8.9 MB): `python3 pipeline/make_camera_trajectory.py --n 131072 --out data/camera_131072`
  (site_config name `CAMERA`). Checked 22.10.2026: the regenerated file equals the one used in the authors' runs array by array
  (t_step, theta, phi, cam, rotmat — `np.array_equal` on all five).
