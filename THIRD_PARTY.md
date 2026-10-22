# Third-party programs and data

This repository contains only our own code. The pipeline calls the programs and reads the data below; they are **not included**.
Get each one from its authors, then point the pipeline at it (`pipeline/site_config.py`: put it at the default location, or set
the environment variable `KNOT_PIPELINE_<NAME>`). `python3 pipeline/site_config.py` lists every location the pipeline will use.

Libraries installed with pip/conda (SnapPy, spherogram, SageMath, Regina, pyknotid, knotpy, NumPy, SciPy, SymPy, pandas, CuPy, pyarrow …)
are listed per conda environment in `environments/`, not here.

| What | Authors | Where to get it | Used by | site_config name |
|---|---|---|---|---|
| **Theta** (`Theta.sage`): the two-variable invariant Θ = (Δ, θ) | Dror Bar-Natan, Roland van der Veen | Sage implementation: https://www.rolandvdv.nl/Theta/ · paper and related files: http://drorbn.net/Theta · paper arXiv:2509.18456 | `theta_crosscheck.py` (stage 8a) | `THETA_DIR` (folder containing `Theta.sage`) |
| **Theta, numerical** (`Theta_numeric.sage`) | our modification of `Theta.sage`, following the paper's Implementation section on numerical specialisation (T1 = 22/7, T2 = 21/13) | not included yet: it is derived from the authors' code, whose licence we have not seen stated | `theta_numeric_crosscheck.py` (stage 8a for large diagrams) | `THETA_DIR` |
| **KnotJob** (`KnotJob.jar`): Khovanov homology, s-invariants | Dirk Schütz | https://www.maths.dur.ac.uk/users/dirk.schuetz/knotjob.html (GNU GPL v3) | `pflicht_invarianten.py`, `pflicht_invarianten_extra.py` | `KNOTJOB_JAR`, `JAVA_BIN` |
| **KnotInfo** index (`index_c_det_sig.json`, `records.json`) | KnotInfo (C. Livingston, A. H. Moore) via the `database_knotinfo` package | built locally from `database_knotinfo` (recipe: `data/README.md`, to come) | `pipeline_crosscheck.py` — **after** identification only, never during it (blind test) | `KNOTINFO_INDEX`, `KNOTINFO_RECORDS` |
| **Fremlin** cache (`fremlin_cache.json`) | Fremlin's knot pages | built locally by parsing his pages (recipe: `data/README.md`, to come) | `pipeline_crosscheck.py` — after identification only | `FREMLIN_CACHE` |
| **Blender** | Blender Foundation | https://www.blender.org/ | `sort_and_render.py` (renders, optional) | `BLENDER`, `RENDER_TEMPLATE` |

Our own files that are needed but generated rather than shipped:
- `camera_131072.npz` (the camera flight, 8.9 MB): `python3 pipeline/make_camera_trajectory.py --n 131072 --out data/camera_131072`
  (site_config name `CAMERA`). Checked 22.10.2026: the regenerated file equals the one used in the authors' runs array by array
  (t_step, theta, phi, cam, rotmat — `np.array_equal` on all five).
