# data/ — how to build the data the pipeline reads

Nothing here is shipped as data: the tables belong to their authors, and the camera flight is regenerated. Each recipe writes to the
location that `pipeline/site_config.py` expects in its portable profile. KnotInfo and Fremlin are read **only** by the cross-check
stage 8b, after the knot has been identified (blind test).

## Camera flight (`data/camera_131072.npz`, needed by the scanner)
```
python3 pipeline/make_camera_trajectory.py --n 131072 --out data/camera_131072
```
Reproduction check, 22.10.2026: the regenerated file equals the one used in the authors' runs, array by array.

## KnotInfo index (`data/knotinfo_index/records.json`, `index_c_det_sig.json`)
Source: KnotInfo (C. Livingston, A. H. Moore), through the Python package `database_knotinfo`.
```
conda run -n sage_stable python data/build_knotinfo_index.py
```
Reproduction check, 22.10.2026, with `database_knotinfo` 2026.10.5: both files are byte-identical to the index used in the
authors' runs (built 14.05.2026): 12 966 knot records, 1 792 (c, det, signature) classes. The script differs from the authors'
original only in its default output folder (proof: abstract syntax trees compared, the default argument is the only difference).

## Fremlin cache (`data/fremlin_cache/fremlin_cache.json`)
Source: David Fremlin, "Knots and their symmetries", https://david.fremlin.de/knots/ (pages not redistributed here).
```
python3 data/build_fremlin_cache.py --fetch      # downloads the index page and the 38 pages it links, 1 s apart
```
The server refuses Python's default user agent, so the script sends its own (`knot-identification-pipeline/1.0 …`).

Reproduction check, 22.10.2026:
- On the pages saved by the authors on 11.05.2026 the output is byte-identical to the authors' cache (37 knot records).
- On the live pages of 22.10.2026 exactly one record differs: knot 15331, where Fremlin added a note on 24.5.26 (its date field
  changes from 21.4.15 to 24.5.26). That note follows the Jones polynomial directly and made the authors' original pattern lose
  the polynomial; the pattern was widened for this release (marked in the code), which leaves the output on the old pages unchanged.
