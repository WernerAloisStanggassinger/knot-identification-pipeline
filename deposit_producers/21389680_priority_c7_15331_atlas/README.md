# Producers for Zenodo 10.5281/zenodo.21389680 — "Priority deposit: c7 inequality at 14/15 crossings, Knot 15331 bounds, and a blind Farey/precessing-ellipse knot atlas"

Section 3 (Hizume atlas, deposit `PRIORITAET_DEPOSIT_20260715/HOCHLADEN_v2/README.md`:96-104). The three scripts are released as `pipeline/precellipse_machine.py`, `pipeline/precellipse_atlas.py` and `pipeline/build_precellipse_atlas_md.py`.

- `pipeline/precellipse_machine.py` — precessing-ellipse curve generator, FFT to `.fseries`, blind measurement via the pipeline's `polygon_resolve` (imported from `freeze_manifest/pipeline_scripts`, plus `maschinen_ideen/01_knotting_oscillations/generate_oscillation_knots.py`).
- `pipeline/precellipse_atlas.py` — sweeps the Farey family, identifies each knot blind, cross-checks det / c / Alexander degree against the closed torus formulas, computes Jones; writes `erzeugte_knoten/precellipse_atlas.json`.
- `pipeline/build_precellipse_atlas_md.py` — turns that JSON into `precellipse_atlas.md`. The deposited `hizume_precellipse_atlas.md` is byte-identical to `erzeugte_knoten/precellipse_atlas.md` (checked with `cmp` on 2026-10-31).

Open: the producer of `hizume_precellipse_atlas_braids.zip` was not found. Section 2 (Knot 15331, `knot15331_sharper_bounds.txt`:11-14, "multi-frame Reidemeister sweep over 5000 generic projections") names no script. The multi-frame sweep of `pipeline/aggressive_simplify.py` (already public) is a candidate, but no receipt ties it to the 5000/500 numbers, so it is listed as open, not as the producer. Section 1 (c7 inequality) is not in this inventory row.

## Where the files are in this repository
The producers of this deposit are released elsewhere in this repository:
- `pipeline/precellipse_machine.py`, `pipeline/precellipse_atlas.py`, `pipeline/build_precellipse_atlas_md.py`
- `machines/01_knotting_oscillations/generate_oscillation_knots.py`
