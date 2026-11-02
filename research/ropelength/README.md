# Ropelength: checker, tightener, smooth upper bound

Tools for tight knots: the polygonal ropelength of a closed polygon, a GPU tightener that shortens a polygon at thickness 1,
and an upper bound on the **smooth** ropelength from a polygon. Convention: unit tube **radius** (as in ridgerunner, KnotInfo
and Tao's list of constants, problem 22a).

**Why we built them.** On the 432-crossing torus knot T(24,19), ridgerunner 2.2.2 (Cantarella and Piatek; the method of Ashton–Cantarella–Piatek–Rawdon) went in
36,000 seconds (49,469 steps) from polygonal ropelength 2674.23 to 1927.16 in our run of 31 October 2026. We wanted an
independent check of every number, a tightener of our own on the GPU, and the smooth bound that Tao's constant is about.
Our first tightener, a SONO-style shrink-and-push scheme (`attempts/tight_gpu.py`, kept as a record), stalled at 38.06 on a
loose trefoil (`calibration/T2_3.vect`, 400 vertices, resampled to 165), so `tight_pen.py` minimises the length directly. From
that loose trefoil `tight_pen.py --opt gd` ended at 39.15 (both constraints active at its end); from the start polygons below
it reaches the values in the table. The loose start was not tried with `--opt lbfgs`.

## Polygonal and smooth ropelength
The ropelength of a knot type is defined for smooth curves (Cantarella–Kusner–Sullivan). Tightening codes work with polygons
and Rawdon's polygonal thickness, which is not the same number. Ashton, Cantarella, Piatek and Rawdon (arXiv:1002.1723, §5)
turn a polygon into a smooth upper bound: round every corner with the circle arc of radius MinRad(v) tangent to both edges,
then bound the thickness of the rounded curve from below. `roundout_check.py` follows it — rounding and exact length as there; for the distance bound it uses its own sample-pair
version instead of their subsquare elimination:
- the rounded curve is C^{1,1}, its smallest radius of curvature is r_min = min MinRad(v), its length is exact:
  Len(polygon) − Σ (2 r_v tan(θ_v/2) − r_v θ_v);
- thickness = min(r_min, half the smallest doubly-critical self-distance) (Theorem 1.3 there, Litherland–Simon–Durumeric–Rawdon);
- a doubly-critical pair is joined by arcs of total curvature at least π (the Fenchel argument in the proof of Lemma 2.7 there),
  so its arc separation is at least π r_min; the smallest self-distance over sample pairs with at least that separation
  (sample spacing ε), minus (1 + 1/r_min) ε² (Proposition 5.1 there), is a lower bound for it;
- smooth ropelength ≤ rounded length / that thickness. (Proposition 5.1 needs the distances to exceed 1/2; the script asserts
  that the sampled minimum minus 2ε does.)

This is floating-point geometry (float64 throughout, distances on the GPU with a margin of 10⁻⁹), not interval arithmetic.

## Files
| File | What it does |
|---|---|
| `ropelength_check.py` | polygonal ropelength (Rawdon): thickness = min(min MinRad, dcsd/2); pure NumPy |
| `tight_pen.py` | the tightener: minimises length + μ·(penalty on non-local segment pairs closer than 2 + penalty on turning angles above the MinRad ≥ 1 limit), μ raised stage by stage; SciPy L-BFGS-B on the exact gradient (`--opt lbfgs`) or backtracking gradient descent with the tangential part of the gradient removed (`--opt gd`); periodic spline re-equilateralisation in both; `--gradcheck` compares the analytic gradient with finite differences. Writes the best polygon and a JSON log line by line |
| `roundout_check.py` | the smooth upper bound described above (CuPy) |
| `polygon_to_pd.py` | PD code of a space polygon (random projection) and its Alexander polynomial (spherogram) — to check that a run kept the knot type |
| `check_torus_alexander.py` | compares that Alexander polynomial with the torus-knot formula |
| `build_torus_knot.py`, `resample_vect.py` | a start polygon for T(p,q) on a standard torus; resampling to N equally spaced vertices (linear or periodic spline) |
| `attempts/tight_gpu.py` | the SONO-style first attempt (record of the failed approach) |
| `calibration/` | start polygons: `3_1_n165.vect`, `10_124_s360.vect` (made from ridgerunner's example files `3.1.vect` and `10.124.vect`, GPL-3: tightened by ridgerunner on their 47 and 102 vertices to polygonal ropelength 33.1499 and 72.1369, then resampled to 165 vertices (linear) and 360 vertices (periodic spline), which raises the polygonal ropelength to 63.61 and 74.77); `circle_R3.vect` (regular 2000-gon on a circle of radius 3); `T2_3.vect` (a loose trefoil) |

## Run
```
python -I ropelength_check.py calibration/circle_R3.vect                                        # any env with NumPy
python -I tight_pen.py calibration/3_1_n165.vect out.vect --opt lbfgs --minutes 20 --stages 16   # gnottheorie_env (CuPy, SciPy)
python -I roundout_check.py out.vect                                                             # gnottheorie_env
python -I polygon_to_pd.py out.vect                                                              # sage_stable (spherogram)
```

## Checked on 2 November 2026, with the files in this folder
Released files, `tight_pen.py --opt lbfgs --minutes 20 --stages 16` from the start polygons in `calibration/` (the circle is only
checked, not tightened). Every number of ours comes from `ropelength_check.py` / `roundout_check.py`; the ridgerunner values are
ridgerunner's own output; the Alexander polynomial of each result (`polygon_to_pd.py`: t²−t+1 resp. t⁸−t⁷+t⁵−t⁴+t³−t+1) is that of
3_1 resp. 10_124 — an invariant check, not a proof of the knot type:

| | polygonal, ours | polygonal, ridgerunner (same start, our run of 31 October) | smooth upper bound, ours | ACPR's published value (arXiv:1002.1723, table) |
|---|---|---|---|---|
| 3_1, 165 vertices | 32.7786 | 32.8033 | 32.7583 | 32.7436 |
| 10_124, 360 vertices | 71.1371 | 71.1360 | 71.1012 | 71.0739 |
| circle (regular 2000-gon, radius 3) | 6.283191 | — | 6.283317 | 2π = 6.283185 (exact) |

`tight_pen.py --gradcheck`: analytic gradient vs. central finite differences on 40 coordinates, worst relative error 1.24·10⁻⁶
(both penalty terms active). The run is time-capped, so its end point depends on the machine; the authors'
earlier L-BFGS runs (15 resp. 20 minutes) reached 32.7832 / 71.1323 (smooth 32.7574 / 71.1011). Our smooth bounds and ACPR's
values are upper bounds for the minimal smooth ropelength (the polygonal values are not); for the circle, 2π is exact. Our logs of these checks and the two
tightened polygons are in `calibration/logs/`.

## How these files relate to the authors' originals
All got the licence header of this repository. `polygon_to_pd.py`: its docstring says "check" instead of "certify" the knot
type. `check_torus_alexander.py` is the authors' `check_alex.sage.py` with a usage line added. `attempts/tight_gpu.py` imports
`ropelength_check.py` from the folder above. Nothing else changed in the scripts. The logs in `calibration/logs/` are our outputs
with the absolute paths shortened to file names (their "certify with ropelength_check.py" is the scripts' own message).
