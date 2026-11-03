# Hyperbolic torsion polynomial against genus and fiberedness (Dunfield–Friedl–Jackson Conjecture 1.7)

Dunfield, Friedl and Jackson (arXiv:1108.3045v3, p. 5) conjecture that for a hyperbolic knot K the hyperbolic torsion polynomial T_K
determines the genus, and that K is fibered if and only if T_K is monic. In the form of Agol–Dunfield (arXiv:1501.02136v3,
Conjecture 1.1): deg T_K = 4g − 2, where deg is the span of the symmetric Laurent polynomial. DFJ checked it for the hyperbolic
knots up to 15 crossings. These scripts compare T_K with the genus and fiberedness of every hyperbolic knot of a Burton census file.

**Why we built it.** SnapPy computes T_K (`Manifold.hyperbolic_torsion`) and knot Floer homology gives genus and fiberedness
exactly, but a census-wide comparison needs: a numerically safe reading of degree and leading coefficient, escalation in precision
where a reading is uncertain, retries when a triangulation is not geometric, and a resumable run over files of up to hundreds of
millions of knots. That is what
this folder adds; T_K itself and the genus and fiberedness are computed by SnapPy and by knot Floer homology.

## The comparison (`dfj_core.py`)
- **T_K**: SnapPy `hyperbolic_torsion(bits_prec=500)` (DFJ's conventions), on a triangulation with all tetrahedra positively oriented;
  up to 6 tries with `randomize()`, otherwise the knot is an ERROR, never "agree". SnapPy sets coefficients with |c| < 2^(−bits/3) to zero.
- **Genus and fibered**: knot Floer homology (SnapPy's `knot_floer_homology`, Szabó's program) on the same diagram, independent of T_K.
- **Degree**: deg T_K = 4g − 2.
- **Monic**: the leading coefficient lc equals +1 as a complex number. |lc − 1| < 2^−100 is "monic", > 2^−30 is "not monic", anything
  in between is "uncertain". The same thresholds decide whether T_K is real (largest imaginary part divided by the larger of 1 and the
  largest coefficient).
- **Status**: `agree` (degree right and monic ⟺ fibered), `MISMATCH`, `uncertain`, or `ERROR`. A result that is not `agree`, or whose
  real-ness is uncertain, is recomputed on a re-randomized triangulation at 1000 and then 2000 bits. An `ERROR` (no geometric
  triangulation after 6 tries, or SnapPy failing) is not escalated; `check_census.py --retry-errors` computes it again at 500 bits.
- **Why the degree reading is safe**: DFJ's Theorem 1.5 gives deg T_K ≤ 4g − 2. A spurious non-zero top coefficient would show up as
  deg > 4g − 2, a true top coefficient below the cleaning threshold as deg < 4g − 2; both are flagged and escalated.
- Where T_K is not certainly non-real, SnapPy's `symmetry_group().is_amphicheiral()` is recorded. That is numerical, not certified.
  (DFJ also observe that up to 15 crossings T_K is real exactly for the amphichiral knots; amphichiral ⟹ real is their Theorem 1.2(e).)

## Files
| File | What it does |
|---|---|
| `dfj_core.py` | the comparison above: T_K with retries, reading of degree, leading coefficient and real-ness, verdict, escalation in precision |
| `check_census.py` | every knot of one Burton census file: knot Floer homology, then `dfj_core.check`; one JSON line per knot written on arrival, resumable, shardable over processes; non-hyperbolic census names are written as SKIP |
| `calibrate.py` | calibration inside DFJ's range and counter-cases (below) |
| `two_precision.py` | stability check: recomputes T_K at 1000 bits on a re-randomized triangulation for N random knots and the 10 numerically tightest non-monic ones (smallest \|lc\| or \|lc − 1\|) of a `check_census.py` output and compares degree, classes and lc |

## Run
```
export BURTON_CENSUS_DIR=/path/to/burton/raw                                  # Burton's census files (see THIRD_PARTY.md)
python -I calibrate.py calib.jsonl 400                                        # SageMath's Python
python -I check_census.py $BURTON_CENSUS_DIR/16n-hyp.csv.bz2 out_0.jsonl 0 10   # worker 0 of 10; one out file per worker
python -I two_precision.py $BURTON_CENSUS_DIR/16n-hyp.csv.bz2 out_0.jsonl 100 two_precision.jsonl
```
`--retry-errors` computes the ERROR rows of an out file again. `check_census.py` and `two_precision.py` limit their address space to
`MEM_GB` (default 6) GB; `calibrate.py` to 6 GB.

## Calibration and counter-cases (`calibrate.py`)
- 19 named knots up to 14 crossings (among them 4_1, 5_2, 8_17, 8_20, 9_42, the Kinoshita–Terasaka knot K11n42 and the Conway knot
  K11n34, K12n242, K14a506) should agree. The script writes one line per case and asserts nothing; the expectations below are read
  from its output.
- A random sample from DFJ's range (Burton census, seed 20261016): 400 from 15n, 100 from 15a, 100 from 14n, 50 from 13n, 50 from 12n.
- Counter-cases that should **not** agree: a wrong genus (4_1 with g + 1, K11n42 with g − 1) and a flipped fibered flag (5_2, 8_20) must
  give MISMATCH; the torus knots T(2,5), T(2,7), T(3,5) and the non-hyperbolic knots of the 16n census (one torus, ten satellites)
  must give ERROR (no geometric triangulation).

## Checked on 3 November 2026
SageMath 10.9, SnapPy 3.3.2.
- `calibrate.py`: 19 of 19 named knots and 700 of 700 sample knots agree; the 4 wrong-input counter-cases give MISMATCH and the 14
  non-hyperbolic knots give ERROR.
- `check_census.py` on `10n-hyp.csv.bz2`: 41 of 41 knots agree.
- `check_census.py` on 21 knots of `16n-hyp.csv.bz2` (worker 0 of 50,000): degree, leading coefficient, monic and real class, genus,
  fibered flag and status are the same as in the authors' earlier run of the same comparison on these knots.

## How these files relate to the authors' originals
`dfj_core.py` is the authors' file of 16 October 2026 with the licence header, and its first comment line (which named an internal
occasion) replaced; the code is unchanged. `calibrate.py` is the authors' file with the licence header, its first comment line
reworded, a usage line, and the census folder and the module folder found from `BURTON_CENSUS_DIR` and the script's own folder
instead of absolute paths. `check_census.py` and `two_precision.py` replace the authors' run scripts, which read diagrams, genus and
fiberedness from an earlier run of theirs; `check_census.py` computes genus and fiberedness itself, from knot Floer homology on the
census diagram; `two_precision.py` picks its 10 tightest cases among the non-monic results only (the authors' script, by a slip in
its sort key, picked monic ones).
