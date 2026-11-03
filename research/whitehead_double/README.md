# Whitehead-double HOMFLY

For a knot K, the HOMFLY polynomial of its untwisted Whitehead doubles (positive and negative clasp), and a comparison of two knots
by it that ignores mirror images. Two knots whose volume, Alexander, Jones and HOMFLY polynomials agree can sometimes be told apart
this way. Stoimenow's tables of mutant knots (https://stoimenov.net/stoimeno/homepage/ptab/index.html) use this invariant, with group-theoretic tests, to separate pairs with
equal invariants that are not listed as mutants.

**Why we built it.** No tool we have does this in one step. Regina builds the double (`Link.whiteheadDouble`) and SageMath computes
the HOMFLY polynomial, but the comparison of whole groups of census knots, both clasps, against both mirror images, has to be
scripted. We needed it for pairs of census knots with equal volume, Alexander and Jones polynomial at different crossing numbers.

## What a result means
- **Different** values for a clasp, against B and against mirror(B), prove that A and B are different knots up to mirror image.
- **Equal** values prove nothing. The Conway mutant pair among the counter-cases below gives equal values. Whether equal values hold
  for every mutant pair is not claimed here, and no theorem saying so is cited; the pair is a regression check (see below).

## Files
| File | What it does |
|---|---|
| `wd_homfly.py` | builds the doubles with Regina (untwisted, both clasps, for the knot and its mirror image), computes their HOMFLY polynomials with SageMath (libhomfly), in pairs mode writes one JSON line per double the moment it is computed, then prints one verdict line per pair |

## Run
```
export BURTON_CENSUS_DIR=/path/to/burton/raw        # only for census names; Burton's <c><a|n>-<hyp|satellite|torus>.csv.bz2 (see THIRD_PARTY.md)
python -I wd_homfly.py counter 4                     # counter-cases, 4 workers (SageMath's Python)
python -I wd_homfly.py pairs.tsv out.jsonl 8         # pairs.tsv: two names per line
```
Names are census names (`18nh_00010932`) or names from SnapPy's link tables (`K14n22185`). A restart with the same `out.jsonl` skips the
doubles already on disk. The time per double depends on the knot and on the load of the machine more than on the crossing number:
in our runs (1 to 3 November 2026, each double in its own worker process) from about 30 s for some 14-crossing knots to about 23 minutes for some
15- to 17-crossing knots.

## Counter-cases (`counter`)
| Pair | Expected | Why |
|---|---|---|
| K14n22185 / K14n22589 | different, both clasps | Stoimenow's pair (14_41721, 14_42125) in his numbering of all 14-crossing knots (= KnotScape 14n 22185 / 22589, after the 19,536 alternating ones): equal volume, Alexander and Jones polynomial; his page lists it as separated by the Whitehead-double skein polynomial |
| 15ah_27440 / 15ah_27441 | equal, both clasps | Regression check, not a consequence of a cited theorem: these Conway mutants (explicit mutation; = KnotScape 15a_13571 / 15a_13629) gave equal values for both clasps on 3 November 2026, so a difference means the computation has changed; their census DT codes are written into the script |

The script stops with exit code 1 if either expectation fails. The first rests on Stoimenow's published separation, the second only on our own earlier result.

## Checked on 3 November 2026
`python -I wd_homfly.py counter 4` (SageMath 10.9, Regina 7.4, SnapPy 3.3.2): K14n22185 / K14n22589 different for both clasps,
15ah_27440 / 15ah_27441 equal for both clasps, `COUNTER-CASES OK`. A pairs file with census names (`10nh_01 10nh_02`, `10nh_01 K14n22185`)
gives "different knots" for both pairs; run again on the same `out.jsonl`, it computes nothing (`12 doubles, 12 on disk, 0 to compute`).

## How this file relates to the authors' originals
Made from the authors' `wd_check_groups_incremental.py` (1 November 2026), which had the groups and the calibration pairs written into
the code and read the 14-crossing pair from a decoded KnotScape table. Here the pairs come from a file, the 14-crossing pair is taken
from SnapPy's tables (SnapPy identifies the authors' decoded KnotScape knots 14n 22185 and 22589 as K14n22185 and K14n22589), the
15-crossing pair is written in as DT codes, and the census folder is read from `BURTON_CENSUS_DIR`. The construction of the double and
of its HOMFLY polynomial is unchanged.
