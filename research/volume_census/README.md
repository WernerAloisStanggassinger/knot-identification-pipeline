# Volume census

SnapPy's hyperbolic volume for every knot of one of Burton's census files, from the census DT code. Resumable, written on arrival.

**Why we built it.** Burton's census files give name, `knot_sig` and DT code, no volumes. To compare knots across crossing
numbers by volume (a mutation does not change the volume) we needed the volume of every knot.

## Files
| File | What it does |
|---|---|
| `volume_census.py` | one census file → `name<TAB>volume<TAB>flag<TAB>attempts` per knot. Up to 8 attempts per knot (re-triangulating) until all tetrahedra are positively oriented; the flag says how it ended (`pos`, another solution type, or `ERR`); 5 s per attempt. Rows are written in chunks (`#chunk <i> <n>` lines), and a restart skips the chunks already written. Checks 4_1 and 5_2 against their known volumes first |
| `calibrate_against_knotscape.py` | compares the multiset of volumes with a second list (name, volume), e.g. volumes computed from another knot table — names differ between tables, so only the multisets are compared — at 6, 8 and 9 digits; rows of the second list with volume ≤ 0.5 are dropped as non-hyperbolic |

These volumes are floating-point values from SnapPy's shape solution, not verified (interval) values.

## Run
```
python -I volume_census.py /path/to/burton/raw/15n-hyp.csv.bz2 vol_15n.tsv 4        # sage_stable (SnapPy)
python3 -I calibrate_against_knotscape.py vol_15n.tsv knotscape_15n.tsv
```

## Checked
- 2 November 2026, authors' run of 1 November against volumes computed on 30 October 2026, also with SnapPy, from KnotScape's
  table of 15-crossing non-alternating knots (a different table of diagrams; the same program):
  168,023 volumes, all with flag `pos`. At 6 digits every one of them occurs in the KnotScape-based list with the same
  multiplicity; that list has 4 further values above 0.5 (3.663938, 3.663981, 5.333573, 5.333580). These lie
  close to (agreeing to 4 to 6 digits with) SnapPy's numerical values for four of Burton's six 15-crossing satellite knots (15ns_2, 15ns_3, 15ns_5, 15ns_6:
  5.333568, 3.663891, 5.333538, 3.663863). Satellite knots are not hyperbolic and are not in the file of hyperbolic knots. At 8 and 9 digits the multisets differ also by rounding at the
  last digit.
  The KnotScape-based list and the script that made it are not part of this release, so this check cannot be rerun from
  the repository alone.
- 2 November 2026, the released `volume_census.py` on `15n-hyp.csv.bz2` (134 s, 4 workers) against the authors' run of 1 November:
  same 168,023 names; volumes equal up to a relative difference of at most 9.4·10⁻¹⁵ (last-digit rounding; 9,970 rows differ
  there). The re-triangulation between attempts is random and not seeded, so the number of attempts can differ between runs:
  one knot, 15nh_000494, needed 5 attempts in the authors' run and ended after 8 with `contains_negatively_oriented_tetrahedra`
  in this one, with the same volume to all printed digits (12.2211259055064).

## How these files relate to the authors' originals
Both got the licence header of this repository; their docstrings were reworded (no occasion, no internal file names, and the
calibration says "is compared with" instead of "must equal" and describes the second list generically instead of naming
the authors' file). Nothing else changed.
