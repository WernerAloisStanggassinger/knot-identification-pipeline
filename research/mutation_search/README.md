# Conway mutation search between two census knots

Given two knots of Burton's census by name, search for a Conway mutation that turns one into the other. Mirror-blind: SnapPy's
`is_isometric_to` answers True only with an explicit isometry (an isomorphism of canonical triangulations), so True proves that
the exteriors are homeomorphic and hence, by Gordon–Luecke, that the knots agree up to mirror image; False is not a proof. A hit
proves: the target is a Conway mutant of the source **or of its mirror image**. An orientation-sensitive invariant can tell the two
apart only where it differs between the target and its mirror (for an amphichiral target there is nothing to tell).

**Why we built it.** We needed, for given pairs of census knots, a direct search for a Conway mutation between them. The hyperbolic
volume and the Alexander, Jones and HOMFLY polynomials do not change under mutation, so equal values of these leave the
question open.

## Files
| File | What it does |
|---|---|
| `mutation_search_census.py` | exhaustive search on the census diagrams of both knots, in both directions: every set of crossings bounded by four edges that forms a tangle (a disc meeting the diagram in four points), every reconnection that is a mutation, then an isometry test of the result against the other knot |
| `mutation_search_diagrams.py` | the same mutation construction on many further, non-minimal diagrams of the source knot (from `spherogram`'s `backtrack()` and from SnapPy's `exterior_to_link()`), within a time budget; one JSON line per diagram, written on arrival |

## Run
```
export BURTON_CENSUS_DIR=/path/to/burton/raw                 # Burton's <c>n-hyp.csv.bz2 files (see THIRD_PARTY.md)
python -I mutation_search_census.py 15nh_048139 15nh_048140                                   # sage_stable
python -I mutation_search_diagrams.py 15nh_048139 15nh_048140 out.jsonl 2 19 600 16           # SRC DST OUT SEED MAXCROSS SECONDS [MINCROSS]
```
Both scripts take names of non-alternating hyperbolic census knots (`<c>nh_…`). Time limits: `mutation_search_census.py` stops
after 3000 s, `mutation_search_diagrams.py` after SECONDS + 900 s (SIGALRM); the census search is exhaustive only if it finishes.

## Checked on 2 November 2026, with the files in this folder
- `mutation_search_census.py 15nh_048139 15nh_048140` finds a mutation in both directions: 15nh_048140 is a mutant of
  15nh_048139 or of its mirror image. The authors' `mutation_check_15.py`, from which the search functions are copied, found the
  mutation 15nh_048139 → 15nh_048140 on 1 November 2026.
- `mutation_search_diagrams.py 15nh_048139 15nh_048140 out.jsonl 2 19 600 16` (only diagrams with at least 16 crossings, so not
  the census diagram) finds the mutation; the found record (tangle, boundary edges, diagram) is identical to the authors' run of
  the same command on 2 November 2026.

A search that finds nothing proves nothing: it covers the diagrams it tried, not all diagrams.

## How these files relate to the authors' originals
Both got the licence header of this repository; the census folder is read from `BURTON_CENSUS_DIR` instead of an absolute path;
the header comments were reworded (no internal occasions, results or folders; the provenance and calibration of the copied
functions are stated as in this README), and the sentence of `mutation_search_diagrams.py` on what a hit proves now says "a mutant
of the source or of its mirror image". Nothing else changed.
