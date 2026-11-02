# Jones polynomial at a point, as a tensor network on the PD code (CPU reference + GPU batch)

The Jones polynomial V(t) at a point, exactly modulo a prime, computed directly from a planar diagram (PD code) as a
tensor-network contraction of the Kauffman bracket. There is a CPU reference, and a GPU version that contracts thousands of knots
at once. It needs no braid and no height function.

**What it is for, honestly.** This is an independent second implementation of the Jones polynomial, a cross-check, and **not a
faster one**. Regina's `Link.jones()` gives the exact polynomial in about 0.20 ms per 15-crossing knot and 0.36 ms per 19-crossing knot on
one CPU core, including `Link.fromDT` (0.02 + 0.18 ms and 0.03 + 0.34 ms in `logs/timing_20261102.log`, first 1,000 rows of each
census file). This route needs 2.4 resp. 3.2 ms per knot of CPU preparation (spherogram's DT→PD conversion 1.7 resp. 2.2 ms,
`prepare_fast` for two points 0.7 resp. 0.9 ms) and gives one value of V mod p per point.

Two things led us here:
- An earlier attempt of ours, through braid words and 2^(strands) transfer matrices, slowed to 34 knots/s at 10 strands and 2 knots/s
  at 12, and from 14 strands on its 2^n × 2^n state no longer fit into memory (measured on a fully loaded machine; that
  experiment's logs are not part of this folder).
- We wanted a version that works directly on the PD code.

V at a point t and at 1/t gives a chirality test: the mirror image of K has Jones polynomial V_K(1/t). So V_K(t) ≠ V_K(1/t) shows
that K is chiral, provided the computation is correct.

## The construction (`jones_tn.py`)
- **Edges.** Every edge e of the diagram carries an index s_e ∈ {+1, −1}: the direction in which a state loop runs along e,
  relative to e's orientation.
- **Crossings.** Every crossing tensor is A·(A-smoothing) + A⁻¹·(B-smoothing). A smoothing joins two pairs of neighbouring ports,
  i.e. two corners. A corner is allowed only if the loop enters through one port and leaves through the other. It gets the factor x
  for a left turn and x⁻¹ for a right turn.
- **Edge turnings.** Every edge gets the factor x^(s_e·τ_e). The integers τ_e (turning of the edge, in units of π/2) solve the face
  equations of the planar 4-valent graph: corner turns plus edge turns equal +4 around every bounded face and −4 around the outer
  face. They come from a spanning tree of the dual graph, which gives an integer solution.
- **Loop value.** The claim we use: every simple closed loop of a Kauffman state then has total turning ±4 units. The reason is
  that summing the face equations over the faces on one side of the loop cancels the interior edge turns. Its two orientations
  then give x⁴ + x⁻⁴, which is d = −A² − A⁻² when x⁴ = −A².
- **Result.** The full contraction is ⟨K⟩ with ⟨O⟩ = d, and V(t) = (−A³)^(−w) ⟨K⟩/d with t = A⁻⁴ (w = writhe).
- **Status.** We did not check the literature for this formulation and claim no novelty for it. The argument above
  is not written out as a full proof; the calibration below checks the results.
- **Arithmetic.** Everything is in GF(p) with p = 536870909, which is ≡ 1 mod 4, so i = √−1 exists. Take A = i x², so t = x⁻⁸.
  Every product is reduced mod p at once, and a CPU tensordot sums at most 16 products below 2⁵⁸. `jones_tn.jones_at` and
  `jones_gpu.jones_from_bracket` stop if d ≡ 0 at the chosen point.

Three conventions are fixed by calibration against SageMath, not assumed: which of the two smoothings is "A", the sign of the
writhe, and the orientation of t.

## Files
| File | What it does |
|---|---|
| `jones_tn.py` | The construction above. `jones_at(pd, x)` is the CPU reference (NumPy, contraction in PD order). `crossing_exponents` and `eval_tensors` give the same tensors in an x-independent monomial form (fast path). |
| `jones_gpu.py` | Batched GPU contraction (CuPy, int64 mod p). Per knot it picks a crossing order that keeps the frontier of open edges small (greedy) and assigns a bit slot to each open edge. Per batch it updates a state of size 2^F per knot, crossing by crossing, over the 16 leg assignments. `prepare_fast` builds one plan for several points. Refuses diagrams with a kink, an edge with both ends at one crossing; the census counts these as errors (0 at 15 crossings). Kinks are not covered by the shipped tests. |
| `jones_census.py` | V(t) and V(1/t) for every knot of Burton's census files. Knots are bucketed by crossing number and frontier width under a memory latch (`JONES_MAXCELLS`). Output is written on arrival: name, crossings, F, V(t), V(1/t). Counter-cases first: 3_1 must differ and 4_1 must agree. There is no reference check inside a run. |
| `calibrate_jones_tn.py` | Tests all 8 combinations of the open conventions against Sage's Jones polynomial (spherogram, `new_convention=False`). Uses 14 named knots plus 40 knots drawn at random (seed 7) from the first 5,000 rows of `14n-hyp.csv.bz2`, each at 3 random points. |
| `test_gpu_vs_ref.py` | Checks that the GPU batch equals the CPU reference knot by knot, at x and at 1/x. Samples from the first 20,000 rows of the given census file. |
| `timing.py` | the per-knot times quoted above |
| `extra_countercases.py` | Compares the point test with Sage's exact test J(q) = J(1/q) on knots whose answer is known. |

## Run
```
export BURTON_CENSUS_DIR=/path/to/burton/raw
python -I calibrate_jones_tn.py 40                                     # Python with SageMath + spherogram
python -I extra_countercases.py                                        # Python with SageMath + spherogram
python -I test_gpu_vs_ref.py $BURTON_CENSUS_DIR/14n-hyp.csv.bz2 200    # Python with CuPy + spherogram
python -I jones_census.py 4 out15 $BURTON_CENSUS_DIR/15n-hyp.csv.bz2   # Python with CuPy + spherogram; 4 = CPU workers
bzcat out15_jones.tsv.bz2 | awk -F'\t' '$4==$5{print $1}' | sort      # the knots symmetric at the point
```

## Checked on 2 November 2026 (logs in `logs/`)
- **Calibration.** Of the 8 combinations, exactly the two that are mirror images of each other agree with Sage on all 54 knots at all
  3 points. Swapping the smoothings, the writhe sign and t ↔ 1/t together describes the same evaluation. Every other combination
  fails. Fixed: A-smoothing = slots (0,1)(2,3), writhe sign +1, t = x⁻⁸ corresponds to Sage's q.
- **GPU = reference.** Both the plain path (`prepare`) and the fast path used by the census (`prepare_fast`), each at x and at 1/x,
  against the CPU reference: 200 knots with 14 crossings (800 checks, frontier width up to 12) and 300 knots with 19 crossings
  (1,200 checks, frontier width up to 14), 0 mismatches (`logs/test_gpu_vs_ref_20261102.log`, `logs/test_gpu_vs_ref_19n_20261102.log`).
  Inside a census run there is no reference check.
- **15-crossing census.** `jones_census.py` with 4 workers on `15n-hyp.csv.bz2` (the final code of this folder), 168,023 knots,
  136 s, frontier width up to 14 (`logs/nh15final.log`, `logs/nh15final_stats.json`).
  - 167,525 knots have V(t) ≠ V(1/t), which makes them chiral if the kernel is correct.
  - 498 knots are symmetric at the point (`logs/15n_hyp_jones_symmetric_at_point.txt`). They include 15nh_008268, the only
    amphichiral prime knot with 15 crossings (OEIS A052401(15) = 1; certified amphichiral in our θ deposit, DOI 10.5281/zenodo.23856098).
  - Two earlier runs the same day, with the first, slower tensor builder and with an intermediate version, gave the same set of
    rows (their logs are not shipped).
- **Further counter-cases.** 9_42, 10_48, 10_71 (chiral, Jones-symmetric) and 8_17, 4_1, 6_3 (amphichiral) are symmetric at the
  point, and 3_1 is not. In every case this
  agrees with Sage's exact test J(q) = J(1/q) (`logs/extra_countercases_20261102.log`).

"Symmetric at the point" means only V(t) = V(1/t) mod p at one t. It does not mean V(t) = V(1/t) as polynomials, and it does not mean
amphichiral. The Jones polynomial misses the chirality of many chiral knots, 9_42 being the classical example. The HOMFLY
polynomial, of which Jones is a specialisation, is at least as strong.

## How these files relate to the authors' originals
They are the authors' files with this repository's licence header; the logs are the authors' logs.
