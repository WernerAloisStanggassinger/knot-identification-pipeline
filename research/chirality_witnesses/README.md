# Chirality witnesses and a chirality / amphichirality certificate

Tools that decide, or try to decide, whether a knot is chiral, i.e. not isotopic to its mirror image. We built them for knots on
which our θ census at a point gives θ = 0: there θ cannot tell a knot from its mirror image, and no single invariant we know of
separates every chiral knot from its mirror image (SnapPy's numerical symmetry group decides many hyperbolic cases, but not
certified).

There are two kinds of tool:
- **Witnesses** compute an invariant whose behaviour under mirroring is known: either its values on K and on the mirror image mK
  are compared (Σ3..Σ6, HOMFLY, knot Floer, Khovanov), or K's value is tested against what amphichirality would force
  (Σ2, Tristram–Levine, Chern–Simons, writhe), or a diagram criterion is tested (flypes). A difference proves K chiral; agreement
  proves nothing.
- **The certificate** tries to prove either answer for hyperbolic knots.

Every witness module except `khovanov_witness.py` has a `counter` mode with counter-cases (knots whose answer is known); for Khovanov
the counter-case PDs come from `make_kh_pd.py counter` (see Run). Run the counter-cases first: the name-list mode does not run
them by itself.

## Files
| File | Invariant / method | Rigour of a "CHIRAL" answer |
|---|---|---|
| `certify_amphichiral_chiral.py` | **Amphichiral:** an orientation-reversing automorphism of SnapPea's proto-canonical triangulation, plus an exact homeomorphism to the exterior (Regina: simplify both, then a common isomorphism signature, searching the exterior's simplified triangulation at retriangulation height 1, then 2, and the other at height 1 — height = extra tetrahedra over that triangulation; the log field `homeo_pachner_distance_bound` records 0, the first height, or that height + 1, not a number of moves); Gordon–Luecke. **Chiral:** SnapPy's verified canonical retriangulation (interval arithmetic) has no orientation-reversing automorphism; Mostow. If that verified retriangulation does have one: **amphichiral** (it triangulates the exterior itself). Otherwise **OPEN**. | exact / verified |
| `sigma2_linking_form.py` | Linking form of the double branched cover Σ2 (from the Seifert matrix, Sage), multiset of q(x) = λ(x,x) compared with its negative | exact |
| `cyclic_covers.py` | Σ3..Σ6: n-fold cyclic branched covers (SnapPy orbifold cover, Regina homological data: Kawauchi–Kojima invariants of the torsion linking form), K vs mK. Where H1(Σn(K)) is finite, it is checked against the Seifert-matrix formula (not for mK, and not when H1 has positive rank). | exact (Regina) |
| `tristram_levine.py` | Tristram–Levine signatures at the midpoints between unit-circle roots of Δ (mpmath, 60 digits); a value is used only if the smallest eigenvalue modulus exceeds 1e-20 | numerical |
| `chern_simons_witness.py` | Chern–Simons invariant (SnapPy `complex_volume(verified_modulo_2_torsion=True)`); amphichiral forces CS ∈ {0, π²/4} mod π²/2 | verified (interval) |
| `homfly_witness.py` | HOMFLY polynomial of K vs mK (Sage, libhomfly) | exact |
| `hfk_witness.py` | Knot Floer homology of K vs mK (SnapPy, Szabó's program): ranks, τ, ε | exact (program output) |
| `khovanov_witness.py` + `make_kh_pd.py` | Khovanov homology over Z (khoca, Lewark), also sl(3) over Z and Bar-Natan's deformation (`KH_ALG`), K vs mK | exact (program output) |
| `writhe_alternating.py` | Writhe of a reduced alternating diagram, which is a knot invariant (Thistlethwaite / Murasugi). The diagram must have no nugatory crossing; this is checked on the Gauss word | exact; reduced alternating diagrams only |
| `flype_signed_graph.py` | Signed interlace graph of a reduced alternating diagram vs the sign-swapped graph (Sage and networkx, both must agree). Rests on the flype lemma for the interlace graph in our A343358 deposit (DOI 10.5281/zenodo.23855791), on the argument in the module that flypes carry the crossing signs along (checked only by the counter-cases), and on Menasco–Thistlethwaite (prime knots) | combinatorial; reduced alternating diagrams of prime knots |
| `flype_orbit.py` | Is the mirror diagram in the complete flype orbit of the diagram (Menasco–Thistlethwaite, prime knots)? Flypes as rewrites of the signed Gauss word. A CHIRAL answer is only as sound as the completeness of this flype generator; that is tested only as a necessary condition, by the 20 amphichiral counter-cases (7 of them have a flype orbit larger than one diagram), and the field `complete` means only that the orbit search hit no size cap; "mirror in orbit" is not used as a proof | combinatorial; reduced alternating diagrams of prime knots |
| `census_codes.py` | Shared driver: DT codes from Burton's census files, counter-case mode, name-list mode; `is_reduced` (no nugatory crossing, from the Gauss word) | — |

Each witness's comments give the argument why "differs" implies "chiral". Several witnesses rest on a program's output (HOMFLY,
HFK, Khovanov): a "CHIRAL" from them is exactly as reliable as that program.

## Run
```
export BURTON_CENSUS_DIR=/path/to/burton/raw          # Burton's <c>{a,n}-{hyp,satellite,torus}.csv.bz2
python -I hfk_witness.py counter                      # counter-cases; exit code 1 if one fails
python -I hfk_witness.py names.txt out.jsonl          # one JSON line per census name, written on arrival, resumable
python -I make_kh_pd.py names.txt pds.json            # Khovanov: PDs in a SnapPy environment ...
KH_PD=pds.json python -I khovanov_witness.py NAME     # ... then khoca in its own environment (one knot per call)
python -I make_kh_pd.py counter kh_counter.json       # Khovanov counter-cases: 3_1, 4_1, 9_42, 10_125
```
Needed (SageMath, SnapPy and Regina importable in one Python where several are named): SnapPy for every module except
`khovanov_witness.py` (khoca only; `hfk_witness.py` uses SnapPy's knot Floer homology add-on); SageMath for `sigma2_linking_form.py`,
`cyclic_covers.py`, `tristram_levine.py` (also mpmath), `chern_simons_witness.py`, `homfly_witness.py`, `flype_signed_graph.py` (also
networkx) and `certify_amphichiral_chiral.py` (its CHIRAL branch uses SnapPy's verified computation, which needs Sage; without Sage that
branch returns OPEN);
Regina (7.4 used) for `cyclic_covers.py`, `flype_orbit.py`, `certify_amphichiral_chiral.py` and the square-knot counter-case of
`writhe_alternating.py`. The flype witnesses assume a prime knot (true for census knots); for other diagrams primality is the
caller's responsibility. Versions used here: SageMath 10.9, SnapPy 3.3.2 (knot_floer_homology 1.2.2), Regina 7.4, networkx 3.7,
mpmath 1.4.1, khoca 1.5 (pip). The sl(3) and Bar-Natan options of `khovanov_witness.py` (`KH_ALG`) are not covered by the logs here.

## Checked on 2 November 2026, with the files in this folder (logs in `logs/`)
- **Counter-cases.** Every module with a `counter` mode passes its own (`logs/counter_*.log`, re-run on the final files of this
  folder). The Khovanov log has no automatic pass/fail line; its four verdicts (3_1, 9_42, 10_125 differ, 4_1 equal) were checked
  by reading it.
  - Each module has its own list (see the logs): chiral knots it must detect, among 3_1, 5_2, 9_42, 10_125, 10_1, 10_48,
    10_71, 10_104 (`cyclic_covers.py` only reports its chiral cases; the certificate must give CHIRAL or OPEN, never
    AMPHICHIRAL, for 5_2 and 15nh_048139); amphichiral knots it must not, among 4_1, 6_3, 8_3, 8_9, 8_12, 8_17, 8_18, 10_17, 10_99 and the 20 amphichiral
    alternating knots in `flype_orbit.py`; and known blind spots that must stay blind (the HOMFLY polynomial of 9_42 equals
    that of its mirror image).
  - `writhe_alternating.py` additionally rejects a 7-crossing alternating but not reduced diagram of the amphichiral square knot
    3_1 # m3_1 (writhe 1), which a weaker reducedness test let through.
  - The certificate gives 4_1 AMPHICHIRAL. It does not give AMPHICHIRAL for 5_2 or for 15nh_048139; for 15nh_048139 it answers
    CHIRAL.
- **The three 15-crossing census knots with θ = 0** (`logs/c15_*`, names in `logs/theta_zeros_c15.txt`). Our θ deposit (DOI 10.5281/zenodo.23856098) shows that
  15nh_008268 is amphichiral and that 15nh_048139 and 15nh_048140 are chiral.
  - The certificate reproduces all three verdicts.
  - The Chern–Simons witness detects the two chiral knots.
  - The linking forms of Σ2..Σ6, the Tristram–Levine signatures, HOMFLY, knot Floer homology and Khovanov homology over Z are all
    blind on these three knots.
  - The writhe witness does not apply, since these are not alternating diagrams.

- **Name-list mode on alternating knots** (`logs/list_*.jsonl`, names in `logs/alt_test.txt`): the writhe and both flype
  witnesses call 15ah_00001 and 15ah_00002 chiral (SnapPy's numerical symmetry group agrees: `snappy_symmetry_check.py`, `logs/snappy_symmetry_alt_test.log`) and give no verdict for the
  non-alternating 15nh_008268.

## How these files relate to the authors' originals
Each witness module is the authors' script up to its `if __name__ == '__main__':`, verbatim apart from this repository's licence
header and a few comment lines: those named the authors' internal result files and run notes, and `flype_*.py` now cite the
A343358 deposit by DOI. Some of these comment edits also update content: the counter-case lists in the headers of
`sigma2_linking_form.py` (adds 8_18, which the original also ran; drops the note on the authors' calibration set) and
`cyclic_covers.py` (adds 9_42, 10_125; drops the note on the authors' calibration set), the construction note of
`cyclic_covers.py` (now names the orbifold route the code uses; an unlogged history note was dropped), an unlogged timing remark
in `khovanov_witness.py`,
and its remark on mutation (now: Viro for n = 2 only), the header of `homfly_witness.py` (dropped the authors' knot sets), the
"why this route", height and counter-case notes of the certificate (its unlogged shape-field remark is now marked as such), `flype_signed_graph.py` ("CALIBRATION" block became counter-cases; the
count of measured pairs of the lemma was dropped in favour of the DOI) and `flype_orbit.py` ("calibration" now reads
"counter-cases"). The other comments are
unchanged and still carry the authors' dates and some run notes.

The original main parts read the authors' run-specific result files. They are replaced by the shared driver in `census_codes.py`,
which runs the originals' counter-cases (with the same expectations) and takes any list of census names. 

Further differences:
- `flype_signed_graph.py` and `flype_orbit.py` get a verdict wrapper in their main part that also checks reducedness.
- `tristram_levine.py` and `homfly_witness.py` get a `test()` wrapper in their main part, with the same verdict rule as the
  original.
- `snappy_symmetry_check.py` is new (a cross-check, not a witness).
- `writhe_alternating.py` is new: a module version of the authors' run script `writhe_alt.py`, with the reducedness test replaced
  (the original used SnapPy `simplify('basic')`, which misses some nugatory crossings — the square-knot counter-case).
- `khovanov_witness.py` and `make_kh_pd.py` are the authors' scripts with the input files generalised (PD file and name list)
  and the comments updated.

Calibrations and results of the authors' own runs on larger knot sets are not part of this release.
