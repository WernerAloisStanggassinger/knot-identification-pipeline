# θ census on the GPU

Bar-Natan–van der Veen θ at the point (T1, T2) = (22/7, 21/13), for every knot of Burton's census files and for its mirror
image, computed exactly **modulo a prime** on the GPU. It is a filter: it proves θ ≠ 0 for almost every knot and hands the few
remaining ones to an exact computation over ℚ.

**Why we built it.** The established route is the Sage program of Bar-Natan and van der Veen (in our numerical version,
`Theta_numeric.sage`), exact over ℚ. On the 253,293 census knots with 15 crossings, our exact CPU run took about 10 minutes
(a run note; program, machine load and number of workers were not recorded, so this is no like-for-like comparison), this
filter 66 seconds including the CPU conversion of the diagrams (31 October 2026; 104 seconds in the re-run of 2 November 2026). The census files with 16 to 19 crossings
are much larger: 19 crossings = 294,130,458 knots, which the filter went through in 92,231 seconds on one RTX 3070.

## Why the filter misses nothing
`theta_roh` (θ_raw) is the quantity `Theta_numeric.sage` returns: θ_raw = θ / (Δ(T1) Δ(T2) Δ(T1T2)), the paper's θ₀
(arXiv:2509.18456v4, eq. (6)); `calibrate_theta_modp.py` checks the kernel against it. Δ has no zero at 22/7, 21/13 or 66/13
(argument in our θ deposit, see below), so θ = 0 at the point exactly when θ_raw = 0, and the matrices A(T) are invertible over ℚ.
Reducing modulo a prime p is a ring homomorphism as long as p divides no denominator; if some A(T) becomes singular modulo p, the
kernel returns `None` instead of a value. So θ = 0 implies θ_raw ≡ 0 mod p or `None`. A knot is flagged when θ_raw(K) ≡ 0, when either value (K or its mirror) is
`None`, when θ_raw(K) + θ_raw(mirror K) ≢ 0 mod p (the conjectured sign rule θ(mirror K) = −θ(K), Conjecture 21 of arXiv:2509.18456v4,
fails at the point; this also catches a mirror with θ_raw ≡ 0 while θ_raw(K) ≢ 0), or when its diagram could not be converted.
Hence, **given that the kernel and the diagram conversion are correct**, every exact θ-zero is among the flagged knots, and a
knot that is not flagged has θ ≠ 0. The modular step is a proof; the kernel and the conversion are tested, not proven: the
kernel against Sage over ℚ on 400 knots with 3 to 11 crossings, the whole chain against the exact zero set at 15 crossings
(below). p = 2147483629.

The flagged knots are then computed exactly over ℚ with our numerical version `Theta_numeric.sage` of Bar-Natan and van der
Veen's program (not included, see `THIRD_PARTY.md`).

## Files
| File | What it does |
|---|---|
| `theta_modp.py` | θ_raw at a point over GF(p): `ref()` pure Python (one knot, reference), `gpu_batch()` CuPy (a batch of knots with the same number of crossings) |
| `pd_to_theta_input.py` | PD code → the input `[Cs, phi]` of Bar-Natan and van der Veen's program; self-test on the trefoil and on a kinked PD, which must be refused |
| `theta_gpu_census.py` | the census driver: DT code → SnapPy link → PD → `[Cs, phi]` for K and its mirror on a CPU pool, then the GPU kernel in batches of 20,000 knots per crossing number; writes `<prefix>_candidates.tsv` and `<prefix>_stats.json`; with `--values` (since 3 November 2026) also `<prefix>_values.tsv.bz2` with θ_raw(K) and θ_raw(mirror K) mod p (or `None`) of every knot whose diagram could be converted (conversion errors appear only in the candidates file), the input of `research/theta_collision` |
| `certify_zeros.py` | for θ-zeros of **hyperbolic** knots, an attempt to certify (it can end in `NOT CERTIFIED`): verified oriented isometry signatures (SnapPy, `verified=True`, needs SageMath) of the exterior and of its mirror. Equal ⟹ amphichiral (an orientation-reversing isometry; Gordon–Luecke). Different ⟹ chiral (Mostow–Prasad: every self-homeomorphism of a hyperbolic exterior is homotopic to an isometry). Counter-cases 4_1 and 5_2 run first |
| `calibrate_theta_modp.py` + `ref_sage_22_7_21_13.json` | calibration: 400 knots with θ_raw over ℚ from Sage; reference and GPU must agree mod p |

## Run
```
python -I calibrate_theta_modp.py                                          # gnottheorie_env (CuPy)
python -I theta_gpu_census.py 8 out_c15 /path/to/burton/raw/15*.csv.bz2     # gnottheorie_env; env THETA_GPU_MEM_GB (default 6)
python -I theta_gpu_census.py --values 8 out_c15 /path/to/burton/raw/15*.csv.bz2   # same, plus out_c15_values.tsv.bz2
python -I certify_zeros.py zeros.tsv /path/to/burton/raw cert.jsonl 3 '15n-hyp.csv.bz2'   # sage_stable
```
Burton's census files (`<c>{a,n}-{hyp,torus,satellite}.csv.bz2`, columns `name`, …, `dt_code`) are not included; see
Benjamin A. Burton, *The next 350 million knots* (SoCG 2020). `zeros.tsv` is written from the exact recomputation of the flagged
knots: one row `name<TAB>crossings<TAB>θ` per knot, and `certify_zeros.py` takes the rows whose third column is exactly `0`
(the candidates file itself has the word `zero` there and cannot be passed directly).

## Checked on 2 November 2026, with the files of that release
- `calibrate_theta_modp.py`: reference vs Sage over ℚ 400/400 equal, GPU vs Sage 400/400 equal. Counter-case (T1 + 1 instead
  of T1): of the 380 knots with θ_raw ≠ 0 at the point, 0 agree; the other 20 have θ_raw = 0 at the point.
- `theta_gpu_census.py` on all five 15-crossing census files: 253,293 knots, 3 flagged (15nh_008268, 15nh_048139,
  15nh_048140, all θ_raw ≡ 0), no sign-rule miss, no `None`, no conversion error. Candidates file and statistics identical to
  the authors' run of 31 October 2026 (apart from the running time, 104 s against 66 s). The three are exactly the θ-zeros at
  15 crossings found by the exact evaluation at the point in our deposit *θ detects chirality for all prime knots up to 14 crossings — and
  first fails at 15*, Zenodo, DOI 10.5281/zenodo.23856098.
- `certify_zeros.py` on these three: 15nh_008268 amphichiral, 15nh_048139 and 15nh_048140 chiral, as in that deposit.

## Diagrams with a Reidemeister-I kink: a fault in our kernel, found and fixed (3 November 2026)
On diagrams with a Reidemeister-I kink (a PD crossing whose 4-tuple repeats a label), `theta_modp.py` gave a wrong θ_raw: of 160
enlargements of 6_2, K8n1, 4_1 and 7_4 by SnapPy's `backtrack()` (2 to 10 steps, seed 7), 99 were wrong. **Cause** (found by an
outside reviewer on 3 November 2026 and reproduced by us the same day): `_A()` and the GPU `build()` *set* the three matrix entries
of a crossing, `A[i][i+1]`, `A[i][j+1]`, `A[j][j+1]`, with `=` instead of *adding* them to the identity matrix as Bar-Natan and van
der Veen's `Theta.m` does (`A[[{i,j},{i+1,j+1}]] += …`). Every half-strand is the incoming strand of exactly one crossing, so the
only entry two terms can share is a diagonal one, and that happens exactly when i = j + 1 at a crossing, which only a kink produces;
there the 1 of the identity was overwritten. All 99 wrong diagrams have such a crossing, none of the 61 others is wrong. Since
3 November 2026 both functions add (`+=`); with that, all 160 diagrams give the θ_raw of the minimal diagram, on the CPU reference
and on the GPU, and `calibrate_theta_modp.py` still gives 400/400 on both. The conversion `pd_to_theta_input.py` was right all
along. (An earlier version of this paragraph blamed it: the Alexander polynomial we had computed from the crossing list as a check
used the same assigning matrix.)
**Published results are not affected:** a crossing with i = j + 1 needs a kink, and a minimal diagram has none. Checked directly:
none of the 313,230 converted census diagrams up to 15 crossings and none of the 839 PD codes of the authors' corpus (834 of its 843 records — the other
9 are unknots without a diagram — and its 5 validation knots) has such a crossing (3 November 2026). `pd_to_cs_phi` still refuses a diagram with a kink (`ValueError`; `theta_gpu_census.py` reports it as
a conversion error); the refusal is now a safety net for any kernel that still assigns. The test scripts for this paragraph are not
part of the repository.

## How these files relate to the authors' originals
All five scripts got the licence header of this repository (in `pd_to_theta_input.py` it replaces the deposit's header, which
pointed to `LICENSE-CODE.txt`). In the release of 2 November 2026, `theta_modp.py` and `pd_to_theta_input.py` were proven
EQUIVALENT to the originals by `tools/strict_equivalence.py` (only comments, docstrings and message texts changed). The other
three changed as follows in that release, and in nothing else:
- `theta_gpu_census.py`: imports its two modules from this folder instead of from two absolute paths; one docstring line names
  the calibration script instead of the authors' folder.
- `certify_zeros.py`: made from the c = 19 copy of the authors' script (the c = 17, 18 and 19 copies differ only in a first
  comment line and in the fixed name of the candidates file). The fixed file name became an optional file pattern argument; the
  first comment line was dropped and the header comment and usage line were rewritten for this release.
- `calibrate_theta_modp.py`: finds its own folder (`import os`, path of the reference file) instead of the parent folder; its
  messages and two variable names were translated; its counter-case line now also counts the knots with θ_raw ≠ 0 separately.

Since 3 November 2026 three files differ in code as well: `pd_to_theta_input.py` has the kink guard described above (a check and a
`raise` at the start of `pd_to_cs_phi`, and a check in its self-test that a kinked PD is refused; the same in
`pipeline/pd_to_theta_input.py`), `theta_gpu_census.py` has the `--values` option (argument parsing and the writer of the values
file), and `theta_modp.py` adds the matrix entries instead of setting them (`_A()` and the GPU `build()`, described above).