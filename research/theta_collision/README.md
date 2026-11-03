# θ on volume-equal knots

Which pairs of census knots with the same hyperbolic volume does the Bar-Natan–van der Veen invariant θ tell apart?

Mutation preserves the hyperbolic volume (D. Ruberman, *Mutation and volumes of knots in S³*, Invent. Math. 90 (1987) 189–215), so every
pair of hyperbolic mutant knots has equal volume. This tool takes θ at the point (T1, T2) = (22/7, 21/13), computed modulo a
prime by `research/theta_census_gpu` for every knot whose diagram it could convert, and the hyperbolic volumes from `research/volume_census`. It forms groups of
knots whose floating-point volumes agree within a tolerance, and for every pair in a group it records whether the θ values prove
the two knots different.

**Why we built it.** We know of no program that combines θ over whole census files with volumes over whole census files, and
comparing pairs by hand does not scale: at 15 crossings there are 36,194 pairs within the tolerance among the non-alternating
hyperbolic knots.

## What the answer means
- The tool compares θ_raw = θ / (Δ(T1) Δ(T2) Δ(T1T2)) (see `research/theta_census_gpu`), a knot invariant, for a knot A and its
  mirror image mA, and for B and mB. The pair counts as **separated** when the values show both A ≠ B and A ≠ mB:
  (θ_raw(A) ≢ θ_raw(B) or θ_raw(mA) ≢ θ_raw(mB)) and (θ_raw(A) ≢ θ_raw(mB) or θ_raw(mA) ≢ θ_raw(B)), modulo p.
- A difference modulo p proves that θ_raw differs, because reduction modulo p is a ring homomorphism, **given that the θ kernel and
  the diagram conversion are correct** (tested, not proven; see `research/theta_census_gpu`). Then A ≠ B and A ≠ mB as
  **oriented** knots: each census diagram fixes one orientation, and ruling out A = −B as well would need θ(−K) = θ(K), which
  Bar-Natan and van der Veen state as Conjecture 22 of arXiv:2509.18456v4. (Census knots are distinct anyway; what the separation adds is that θ_raw — θ itself where the
  Alexander polynomials agree — tells them apart as oriented knots.)
- θ_raw contains the Alexander polynomial. For pairs with equal Alexander polynomial, such as mutants, a difference in θ_raw is a
  difference in θ itself. For other pairs it may come from Δ alone.
- Agreement modulo p means only "equal at this point, modulo p". It is no proof that θ is equal, and no proof that the knots are
  mutants.
- A pair within the tolerance is not necessarily a mutant pair, or even a pair of equal volume: the volumes are SnapPy numbers,
  not verified, and the volume tables store them with 15 significant digits (about 13 decimals for volumes between 10 and 100).
- Census names ignore mirror images, so the comparison is mirror-blind.

## The tolerance
We use 1e-9. It is a choice. At 15 crossings the gaps between neighbouring volumes in sorted order are as below. Gaps under about
1e-13 are at the printing resolution of the tables: most of them are exactly 0, the rest are one unit in the last printed digit,
and whether such a gap lands in the first or the second row is floating-point rounding. The table counts only gaps below 1e-8, and it was measured on our volume run of 1 November 2026; the
released `volume_census.py` is not seeded, and in its re-run the last printed digit differed on some rows (see
`research/volume_census/README.md`), so the two lowest rows need not reproduce exactly; such a re-run can also change which rows
carry flag `pos` (one 15-crossing knot lost it there), and with it the knot count.

| gap | count |
|---|---|
| below 1e-13 | 26,229 |
| 1e-13 to 1e-12 | 232 |
| 1e-12 to 1e-11 | 0 |
| 1e-11 to 1e-10 | 1 |
| 1e-10 to 1e-9 | 1 |
| 1e-9 to 1e-8 | 19 |

The gap of 4.4e-11 lies between 15nh_065614 and 15nh_065615; 15nh_065613 agrees with 15nh_065614 in every printed digit, so
the three form one group. The gap of 1.68e-10 lies between 15nh_154383 and 15nh_154384. θ separates all four pairs in these two
groups, so all five knots are different. A tolerance near 1e-11 would drop three of the four pairs; 15nh_065613 / 15nh_065614
would stay. We keep 1e-9 because a group that is too large costs only some
extra pairs, which θ (or another invariant) then separates, while a group that is too small can lose a real pair.

## Files
| File | What it does |
|---|---|
| `collide.py` | Groups by volume and counts pairs, as described below. `--selftest` runs 22 checks on synthetic data. |
| `../theta_census_gpu/theta_gpu_census.py --values` | Writes `<prefix>_values.tsv.bz2`: name, diagram crossings, θ_raw(K) and θ_raw(mirror K) modulo p (or `None`), for every knot whose diagram could be converted. Conversion errors appear only in the candidates file. |
| `../volume_census/volume_census.py` | Writes a volume table (name, volume, flag, attempts). Its output file name is an argument; we name them `vol_<c>n.tsv`. |

How `collide.py` works:
- It sorts the knots by volume. Neighbours closer than the tolerance join one group, so groups can chain.
- Every pair inside a group is counted. A name that appears in several volume tables is used once.
- Volume rows (flag `pos`) whose name has no θ value are left out and counted as `volume_rows_without_theta`;
  `knots_with_volume` counts the knots that have both.
- It writes `<out>_stats.json` (totals, and a breakdown by the unordered pair of crossing numbers) and `<out>_not_separated.jsonl`
  (the pairs θ does not separate, each with `theta_zero_both`, which looks at θ_raw(K) of both knots).

## Run
```
python -I collide.py --selftest
python -I ../theta_census_gpu/theta_gpu_census.py --values 8 c15 /path/to/burton/raw/15*.csv.bz2      # gnottheorie_env (CuPy)
python -I collide.py c15_values.tsv.bz2 1e-9 coll15 vol_15n.tsv
```
Several values files can be given as one comma-separated argument. The volume tables of several crossing numbers can be given
together, and then pairs across crossing numbers are counted too. Only volume rows with flag `pos` are used: a solution with all
tetrahedra positively oriented. Burton's census files and our volume tables are not included; see the README of
`research/theta_census_gpu` for the census files.

## Checked on 3 November 2026, with the programs in this folder
The inputs (census files, volume table, values file) are not in the repository.
- `collide.py --selftest`: 22 of 22 checks pass. They cover the counts (groups, separated, not separated, `none`, knots with
  θ, knots with a volume, volume rows without θ); each kind of pair (equal θ, θ equal to the other knot's mirror, θ = 0 on both
  sides, θ = 0 on one side only, different θ, and the two half-matches that still prove A ≠ B and A ≠ mB); chaining (three knots
  whose ends are 1.2e-9 apart form one group), a gap of 2e-9 and a gap exactly equal to the tolerance (not grouped); the key of
  the crossing numbers (smaller first, numeric: `9-10`, `17-18`); a row without flag `pos`; a name in two volume tables and in two
  values files (used once; conflicting θ values of a knot in a volume group raise); the built-in counter-case, applied and passing, and a failing one, which
  raises.
- `theta_gpu_census.py --values` (the released file, with the kink guard) on the five 15-crossing census files: 253,293 knots
  in 76 s on one RTX 3070 (run on the final files, 3 November 2026), no conversion error.
  - The candidates file is identical to the one of our census run of 31 October 2026 (the authors' script, without `--values`):
    the same 3 θ-zeros 15nh_008268, 15nh_048139 and 15nh_048140.
  - The 168,023 values of the 15nh knots are identical to those of our earlier census run (2 November 2026).
- `collide.py` on these values with our volume table of the non-alternating hyperbolic knots with 15 crossings:
  - 168,023 knots with a volume, 20,678 groups, 36,194 pairs within the tolerance.
  - θ separates 31,308 pairs (86.5 %) and does not separate 4,886.
  - The built-in counter-case passes: 15nh_048139 and 15nh_048140 (θ = 0, equal volume) are reported as not separated.

Runs on census files with more crossings are not part of this release.
