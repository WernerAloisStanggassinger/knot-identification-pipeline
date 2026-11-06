# Kauffman polynomial F(a, z)

The two-variable Kauffman polynomial F(a, z) of a knot or link, computed from a PD code in pure Python with exact integers.

**Why we built it.** We needed F for 12- and 13-crossing knots, and in our own hands, not only as a table entry. Our earlier
attempts did not reach that far (one stopped at braid index 4), and a table value cannot be checked by the table itself.

## Algorithm
Kauffman's descending-diagram recursion (L. H. Kauffman, An invariant of regular isotopy, Trans. Amer. Math. Soc. 318 (1990) 417–471),
the same scheme as `Kauffman[]` in the Mathematica package KnotTheory`:
- L is the regular-isotopy invariant with L(unknot) = 1, L(positive kink) = a·L and L(X+) + L(X−) = z·(L(X0) + L(X∞)).
- A crossing is *good* if the traversal (base point, direction and order of the components fixed) meets it first on the over strand.
  A diagram without bad crossings is descending, i.e. an unlink with c components and writhe w: L = a^w · μ^(c−1), μ = (a + 1/a)/z − 1.
- Otherwise the first bad crossing is resolved by L(D) = −L(D switched) + z·(L(D0) + L(D∞)); the switched diagram has one bad crossing
  fewer, the two smoothings one crossing fewer.
- F(a, z) = a^(−w(D)) · L(D).

The cost grows like 3^n in the number n of crossings (about 5·10^5 leaves at n = 12); intermediate diagrams are memoised on their
canonical PD code.

## Conventions
PD codes as in KnotInfo and KnotTheory`: X[i, j, k, l], i the incoming under strand, labels counterclockwise. The orientation of the
over strand is read from the traversal. The sign of a crossing follows the right-hand rule. The mirror convention a ↔ 1/a is fixed by
the calibration below.

## What has been checked
- `--selftest` (no network; two reference values, F of 3_1 and 4_1, are copied from KnotInfo into the script), nine checks:
  the unknot with one positive kink gives F = 1; the trefoil gives KnotInfo's F (this fixes the mirror convention) and F ≠ its mirror
  (chiral); the figure-eight knot gives KnotInfo's F and F = its mirror (amphichiral); the trefoil with one crossing switched (an
  unknot diagram) gives F = 1 and, as counter-case, F ≠ F(trefoil); the parser reads KnotInfo's string for 3_1; the a-span of F(trefoil) is 3.
- `--calibrate 10`: F of every prime knot with at most 10 crossings (249 knots) from its KnotInfo PD code, compared with KnotInfo's
  stored F — 249 of 249 equal, 0 different (run 6 November 2026, 131 s; needs the Python package `database_knotinfo`).

## Run
```
python3 kauffman_f.py --pd "[[1,5,2,4],[3,1,4,6],[5,3,6,2]]"   # trefoil: prints F and its a-span
python3 kauffman_f.py --selftest
python3 kauffman_f.py --calibrate 10                         # needs database_knotinfo (e.g. in a SageMath environment)
python3 kauffman_f.py --knotinfo 3_1 8_19                     # F from KnotInfo's PD code, compared with KnotInfo's F
```
Polynomials are dictionaries {(i, j): coefficient}, meaning coefficient · a^i · z^j.

## Licence
GPL-3.0-or-later, like the rest of this repository.
