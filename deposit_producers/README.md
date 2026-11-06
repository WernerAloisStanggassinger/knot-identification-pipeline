# Deposit producers

These are the scripts that computed the numbers in our Zenodo deposits, where the deposits themselves do not include the code.
There is one folder per deposit, named by its record number. Each folder has its own README saying which script produced which
quantity, with the receipt line where a receipt names the script, and, where it applies, what is not shipped.

| Folder | Zenodo DOI | Content |
|---|---|---|
| `21389680_priority_c7_15331_atlas` | 10.5281/zenodo.21389680 | Points only to producers released elsewhere: `pipeline/precellipse_machine.py`, `pipeline/precellipse_atlas.py`, `pipeline/build_precellipse_atlas_md.py` and machine 01 |
| `21590760_amphicheiral_338` | 10.5281/zenodo.21590760 | 338 positive amphicheiral knots: curve generator and sieve, M/E solver, deposit builder, count checks |
| `21933898_u_lower_bounds` | 10.5281/zenodo.21933898 | Lower bounds for the unknotting number over a corpus of built curves (plus the witness search for upper bounds) |
| `21946255_theta_mutant_pairs` | 10.5281/zenodo.21946255 | θ on constructed Conway-mutant pairs |
| `22893039_rs_12n706` | 10.5281/zenodo.22893039 | Symmetric ribbon numbers: r_s(12n706) = 4. The 12n504 value of 5 was withdrawn on 17 October 2026; 12n504 is only a candidate. |
| `22975395_ak_volume_conjecture` | 10.5281/zenodo.22975395 | Deposit builder and FAMED calibration for the Andersen–Kashaev volume conjecture on 220 knots |

The code is the code as it was run. Only the comments and messages are translated into English, and every translated file was
checked to change nothing else (`tools/strict_equivalence.py`; per-file results in `TRANSLATION_CHECK.md`). Of the 56 scripts
here, 36 went through that check (the 37th checked file is held back). The other 20 were already in English and are
byte-identical to the files as run. One of them, `21933898_u_lower_bounds/myf_ebenen/nach_lauf2_zweite_schicht.sh`, had a
German log message in line 21; that message was translated on 6 November 2026, and nothing else in the file changed. The paths inside
the scripts are the authors' paths of 2026. The scripts are receipts, not maintained tools.
