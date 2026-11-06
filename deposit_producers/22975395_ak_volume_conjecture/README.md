# Producers for Zenodo 10.5281/zenodo.22975395 — "The Andersen–Kashaev volume conjecture for 220 knots beyond the 14-crossing census of Ben Aribi, Guilloux and Wong"

The FAMED checker itself is public in the deposit (`reproduce.py --full`). These translated copies are the parts that are not; hard-coded paths are unchanged.

- `deposit_builder/bau_deposit_ak.py` — builds the deposit from `famed/deposit_frisch.jsonl`: asserts FAMED/geometric/n1=0/H1=Z, fresh HT-census exclusion (c_lb >= 15 or `identify` on geometric retriangulations, positive controls K14n1207, K12n706), isometry duplicates and mirror pairs; cuts the checking code into `reproduce.py`. Receipt: `PRIORITAET_DEPOSIT_20260926_AK/BELEG_zahlen_deposit_ak_20260926.md`:8-10.
- `famed/identify_eichung.py` — calibration of `HTLinkExteriors.identify` on 300 census knots. Receipt: BELEG:11.
- `famed/gegenfall_4.py` — counter-case: 15/15 perturbations rejected by condition (4). Receipt: BELEG:15. It `exec`s `famed_unabhaengig.py` (the public checker) from the working directory, and `bau_deposit_ak.py` cuts `mit_stoerung` out of it by AST.

Not shipped: the producer of `deposit_frisch.jsonl` (an input).

## Where the files are in this repository
The scripts are in this folder, as they were run (paths inside them are the authors' paths of 2026, kept so the code is
the code that produced the deposit; they are receipts, not maintained tools):
- `deposit_builder/bau_deposit_ak.py`
- `famed/gegenfall_4.py`
- `famed/identify_eichung.py`
