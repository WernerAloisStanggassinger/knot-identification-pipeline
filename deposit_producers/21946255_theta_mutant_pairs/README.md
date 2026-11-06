# Producers for Zenodo 10.5281/zenodo.21946255 — "Theta on constructed Conway-mutant pairs: thirty separated, one not"

Translated copies of the scripts that produced the deposited pairs; hard-coded paths are unchanged.

- `mutanten_stufe1.py` — builds the pretzel mutant pairs P(2,b,c,d) vs P(2,c,b,d), calibrates on P(2,3,5,7)=247 and P(2,3,5,3)=123, and sieves (hyperbolic, isometry_signature different, Alexander/Jones/det equal); writes `mutanten_kandidaten.json`. Receipt: `PRIORITAET_DEPOSIT_20260814/BELEG_zahlen_deposit_theta.md`:10 (31 parameter sets), :20-23 (calibration, faulty first build det 17).
- `invarianz_probe.py` — theta invariance probe: knots built under two parameter sets give the same theta value (28/28, 9 knots). Receipt: BELEG:24-25.

- `mutanten_stufe2_theta.py` — theta at three points on the de-duplicated pairs; writes `ergebnis_mutanten_theta.json` (BELEG:18, 30/30 separated at points).
- `nachziehen_homfly_symbolisch.py` — HOMFLY via Sage libhomfly (BELEG:16-17) and symbolic theta comparison; writes `ergebnis_nachziehen.json`, verdict `echte_trennung` (BELEG:19).

Producer links for the last two: found by grep for the result file names (`json.dump` targets), not stated in a receipt. Third-party input: Bar-Natan/van der Veen `Theta.sage` / `Theta_numeric.sage`, and `pd_to_theta_input.py` from the pipeline (hard-coded path `freeze_manifest/pipeline_scripts`).

## Where the files are in this repository
The scripts are in this folder, as they were run (paths inside them are the authors' paths of 2026, kept so the code is
the code that produced the deposit; they are receipts, not maintained tools):
- `invarianz_probe.py`
- `mutanten_stufe1.py`
- `mutanten_stufe2_theta.py`
- `nachziehen_homfly_symbolisch.py`
