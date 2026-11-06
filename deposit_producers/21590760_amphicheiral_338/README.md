# Producers for Zenodo 10.5281/zenodo.21590760
Deposit: "Priority deposit v5: 338 positive amphicheiral hyperbolic knots whose Conway polynomial does not split -- Part III, order 16 and a discrepancy in the rotation parameter" (concept DOI, v5; source `claude_setup_v2/ZENODO_DOIS.md:22`). These scripts produced results reported there but were not shipped with it (`HOCHLADEN_v5/README.md:28-32`, "not included"; `:221-223`, "the solver is not shipped").
- `sieb/forge_c.py` (curve generator with Conant's rotation parameter c), `sieb/breit_c.py` (sieve, deadline harness), `sieb/alex_roh_voll.py` (Alexander polynomial by point evaluation on the raw diagram), with helpers `sieb/breit.py`, `sieb/stufe_a_pd.py`, `sieb/alex_punkt.py`: generation and sieving of the candidates; deadline stages 90/600/7200 s (`PRIORITAET_DEPOSIT_20260730/BELEG_zahlen_deposit_v3.md:209`), 341 calibration cases of the interpolated Alexander polynomial (`BELEG_zahlen_deposit_v3.md:221`).
- `me_lauf_a3.py` + `me_rechnen.py`: the Corollary 7.3 M/E solver ("identity solvable 338/338", `HOCHLADEN_v5/README.md:157`; method `HOCHLADEN_v5/verification_log.md:108-114`). It writes `me_ergebnis_a3.json` (`me_lauf_a3.py:30`, field `loesung_existiert`), which `deposit_builder/bau_v3.py:41,66` turns into the record flag `me_loesbar`. The 23 knots added in v4 were run with a byte-identical copy in `konversationen/kawauchi_20260725/lc_strasse_20261019/`.
- `deposit_builder/bau_v3.py`: v3 data builder, 60 randomized `simplify()` runs per record (`BELEG_zahlen_deposit_v3.md:113,214`). `deposit_builder/wache_volumen.py`: the volume guard, sample of 25 (`BELEG_zahlen_deposit_v3.md:216`).
- `conant_rechnung/fristzaehlung.py`: count of candidates that were never evaluated (2586, `BELEG_zahlen_deposit_v3.md:146`).
- `conant_rechnung/turm.py` + `turm_a4_keyed.py` (run with `A_MIN=2`), with helpers `turm_lauf.py` and `mod4_test.py`: the Theorem 7.1 quotient tower (P1/P2/P3 315/315, `konversationen/conant/BELEG_zahlen_nachfassen_20260729.md:147-158`). `turm_lauf.py` on its own is the earlier 250-knot run (`konversationen/conant/rechnung/BEFUND_lc_mod4_mechanismus_20260801.md:7`). The deposit itself mentions the tower only in `HOCHLADEN_v5/NOTE.md:112`.
Not here: the second Alexander computation, the symmetry-type reading and the splitting decision. No receipt names their scripts. The v4/v5 data path (`V4_338_20261019/`, `lc_strasse_20261019/`) is out of scope.
Hard-coded absolute paths (`/media/squeezy/...`, `/home/squeezy/...`) have been left as they are: this is the code as it was run.

## Where the files are in this repository
The scripts are in this folder, as they were run (paths inside them are the authors' paths of 2026, kept so the code is
the code that produced the deposit; they are receipts, not maintained tools):
- `conant_rechnung/fristzaehlung.py`
- `conant_rechnung/mod4_test.py`
- `conant_rechnung/turm.py`
- `conant_rechnung/turm_a4_keyed.py`
- `conant_rechnung/turm_lauf.py`
- `deposit_builder/bau_v3.py`
- `deposit_builder/wache_volumen.py`
- `me_lauf_a3.py`
- `me_rechnen.py`
- `sieb/alex_punkt.py`
- `sieb/alex_roh_voll.py`
- `sieb/breit.py`
- `sieb/breit_c.py`
- `sieb/forge_c.py`
- `sieb/stufe_a_pd.py`
