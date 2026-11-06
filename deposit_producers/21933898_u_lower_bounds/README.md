# Producers for Zenodo 10.5281/zenodo.21933898 — "Lower bounds for the unknotting number over a corpus of built curves"

These scripts computed the bound quantities shipped in the deposit (`u_bounds_data.json`); the deposit's own `reproduce.py` checks consistency only (deposit `reproduce.py`:7-18). Translated copies; hard-coded `/media/squeezy/...` paths are unchanged.

- `myf_ebenen/ernte_u_schranken.py` — harvests sigma, s, tau, g4 into u_lower / witness-backed u_obere, with counter-checks P1–P4 (receipt `PRIORITAET_DEPOSIT_20260813/BELEG_zahlen_deposit_u.md`:52-56).
- `myf_ebenen/wendt_sweep.py` (+ helper `wendt_probe.py`) — branched-cover bound mu(H_1(Sigma_2)) = field `branched_cover_mu`; check G1 (BELEG:57).
- `myf_ebenen/gegenprobe_goeritz.py` — independent Goeritz + Smith counter-check of that mu.
- `myf_ebenen/nakanishi.py` — Nakanishi index lower bound = field `nakanishi_index`; checks N1, N3 (BELEG:59-60; deposit README `HOCHLADEN_v3/README.md`:86, 111-112).
- `myf_ebenen/neue_obere_schranken2.py` (lauf2) and `myf_ebenen/zeugen_andere_diagramme.py` (lauf3, launched by `nach_lauf2_zweite_schicht.sh`) — witness search for upper bounds; lauf3 gives `Ndiagramme` 49 703 / `Nwechsel` 5.96 million (BELEG:35-36).
- `myf_ebenen/nachweis_obere_schranken.py` — re-proves exact u on the stored diagram.

Not shipped: `pruef_blatt_zahlen.py` (the sieve `sammeln()` that `bau_deposit_u.py`:39 imports for every deposit number; held back because its receipt template names the person a sheet was handed to), and the inputs these scripts read: `bau_deposit_u.py` (deposit builder), `bilanz.py` → `bilanz_daten.json`, `audit_unknotting_diagramm.py`, `neue_obere_schranken.py` (lauf1), the per-knot records, and `corpus_db.json` as it stood when the deposit was built (a later snapshot ships as `pipeline/corpus_db.json`; see its `meta` for the build date).

## Where the files are in this repository
The scripts are in this folder, as they were run (paths inside them are the authors' paths of 2026, kept so the code is
the code that produced the deposit; they are receipts, not maintained tools):
- `myf_ebenen/ernte_u_schranken.py`
- `myf_ebenen/gegenprobe_goeritz.py`
- `myf_ebenen/nach_lauf2_zweite_schicht.sh`
- `myf_ebenen/nachweis_obere_schranken.py`
- `myf_ebenen/nakanishi.py`
- `myf_ebenen/neue_obere_schranken2.py`
- `myf_ebenen/wendt_probe.py`
- `myf_ebenen/wendt_sweep.py`
- `myf_ebenen/zeugen_andere_diagramme.py`
