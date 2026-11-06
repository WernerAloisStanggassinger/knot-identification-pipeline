# Producers for Zenodo 10.5281/zenodo.22893039
**Symmetric ribbon numbers: r_s(12n706) = 4, a candidate for 12n504, and a computed subset of R^s_5.**
These are the original working files of the search, which the deposit does not include (deposit README, `PRIORITAET_DEPOSIT_20260922/HOCHLADEN_v2/README.md:235`). Receipts are in `maschinen_ideen/RS4_12N706_20260919/`. Comments and messages are translated; the code is unchanged (strict-equivalence gate).
- `enum_bogendiagramme.py`, `pruef_aglet_flaeche.py`, `abschluss.py`, `filter_det7.py`, `dt_bau.py`, `identifiziere_kandidaten.py` (+ `karten_analyse.py`, `gegentest_enum.py`, `det_aus_pd.py`): enumeration calibrated 1/2/8; aglet faces; closure det; det-7 filter; partial knot 5_2. Receipt: STAND.md:37,52,69,87,109.
- `band_rand.py`, `eiche_union.py`, `rs4_voll.py`, `jagd_12n706.py`, `kontrolle.py`: union constructor, R^s_4 = 27 exactly, witness for 12n706, controls 12/12 negative and 5/5 positive. Receipt: STAND.md:129; BELEG_rs_12n706_20260920.md:41-42,65,108-110.
- `jagd_fork.py`, `offene_restliche.py`, `komplement.py`, `fkp_stellen.py`: the remaining Table 3 cases, the twist-circle complement and FKP thresholds (finite search). Receipt: STAND.md:158-160,173.
- `rs5.py`: the computed subset of R^s_5 (145 polynomials, `rs5_ergebnis.json`). Receipt: BELEG_rs5_20260920.md:24.
- `zeugen5.py` (+ `sweep12.sh`, `sweep11_g5.sh`): c=5 witness search; the 12n504 witness (`zeuge_12n504.json`) is only a candidate, because r_s(12n504) = 5 was WITHDRAWN on 17.10.2026. Receipt: BELEG_rs_12n504_20260920.md:1,61-63.
- `zweite_implementierung.py`, `blasen_gegenprobe.py`, `konstruktor2.py`: second implementation of the enumeration and of the constructor, plus the check on the bubble diagrams. Receipt: BELEG_zweite_implementierung_20260922.md:11,77,93,112-115,126.
Note: most files still contain the absolute path `/media/squeezy/kauffman_storage/maschinen_ideen/RS4_12N706_20260919` (as `HIER`/`sys.path`); `rs5_ergebnis.json` and `paper_tabelle1.txt` are inputs that are not in this folder.

## Where the files are in this repository
The scripts are in this folder, as they were run (paths inside them are the authors' paths of 2026, kept so the code is
the code that produced the deposit; they are receipts, not maintained tools):
- `abschluss.py`
- `band_rand.py`
- `blasen_gegenprobe.py`
- `det_aus_pd.py`
- `dt_bau.py`
- `eiche_union.py`
- `enum_bogendiagramme.py`
- `filter_det7.py`
- `fkp_stellen.py`
- `gegentest_enum.py`
- `identifiziere_kandidaten.py`
- `jagd_12n706.py`
- `jagd_fork.py`
- `karten_analyse.py`
- `komplement.py`
- `konstruktor2.py`
- `kontrolle.py`
- `offene_restliche.py`
- `pruef_aglet_flaeche.py`
- `rs4_voll.py`
- `rs5.py`
- `sweep11_g5.sh`
- `sweep12.sh`
- `zeugen5.py`
- `zweite_implementierung.py`
