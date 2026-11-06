# Translation check of the deposit producers (5 November 2026)

Command per file: `python3 -I tools/strict_equivalence.py ORIGINAL TRANSLATED` → must print EQUIVALENT.
Result: 37 of 37 EQUIVALENT (36 of these files ship; one is held back, see the 21933898 README) (run by the translating agent, then re-run independently on the whole list: eq=37 bad=0). 3 files changed, 34 were already English. Per-file report of the translating agent follows.

# Translation wave B, deposit producers — 2026-11-05

Original tree: `<authors-open_source>/wave_b/port_20261105_ORIG`
Target tree: `<authors-open_source>/wave_b/port_20261105`
File list: `TRANSLATE_LIST_PRODUCERS_20261105.txt` (37 .py files, no .sh).
Gate per file: `strict_equivalence.py ORIGINAL TARGET` must print EQUIVALENT, plus `ast.parse` syntax check (no __pycache__ written).
Hard-coded /media/squeezy paths are kept "as run" and only listed.

## PASS — `deposit_producers/21590760_amphicheiral_338/conant_rechnung/turm_lauf.py`
- gate: EQUIVALENT; syntax: OK; edited: yes (docstring line 7 was a German duplicate of the English line 8; German line removed)
- German data strings (kept): dict keys/values `knoten`, `erzeugung`, `achsen_verschlingung_ell`, `kein_diagramm`, `nu1_unknoten`, `fehler`; status values `stufe{lvl}:konversion`, `stufe{lvl}:uneinig` (L73-78); output file `ergebnis_turm_20260729.json`
- /media/squeezy paths: none (paths via `site_config.ROOT`)
- own modules NOT in port_20261105: site_config (present: turm, mod4_test; forge_c loaded by path from ROOT/konversationen/kawauchi_20260725/sieb)
- private data / credentials: none

## PASS — `deposit_producers/21590760_amphicheiral_338/deposit_builder/bau_v3.py`
- gate: EQUIVALENT; syntax: OK; edited: yes (docstring L21 garbled sentence fixed; rest was already English)
- German data strings (kept): L202 `'fuer alle Records das beste gefundene Diagramm; Alexander gegen das v2-Diagramm geprueft'` (meta value written into the published JSON); German dict keys throughout (`knoten`, `erzeugung`, `kreuzungs_obergrenze`, `pruefungen`, `aussage`, `aenderung_gegen_v2`, ...); `quelle_fehlt` counter key; output file `counterexamples_data_ENTWURF.json`
- /media/squeezy paths: none; L212 print message contains `/home/squeezy/miniforge3/envs/gnottheorie_env/bin/python` (local interpreter path, kept)
- own modules NOT in port_20261105: site_config (third-party: regina, snappy)
- private data / credentials: none

## PASS — `deposit_producers/21590760_amphicheiral_338/deposit_builder/wache_volumen.py`
- gate: EQUIVALENT; syntax: OK; edited: no (comments, docstring and messages already English)
- German data strings (kept): dict keys `knoten`, `volumen`, `pd_code`, `kreuzungs_obergrenze_v2`; CLI flag `--selbsttest`; file names `counterexamples_data_ENTWURF.json`
- /media/squeezy paths: none (via `site_config.ROOT`)
- own modules NOT in port_20261105: site_config (third-party: snappy, spherogram)
- private data / credentials: none

## PASS — `deposit_producers/21590760_amphicheiral_338/me_lauf_a3.py`
- gate: EQUIVALENT; syntax: OK; edited: no (already English)
- German data strings (kept): L209 group labels `'a=2 (k=2, Ordnung 4)'`, `'a=3 (k=4, Ordnung 8)'`, `'a=4 (k=8, Ordnung 16)'` (dict values, glossed by a comment); status values `FEHLER ...`, `leer`, `lesefehler ...`; result keys `vorhersage_E0`, `treffer`, `loesung_existiert`, `E0_passt`, `sek`
- /media/squeezy paths: none
- own modules NOT in port_20261105: none (me_rechnen.py and sieb/breit.py loaded by path, both present)
- private data / credentials: none

## PASS — `deposit_producers/21590760_amphicheiral_338/sieb/alex_roh_voll.py`
- gate: EQUIVALENT; syntax: OK; edited: no (already English)
- German data strings (kept): verdict labels `BYTE-GLEICH`, `WEICHT` (printed with English gloss at L425/L430), reject reason codes raised as `ValueError('nicht_palindrom' / 'leiter_erschoepft' / 'delta1' / 'null' / 'verify')` that end up as `rej:` values in records; keys `urteil`, `gesamt_s`, `alex_neu`, `confirm_4320_punkt`; CLI flags `--selbsttest`, `--eiche`, `--sieb`
- /media/squeezy paths: none
- own modules NOT in port_20261105: none (breit_c.py, alex_punkt.py loaded by path, both present; third-party: numpy, regina)
- private data / credentials: none

## PASS — `deposit_producers/21590760_amphicheiral_338/sieb/breit_c.py`
- gate: EQUIVALENT; syntax: OK; edited: no (already English)
- German data strings (kept): `ORDER=verschraenkt` env value; `rej='frist'` value
- /media/squeezy paths: none
- own modules NOT in port_20261105: none (breit.py, forge_c.py loaded by path, both present; third-party: numpy, regina)
- private data / credentials: none

## PASS — `deposit_producers/21590760_amphicheiral_338/sieb/stufe_a_pd.py`
- gate: EQUIVALENT; syntax: OK; edited: no (already English)
- German data strings (kept): L82 reject value `'sigma nicht bitgenau'`
- /media/squeezy paths: none
- own modules NOT in port_20261105: site_config; via sys.path L19-20 (ROOT/maschinen_ideen, ROOT/freeze_manifest) — machine_common present (third-party: numpy, regina)
- private data / credentials: none

## PASS — `deposit_producers/21933898_u_lower_bounds/myf_ebenen/nachweis_obere_schranken.py`
- gate: EQUIVALENT; syntax: OK; edited: no (already English)
- German data strings (kept): verdict values L81 `'NICHT-GEPRUEFT (zu viele Kombinationen)'`, L87 `'AM RECORD NACHWEISBAR'` / `'aus dem Record nicht reproduzierbar'` (compared again at L95-97); keys `urteil`, `gefunden`, `sekunden`, `zeilen`, `ist_vertreter`
- /media/squeezy paths: none; docstring L29-30 usage line contains `/home/squeezy/miniforge3/envs/sage_stable/bin` (kept)
- own modules NOT in port_20261105: site_config (third-party: spherogram)
- private data / credentials: none

## PASS — `deposit_producers/21933898_u_lower_bounds/myf_ebenen/nakanishi.py`
- gate: EQUIVALENT; syntax: OK; edited: no (already English)
- German data strings (kept): keys `zeilen`, `fehler`, `treffer`, `u_obere`, `korang_bei_alpha_1`; CLI flag `--selbsttest`
- /media/squeezy paths: none
- own modules NOT in port_20261105: site_config (third-party: numpy, spherogram, sympy)
- private data / credentials: none

## PASS — `deposit_producers/21933898_u_lower_bounds/myf_ebenen/neue_obere_schranken2.py`
- gate: EQUIVALENT; syntax: OK; edited: no (already English)
- German data strings (kept): L79 `'Gruppe des Aussenraums = Z'`, L92 `f"backtrack(steps={s}) auf 0 Kreuzungen"`, L156 `f"erschoepfend: KEIN Zeuge bei k={k} im gespeicherten Diagramm"`, L158 `f"kein Zeuge in {BUDGET:.0f}s ({geprueft} von {anzahl}) — offen"`; verdict codes `'zeuge'`, `'widerlegt'`; dir names `lauf2`
- /media/squeezy paths: none; docstring L31-32 `/home/squeezy/miniforge3/envs/sage_stable/bin` (kept)
- own modules NOT in port_20261105: none (third-party: snappy, spherogram)
- private data / credentials: none

## PASS — `deposit_producers/21933898_u_lower_bounds/myf_ebenen/wendt_probe.py`
- gate: EQUIVALENT; syntax: OK; edited: no (already English)
- German data strings (kept): none of note (keys only)
- /media/squeezy paths: none
- own modules NOT in port_20261105: site_config (third-party: spherogram, sympy)
- private data / credentials: none

## PASS — `deposit_producers/21933898_u_lower_bounds/myf_ebenen/zeugen_andere_diagramme.py`
- gate: EQUIVALENT; syntax: OK; edited: no (already English)
- German data strings (kept): L150 `urteil=f"kein Zeuge in {diagramme} anderen Diagrammen ({geprueft} Proben) — offen"`; dir names `lauf2`, `lauf3`; key `zeilen`
- /media/squeezy paths: none
- own modules NOT in port_20261105: none (neue_obere_schranken2 present; third-party: snappy, spherogram)
- private data / credentials: none


## PASS — `deposit_producers/21946255_theta_mutant_pairs/mutanten_stufe1.py`
- gate: EQUIVALENT; syntax: OK; edited: no (already English)
- German data strings (kept): L117 `'kein Knoten (Verschlingung) — verworfen'`, L119 `'keine geometrische Loesung — verworfen'`, L121 `'SELBER Knoten (isosig gleich) — kein Prueffall'` (verdict values)
- /media/squeezy paths: none
- own modules NOT in port_20261105: none (third-party: snappy, spherogram)
- private data / credentials: none

## PASS — `deposit_producers/21946255_theta_mutant_pairs/nachziehen_homfly_symbolisch.py`
- gate: EQUIVALENT; syntax: OK; edited: no (already English)
- German data strings (kept): L77 `'Einheiten-Test nicht durchfuehrbar: %s'`, L80 `'Quotient kein Monom'` (result values)
- /media/squeezy paths: none
- own modules NOT in port_20261105: site_config, pd_to_theta_input (via sys.path L37 ROOT/freeze_manifest/pipeline_scripts) (third-party: sage)
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/abschluss.py`
- gate: EQUIVALENT; syntax: OK; edited: no (already English)
- German data strings (kept): none
- /media/squeezy paths: none
- own modules NOT in port_20261105: site_config (pruef_aglet_flaeche present; sys.path L43 ROOT/maschinen_ideen/RS4_12N706_20260919)
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/band_rand.py`
- gate: EQUIVALENT; syntax: OK; edited: no
- German data strings (kept): none
- /media/squeezy paths: none
- own modules NOT in port_20261105: site_config (abschluss, dt_bau, pruef_aglet_flaeche present)
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/blasen_gegenprobe.py`
- gate: EQUIVALENT; syntax: OK; edited: no (already English)
- German data strings (kept): L105 return value `f"nicht bestimmt ({e})"`
- /media/squeezy paths: none
- own modules NOT in port_20261105: none (abschluss, band_rand, eiche_union present; third-party: sympy)
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/det_aus_pd.py`
- gate: EQUIVALENT; syntax: OK; edited: no
- German data strings (kept): none
- /media/squeezy paths: none
- own modules NOT in port_20261105: none
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/dt_bau.py`
- gate: EQUIVALENT; syntax: OK; edited: no
- German data strings (kept): L98-99 labels `'nicht-alt A'`, `'nicht-alt B'`
- /media/squeezy paths: none
- own modules NOT in port_20261105: site_config (sys.path L30 ROOT/maschinen_ideen/RS4_12N706_20260919; third-party: spherogram)
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/eiche_union.py`
- gate: EQUIVALENT; syntax: OK; edited: no
- German data strings (kept): none
- /media/squeezy paths: none
- own modules NOT in port_20261105: site_config (abschluss, band_rand present)
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/enum_bogendiagramme.py`
- gate: EQUIVALENT; syntax: OK; edited: no
- German data strings (kept): none
- /media/squeezy paths: none
- own modules NOT in port_20261105: none
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/filter_det7.py`
- gate: EQUIVALENT; syntax: OK; edited: no
- German data strings (kept): none
- /media/squeezy paths: none
- own modules NOT in port_20261105: site_config (abschluss, pruef_aglet_flaeche present; sys.path L19 ROOT/maschinen_ideen/RS4_12N706_20260919)
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/fkp_stellen.py`
- gate: EQUIVALENT; syntax: OK; edited: no — an edit was tried and reverted: the embedded child-code template `KIND` (L63-104) contains the German comment L72 `# erst damit funktioniert spherograms exterior()` ("only with this import does spherogram's exterior() work"). The template holds `%(hier)r` placeholders, so the gate cannot parse it as code and compares it byte-for-byte; translating the comment made the gate print DIFFERENT. File copied back from ORIGINAL (cmp identical), gate EQUIVALENT.
- German data strings (kept): L72 comment inside the `KIND` template (above); status values in the template `'kein_dt'`, `'nicht_hyperbolisch'`; keys `hyperbolisch`, `vol_M_intervall`, `fehler`
- /media/squeezy paths: none
- own modules NOT in port_20261105: none (band_rand, eiche_union, pruef_aglet_flaeche, komplement present; third-party: snappy, spherogram)
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/gegentest_enum.py`
- gate: EQUIVALENT; syntax: OK; edited: no
- German data strings (kept): none
- /media/squeezy paths: none
- own modules NOT in port_20261105: site_config (enum_bogendiagramme present; third-party: networkx)
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/identifiziere_kandidaten.py`
- gate: EQUIVALENT; syntax: OK; edited: no
- German data strings (kept): none
- /media/squeezy paths: none; docstring L18-19 `/home/squeezy/miniforge3/envs/sage_stable/bin` (kept)
- own modules NOT in port_20261105: site_config (abschluss, dt_bau, pruef_aglet_flaeche present)
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/jagd_12n706.py`
- gate: EQUIVALENT; syntax: OK; edited: no
- German data strings (kept): none
- /media/squeezy paths: none
- own modules NOT in port_20261105: site_config (abschluss, band_rand, dt_bau, pruef_aglet_flaeche present; third-party: snappy, spherogram)
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/jagd_fork.py`
- gate: EQUIVALENT; syntax: OK; edited: no
- German data strings (kept): status values `'kein_dt'`, `'nicht_hyperbolisch'`; keys `aufgaben`, `erledigt`, `treffer`, `stufen`, `wort`
- /media/squeezy paths: none
- own modules NOT in port_20261105: site_config (band_rand, dt_bau, eiche_union, pruef_aglet_flaeche present; third-party: snappy, spherogram)
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/karten_analyse.py`
- gate: EQUIVALENT; syntax: OK; edited: no
- German data strings (kept): none
- /media/squeezy paths: none
- own modules NOT in port_20261105: site_config (enum_bogendiagramme present; sys.path L19 ROOT/maschinen_ideen/RS4_12N706_20260919)
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/komplement.py`
- gate: EQUIVALENT; syntax: OK; edited: no
- German data strings (kept): none
- /media/squeezy paths: none
- own modules NOT in port_20261105: site_config (band_rand, pruef_aglet_flaeche present)
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/konstruktor2.py`
- gate: EQUIVALENT; syntax: OK; edited: no (already English; the docstring quotes Werner's German instruction verbatim with an English rendering, kept on purpose)
- German data strings (kept): `'ALLE_KARTEN'`, `'KEIN_LINK'`, `'FEHLER:'`, `'kein DT'`; CLI word `eichen`, `--laut`
- /media/squeezy paths: none; docstring L35-36 `/home/squeezy/miniforge3/envs/sage_stable/bin` (kept)
- own modules NOT in port_20261105: none (zweite_implementierung present; third-party: networkx, sage, spherogram)
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/kontrolle.py`
- gate: EQUIVALENT; syntax: OK; edited: no
- German data strings (kept): none
- /media/squeezy paths: none
- own modules NOT in port_20261105: site_config (abschluss, band_rand, dt_bau, eiche_union, pruef_aglet_flaeche present; third-party: snappy, spherogram)
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/offene_restliche.py`
- gate: EQUIVALENT; syntax: OK; edited: no
- German data strings (kept): none
- /media/squeezy paths: none
- own modules NOT in port_20261105: site_config (band_rand, dt_bau, eiche_union, pruef_aglet_flaeche present; third-party: snappy, spherogram)
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/pruef_aglet_flaeche.py`
- gate: EQUIVALENT; syntax: OK; edited: no
- German data strings (kept): control-case names L121-124 `"leeres Wort (eingebetteter Bogen)"`, `"'11'   R1-Schlinge, beide Enden aussen"`, `"'1122' zwei R1-Schlingen"`, `"'123123' = Kleeblatt aufgeschnitten"` (tuple data, printed later)
- /media/squeezy paths: none
- own modules NOT in port_20261105: none
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/rs4_voll.py`
- gate: EQUIVALENT; syntax: OK; edited: no
- German data strings (kept): none
- /media/squeezy paths: none
- own modules NOT in port_20261105: site_config (band_rand, eiche_union, pruef_aglet_flaeche present)
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/rs5.py`
- gate: EQUIVALENT; syntax: OK; edited: no
- German data strings (kept): JSON keys (`aufgaben`, ...)
- /media/squeezy paths: none
- own modules NOT in port_20261105: site_config (abschluss, band_rand, eiche_union, enum_bogendiagramme present)
- private data / credentials: none

## PASS — `deposit_producers/22893039_rs_12n706/zeugen5.py`
- gate: EQUIVALENT; syntax: OK; edited: yes — print messages L109-110 (`ZIEL ... Delta-Stellen ... (nur NEU gegenueber ...)` -> `TARGET ... Delta sites ... (only NEW compared to ...)`), L125 (`Aufgaben` -> `tasks`), L152 (`ZEUGE` -> `WITNESS`), L164-165 (`FERTIG ... Aufgaben · Zeugen ... · tot ... · Stufen` -> `DONE ... tasks · witnesses ... · dead ... · stages`). Note: run logs in the deposit, if any, still show the German wording.
- German data strings (kept): status values `'kein_knoten'`, `'nicht_hyperbolisch'`; state keys `ziel`, `grenze`, `aufgaben`, `zeugen`, `stufen`, `erledigt`, `tot`, `wort`
- /media/squeezy paths: none
- own modules NOT in port_20261105: site_config (abschluss, band_rand, dt_bau, eiche_union, enum_bogendiagramme, pruef_aglet_flaeche present; third-party: snappy, spherogram)
- private data / credentials: none

## Totals
- Files: 37 .py (no .sh in this list). PASS 37 · RESTORED (failed gate) 0.
- Edited: 3 (turm_lauf.py: German duplicate docstring line removed; bau_v3.py: one garbled docstring sentence fixed; zeugen5.py: 6 German print lines translated). 34 files were already English in comments/docstrings/messages and are byte-identical to ORIGINAL (fkp_stellen.py: an edit was reverted, see its entry).
- German left on purpose: data strings (dict keys/values, verdict and status codes, generated receipt Markdown) in about 20 files, listed per file; one German comment inside the byte-compared embedded template of fkp_stellen.py (L72).
- Hard-coded /media/squeezy paths: none in these 37 files (paths go through `site_config.ROOT`). `/home/squeezy/miniforge3/...` interpreter paths appear in docstrings/messages of bau_v3.py, nachweis_obere_schranken.py, neue_obere_schranken2.py, identifiziere_kandidaten.py, konstruktor2.py (kept).
- Own module missing from port_20261105: `site_config` (imported by most files), `pd_to_theta_input` (nachziehen_homfly_symbolisch.py). All other own modules are present.

## Addendum 6 November 2026 (cold read 2)
The folder holds 56 scripts (53 .py, 3 .sh): 36 from the list above, plus 20 that were not on the list because they were already
in English. All 20 were cmp-identical to `port_20261105_ORIG` on 6 November 2026. In one of them,
`21933898_u_lower_bounds/myf_ebenen/nach_lauf2_zweite_schicht.sh`, line 21 had a German log message. It was translated, the
diff against the original is that one line, and `bash -n` passes. Nothing reads that message.
