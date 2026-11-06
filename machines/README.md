# Knot machines

These folders hold the generators that built the machine knots of the corpus, plus the scripts that compared the corpus with the
census tables and checked knot identities. Each numbered folder is one machine (`01_knotting_oscillations` …
`259_over_knots`): a construction idea (oscillations, tangles, cables, braids, foldings, …) turned into closed curves, written as
Fourier series (`.fseries`). Which knot a curve is gets decided afterwards, by the identification pipeline (`pipeline/`), never
by the machine itself.

## What is here

- **77 numbered machines.**
  - 72 of them have a `generate_*.py`. Most of these have a `demo` subcommand, which is the smallest run; some also have
    production subcommands (`make`, `plan`, `scan`, …).
  - 62 machines import the shared helpers in `machine_common.py`: the embedding report and polyline → `.fseries`.
  - Only machine 81 uses `machine_topo_gate.py`.
- **`the_knotfamily_machine/`** holds 34 family generators (`make_family*.py`) and `make_connectsum_ring.py`.
  - Its `fseries_io.py` is used throughout: by `machine_common.py`, and directly by 9 machines. Machines import it from this
    sibling folder, so the folder has to stay here.
  - `knot_machine.py` is used only inside the folder, and by `pipeline/r3_polygon_regina_det.py`.
  - `double3_module_extracted.npy` comes from the sculpture Double-3 and carries its own © notice.
- **`corpus_vs_burton/` and `corpus_vs_maguire/`** compare the corpus with Burton's census of prime knots up to 19 crossings and
  with Maguire's data. They need those third-party data sets (`external_tools/burton_census/raw`, `maguire_data/`), which are not
  in this repository.
- **`trainstation_sundown_isometry_20260917/`, `iso_audit_distinct_20260917/` and `KNOWN_DOUBLE_COUNT_20260919/`** are
  identity checks: isometry of the knot complements, and normalisation of the Alexander polynomial for the distinctness key.
- **`_receipts_20261105/`** holds the logs of the test runs described below and the table of sibling references.

## Running a machine

```
cd machines/13_shrinking_into_itself
python generate_selfsimilar_knots.py demo
```

The numbered machines write into their own folder: 72 scripts write to `demo_smoke/`, 17 to `out_fseries/` and one to
`new_knots/`.

Some scripts read `pipeline/site_config.py`, so they need `pipeline/` on `PYTHONPATH`. These are the scripts in 6 machine
folders, plus all `make_family*.py` and `make_connectsum_ring.py`. In the default (portable) profile, `site_config` puts the
knots they generate under `output/`.

Requirements:

- Every machine except 19 needs `numpy`. Machine 19 uses only the Python standard library.
- Machine 258 also needs `scipy`, and machine 247 needs `sympy`.
- The topology gate of machine 81 needs `regina` and `pyknotid`.
- The end-to-end script of machine 191 (`ondemand_e2e_2026-06-29/ondemand_e2e.py`) needs `snappy`.
- The census and identity checks need `snappy` and `spherogram`. The comparisons with Burton's census also need `regina`
  (5 scripts in `corpus_vs_burton/work/`, `work_c19_20260904/` and `work_detjagd_20260907/`). `nitro_isometry.py` also needs `sympy` for its last step
  (determinant and symmetry summary), so run it in an environment that has all three, for example `sage_stable`.

Keep `machines/` together. 53 of the 77 machines import from, or read files of, sibling folders: for example,
`generate_braid_closures.py` comes from machine 07 and `generate_writhe_pumped.py` from machine 27. A machine copied out on its
own will not start. The full table is `_receipts_20261105/sibling_references_20261106.txt`.

`output/bekannt/` holds 222 machine curves of known table knots. Machines 47 and 49 read them as their input stock.

## State on 6 November 2026

Every generator script, 121 in all, was run once in a fresh copy of this repository: with `demo` where it exists, otherwise with
its defaults, each capped at 180 s. The run used the environment of the first test on 5 November for each script, and `base`
for 32, 47 and 80. Any script that failed was run a second time after all the others. The results, one line per run, are in
`_receipts_20261105/smoke_all_post_port_20261106.jsonl`. The first test, `smoke_run_118.json`, was made in a copy that still
lacked parts of this repository and that held older demo outputs; it is kept as a record, under the folder names of its day.

Of the 121 scripts:

- **116 finished without error.** Two of them need care:
  - `114_cabling_the_ghosts` reads the demo output of machine 21, so run machine 21's demo first. In the first pass it came
    before 21 and stopped (`21er-Geist fehlt`, see German strings below); in the second pass it was green.
  - `64_counting_the_sticks` needs about 150–220 s. It hit the 180 s cap on a fully loaded machine, and run alone without
    the cap it finished green in 218 s (`_receipts_20261105/smoke_64_alone_20261106.log`).
- **1 was still computing at the cap:** `259_over_knots/rebuild.py`, as in the first test.
- **3 demos fail:**
  - `81_clasping_the_double` stops at its topology-gate check (`assert ok3`).
  - `47_breeding_for_strangeness` stops at its embedding gate (`assert ok4`).
  - `111_reversing_the_arrow` builds from the demo outputs of machines 47, 84 and 80 and needs at least two of them. It stops
    at the first source, machine 47, which has produced no output.

  For 81 and 47, the original code in the authors' tree, run in today's environments, stops at the same assertion. The logs
  are `_receipts_20261105/orig_run_81_demo_20261105.log` and `orig_run_47_demo_20261105.log`. The cause is therefore not
  this release. Machine 111 was not run with its original code.
- **1 needs arguments:** `191_twisted_torus_hybrid/ondemand_e2e_2026-06-29/braid_to_fseries.py` takes them on the command
  line. Without them it stops with `ValueError`. With them, for example `braid_to_fseries.py 2 3 2 1 64 out.fseries`, it
  works.

The three identity checks were run on the same tree: `isometry.py` (Trainstation = Sundown = Burton's 18ah_5090151),
`three_sculptures.py` (Double-3 = Always a Pleasure = Proxima B = 12a_1019) and `nitro_isometry.py` (Nitro = kf191, kf201, kf202). Each one
found its identities, and each counter-case came out different (`_receipts_20261105/identity_checks_post_port_20261106.log`).
The first `nitro_isometry.py` run stopped after its result at `import sympy`, because that environment lacked sympy. In
`sage_stable` it ran through.

Machine 80's demo used to shuffle with numpy's global, unseeded random generator, and its check (1) failed in 4 of 16 runs
(`_receipts_20261105/rerun_demos_20261106.log`). On 6 November 2026 the shuffle was moved to the machine's own seeded
generator (`generate_hp_folding_knots.py`, line 52). Since then the demo was green in 8 of 8 runs, with byte-identical output in
both environments. Curves made before that date cannot be reproduced bit for bit from their seed.

## German strings

Comments, docstrings and messages are translated into English. Some data values stay German, for example rejection reasons,
verdict labels and dict keys such as `"rejected": "embedding gate (Voll-Auflösung)"` or `urteil = "EINGEBETTET"`. They are
written into result files that the rest of the pipeline and the published corpus read, so translating them would change the
data. The same holds for a few short labels inside f-strings.
Two assertion messages also stay German, because the equivalence gate compares them as code: `21er-Geist fehlt` in machine 114
("the ghost knot of machine 21 is missing": run machine 21's demo first) and `Richtungs-Tripel muss teilerfremd sein` in machine 57
("the direction triple must be coprime").

## Names
Folder and file names are English. On 6 November 2026 these were renamed, and every reference to them inside `machines/` was
changed with them:

| now | in the authors' tree |
|---|---|
| `258_double_3_lattice_module/` (`build_lattice.py`, `embedding.py`) | `258_double_3_gitter_modul/` (`gitter_bauen.py`, `einbettung.py`) |
| `259_over_knots/` (`build.py`, `rebuild.py`) | `259_ueberknoten/` (`bauen.py`, `neubau.py`) |
| `trainstation_sundown_isometry_20260917/` (`isometry.py`, `nitro_isometry.py`, `three_sculptures.py`) | `trainstation_sundown_isometrie_20260917/` (`isometrie.py`, `nitro_isometrie.py`, `drei_skulpturen.py`) |
| `KNOWN_DOUBLE_COUNT_20260919/` (`alex_normal_form_countercase.py`, `alex_sign_sweep.py`) | `KNOWN_DOPPELZAEHLUNG_20260919/` (`alex_normalform_gegenfall.py`, `alex_vorzeichen_sweep.py`) |
| `iso_audit_distinct_20260917/` | `iso_audit_distinkt_20260917/` |
| `32_effects_units_machine.txt`, `47_evolution_machine.txt`, `80_chain_folding.txt` | `32_effektgeraete_maschine.txt`, `47_evolutions_maschine.txt`, `80_ketten_faltung.txt` |

Published data (`corpus_db.json`, `proven_distinctions.json`) and older receipts give paths in the authors' tree, so they use
the right-hand names. A few names were kept on purpose:
- the receipt files in `corpus_vs_burton/` (`beleg_*.json`, `burton_beleg_20260904.json`, `work_detjagd_20260907/`), because
  published data names them;
- output names that scripts create at run time (`kandidaten/`);
- the file names in `deposit_producers/`, which are receipts of the runs as they were made.

## Licence

- Code: GPL-3.0-or-later.
- Curves the machines produce: CC BY-NC 4.0.
- Inputs derived from Werner Alois Stanggassinger's sculptures: © all rights reserved. Downloading them and computing with
  them for study is allowed; see the NOTICE next to each file.
