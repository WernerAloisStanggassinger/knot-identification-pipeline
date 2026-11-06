#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger. With Claude Code (Anthropic).
"""
theta_pointeval_corpus.py — θ@(22/7,21/13) POINT EVALUATION over all machine-
knots. NO cliff: Theta_numeric (three rational matrix inversions, paper "up to
c=600") gets EVERY boss in seconds where the full symbol times out (Disc #246
addendum, verified 2026-06-02: kf372 c32 in 0.05s, poly_c09 c144 in 3.9s).

Step 1 (here, base): harvest (= same PD sources as theta_over_corpus:
stored_simplified or denser simppd for the n20 family) -> (Cs,phi) literals.
Step 2 (sage_stable): ONE Sage session loads Theta_numeric, evaluates all.
Pure diagnosis. corpus_db UNTOUCHED.
"""
import site_config
import sys, json, subprocess
from pathlib import Path

ROOT = Path(site_config.ROOT)
sys.path.insert(0, str(site_config.TOOLS_DIR))
sys.path.insert(0, str(site_config.PIPELINE_DIR))
import theta_over_corpus as toc                                  # noqa
from pd_to_theta_input import pd_to_cs_phi, cs_phi_to_sage_literal  # noqa

LITS = ROOT / "claude_setup/tmp_durable/_pointeval_literals.json"
OUT  = ROOT / "claude_setup/tmp_durable/theta_pointeval_RESULTS.json"
SAGE = Path(site_config.TOOLS_DIR) / "theta_pointeval_run.sage"


def main():
    reps = toc.load_reps()
    print(f"== harvest {len(reps)} machine reps ==", flush=True)
    simppd = json.loads(toc.SIMPPD.read_text())
    jobs = toc.harvest(reps, toc.index_complete(), simppd)
    lits = {}
    for j in jobs:
        stem = j["stem"]
        if "pd" not in j:
            lits[stem] = {"error": j.get("error", "kein PD")}
            continue
        try:
            cs, phi = pd_to_cs_phi([list(t) for t in j["pd"]])
            lits[stem] = {"lit": cs_phi_to_sage_literal(cs, phi),
                          "c": j.get("simp_c"), "src": j.get("source")}
        except Exception as e:
            lits[stem] = {"error": f"pd_to_cs_phi: {type(e).__name__}: {e}"}
    LITS.write_text(json.dumps(lits, indent=2))
    ok_lit = sum(1 for v in lits.values() if "lit" in v)
    print(f"== {ok_lit}/{len(lits)} literals -> sage point-eval ==", flush=True)

    conda = site_config.CONDA_BIN
    p = subprocess.run([conda, "run", "-n", "sage_stable", "--no-capture-output",
                        "sage", str(SAGE)], text=True)
    print(f"sage exit={p.returncode}")
    if OUT.exists():
        res = json.loads(OUT.read_text())
        okv = sum(1 for v in res.values() if v.get("ok"))
        z = sum(1 for v in res.values() if v.get("ok") and v.get("theta_zero"))
        print(f"== RESULT: {okv}/{len(res)} point-eval ok | {z} θ=0 (amphi) ==")


if __name__ == "__main__":
    main()
