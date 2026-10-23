#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""theta_numeric_crosscheck.py — stage 8a, numeric variant.

Background: theta_crosscheck.py with Theta.sage (symbolic ℚ(t₁,t₂) inversion)
runs into the timeout at c ≳ 35. Paper workaround (arXiv 2509.18456,
Implementation.tex l.78): rational specialisation t₁=22/7, t₂=21/13. Three
separate ℚ matrix inversions instead of one symbolic inversion over Q(t1,t2).
Computable up to c≈600 within seconds.

Usage:
    python theta_numeric_crosscheck.py <pflicht.json> [--output X.json]
                                       [--t1 22/7] [--t2 21/13] [--timeout 600]

Output: {name, pd_crossings, delta1_at_t1, theta_full_at, theta_normalized,
         t1, t2, duration_s, expo, error, raw_stdout}
"""
from __future__ import annotations
import argparse
import site_config
import ast
import json
import os
import re
import subprocess
import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR))
from pd_to_theta_input import pd_to_cs_phi  # noqa: E402

THETA_DIR = site_config.THETA_DIR
SAGE_BIN = site_config.SAGE_BIN


def _parse_pd(mandatory: dict) -> list:
    pd_list_str = mandatory.get("pd_list_str", "")
    if pd_list_str:
        try:
            pd_list = ast.literal_eval(pd_list_str)
            return [[int(x) for x in tup] for tup in pd_list]
        except Exception:
            pass
    diag = mandatory.get("mandatory_list", {}).get("DIAGRAM", {})
    pd_code = diag.get("pd_code")
    if pd_code:
        crossings = re.findall(r"X\s*\[\s*([\d\s,]+?)\s*\]", str(pd_code))
        if crossings:
            out = []
            for c in crossings:
                parts = [int(x.strip()) for x in c.split(",") if x.strip()]
                if len(parts) == 4:
                    out.append(parts)
            if out:
                return out
    raise ValueError("no PD code found in pflicht JSON")


def _format_cs_phi_for_sage(cs_rows, phi):
    cs_str = "matrix(ZZ, [" + ",".join(
        f"[{int(s)},{int(i)},{int(j)}]" for (s, i, j) in cs_rows) + "])"
    phi_str = "[" + ",".join(str(int(p)) for p in phi) + "]"
    return cs_str, phi_str


SAGE_TEMPLATE = r"""
import sys, time, json
load("{THETA_DIR}/Theta_numeric.sage")
Cs = {CS}
phi = {PHI}
t1 = QQ({T1_NUM})/{T1_DEN}
t2 = QQ({T2_NUM})/{T2_DEN}
t0 = time.time()
try:
    D1, theta, theta_total, expo = Theta_numeric([Cs, phi], t1=t1, t2=t2)
    dt = time.time() - t0
    out = dict(
        ok=True,
        duration_s=float(dt),
        t1=str(t1),
        t2=str(t2),
        expo=str(expo),
        delta1=str(D1),
        theta=str(theta),
        theta_total=str(theta_total),
    )
except Exception as e:
    dt = time.time() - t0
    out = dict(ok=False, error=repr(e), duration_s=float(dt))
print("THETA_NUMERIC_JSON_START")
print(json.dumps(out))
print("THETA_NUMERIC_JSON_END")
"""


def run_numeric_theta(cs_rows, phi, t1_str="22/7", t2_str="21/13",
                      timeout=600):
    cs_str, phi_str = _format_cs_phi_for_sage(cs_rows, phi)

    def parse_frac(s):
        if "/" in s:
            a, b = s.split("/")
            return int(a), int(b)
        return int(s), 1

    t1_num, t1_den = parse_frac(t1_str)
    t2_num, t2_den = parse_frac(t2_str)
    sage_src = SAGE_TEMPLATE.format(
        THETA_DIR=THETA_DIR, CS=cs_str, PHI=phi_str,
        T1_NUM=t1_num, T1_DEN=t1_den, T2_NUM=t2_num, T2_DEN=t2_den,
    )
    p = subprocess.run(
        [SAGE_BIN, "-c", sage_src],
        capture_output=True, text=True, timeout=timeout,
        cwd=THETA_DIR,
    )
    stdout = p.stdout
    stderr = p.stderr
    m = re.search(
        r"THETA_NUMERIC_JSON_START\s*(.*?)\s*THETA_NUMERIC_JSON_END",
        stdout, flags=re.S)
    if not m:
        return {"ok": False, "error": "no JSON marker in sage stdout",
                "stdout_tail": stdout[-2000:], "stderr_tail": stderr[-2000:]}
    return json.loads(m.group(1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mandatory_json", help="Path to the mandatory_invariants[_complete].json")
    ap.add_argument("--output", default=None)
    ap.add_argument("--t1", default="22/7")
    ap.add_argument("--t2", default="21/13")
    ap.add_argument("--timeout", type=int, default=600)
    args = ap.parse_args()

    mandatory = json.loads(Path(args.mandatory_json).read_text())
    name = mandatory.get("name") or Path(args.mandatory_json).stem.split("_")[0]
    pd_list = _parse_pd(mandatory)
    cs_rows, phi = pd_to_cs_phi(pd_list)

    print(f"[theta_numeric] {name}: c={len(pd_list)} PD-crossings")
    print(f"[theta_numeric]   Cs[0:3]={cs_rows[:3]}..., "
          f"phi-length={len(phi)}, sum(phi)={sum(phi)}, "
          f"sum(signs)={sum(c[0] for c in cs_rows)}")
    print(f"[theta_numeric]   t1={args.t1}  t2={args.t2}  timeout={args.timeout}s")

    res = run_numeric_theta(cs_rows, phi, args.t1, args.t2, args.timeout)

    out_json = args.output or args.mandatory_json.replace(".json", "_theta_numeric.json")
    # res may carry its own 't1'/'t2' keys — avoid the clash
    res_clean = {k: v for k, v in res.items() if k not in ("t1", "t2")}
    record = dict(
        name=name,
        source_mandatory_json=str(args.mandatory_json),
        pd_crossings=len(pd_list),
        pd_input=pd_list,
        cs_rows=[list(map(int, r)) for r in cs_rows],
        phi=list(map(int, phi)),
        sum_phi=int(sum(phi)),
        sum_signs=int(sum(c[0] for c in cs_rows)),
        t1=args.t1,
        t2=args.t2,
        t1_from_sage=res.get("t1"),
        t2_from_sage=res.get("t2"),
        **res_clean,
    )
    Path(out_json).write_text(json.dumps(record, indent=2))

    if res.get("ok"):
        print(f"[theta_numeric] OK in {res['duration_s']:.2f}s")
        print(f"[theta_numeric]   Δ(t1={args.t1}) = {res['delta1']}")
        # theta values are often very long; truncate
        th = res.get("theta_total", "")
        if len(th) > 200:
            th = th[:200] + f"... ({len(th)} chars)"
        print(f"[theta_numeric]   theta_total      = {th}")
        print(f"[theta_numeric]   expo             = {res.get('expo')}")
    else:
        print(f"[theta_numeric] ERROR: {res.get('error')}")
        if "stderr_tail" in res:
            print(res["stderr_tail"][-500:])
    print(f"[theta_numeric] → {out_json}")


if __name__ == "__main__":
    main()
