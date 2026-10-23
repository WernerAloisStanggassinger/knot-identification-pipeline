#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""theta_crosscheck.py  —  Stage 8a: Bar-Natan / van der Veen Theta invariant
                                      as a 2nd independent identification path.

Background (Discovery #133, roadmap step 4 "Dror's QR-code method"):
- Theta is a 2-variable polynomial (t1, t2) that separates
  Conway/Kinoshita-Terasaka pairs on certain knots (same Alexander, Khovanov,
  HOMFLY, hyperbolic volume — Theta DISTINGUISHES them).
- Bar-Natan's tool is used here as a SECOND, fully independent path
  for pipeline identification. NO table lookup; Theta is computed
  reproducibly from the PD code.
- For corpus knots with a Theta result we can later (roadmap
  step 4b/web) match against a small Theta DB.

Workflow (stage 8a):
  1. Load the pflicht JSON, extract the PD code from mandatory_list.DIAGRAMM.pd_code
     (or the top-level pd_list_str).
  2. Convert PD -> (Cs, phi) via pd_to_theta_input.pd_to_cs_phi.
  3. Call sage in a sage_stable subprocess, which loads Theta.sage and computes
     theta = Theta([Cs, phi]) as well as the separately extracted Alexander
     Delta1.
  4. Output: <knot>_theta.json with both polynomials as strings, plus an
     Alexander consistency check (Alexander from the pipeline pflicht vs.
     Theta Delta1).

Usage:
  python theta_crosscheck.py <knot>_pflicht_invarianten_complete.json [--output X.json]

ENV gating (in pipeline_run.sh): run only if RUN_THETA=1.
"""
from __future__ import annotations
import argparse
import site_config
import ast
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR))

from pd_to_theta_input import pd_to_cs_phi, cs_phi_to_sage_literal  # noqa: E402

THETA_DIR = site_config.THETA_DIR


def _parse_pd(mandatory: dict) -> list:
    """Extract the PD list as a list of 4-tuples (int) from a pflicht JSON.

    Searches in this order:
      1. pd_list_str (top level, e.g. '[(1, 5, 2, 4), (3, 1, 4, 6), ...]')
      2. mandatory_list.DIAGRAMM.pd_code as a string 'PD[X[1, 5, 2, 4], X[3,...], ...]'
    """
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
        # Parse "PD[X[1, 5, 2, 4], X[3, 1, 4, 6], X[5, 3, 6, 2]]"
        crossings = re.findall(r"X\s*\[\s*([\d\s,]+?)\s*\]", str(pd_code))
        if crossings:
            out = []
            for c in crossings:
                parts = [int(x.strip()) for x in c.split(",") if x.strip()]
                if len(parts) == 4:
                    out.append(parts)
            if out:
                return out
    raise ValueError("No PD code found in pflicht JSON (neither pd_list_str nor pd_code).")


SAGE_TEMPLATE = r"""
import os, json, sys, signal
THETA_DIR = "__THETA_DIR__"
load(os.path.join(THETA_DIR, "Theta.sage"))

cs_rows = __CS_ROWS__
phi = __PHI__

n = len(cs_rows)
if int(n) == 0:
    # Unknot special case: no PD, no Theta to compute.
    print(json.dumps({
        "theta_full_str": "0",
        "delta1_str":     "1",
        "theta_norm_str": "0",
        "n_crossings":    int(0),
        "ok": True,
        "note": "Unknot (PD empty) - theta=0 by convention",
    }))
    sys.exit(0)

Cs_mat = matrix(cs_rows)

# Sub-timeouts via SIGALRM: for certain PD configurations pd_to_cs_phi can
# produce a Cs vector in an order for which Theta() blows up polynomially in the
# inner double loop (>20min observed on
# kf13_shear+1.0). We guard Theta() itself and theta_norm separately.
def _alrm(*_):
    raise TimeoutError("theta-substep timeout")
signal.signal(signal.SIGALRM, _alrm)

theta_poly = None
theta_full_str = "<not_computed>"
theta_note = ""
try:
    signal.alarm(__THETA_TIMEOUT__)
    theta_poly = Theta([Cs_mat, phi])
    signal.alarm(0)
    theta_full_str = str(theta_poly)
except TimeoutError:
    signal.alarm(0)
    theta_full_str = "<timeout_theta_main>"
    theta_note = "Theta() sub-timeout: pathologische Cs-Reihenfolge"
    print(json.dumps({
        "theta_full_str": theta_full_str,
        "delta1_str":     "<skipped>",
        "theta_norm_str": "<skipped>",
        "theta_note":     theta_note,
        "n_crossings":    int(n),
        "ok": False,
        "error":          "theta_main_timeout",
    }))
    sys.exit(0)

# Extract Alexander Delta1 separately (cleaner than the composite Theta).
R = PolynomialRing(QQ, 't1,t2'); Fr = R.fraction_field()
t1, t2 = R.gens()
A1 = identity_matrix(Fr, 2*n+1)
for c in Cs_mat:
    A1[c[1], c[1]+1] = -t1**c[0]
    A1[c[1], c[2]+1] = t1**c[0] - 1
    A1[c[2], c[2]+1] = -1
Delta1 = t1**(QQ((-sum(phi) - sum(Cs_mat[:,0])[0])/2)) * det(A1)

# theta_norm computation: 2-variable polynomial division in the fraction field;
# can get expensive for large polynomials — own sub-timeout.
theta_norm_str = "<not_computed>"
theta_norm_note = ""
try:
    signal.alarm(60)
    Delta2 = Delta1(t1=t2)
    Delta3 = Delta1(t1=t1*t2)
    theta_norm = theta_poly / (Delta1 * Delta2 * Delta3)
    theta_norm_str = str(theta_norm)
    signal.alarm(0)
except TimeoutError:
    theta_norm_str = "<timeout_60s>"
    theta_norm_note = "theta_norm sub-timeout in 2-var frac-field division"
except Exception as e:
    theta_norm_str = "<error>"
    theta_norm_note = f"theta_norm error: {e}"
finally:
    signal.alarm(0)

out = {
    "theta_full_str":   theta_full_str,
    "delta1_str":       str(Delta1),
    "theta_norm_str":   theta_norm_str,
    "theta_norm_note":  theta_norm_note,
    "theta_note":       theta_note,
    "n_crossings":      int(n),
    "ok": True,
}
print(json.dumps(out))
"""


def run_theta_sage(cs_rows, phi, env_name: str = "sage_stable",
                   timeout_s: int = 600) -> dict:
    """Call Theta.sage in a sage_stable subprocess, fetch the polynomial strings.

    sage_stable has sage 9.x and contains the knot-theory/sage standard packages.
    """
    conda = shutil.which("conda") or site_config.CONDA_BIN
    if not Path(conda).exists():
        return {"ok": False, "error": f"conda not found at {conda}"}

    cs_rows_lit = "[" + ", ".join(
        "[{}, {}, {}]".format(int(s), int(i), int(j)) for s, i, j in cs_rows
    ) + "]"
    phi_lit = "[" + ", ".join(str(int(x)) for x in phi) + "]"
    # Theta() sub-timeout = ~80% of the wrapper timeout (the rest for Delta1 + IO)
    theta_main_timeout = max(30, int(timeout_s * 0.8))
    code = (SAGE_TEMPLATE
            .replace("__THETA_DIR__", THETA_DIR)
            .replace("__CS_ROWS__", cs_rows_lit)
            .replace("__PHI__", phi_lit)
            .replace("__THETA_TIMEOUT__", str(theta_main_timeout)))

    # Write the Sage script to a temp file (Sage has trouble with -c stdin and
    # multi-line scripts via sage -c). We write to a temp .sage:
    import tempfile
    with tempfile.NamedTemporaryFile(
        "w", suffix=".sage", delete=False, prefix="theta_run_"
    ) as fh:
        fh.write(code)
        sage_path = fh.name

    # Process-group isolation: conda run starts sage as a child; on a subprocess
    # timeout we must kill the whole group, otherwise Sage zombies live on
    # at 100% CPU (observed on kf13: >20min CPU hang after the timeout).
    import os as _os
    try:
        proc = subprocess.Popen(
            [conda, "run", "-n", env_name, "sage", sage_path],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            preexec_fn=_os.setsid,
        )
        try:
            stdout, stderr = proc.communicate(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            try:
                _os.killpg(proc.pid, 9)  # SIGKILL the whole process group
            except ProcessLookupError:
                pass
            proc.wait()
            return {"ok": False, "error": f"timeout after {timeout_s}s"}

        class _PR:  # minimal adapter
            pass
        proc_result = _PR()
        proc_result.stdout = stdout
        proc_result.stderr = stderr
        proc_result.returncode = proc.returncode
        proc = proc_result
    finally:
        try:
            os.unlink(sage_path)
            # Sage additionally creates a .sage.py file next to it:
            py_path = sage_path + ".py"
            if os.path.exists(py_path):
                os.unlink(py_path)
        except OSError:
            pass

    stdout = (proc.stdout or "").strip()
    if not stdout:
        return {
            "ok": False,
            "error": "no stdout",
            "stderr_tail": (proc.stderr or "")[-800:],
        }
    payload = None
    for line in reversed(stdout.splitlines()):
        line = line.strip()
        if line.startswith("{") and line.endswith("}"):
            try:
                payload = json.loads(line)
                break
            except json.JSONDecodeError:
                continue
    if payload is None:
        return {
            "ok": False,
            "error": "no JSON in sage stdout",
            "stdout_tail": stdout[-800:],
            "stderr_tail": (proc.stderr or "")[-800:],
        }
    return payload


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mandatory_json", type=Path)
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--sage-env", default="sage_stable")
    ap.add_argument("--timeout", type=int, default=600,
                    help="Timeout for the Theta subprocess in seconds")
    args = ap.parse_args()

    if not args.mandatory_json.exists():
        raise SystemExit(f"Pflicht JSON missing: {args.mandatory_json}")
    mandatory = json.loads(args.mandatory_json.read_text())
    name = mandatory.get("name", args.mandatory_json.stem)

    # 1. Extract PD
    pd_list = _parse_pd(mandatory)
    n_cross = len(pd_list)
    print(f"[theta_crosscheck] {name}: c={n_cross} PD-crossings")

    # 2. PD -> (Cs, phi)
    cs_rows, phi = pd_to_cs_phi(pd_list)
    print(f"[theta_crosscheck]   Cs={cs_rows[:3]}{'...' if len(cs_rows)>3 else ''}, "
          f"phi=[{phi[0]},...,{phi[-1]}] (length {len(phi)})")

    # 3. Theta via Sage
    result = run_theta_sage(cs_rows, phi, env_name=args.sage_env,
                            timeout_s=args.timeout)

    # Pipeline doctrine (Discovery #119): on failure we write JSON
    # WITH not_computable markers instead of exiting. That way the knot stays
    # documented in the Theta corpus (failed but recorded), and the
    # corpus table records it as failed.
    if not result.get("ok"):
        print(f"[theta_crosscheck] ERROR: {result.get('error')}")
        if "stderr_tail" in result:
            print("stderr_tail:", result["stderr_tail"][:300])
        result.setdefault("theta_full_str", "<not_computable>")
        result.setdefault("delta1_str", "<not_computable>")
        result.setdefault("theta_norm_str", "<not_computable>")

    # 4. Alexander consistency with the pipeline pflicht (optional, informative).
    alex_pipeline = mandatory.get("mandatory_list", {}).get("ALGEBRA", {}).get(
        "alexander_polynomial", ""
    )
    alex_theta = result.get("delta1_str", "")

    # 4b. HARD det gate (Disc #533): |Delta1(-1)| MUST == the pflicht quartet det.
    #     The string Alexander check (below) was a DEAD GATE: variable/sign/
    #     mirror conventions -> always False -> defused as "indicative only". As a result
    #     stale/wrong-knot theta slipped through (trainstation det1215!=3969, double-6
    #     det49!=483428169). det is mirror/convention invariant = a usable hard bound.
    mandatory_det = _mandatory_determinant(mandatory)
    theta_det = _theta_determinant(result.get("delta1_str", ""))
    det_gate_ok = (
        isinstance(mandatory_det, int) and isinstance(theta_det, int)
        and mandatory_det == theta_det
    )

    out = {
        "name": name,
        "source_mandatory_json": str(args.mandatory_json),
        "pd_input": pd_list,
        "theta_input": {
            "cs_rows": [list(r) for r in cs_rows],
            "phi": list(phi),
            "n_crossings": n_cross,
        },
        "theta_output": {
            "theta_full":  result["theta_full_str"],
            "delta1":      result["delta1_str"],
            "theta_norm":  result["theta_norm_str"],
            "theta_note":  result.get("theta_note", ""),
            "ok":          bool(result.get("ok")),
            "error":       result.get("error"),
        },
        "alexander_consistency": {
            "pipeline_alexander": alex_pipeline,
            "theta_delta1":       alex_theta,
            "note": ("variable/sign/mirror conventions can "
                     "differ; a hard string match is only indicative."),
            "string_equal_after_normalize": _alex_match(alex_pipeline, alex_theta),
        },
        "determinant_gate": {
            "mandatory_determinant": mandatory_det,
            "theta_determinant": theta_det,
            "match": det_gate_ok,
            "note": ("HARD (Disc #533): |Delta1(-1)| == pflicht quartet det. "
                     "Mismatch = stale/wrong-knot theta -> exit 4, do not trust. "
                     "theta_det=None (empty/<...> delta1) = gate not applicable."),
        },
        "tool_version": "Bar-Natan/van der Veen Theta.sage (external_tools/bar_natan_theta)",
    }

    out_path = args.output or args.mandatory_json.with_name(
        args.mandatory_json.stem.replace("_mandatory_invariants_complete", "")
        + "_theta.json"
    )
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False))
    print(f"[theta_crosscheck] {name}: theta file -> {out_path}")
    if result.get("ok"):
        print(f"[theta_crosscheck]   Delta1 = {result['delta1_str']}")
        # det-gate verdict: the computed theta MUST hit the pflicht det, otherwise
        # it ran on the wrong knot (Disc #533) -> loud + exit != 0.
        if theta_det is not None and mandatory_det is not None and not det_gate_ok:
            print(f"[theta_crosscheck] 🔴 DET-GATE FAIL {name}: theta_det={theta_det} "
                  f"!= pflicht_det={mandatory_det} -> stale/wrong-knot theta, do NOT trust!")
            sys.exit(4)  # exit 4 = det-gate breach (distinguishable from exit 3 = compute fail)
    else:
        sys.exit(3)  # exit 3 for sweep tracking, JSON has been written


def _alex_match(pipeline_alex: str, theta_delta1: str) -> bool:
    """Very coarse consistency check: parse both polynomials into coefficient tuples
    (shift-normalised) via pipeline_crosscheck.poly_canonical, compare."""
    try:
        from pipeline_crosscheck import poly_canonical, jones_canonical_with_mirror
        p = poly_canonical(pipeline_alex)
        k_can, k_mir = jones_canonical_with_mirror(theta_delta1)
        if p is None or k_can is None:
            return False
        if p == k_can or p == k_mir:
            return True
        # Negated coefficient match (Alexander +/- sign)
        neg = tuple(sorted((e, -c) for (e, c) in p))
        return neg == k_can or neg == k_mir
    except Exception:
        return False


def _mandatory_determinant(mandatory: dict):
    """Quartet det from pflicht (konsistenz_quartett preferred, otherwise ALGEBRA)."""
    dv = mandatory.get("consistency_quartet", {}).get("det_value")
    if isinstance(dv, int):
        return dv
    d = mandatory.get("mandatory_list", {}).get("ALGEBRA", {}).get("determinant")
    return d if isinstance(d, int) else None


def _theta_determinant(delta1_str: str):
    """|Delta1(-1)| as int; None if not evaluable (empty / <...> marker).

    A genuine knot has Delta1(-1) = odd, != 0. delta1='1'/'0' (unknot-like,
    trivial) evaluates normally and FAILS against a nontrivial pflicht det."""
    s = str(delta1_str).strip()
    if not s or s.startswith("<"):
        return None
    try:
        import sympy as sp
        t1, t2 = sp.symbols("t1 t2")
        e = sp.sympify(s, locals={"t1": t1, "t2": t2})
        for sym in e.free_symbols:
            e = e.subs(sym, -1)
        return abs(int(sp.nsimplify(e)))
    except Exception:
        return None


if __name__ == "__main__":
    main()
