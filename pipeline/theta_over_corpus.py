#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
"""
theta_over_corpus.py — Dror's machine (Bar-Natan/van-der-Veen theta) over the
WHOLE machine corpus. First time comprehensive (2026-06-02, Werner "wir lassen
ueber alles laufen" (EN: we let it run over everything)). Diagnosis/enrichment ONLY — corpus_db UNTOUCHED.

DOCTRINE ANCHORING:
 - theta = cluster-hash BONUS, NOT proof of distinctness (Disc #246).
 - PD source = corpus-certified stored PD from _complete.json,
   regina-FURTHER-simplified (Disc #246: simplification is the lever, not
   timeout). Bypasses the run_chain min-n frame (Disc #243/244 = mis-resolving).
 - poly_cN/kf FAMILY (11): on-disk PD is n20-tainted (Disc #247) -> source =
   dense validated simppd (r3_dense_simppd.json); finished theta adopted from
   r3_dense_theta_DONE.json.
 - own_knots (sculptures): EXCLUDED (Werner: viewing problems, separate).
 - Full throttle: c<=30 fully parallel, c>30 throttled (OOM path wave29 = the ONE
   documented crash -> only this path reduced, Disc/Hard-Rule 6).
"""
import site_config
import sys, json, os, re, glob, argparse, subprocess, resource
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = Path(site_config.ROOT)
sys.path.insert(0, str(site_config.PIPELINE_DIR))
from pd_to_theta_input import pd_to_cs_phi              # noqa
from theta_crosscheck import run_theta_sage             # noqa

DB     = Path(site_config.PIPELINE_DIR) / "corpus_db.json"
SIMPPD = ROOT / "claude_setup/tmp_durable/r3_dense_simppd.json"
DENSE_THETA = ROOT / "claude_setup/tmp_durable/r3_dense_theta_DONE.json"
OUT    = ROOT / "claude_setup/tmp_durable/theta_over_corpus_RESULTS.json"
# ANTI-LOOP (Werner "wir duerfen nicht schleifen" (EN: we must not loop)): persistent "already
# tried" marker, written BEFORE every theta computation. If the box crashes (OOM
# at c144/149), --resume sees "attempted" and skips -> no reloop.
ATTEMPTED = ROOT / "claude_setup/tmp_durable/theta_over_corpus_ATTEMPTED.json"
import threading
_ATT_LOCK = threading.Lock()

# dense family (n20-tainted on-disk -> use dense source)
DENSE_STEMS = ["poly_c04", "poly_c05", "poly_c09", "poly_s08",
               "kf86", "kf106", "kf134", "kf285", "kf295", "kf372", "kf373"]

RLIM_GB = 26  # per regina-simplify subprocess (32GB box, OOM protection only here)


def log(m): print(m, flush=True)


def has(r, k):
    v = r.get(k); return v not in (None, "", [], {}, "null")


def load_reps():
    db = json.loads(DB.read_text())
    recs = db["knots"]; recs = list(recs.values()) if isinstance(recs, dict) else recs
    groups = {}
    for r in recs:
        did = r.get("distinct_knot_id") or r.get("folder")
        groups.setdefault(did, []).append(r)
    reps = []
    for did, g in groups.items():
        rep = next((r for r in g if r.get("is_distinct_representative")), g[0])
        reps.append(rep)
    # ONLY machine, without existing theta_full
    return [r for r in reps if r.get("category") == "machine" and not has(r, "theta_full")]


def index_complete():
    idx = {}
    for sub in ("erzeugte_knoten", "new_found_knots", "the_knotfamily_machine"):
        for f in glob.glob(str(ROOT / sub / "**/*_complete.json"), recursive=True):
            idx.setdefault(os.path.basename(os.path.dirname(f)), []).append(f)
    return idx


def parse_pd(s):
    """'[(2,12,3,11),...]' or 'PD[X[..],..]' -> list of 4-tuples (int)."""
    tuples = re.findall(r'[\(\[]\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*[\)\]]', str(s))
    return [tuple(int(x) for x in t) for t in tuples]


def _rlimit():
    b = RLIM_GB * 1024**3
    try: resource.setrlimit(resource.RLIMIT_AS, (b, b))
    except Exception: pass


def regina_simplify(pd, timeout=600):
    """regina fromPD -> simplify -> reduced PD (list of 4-tuples). base env."""
    script = (
        "import regina, re, sys, json\n"
        f"pd={[list(t) for t in pd]}\n"
        "L=regina.Link.fromPD(pd)\n"
        "n0=int(L.size())\n"
        "L.simplify()\n"
        "s=str(L.pd())\n"
        "ts=re.findall(r'X\\[\\s*(\\d+)\\s*,\\s*(\\d+)\\s*,\\s*(\\d+)\\s*,\\s*(\\d+)\\s*\\]',s)\n"
        "print('REG'+json.dumps({'raw':n0,'simp':int(L.size()),"
        "'pd':[[int(a),int(b),int(c),int(d)] for a,b,c,d in ts]}))\n"
    )
    try:
        p = subprocess.run(
            ["conda", "run", "-n", "base", "--no-capture-output", "python", "-c", script],
            capture_output=True, text=True, timeout=timeout, preexec_fn=_rlimit)
    except subprocess.TimeoutExpired:
        return {"error": "regina-simplify-timeout"}
    for line in p.stdout.splitlines():
        if line.startswith("REG"):
            return json.loads(line[3:])
    return {"error": f"regina-fail: {(p.stderr or '')[-160:]}"}


def harvest(reps, comp_idx, simppd):
    """Per rep: (stem, source, pd, raw_c, simp_c) or preliminary result."""
    jobs = []
    for r in reps:
        folder = str(r.get("folder"))
        base = os.path.basename(folder)
        stem = base
        dense_key = next((d for d in DENSE_STEMS
                          if base == d or base.startswith(d + "_") or ("__" + d) in base
                          or base.split("__")[-1].startswith(d + "_")
                          or base.split("__")[-1] == d), None)
        if dense_key:
            # dense source (avoid n20-tainted on-disk)
            sk = next((k for k in simppd if k == dense_key or k.startswith(dense_key + "_")), None)
            if sk and "pd" in simppd[sk]:
                pd = [tuple(t) for t in simppd[sk]["pd"]]
                jobs.append({"stem": stem, "dense_key": dense_key, "source": "dense_simppd",
                             "pd": pd, "raw_c": len(pd), "simp_c": len(pd)})
                continue
            jobs.append({"stem": stem, "dense_key": dense_key, "source": "dense_missing",
                         "error": "no dense simppd"})
            continue
        # normal knot: stored PD -> regina-simplify
        cands = comp_idx.get(base) or comp_idx.get(folder)
        if not cands:
            jobs.append({"stem": stem, "source": "none", "error": "no _complete.json"})
            continue
        d = json.loads(Path(cands[0]).read_text())
        pd_raw = parse_pd(d.get("pd_list_str") or d.get("pd") or "")
        if not pd_raw:
            jobs.append({"stem": stem, "source": "stored", "error": "no PD in _complete.json"})
            continue
        rs = regina_simplify(pd_raw)
        if "error" in rs:
            jobs.append({"stem": stem, "source": "stored", "error": rs["error"]})
            continue
        jobs.append({"stem": stem, "source": "stored_simplified", "pd": rs["pd"],
                     "raw_c": rs["raw"], "simp_c": rs["simp"]})
    return jobs


def _mark_attempted(stem):
    """Persistent 'already tried' BEFORE the computation -> crash-safe against reloop."""
    with _ATT_LOCK:
        att = json.loads(ATTEMPTED.read_text()) if ATTEMPTED.exists() else []
        if stem not in att:
            att.append(stem); ATTEMPTED.write_text(json.dumps(att, indent=2))


def theta_one(job, timeout):
    _mark_attempted(job["stem"])   # write the marker first, then compute (OOM-crash-safe)
    try:
        cs, phi = pd_to_cs_phi([list(t) for t in job["pd"]])
    except Exception as e:
        return {"error": f"pd_to_cs_phi: {type(e).__name__}: {e}"}
    res = run_theta_sage(cs, phi, env_name="sage_stable", timeout_s=timeout)
    return {"ok": bool(res.get("ok")), "theta_full_str": res.get("theta_full_str"),
            "theta_delta1": res.get("theta_delta1"), "error": res.get("error")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--timeout", type=int, default=3600, help="per-knot theta timeout")
    ap.add_argument("--cliff", type=int, default=30, help="simp_c <= cliff = Tier-1 (parallel)")
    ap.add_argument("--workers-low", type=int, default=12)
    ap.add_argument("--workers-high", type=int, default=3)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--retry", action="store_true",
                    help="also tackle tried-without-result (crash/timeout) AGAIN "
                         "(default: NOT -> we must not loop)")
    ap.add_argument("--harvest-only", action="store_true")
    args = ap.parse_args()

    out = json.loads(OUT.read_text()) if (args.resume and OUT.exists()) else {}
    # Attempted set: everything that ever went into the computation (even if the box
    # died in the process and NO out entry arose). --resume skips these -> no reloop.
    attempted = set(json.loads(ATTEMPTED.read_text())) if (args.resume and ATTEMPTED.exists()) else set()
    if not args.resume:
        try: ATTEMPTED.unlink()
        except FileNotFoundError: pass
    simppd = json.loads(SIMPPD.read_text())
    dense_theta = json.loads(DENSE_THETA.read_text()) if DENSE_THETA.exists() else {}

    reps = load_reps()
    log(f"== {len(reps)} machine distinct reps without theta_full ==")
    jobs = harvest(reps, index_complete(), simppd)

    # adopt finished dense theta
    runnable = []
    for j in jobs:
        stem = j["stem"]
        if stem in out and out[stem].get("ok"):   # resume: already ok
            continue
        # ANTI-LOOP: already tried (also crash/timeout without out entry) -> not again,
        # except --retry. Prevents reloop on the OOM monsters after a box crash.
        if not args.retry and stem in attempted:
            if stem not in out:
                out[stem] = {"ok": False, "source": j.get("source"),
                             "error": "attempted_no_result (Crash/Timeout, resume-skip)",
                             "simp_c": j.get("simp_c")}
            log(f"  SKIP(attempted) {stem}")
            continue
        if j.get("source", "").startswith("dense") and j.get("dense_key"):
            dk = j["dense_key"]
            dkey = next((k for k in dense_theta if k == dk or k.startswith(dk + "_")), None)
            if dkey and dense_theta[dkey].get("ok"):
                dt = dense_theta[dkey]
                out[stem] = {"ok": True, "source": "dense_theta_DONE", "simp_c": dt.get("crossings"),
                             "theta_full_str": dt.get("theta_full_str"),
                             "theta_delta1": dt.get("theta_delta1")}
                log(f"  adopt dense theta: {stem} (c={dt.get('crossings')})")
                continue
        if "error" in j and "pd" not in j:
            out[stem] = {"ok": False, "source": j.get("source"), "error": j["error"]}
            log(f"  SKIP {stem}: {j['error']}")
            continue
        runnable.append(j)
    OUT.write_text(json.dumps(out, indent=2))

    if args.harvest_only:
        log("harvest-only -> stop")
        for j in sorted(runnable, key=lambda x: x.get("simp_c", 999)):
            log(f"  {j['stem']:42s} simp_c={j.get('simp_c')} src={j['source']}")
        return

    low  = sorted([j for j in runnable if j.get("simp_c", 999) <= args.cliff], key=lambda x: x["simp_c"])
    high = sorted([j for j in runnable if j.get("simp_c", 999) >  args.cliff], key=lambda x: x["simp_c"])
    log(f"== Tier-1 (simp_c<={args.cliff}): {len(low)} | Tier-2 (>{args.cliff}): {len(high)} ==")

    def run_tier(tier, workers, label):
        with ThreadPoolExecutor(max_workers=workers) as ex:
            fut = {ex.submit(theta_one, j, args.timeout): j for j in tier}
            for f in as_completed(fut):
                j = fut[f]; r = f.result()
                r["source"] = j["source"]; r["simp_c"] = j.get("simp_c"); r["raw_c"] = j.get("raw_c")
                out[j["stem"]] = r
                OUT.write_text(json.dumps(out, indent=2))
                tf = (r.get("theta_full_str") or "")[:32]
                log(f"  [{label}] {j['stem']:42s} c={j.get('simp_c')} ok={r.get('ok')} "
                    f"{('θ='+tf) if r.get('ok') else ('err='+str(r.get('error'))[:30])}")

    run_tier(low, args.workers_low, "T1")
    run_tier(high, args.workers_high, "T2")

    okn = sum(1 for v in out.values() if v.get("ok"))
    log(f"\n== DONE: {okn}/{len(out)} theta_ok -> {OUT} ==")


if __name__ == "__main__":
    main()
