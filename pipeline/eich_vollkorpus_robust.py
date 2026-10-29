#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Knot-Identification Pipeline — Werner Alois Stanggassinger & Claude Code (Anthropic).
# ============================================================================
# eich_vollkorpus_robust.py — full-corpus re-calibration on the robust
#   route-C SONO form (sixd_robust, Disc #263), GPU at full steam.
#
# WHY: corpus_db still sits on the old add.at form (order-dependent,
#   not well-defined). Only after this calibration do ALL corpus numbers rest
#   (headline, NN, trainstation≡sundown, S-cluster distinctness) on the
#   well-defined, backend-portable form.
#
# MECHANICS (batch_verify_gpu pattern, source = corpus_db instead of 5 verify knots):
#   - jobs = every record (knots + validation_knots), fseries via pflicht_path.
#   - spawn executor (CUDA cannot fork); every worker: sixd_robust.install('gpu')
#     -> GPU overlap SET + det_accumulate (canonical order) + numpy
#     geometry. OUT -> out2_robust/<stem>_sixd2.json.
#   - Full throttle: worker default fat; GPU time-sliced over the jobs.
#
# ----------------------------------------------------------------------------
# 2026-08-25 — SURVIVING WORKER DEATH (Disc #704, Werner: "you have full authority for this")
# ----------------------------------------------------------------------------
# BEFORE, this ran via multiprocessing.Pool.imap_unordered. If a
# worker dies in the middle of a task, Pool does NOT reissue the task and reports
# NOTHING — the iterator waits forever. Exactly that happened on 24.08.: the run
# stood 20 h 43 min at [407/413] and looked like work (ps showed 45 %
# CPU — that is the lifetime average, not the actual value; measured were
# 0 bytes IO and 0.2 % CPU). 413 issued, 407 returned, 6 silently lost.
#
# On this box worker death is no exception: journalctl -u earlyoom
# counted 651 SIGTERMs in the run window (earlyoom prefers exactly python).
#
# COUNTER-CASE, run in both directions (Python 3.11.15):
#   mp.Pool(1).imap_unordered, worker kills itself at task 2  -> HANGS, rc=124
#   concurrent.futures.ProcessPoolExecutor, the same             -> BrokenProcessPool
#
# HENCE HERE: ProcessPoolExecutor + self-healing. A dead worker tears the
# whole executor down with it (all open futures get BrokenProcessPool) — the
# run therefore does NOT abort, but rebuilds the executor and reissues the
# unfinished tasks. Loud, counted, with an upper limit per task.
# Running on blindly would be the old error in green.
#
# SECOND LESSON from #704, built in here: EVERY progress line carries a
# TIME OF DAY. A progress log without a timestamp cannot show standstill —
# the standstill only became visible via the mtime of the output files.
#
# Call (gnottheorie_env python):
#   python eich_vollkorpus_robust.py [--workers 8] [--min-n 0] [--only s1 s2]
#                                    [--resume] [--max-versuche 3]
# ============================================================================
import argparse, concurrent.futures as cf, json, multiprocessing as mp
import site_config
import os, sys, time
from concurrent.futures.process import BrokenProcessPool

SCANNER_DIR = site_config.SIXD_DIR
CORPUS_DB = "@PIPELINE_DIR@/corpus_db.json".replace('@PIPELINE_DIR@', site_config.PIPELINE_DIR)
OUT = os.path.join(SCANNER_DIR, "out2_robust")
CURVES = os.path.join(OUT, "curves")


def _worker(job):
    stem, fpath, min_n = job
    sys.path.insert(0, SCANNER_DIR)
    import canonical_flow as cf_mod
    import sixd_robust
    cf_mod.OUT_DIR = OUT
    cf_mod.CURVE_DIR = CURVES
    # install patches cf.overlap_pairs (GPU, min_n) + cf.wall_block (robust)
    import sixd_gpu
    cf_mod._orig_overlap = cf_mod.overlap_pairs
    cf_mod.overlap_pairs = sixd_gpu.make_dispatcher(cf_mod._orig_overlap, min_n=min_n)
    cf_mod.wall_block = sixd_robust.wall_block_robust
    gpu = sixd_gpu.available()
    t0 = time.time()
    try:
        line = cf_mod.run_record((stem, fpath))
    except Exception as e:
        import traceback
        line = f"FAIL {stem}: {e}\n{traceback.format_exc()[:400]}"
    return stem, time.time() - t0, gpu, line


def collect_jobs(min_n, only=None):
    db = json.load(open(CORPUS_DB))
    pool = list(db["knots"]) + list(db.get("validation_knots", []))
    jobs, missing = [], []
    for rec in pool:
        stem = rec["stem"]
        if only and stem not in only:
            continue
        folder = os.path.dirname(rec["pflicht_path"])
        fpath = os.path.join(folder, stem + ".fseries")
        if not os.path.isfile(fpath):
            cands = ([f for f in os.listdir(folder)
                      if f.endswith(".fseries") and not f.startswith("_sculpture_")]
                     if os.path.isdir(folder) else [])
            if len(cands) == 1:
                fpath = os.path.join(folder, cands[0])
            else:
                missing.append(stem); continue
        jobs.append((stem, fpath, min_n))
    return jobs, missing


def _uhr():
    return time.strftime("%m-%d %H:%M:%S")


def run_jobs(jobs, workers, max_versuche=3, fn=None):
    """Work off all jobs and survive worker death.

    Returns: (n_ok, n_fail, n_aufgegeben, neustarts).

    A dead worker tears the whole executor down: ALL open futures
    get BrokenProcessPool, even those of tasks that have not even
    run. Hence two rules, both worked out on the counter-case:

    (1) An attempt counts against a task only if it ran ALONE.
        Otherwise a poisonous record punishes its innocent neighbours — in the
        first draft exactly such a healthy task was given up that never
        started (counter-case C, first red, then green).
    (2) If a round makes NO progress, the next one goes into isolation:
        one task per executor. Thereby the poisonous one is found precisely,
        and the healthy ones run through.
    """
    fn = fn or _worker
    ctx = mp.get_context("spawn")
    offen = list(jobs)
    versuche = {j[0]: 0 for j in jobs}
    n_ges = len(jobs)
    n_ok = n_fail = n_auf = neustarts = 0
    fertig = 0
    isoliert = False

    while offen:
        runde, offen = offen, []
        portion = 1 if isoliert else len(runde)
        if isoliert:
            print(f"[{_uhr()}] 🔍 ISOLATION — {len(runde)} task(s) one by one, "
                  f"to find the poisonous record", flush=True)
        fortschritt = 0
        i = 0
        while i < len(runde):
            teil = runde[i:i + portion]
            i += portion
            ex = cf.ProcessPoolExecutor(max_workers=workers, mp_context=ctx)
            futs = {ex.submit(fn, j): j for j in teil}
            gebrochen = []
            try:
                for f in cf.as_completed(futs):
                    job = futs[f]
                    try:
                        stem, dt, gpu, line = f.result()
                    except BrokenProcessPool:
                        # Exactly the spot where mp.Pool silently stood still.
                        gebrochen.append(job)
                        continue
                    except Exception as e:
                        n_fail += 1; fertig += 1; fortschritt += 1
                        print(f"[{_uhr()}] [{fertig}/{n_ges}] FAIL {job[0]}: "
                              f"{type(e).__name__}: {e}", flush=True)
                        continue
                    bad = line.startswith("FAIL")
                    n_fail += bad; n_ok += (not bad); fertig += 1; fortschritt += 1
                    print(f"[{_uhr()}] [{fertig}/{n_ges}] {dt:6.1f}s gpu={gpu} "
                          f"{line}", flush=True)
            finally:
                ex.shutdown(wait=False, cancel_futures=True)

            if gebrochen:
                neustarts += 1
                print(f"[{_uhr()}] ⚠️  WORKER DIED — new executor, "
                      f"{len(gebrochen)} task(s) back in the queue "
                      f"(restart {neustarts}). Cause usually earlyoom: "
                      f"journalctl -u earlyoom --since '-10min'", flush=True)
                for job in gebrochen:
                    if len(teil) == 1:          # rule (1)
                        versuche[job[0]] += 1
                    if versuche[job[0]] >= max_versuche:
                        n_auf += 1
                        print(f"[{_uhr()}] 🔴 GIVEN UP after {max_versuche} "
                              f"solo runs (the worker died every time): "
                              f"{job[0]}", flush=True)
                    else:
                        offen.append(job)

        isoliert = (fortschritt == 0 and len(runde) > 1)   # rule (2)
    return n_ok, n_fail, n_auf, neustarts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--min-n", type=int, default=0)   # 0 = always GPU
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--resume", action="store_true",
                    help="skip out2_robust/<stem>_sixd2.json files that already exist")
    ap.add_argument("--max-versuche", type=int, default=3,
                    help="how often a task is reissued after worker death")
    args = ap.parse_args()
    os.makedirs(CURVES, exist_ok=True)
    jobs, missing = collect_jobs(args.min_n, set(args.only) if args.only else None)
    if args.resume:
        done = {f[:-len("_sixd2.json")] for f in os.listdir(OUT)
                if f.endswith("_sixd2.json")}
        jobs = [j for j in jobs if j[0] not in done]
    print(f"[{_uhr()}] [EICH-ROBUST] {len(jobs)} Records, {args.workers} Worker "
          f"(spawn/PPE), GPU_MIN_N={args.min_n}, missing={len(missing)} -> {OUT}",
          flush=True)
    for m in missing:
        print(f"  MISSING {m}")
    t0 = time.time()
    n_ok, n_fail, n_auf, neustarts = run_jobs(jobs, args.workers, args.max_versuche)
    wall = time.time() - t0
    print(f"\n[{_uhr()}] [EICH-ROBUST] DONE {n_ok} ok, {n_fail} fail, "
          f"{n_auf} given up, {len(missing)} missing, {neustarts} worker restarts  "
          f"WALL {wall/60:.1f} min", flush=True)
    # Bookkeeping that the old run did not have: what went in must come out.
    if n_ok + n_fail + n_auf != len(jobs):
        print(f"🔴 BALANCE DOES NOT ADD UP: {len(jobs)} in, "
              f"{n_ok + n_fail + n_auf} out — check!", flush=True)
        sys.exit(2)
    sys.exit(1 if n_auf else 0)


if __name__ == "__main__":
    main()
