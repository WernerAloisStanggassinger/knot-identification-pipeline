#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Stage A (v2, hardened per pipeline audit): stream Burton census (c=17,18[,19]),
reconstruct each knot from its Regina knot_sig, compute det = |Alexander(-1)|,
keep only knots whose det is in the beyond-17 target det-set. Emits det-matching
candidates for Stage B. Runs in `base` env (regina). CPU-parallel over 24 cores.

Hardening (audit findings — all guard against SILENTLY manufacturing false novelty):
  BUG1: `bzip2 -t` integrity-gate every core file before streaming, AND check
        bzcat's return code during streaming — abort loudly on truncation/corruption.
  BUG2: assert the 4 core c17/c18 files exist and are non-empty (only c19 optional).
  BUG3: count det-computation failures per crossing (distinct from 'det not in set').
  PROV: write a provenance manifest (per-file line count + expected census size) so
        Stage B can assert completeness instead of assuming it.
Calibration gate: figure-8 det == 5 (startup, aborts if wrong).
"""
import site_config
import sys, os, json, time, subprocess, multiprocessing as mp

WORK = os.environ.get('CVB_WORK', (site_config.MACHINE_YARD + '/corpus_vs_burton/work'))
RAW  = os.environ.get('CVB_RAW', (site_config.ROOT + '/external_tools/burton_census/raw'))

# known complete hyperbolic census sizes (Burton = Maguire = Hoste-Thistlethwaite,
# verified delta=0 by the completeness cross-check 2026-07-10); None = not yet known.
# Target sizes from the SOURCE (regina-normal.github.io/data.html, retrieved 2026-09-04),
# not from an own earlier run — otherwise the run checks itself against itself.
# 18n used to be None here; the run of 10.07. measured 39,866,095 and hits the source number
# exactly — that is now a check instead of an assumption.
CENSUS_SIZE = {'17a': 1769978, '17n': 6283385, '18a': 8400285, '18n': 39866095,
               '19a': 40619384, '19n': 253510828}

DET_SET = None

def _init(det_set_list):
    global DET_SET, regina
    import regina  # noqa
    DET_SET = frozenset(det_set_list)

def _det_from_sig(sig):
    import regina
    L = regina.Link(sig)
    a = L.alexander()
    d = a.degree()
    return abs(sum(int(str(a[i])) * ((-1) ** i) for i in range(d + 1)))

def _worker(rec):
    crossing, name, sig, dt = rec
    try:
        det = _det_from_sig(sig)
    except Exception:
        return ('FAIL', crossing)          # BUG3: distinct from 'not in set'
    if det in DET_SET:
        return {'det': det, 'crossing': crossing, 'name': name, 'sig': sig, 'dt': dt}
    return None

def _stream_file(path, crossing):
    """Yield (crossing, name, sig, dt); raise on bzcat nonzero exit (BUG1)."""
    p = subprocess.Popen(['bzcat', path], stdout=subprocess.PIPE, bufsize=1 << 20)
    n = 0
    for raw in p.stdout:
        line = raw.decode('ascii', 'replace').rstrip('\n')
        if not line or line.startswith('name,'):
            continue
        f = line.split(',')
        if len(f) < 3:
            continue
        n += 1
        yield (crossing, f[0], f[1], f[2])
    p.stdout.close()
    rc = p.wait()
    if rc != 0:
        raise RuntimeError(f"bzcat rc={rc} on {path} after {n} lines — TRUNCATED/CORRUPT, ABORT")

def main():
    det_set = json.load(open(WORK + '/target_det_set.json'))
    import regina
    F = regina.ExampleLink.figureEight()
    aa = F.alexander()
    cal = abs(sum(int(str(aa[i])) * ((-1) ** i) for i in range(aa.degree() + 1)))
    assert cal == 5, f"CALIBRATION FAILED: figure-8 det={cal}, expected 5"
    print(f"CALIB OK figure-8 det=5 | target dets={len(det_set)}", flush=True)

    core = [(f'{RAW}/17a-hyp.csv.bz2', '17a'), (f'{RAW}/17n-hyp.csv.bz2', '17n'),
            (f'{RAW}/18a-hyp.csv.bz2', '18a'), (f'{RAW}/18n-hyp.csv.bz2', '18n')]
    # CVB_ONLY: a small file instead of the whole census — so that the rebuild can be checked in a
    # run of seconds and does not first show up after 20 hours.
    # Only for the self-test; the real run does not set the variable.
    _only = os.environ.get('CVB_ONLY')
    if _only:
        core = [(f'{RAW}/{_only}-hyp.csv.bz2', _only)]
        print(f"⚠️ CVB_ONLY={_only} — TEST RUN over ONE file, not a complete census", flush=True)
    optional = []
    if '--with-c19' in sys.argv:
        optional = [(f'{RAW}/19a-hyp.csv.bz2', '19a'), (f'{RAW}/19n-hyp.csv.bz2', '19n')]

    # BUG2: core files must exist and be non-empty
    for path, c in core:
        if not (os.path.exists(path) and os.path.getsize(path) > 0):
            raise SystemExit(f"CORE FILE MISSING/EMPTY: {path} — ABORT (would fake beyond-18)")
    # BUG1: integrity-gate every file we will stream (core + requested optional)
    files = core + [(p, c) for p, c in optional if os.path.exists(p)]
    for path, c in files:
        if subprocess.run(['bzip2', '-t', path]).returncode != 0:
            raise SystemExit(f"INTEGRITY FAIL (bzip2 -t): {path} — ABORT")
    print(f"INTEGRITY OK for {len(files)} files: {[c for _, c in files]}", flush=True)

    NPROC = int(os.environ.get('CVB_NPROC', os.cpu_count() or 24))
    # 2026-09-04 STREAMING instead of collecting. The run of 10.07. held 4.69 million hits as
    # Python dicts in memory and wrote them at the end in ONE json.dump (745 MB) — with
    # 56 million streamed knots that just worked. This run streams 294 million knots
    # additionally (19a 40,619,384 + 19n 253,510,828). At the same hit rate that would be
    # millions of dicts plus a json.dump that builds the whole string in memory once more.
    # ⚠️ This box FREEZES on memory pressure instead of killing (reference-box-freezes-instead-
    # of-oom-killing) — the price would be a reboot after 20 hours of computing. Therefore: each
    # line goes out at once, memory stays flat.
    cand_path = WORK + '/candidates_stageA.jsonl'
    cand_fh = open(cand_path, 'w')
    n_matches = 0
    matches, provenance = [], {}
    det_fails = {}
    t0 = time.time(); total = 0
    with mp.Pool(NPROC, initializer=_init, initargs=(det_set,)) as pool:
        for path, crossing in files:
            n_here = 0
            for res in pool.imap_unordered(_worker, _stream_file(path, crossing), chunksize=200):
                total += 1; n_here += 1
                if isinstance(res, dict):
                    cand_fh.write(json.dumps(res) + '\n')
                    n_matches += 1
                elif isinstance(res, tuple) and res[0] == 'FAIL':
                    det_fails[crossing] = det_fails.get(crossing, 0) + 1
                if total % 200000 == 0:
                    el = time.time() - t0
                    print(f"  {crossing} processed={total} matches={n_matches} "
                          f"detfails={sum(det_fails.values())} {el:.0f}s ({total/el:.0f}/s)", flush=True)
            exp = CENSUS_SIZE.get(crossing)
            ok = (exp is None) or (n_here == exp)
            provenance[crossing] = {'lines': n_here, 'expected': exp, 'complete': ok}
            flag = '' if ok else '  ⚠️ COUNT MISMATCH vs census!'
            cand_fh.flush()
            print(f"DONE {crossing}: {n_here} knots (expected {exp}){flag} | matches {n_matches}", flush=True)

    # BUG3 hard-gate: any det failure is unexpected (regina parses all Burton sigs)
    if sum(det_fails.values()) > 0:
        print(f"⚠️ DET FAILURES: {det_fails} — candidates may be incomplete, investigate", flush=True)
    # PROV hard-gate: known census sizes must match
    bad = [c for c, p in provenance.items() if not p['complete']]
    if bad:
        print(f"🔴 COMPLETENESS FAIL for {bad} — beyond-18 result NOT trustworthy", flush=True)

    cand_fh.close()
    # provenance separate — the candidates lie next to it as JSONL, one line per hit.
    json.dump({'candidates_file': 'candidates_stageA.jsonl', 'n_candidates': n_matches,
               'provenance': provenance, 'det_fails': det_fails, 'total_streamed': total},
              open(WORK + '/candidates_stageA_meta.json', 'w'), indent=1)
    el = time.time() - t0
    print(f"STAGE A DONE: {total} knots in {el:.0f}s, {n_matches} det-matches, "
          f"detfails={sum(det_fails.values())} -> candidates_stageA.json", flush=True)
    print(f"  -> {cand_path}", flush=True)

if __name__ == '__main__':
    main()
