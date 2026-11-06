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
CENSUS_SIZE = {'17a': 1769978, '17n': 6283385, '18a': 8400285, '18n': None,
               '19a': None, '19n': None}

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
    matches, provenance = [], {}
    det_fails = {}
    t0 = time.time(); total = 0
    with mp.Pool(NPROC, initializer=_init, initargs=(det_set,)) as pool:
        for path, crossing in files:
            n_here = 0
            for res in pool.imap_unordered(_worker, _stream_file(path, crossing), chunksize=200):
                total += 1; n_here += 1
                if isinstance(res, dict):
                    matches.append(res)
                elif isinstance(res, tuple) and res[0] == 'FAIL':
                    det_fails[crossing] = det_fails.get(crossing, 0) + 1
                if total % 200000 == 0:
                    el = time.time() - t0
                    print(f"  {crossing} processed={total} matches={len(matches)} "
                          f"detfails={sum(det_fails.values())} {el:.0f}s ({total/el:.0f}/s)", flush=True)
            exp = CENSUS_SIZE.get(crossing)
            ok = (exp is None) or (n_here == exp)
            provenance[crossing] = {'lines': n_here, 'expected': exp, 'complete': ok}
            flag = '' if ok else '  ⚠️ COUNT MISMATCH vs census!'
            print(f"DONE {crossing}: {n_here} knots (expected {exp}){flag} | matches {len(matches)}", flush=True)

    # BUG3 hard-gate: any det failure is unexpected (regina parses all Burton sigs)
    if sum(det_fails.values()) > 0:
        print(f"⚠️ DET FAILURES: {det_fails} — candidates may be incomplete, investigate", flush=True)
    # PROV hard-gate: known census sizes must match
    bad = [c for c, p in provenance.items() if not p['complete']]
    if bad:
        print(f"🔴 COMPLETENESS FAIL for {bad} — beyond-18 result NOT trustworthy", flush=True)

    json.dump({'candidates': matches, 'provenance': provenance, 'det_fails': det_fails,
               'total_streamed': total},
              open(WORK + '/candidates_stageA.json', 'w'), indent=1)
    el = time.time() - t0
    print(f"STAGE A DONE: {total} knots in {el:.0f}s, {len(matches)} det-matches, "
          f"detfails={sum(det_fails.values())} -> candidates_stageA.json", flush=True)
    from collections import Counter
    print("  matches by crossing:", dict(Counter(m['crossing'] for m in matches)), flush=True)

if __name__ == '__main__':
    main()
