#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""check_census.py — DFJ Conjecture 1.7 (dfj_core.py) on every hyperbolic knot of one Burton census file.

For each knot: genus and fiberedness from knot Floer homology (SnapPy's knot_floer_homology, Szabó's program) on the census diagram,
the hyperbolic torsion polynomial T_K from SnapPy (500 bits, escalating to 1000 and 2000 bits when a class is uncertain), the
comparison deg T_K = 4g − 2 and (T_K monic ⟺ fibered); where T_K is not certainly non-real, SnapPy's numerical symmetry group says
whether the knot is amphichiral.

  python -I check_census.py <census.csv.bz2> out.jsonl [WORKER NWORKERS] [--retry-errors]

One JSON line per knot, written the moment it is computed; names already in out.jsonl are skipped (with --retry-errors, rows with
status ERROR are computed again and the new line is appended; the last line of a name counts). WORKER/NWORKERS: this process takes
the knots with index ≡ WORKER mod NWORKERS (run one process per worker, each with its own out file). Names whose census kind is not
'h' (torus, satellite) are written as SKIP. Run with SageMath's Python. Each process limits its address space to MEM_GB (default 6) GB.
"""
import bz2, csv, json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dfj_core as D

DT = lambda s: [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    retry = '--retry-errors' in sys.argv
    census, out_path = args[0], args[1]
    w, nw = (int(args[2]), int(args[3])) if len(args) >= 4 else (0, 1)
    D.set_memlimit(int(os.environ.get('MEM_GB', '6')))
    import snappy
    last = {}
    if os.path.exists(out_path):
        for ln in open(out_path):
            try:
                r = json.loads(ln); last[r['name']] = r.get('status')
            except Exception:  # noqa: BLE001  (truncated last line after a kill)
                pass
    done = {n for n, st in last.items() if not (retry and st == 'ERROR')}
    n_new = 0
    with open(out_path, 'a') as fo:
        for i, row in enumerate(csv.DictReader(bz2.open(census, 'rt'))):
            if i % nw != w: continue
            name = row['name']
            if name in done: continue
            if name.split('_')[0][-1] != 'h':
                fo.write(json.dumps(dict(name=name, status='SKIP_nonhyperbolic_by_census_name')) + '\n'); fo.flush(); continue
            t0 = time.time()
            try:
                L = snappy.Link('DT:' + str([tuple(DT(row['dt_code']))]))
                h = L.knot_floer_homology()
                M = L.exterior()
                r = D.check(M, h['seifert_genus'], h['fibered'])
                if r['real'] != 'nonzero':
                    try:
                        r['amphichiral'] = bool(M.symmetry_group().is_amphicheiral())
                    except Exception as e:  # noqa: BLE001
                        r['amphichiral'] = 'ERR ' + repr(e)[:100]
            except Exception as e:  # noqa: BLE001
                r = dict(status='ERROR', err=repr(e)[:300], sec=round(time.time() - t0, 3))
            r['name'] = name
            fo.write(json.dumps(r) + '\n'); fo.flush(); n_new += 1
    print(f'worker {w}/{nw}: {n_new} knots computed', flush=True)


if __name__ == '__main__':
    main()
