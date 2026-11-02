# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""Volume census of one of Burton's census files: SnapPy volume of every knot, from the census DT code.

usage: python -I volume_census.py <census.csv.bz2> <out.tsv> <workers>
- one reader decompresses the file once; workers get chunks of CHUNK rows (multiprocessing imap)
- resumable: chunks already present in <out.tsv> (lines '#chunk <i> <n>') are skipped
- per knot: up to 8 attempts (randomize) until 'all tetrahedra positively oriented'; the flag column says
  how it ended (pos / other solution type / ERR); 5 s alarm per attempt
- output row: name, volume (repr), flag, attempts
Calibration in the main process before the run: 4_1 and 5_2 volumes."""
import sys, bz2, csv, os, signal, time, json
import multiprocessing as mp

CHUNK = 2000
POS = 'all tetrahedra positively oriented'
dt = lambda s: [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]


class Late(Exception):
    pass


def _alarm(*a):
    raise Late()


def one(code):
    import snappy
    M = snappy.Link('DT:' + str([tuple(dt(code))])).exterior()
    last = None
    for att in range(1, 9):
        signal.alarm(5)
        try:
            st = M.solution_type()
            v = M.volume()
        except Late:
            st, v = 'TIMEOUT', None
        finally:
            signal.alarm(0)
        last = (v, st, att)
        if st == POS:
            return v, 'pos', att
        M.randomize()
    v, st, att = last
    return v, st.replace(' ', '_'), att


def work(job):
    i, rows = job
    signal.signal(signal.SIGALRM, _alarm)
    out = []
    for name, code in rows:
        try:
            v, flag, att = one(code)
        except Exception as e:
            v, flag, att = None, 'ERR_' + type(e).__name__, 0
        out.append(f"{name}\t{v!r}\t{flag}\t{att}")
    return i, out


def chunks(path, done):
    rows, i = [], 0
    for r in csv.DictReader(bz2.open(path, 'rt')):
        rows.append((r['name'], r['dt_code']))
        if len(rows) == CHUNK:
            if i not in done:
                yield i, rows
            rows, i = [], i + 1
    if rows and i not in done:
        yield i, rows


def main():
    path, outp, nw = sys.argv[1], sys.argv[2], int(sys.argv[3])
    import snappy
    for nm, ref in (('4_1', 2.029883212819307), ('5_2', 2.828122088330783)):
        v = snappy.Manifold(nm).volume()
        assert abs(v - ref) < 1e-9, (nm, v)
    print('CAL ok 4_1 5_2', flush=True)
    done = set()
    if os.path.exists(outp):
        for line in open(outp):
            if line.startswith('#chunk '):
                done.add(int(line.split()[1]))
    print(json.dumps({'resume_chunks_done': len(done)}), flush=True)
    t0, n = time.time(), 0
    with open(outp, 'a') as f, mp.Pool(nw) as pool:
        for i, out in pool.imap_unordered(work, chunks(path, done), chunksize=1):
            f.write('\n'.join(out) + '\n' + f'#chunk {i} {len(out)}\n')
            f.flush()
            n += len(out)
            if (n // CHUNK) % 50 == 0:
                print(json.dumps({'knots_this_run': n, 'sec': round(time.time() - t0), 'per_knot_ms_wall': round(1000 * (time.time() - t0) / n, 3)}), flush=True)
    print('DONE', json.dumps({'knots_this_run': n, 'sec': round(time.time() - t0)}), flush=True)


if __name__ == '__main__':
    main()
