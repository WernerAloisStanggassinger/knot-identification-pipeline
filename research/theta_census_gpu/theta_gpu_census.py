# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""θ at (22/7, 21/13) mod p on the GPU for every knot of Burton census files, K and its mirror image.
CPU pool: DT code -> SnapPy link -> PD -> pd_to_cs_phi (same route as the θ deposit's theta_census.py), for L and L.mirror().
GPU: theta_modp.gpu_batch (calibrated 400/400 against Sage over Q, calibrate_theta_modp.py), bundled per crossing number.
Exactness: θ = 0 over Q ⟹ θ_raw ≡ 0 mod p (p divides no denominator, else the kernel returns None = flagged). So every θ-zero is
among the flagged knots; θ_raw ≢ 0 mod p ⟹ θ ≠ 0 (certificate). θ_raw(mK) + θ_raw(K) ≢ 0 mod p ⟹ the conjectured sign rule θ(mK) = −θ(K) (Conjecture 21 of arXiv:2509.18456v4) fails at the point.
Flagged (written to the candidates file, then recomputed exactly on the CPU): zero mod p, mirror-rule miss mod p, None, conversion error.
Usage: python -I theta_gpu_census.py <workers> <out_prefix> <census.csv.bz2> ...   (env THETA_GPU_MEM_GB, default 6)"""
import sys, os, bz2, csv, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))     # theta_modp.py and pd_to_theta_input.py live here
os.environ.setdefault('THETA_GPU_MEM_GB', '6')
from multiprocessing import Pool
def dt(s): return [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]
def init():
    global snappy, pd_to_cs_phi
    import snappy
    from pd_to_theta_input import pd_to_cs_phi
def conv(row):
    name, code = row
    try:
        L = snappy.Link('DT:' + str([tuple(dt(code))]))
        out = []
        for X in (L, L.mirror()):
            pd = [tuple(x + 1 for x in c) for c in X.PD_code()]
            cs, phi = pd_to_cs_phi(pd); out.append(([list(c) for c in cs], list(phi)))
        return (name, out)
    except Exception as e:
        return (name, 'ERR:' + repr(e)[:100])
def rows(files):
    for f in files:
        with bz2.open(f, 'rt') as fh:
            r = csv.reader(fh); next(r)
            for x in r: yield (x[0], x[2])
if __name__ == '__main__':
    from theta_modp import gpu_batch, frac_mod, P_DEFAULT as p
    workers, pre, files = int(sys.argv[1]), sys.argv[2], sys.argv[3:]
    t1, t2 = frac_mod(22, 7, p), frac_mod(21, 13, p)
    stats = dict(knots=0, zero=0, mirror_miss=0, none=0, conv_err=0); t0 = time.time()
    cand = open(pre + '_candidates.tsv', 'w'); buf = {}
    def flush(n):
        items = buf.pop(n)
        rK = gpu_batch([it[1][0] for it in items], t1, t2, p); rM = gpu_batch([it[1][1] for it in items], t1, t2, p)
        for (name, _), a, b in zip(items, rK, rM):
            stats['knots'] += 1; ta, tb = a['theta_roh'], b['theta_roh']; why = []
            if ta is None or tb is None: why.append('none'); stats['none'] += 1
            else:
                if ta == 0: why.append('zero'); stats['zero'] += 1
                if (ta + tb) % p: why.append('mirror_miss'); stats['mirror_miss'] += 1
            if why: cand.write(f'{name}\t{n}\t{",".join(why)}\t{ta}\t{tb}\n')
    with Pool(workers, initializer=init) as pool:
        for name, res in pool.imap_unordered(conv, rows(files), chunksize=500):
            if isinstance(res, str):
                stats['conv_err'] += 1; cand.write(f'{name}\t-1\tconv_err\t{res}\t\n'); continue
            n = len(res[0][0]); buf.setdefault(n, []).append((name, res))
            if len(buf[n]) >= 20000:
                flush(n); cand.flush()
                print(json.dumps(dict(stats, sec=round(time.time() - t0))), flush=True)
        for n in list(buf): flush(n)
    cand.close(); stats['sec'] = round(time.time() - t0); stats['p'] = p
    json.dump(stats, open(pre + '_stats.json', 'w')); print('DONE', json.dumps(stats), flush=True)
