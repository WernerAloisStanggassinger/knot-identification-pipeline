#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""SIEVE c — machine 115 with Conant's rotation parameter c, order 8 (k=4, a=3).

TARGET that the c=1 sieve structurally cannot reach: knots with l = +-3 (mod 8). Only there
is Conant's constant congruence E(0) = (l^2-1)/8 (mod 2) SHARP (E(0) must be odd).
All 84 previous a=3 knots have l = 1 (mod 8) — not by chance: Conant Lemma 4.2 says
c*l = 1 (mod 2^a), and the old machine is the case c=1. For a<=3 that is l = c
(every element of (Z/8)* is its own inverse); FROM a=4 l = c^{-1} holds, not l = c.

  c = 3  ->  l = 3 (mod 8)  ->  E(0) ODD   (sharp)
  c = 5  ->  l = 5 (mod 8)  ->  E(0) ODD   (sharp)
  c = 7  ->  l = 7 (mod 8)  ->  E(0) even     (control: must behave like c=1)

The rotation parameter is called `crot` in the record, NOT `c` — `c` there is the crossing number.
Harness identical to breit3.py (worker recycling, WNOHANG, MemFloor).
"""
import json, os, sys, time, glob, collections
import numpy as np
import importlib.util

HIER = os.path.dirname(os.path.abspath(__file__))
spb = importlib.util.spec_from_file_location('b', os.path.join(HIER, 'breit.py'))
B = importlib.util.module_from_spec(spb)
sys.argv = [sys.argv[0]]
spb.loader.exec_module(B)
spf = importlib.util.spec_from_file_location('f', os.path.join(HIER, 'forge_c.py'))
F = importlib.util.module_from_spec(spf); spf.loader.exec_module(F)
nabla, kawauchi_ok = B.nabla, B.kawauchi_ok
A = B.A
TMO_KAND = int(os.environ.get('TMO_KAND', '300'))


def eval_one_c(k, seed, nxy, nz, dc, crot, NS):
    import regina
    P = F.forge_c(k, seed, nxy, nz, dc, c=crot)
    if F.sigma_c_ok(P, k, crot) > 1e-11:
        return dict(rej='sigma')
    idx = np.linspace(0, len(P), NS, endpoint=False).astype(int)
    sg = A.signed_gauss(P[idx])
    if sg is None:
        return dict(rej='tangenz')
    L = regina.Link.fromSignedGauss(sg); L.simplify()
    if L.size() == 0:
        return dict(rej='unknoten')
    return dict(c=int(L.size()), alex=str(L.alexander()), pd=str(L.pd()), c_raw=len(sg) // 3)


def erledigt():
    """Done = DECIDED. A deadline kill is NOT a decision.

    🔴 Finding 29.07.2026: until here every written row counted as 'computed' —
    also `rej='frist'`. A rerun with a larger deadline would thus have stayed at exactly the
    candidates that it was meant to rescue: the undecided would have stayed
    PERMANENTLY undecided, and the run would have reported 'nothing open'.
    Measured on that day: at crot=13, 94.1 % of the rows were deadline kills.
    Three categories, not two — rejected / passed / undecided.
    """
    d = set()
    offen_frist = 0
    for f in glob.glob(os.path.join(HIER, 'breit_c_out', '*.jsonl')):
        for ln in open(f):
            try: r = json.loads(ln)
            except Exception: continue
            if str(r.get('rej', '')).startswith('frist'):
                offen_frist += 1
                continue
            d.add((r['k'], r['seed'], r['nxy'], r['nz'], r['dc'], r['crot']))
    if offen_frist:
        print(f'  ↻ {offen_frist} deadline kills count as UNDECIDED and are retried '
              f'(TMO_KAND={TMO_KAND}s). With the same deadline they die again — raise the deadline.',
              flush=True)
    return d


def _ein_kandidat(k, seed, nxy, nz, dc, crot):
    r = eval_one_c(k, seed, nxy, nz, dc, crot, 1440)
    if 'rej' in r:
        return r
    try:
        nb = nabla(r['alex'])
        r['kaw'] = None if nb is None else bool(kawauchi_ok(nb))
        r['nabla'] = None if nb is None else str(nb.as_expr())
    except Exception as e:
        r['kaw'] = None; r['nabla_err'] = f'{type(e).__name__}'
    if r.get('kaw') is False:
        try:
            hi = eval_one_c(k, seed, nxy, nz, dc, crot, 4320)
            r['confirm_4320'] = (hi.get('alex') == r['alex'])
            r['pd_hi'] = hi.get('pd'); r['c_hi'] = hi.get('c')
        except Exception:
            r['confirm_4320'] = None
    return r


def worker(jobs, path):
    """Like breit3: throwaway child per candidate, deadline in the parent (Disc #592).
    signal.alarm in the computing process would be useless here — regina.simplify() is a
    long C++ call during which Python delivers no signals."""
    import regina  # noqa: F401
    f = open(path, 'w')
    for (k, seed, nxy, nz, dc, crot) in jobs:
        rd, wr = os.pipe()
        pid = os.fork()
        if pid == 0:
            os.close(rd)
            try:
                r = _ein_kandidat(k, seed, nxy, nz, dc, crot)
            except Exception as e:
                r = dict(rej=f'{type(e).__name__}')
            try:
                os.write(wr, json.dumps(r).encode())
            finally:
                os.close(wr); os._exit(0)
        os.close(wr)
        schluss = time.time() + TMO_KAND
        roh = b''
        os.set_blocking(rd, False)
        while True:
            fertig_pid, _ = os.waitpid(pid, os.WNOHANG)
            try:
                d = os.read(rd, 65536)
                if d: roh += d
            except BlockingIOError:
                pass
            except OSError:
                pass
            if fertig_pid:
                try:
                    while True:
                        d = os.read(rd, 65536)
                        if not d: break
                        roh += d
                except Exception:
                    pass
                break
            if time.time() > schluss:
                try: os.kill(pid, 9)
                except ProcessLookupError: pass
                os.waitpid(pid, 0); roh = b''; break
            time.sleep(0.05)
        os.close(rd)
        try:
            r = json.loads(roh) if roh else dict(rej='frist', frist_s=TMO_KAND)
        except Exception:
            r = dict(rej='frist_teilantwort', frist_s=TMO_KAND)
        f.write(json.dumps(dict(k=k, seed=seed, nxy=nxy, nz=nz, dc=dc, crot=crot, **r)) + '\n')
        f.flush()
    f.close(); os._exit(0)


def mem_mb():
    for ln in open('/proc/meminfo'):
        if ln.startswith('MemAvailable:'):
            return int(ln.split()[1]) // 1024
    return 1 << 30


if __name__ == '__main__':
    NW = int(os.environ.get('NW', '22'))
    CHUNK = int(os.environ.get('CHUNK', '80'))
    MEMFLOOR = int(os.environ.get('MEMFLOOR', '4000'))
    SEEDS = int(os.environ.get('SEEDS', '1500'))
    CROTS = [int(x) for x in os.environ.get('CROTS', '3,5,7').split(',')]
    K = int(os.environ.get('K', '4'))
    alle = [(K, s, nxy, nz, dc, cr) for cr in CROTS for s in range(1, SEEDS + 1)
            for nxy in (2, 3, 4) for nz in (2, 3) for dc in (0.45, 0.55, 0.7)]
    fertig = erledigt()
    jobs = [j for j in alle if j not in fertig]
    # ORDER (29.07., a=4 run): sorting by crot means that crot=3 runs through completely
    # before 5/11/13 begin — an abort then delivers ONLY one parameter.
    # With ORDER=verschraenkt all crot run forward together, every state is evaluable.
    if os.environ.get('ORDER', 'crot') == 'verschraenkt':
        jobs.sort(key=lambda j: (j[1], j[5]))      # seed first -> all crot together
    else:
        jobs.sort(key=lambda j: (j[5], j[1]))      # c=3 first (sharp), then 5, then 7
    cc = collections.Counter(j[5] for j in jobs)
    print(f'=== SIEVE c: k={K} (a=3), {len(jobs)} open, per crot {dict(sorted(cc.items()))}, '
          f'{NW} workers x CHUNK {CHUNK} ===', flush=True)
    if not jobs:
        print('nothing open.'); sys.exit(0)

    os.makedirs(os.path.join(HIER, 'breit_c_out'), exist_ok=True)
    stamp = time.strftime('%Y%m%d_%H%M')
    bloecke = [jobs[i:i + CHUNK] for i in range(0, len(jobs), CHUNK)]
    print(f'{len(bloecke)} blocks of {CHUNK}', flush=True)

    t0 = time.time(); lebt = {}; naechster = 0; nw = NW; letzte = 0.0; fertig_bl = 0
    while naechster < len(bloecke) or lebt:
        while naechster < len(bloecke) and len(lebt) < nw:
            pid = os.fork()
            if pid == 0:
                worker(bloecke[naechster],
                       os.path.join(HIER, 'breit_c_out', f'c{stamp}_b{naechster:05d}.jsonl'))
            lebt[pid] = (naechster, time.time()); naechster += 1
        pid, st = os.waitpid(-1, os.WNOHANG)
        if pid and pid in lebt:
            lebt.pop(pid); fertig_bl += 1; continue
        mv = mem_mb()
        if mv < MEMFLOOR and len(lebt) > 4:
            nw = max(4, len(lebt) - 2)
            opfer = max(lebt, key=lambda p: lebt[p][1]); bl, _ = lebt[opfer]
            print(f'🔴 MemFloor {mv} MB: block {bl} back into the queue, workers -> {nw}', flush=True)
            try: os.kill(opfer, 9)
            except ProcessLookupError: pass
            lebt.pop(opfer, None); bloecke.append(bloecke[bl]); time.sleep(10); continue
        if mv > MEMFLOOR * 2 and nw < NW:
            nw = min(NW, nw + 1)
        if time.time() - letzte > 300:
            letzte = time.time()
            n = sum(1 for f in glob.glob(os.path.join(HIER, 'breit_c_out', f'c{stamp}_*.jsonl'))
                    for _ in open(f))
            print(f'  … {n}/{len(jobs)}  blocks {fertig_bl}/{len(bloecke)}  '
                  f'{time.time()-t0:.0f}s  MemAvail {mv} MB  workers {len(lebt)}/{nw}', flush=True)
        time.sleep(2)
    print(f'ALL DONE {time.time()-t0:.0f}s', flush=True)
