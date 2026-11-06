#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""HFK genus hunt over the open (unresolved_c_le_20) representatives WITHOUT
HFK-exact Seifert genus. Goal: g >= 10 HFK-exact => c >= 2g+1 >= 21
=> promotion to novel_hyperbolic (the genus route in the labeler, Disc #757).

Construction after footgun 13/17/18 (00_START):
  - ONE os.fork per task (no multiprocessing; child dies, memory gone)
  - Memory latch RLIMIT_AS PER CHILD (the box FREEZES on memory shortage,
    Disc #701; HFK dies at the latch as std::bad_alloc/SIGABRT — intended)
  - Parent loop ALWAYS os.waitpid(WNOHANG) + short sleep, deadline kill
    via SIGKILL in the parent (never wait blockingly, footgun 18)
  - Result per child as its own JSON file (artefact truth, no pipe loss)

Counter-case (--selbsttest): 4_1 must deliver g=1; a child with a mini latch
must DIE at the latch (a latch that never bites is decoration).

Reading the results: tools/label_distinct_knots does NOT read these files —
a hit (g>=10) is written into the Pflicht record by a follow-up step
(own, deliberate step; result file != corpus).
"""
import site_config
import ast, json, os, resource, signal, sys, time
from pathlib import Path

ROOT = Path(site_config.ROOT)
DB = Path(site_config.PIPELINE_DIR) / 'corpus_db.json'
OUT = Path(site_config.MACHINE_YARD) / 'corpus_vs_burton/hfk_jagd_20260906'
GB = 1024 ** 3

# (latch, deadline) by c_ub: small knots short+narrow, giants fat+long.
# Budget is PER BOX (31 GB): worst-case sum of the latches < ~22 GB.


def budget(c_ub):
    if c_ub is None or c_ub <= 40:
        return 2 * GB, 3600
    return 4 * GB, 6 * 3600


def targets():
    db = json.load(open(DB))
    out = []
    for r in db['knots']:
        if not (r.get('is_distinct_representative')
                and r.get('realizes_class') == 'unresolved_c_le_20'):
            continue
        pf = r.get('pflicht_path')
        if not pf or not Path(pf).exists():
            continue
        d = json.load(open(pf))
        hfk = (d.get('pflichtliste', {}).get('ALGEBRA', {})
               .get('knot_floer_homology') or {})
        if hfk.get('seifert_genus') is not None:
            continue  # HFK-exact already there -> labeler has decided
        pd = d.get('pd_list_str')
        if not pd:
            continue
        out.append({'did': r.get('distinct_knot_id'), 'stem': r.get('stem'),
                    'c_ub': r.get('c'), 'genus_feld': r.get('genus'),
                    'det': r.get('det'), 'pflicht_path': pf, 'pd': pd})
    return out


def kind(t, riegel, out_file):
    """Runs ONLY in the child. Computes HFK, writes JSON, exit."""
    resource.setrlimit(resource.RLIMIT_AS, (riegel, riegel))
    t0 = time.time()
    res = {'did': t['did'], 'stem': t['stem'], 'c_ub': t['c_ub'],
           'genus_feld': t['genus_feld'], 'riegel_gb': riegel / GB}
    try:
        import snappy
        pd = ast.literal_eval(t['pd'])
        L = snappy.Link(pd)
        h = L.knot_floer_homology()
        res['hfk'] = {
            'seifert_genus': int(h['seifert_genus']),
            'tau': int(h['tau']), 'nu': int(h['nu']),
            'epsilon': int(h['epsilon']), 'fibered': bool(h['fibered']),
            'L_space_knot': bool(h['L_space_knot']),
            'total_rank': int(h['total_rank']), 'modulus': int(h['modulus']),
        }
        res['sekunden'] = round(time.time() - t0, 1)
    except MemoryError as e:
        res['fail'] = f'MemoryError (latch): {e}'
    except Exception as e:
        res['fail'] = f'{type(e).__name__}: {e}'
    json.dump(res, open(out_file, 'w'), indent=1)
    os._exit(0)


def selbsttest():
    OUT.mkdir(parents=True, exist_ok=True)
    # Probe 1: 4_1 delivers g=1
    t = {'did': 'TEST_4_1', 'stem': 'test', 'c_ub': 4, 'genus_feld': None,
         'pd': '[(4,2,5,1),(8,6,1,5),(6,3,7,4),(2,7,3,8)]'}
    f = OUT / '_selbsttest_4_1.json'
    pid = os.fork()
    if pid == 0:
        kind(t, 4 * GB, f)
    os.waitpid(pid, 0)
    g = json.load(open(f)).get('hfk', {}).get('seifert_genus')
    assert g == 1, f'Counter-case 4_1: g={g} != 1'
    print('Probe 1 OK: 4_1 -> g=1')
    # Probe 2: mini latch MUST bite (32 MB is not even enough for snappy)
    f2 = OUT / '_selbsttest_riegel.json'
    if f2.exists():
        f2.unlink()
    pid = os.fork()
    if pid == 0:
        kind(t, 32 * 1024 * 1024, f2)
    _, status = os.waitpid(pid, 0)
    beisst = (not f2.exists()) or ('fail' in json.load(open(f2)))
    assert beisst, 'Latch counter-case: 32 MB child ran THROUGH — latch is decoration'
    print(f'Probe 2 OK: mini latch bites (exit status {status})')
    print('SELFTEST GREEN (2 probes)')


def main():
    if '--selbsttest' in sys.argv:
        selbsttest()
        return
    OUT.mkdir(parents=True, exist_ok=True)
    ts = targets()
    print(f'{len(ts)} targets without HFK-exact genus', flush=True)
    lauf = {}  # pid -> (target, deadline, out_file)
    for t in ts:
        riegel, frist = budget(t['c_ub'])
        safe = t['did'].replace('/', '_').replace('#', '_')
        f = OUT / f'{safe}.json'
        if f.exists() and 'hfk' in json.load(open(f)):
            print(f'  skip (done): {t["did"]}', flush=True)
            continue
        pid = os.fork()
        if pid == 0:
            kind(t, riegel, f)
        lauf[pid] = (t, time.time() + frist, f)
        print(f'  started pid={pid} {t["did"]} (c_ub={t["c_ub"]}, '
              f'{riegel // GB} GB, {frist // 3600} h)', flush=True)
    # Parent loop: WNOHANG ALWAYS (footgun 18), deadline kill via SIGKILL
    while lauf:
        try:
            pid, status = os.waitpid(-1, os.WNOHANG)
        except ChildProcessError:
            break
        if pid and pid in lauf:
            t, _, f = lauf.pop(pid)
            verdict = '?'
            if f.exists():
                r = json.load(open(f))
                g = r.get('hfk', {}).get('seifert_genus')
                verdict = (f'g={g}' + (' >=10 PROMOTION CANDIDATE!' if isinstance(g, int) and g >= 10 else '')
                           ) if g is not None else f'FAIL: {r.get("fail", "?")[:100]}'
            else:
                verdict = f'CHILD DEAD without artefact (status={status}) — latch/signal'
            print(f'  done {t["did"]}: {verdict}', flush=True)
            continue
        now = time.time()
        for pid in list(lauf):
            t, deadline, f = lauf[pid]
            if now > deadline:
                os.kill(pid, signal.SIGKILL)
                os.waitpid(pid, 0)
                lauf.pop(pid)
                json.dump({'did': t['did'], 'fail': 'DEADLINE (SIGKILL in parent)'},
                          open(f, 'w'), indent=1)
                print(f'  DEADLINE KILL {t["did"]}', flush=True)
        time.sleep(5)
    # Balance
    treffer = []
    for f in OUT.glob('*.json'):
        if f.name.startswith('_selbsttest'):
            continue
        r = json.load(open(f))
        g = r.get('hfk', {}).get('seifert_genus')
        if isinstance(g, int) and g >= 10:
            treffer.append((r['did'], g))
    print(f'=== HUNT DONE: {len(treffer)} promotion candidates (g>=10) ===', flush=True)
    for did, g in treffer:
        print(f'  {did}: g={g} => c>={2*g+1}', flush=True)
    print('Result file != corpus — the follow-up entry into Pflicht records is a step of its own.', flush=True)


if __name__ == '__main__':
    main()
