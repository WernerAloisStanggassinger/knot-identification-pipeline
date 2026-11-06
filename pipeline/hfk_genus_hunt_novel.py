#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""HFK genus hunt for the YAMADA route (c >= 2g + n - 1), 09.09.2026.

DIFFERENCE to hfk_genus_hunt.py (which stays untouched):
  * Target group: novel_hyperbolic representatives that have an MFW lower bound n,
    but NO HFK-exact genus. There every g immediately pays into a sharper
    crossing-number bound — and where 2g+n-1 meets the upper bound c_ub,
    the crossing number is EXACT instead of estimated.
  * PARALLELISM CAP. The model forks ALL targets at once; with eight
    targets that works, with 157 it would be 157 x 2 GB latches and the box freezes
    (Disc #701 — it freezes instead of killing). Cap here: concurrently
    running children x latch stays below the ~22 GB budget PER BOX documented in the model.
The proven parts (kind/budget/latch/selftest) are IMPORTED,
not copied — no second place where the latch can rot.
"""
import site_config
import json, os, signal, sys, time, argparse
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hfk_genus_hunt import kind, budget, GB, selbsttest, DB
from yamada_c_bound import collect, je_cluster, bound
from span_lower_bound import laden

OUT = Path((site_config.MACHINE_YARD + '/corpus_vs_burton/hfk_jagd_novel_20260909'))


def ziele(max_c, klassen=('novel_hyperbolic',), min_c=0):
    """Representatives of the chosen classes WITH MFW n and WITHOUT HFK g, small first."""
    ks = laden()
    rows = collect(ks)
    tab, _ = je_cluster(rows)
    by_stem = {r.get('stem'): r for r in ks}
    out = []
    for r in rows:
        if not (r['rep'] and r['rc'] in klassen):
            continue
        t = tab.get(r['did']) or {}
        if t.get('hfk_g') is not None or not isinstance(t.get('bi_l'), int):
            continue
        c_ub = r['c']
        if not isinstance(c_ub, int) or c_ub > max_c or c_ub < min_c:
            continue
        rec = by_stem.get(r['stem']) or {}
        pf = rec.get('pflicht_path')
        if not pf or not Path(pf).exists():
            continue
        try:
            pd = json.load(open(pf)).get('pd_list_str')
        except Exception:
            continue
        if not pd:
            continue
        out.append({'did': r['did'], 'stem': r['stem'], 'c_ub': c_ub,
                    'genus_feld': r['genus'], 'n_mfw': t['bi_l'],
                    'pflicht_path': pf, 'pd': pd})
    out.sort(key=lambda t: t['c_ub'])       # cheap first = early yield
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--selbsttest', action='store_true')
    ap.add_argument('--max-c', type=int, default=40,
                    help='c_ub upper limit; above it HFK dies at the girth bound')
    ap.add_argument('--parallel', type=int, default=10,
                    help='upper limit of concurrent children (second bound; the first is --gb)')
    ap.add_argument('--riegel-gb', type=float, default=0.0,
                    help='latch PER CHILD in GB, overrides budget(). MEASURED 09.09.2026: '
                         'budget() gives c_ub<=40 only 2 GB — widmark_hex24 (c_ub=36) died there '
                         'with SIGABRT (C++ abort), ran through with 12 GB and delivered g=7, '
                         'c>=21. So for c_ub around 36 use 12 instead of 2. budget() itself stays '
                         'untouched (inherited tool, no drift). 0 = budget() '
                         'decides. Without this switch --gb only turns the SUM and '
                         'each child stays at its 2 or 4 GB — a diagnosis ran past '
                         'that on 09.09.2026.')
    ap.add_argument('--gb', type=float, default=20.0,
                    help='Sum of the LATCHES of all running children in GB. The child cap '
                         'alone is not enough: budget() gives large knots 4 GB instead of 2, '
                         '10 children thus demanded 34 GB against 15 GB free (measured 09.09.2026). '
                         'This box FREEZES on memory shortage instead of killing (Disc #701) — '
                         'that is the documented exception to full throttle, so GB counts here.')
    ap.add_argument('--nur-zaehlen', action='store_true')
    ap.add_argument('--min-c', type=int, default=0)
    ap.add_argument('--klassen', default='novel_hyperbolic',
                    help='comma list, e.g. novel_hyperbolic,unresolved_c_le_20')
    a = ap.parse_args()

    if a.selbsttest:
        selbsttest()
        return

    OUT.mkdir(parents=True, exist_ok=True)
    ts = ziele(a.max_c, tuple(a.klassen.split(',')), a.min_c)
    offen = [t for t in ts if not ((OUT / f'{t["did"].replace("/","_")}.json').exists()
                                   and 'hfk' in json.load(open(OUT / f'{t["did"].replace("/","_")}.json')))]
    print(f'{len(ts)} targets (c_ub <= {a.max_c}), of which {len(offen)} open · '
          f'cap: {a.gb:.0f} GB latch sum OR {a.parallel} children — '
          f'whichever applies first', flush=True)
    if a.nur_zaehlen:
        for t in offen[:20]:
            print(f'   c_ub={t["c_ub"]:<4} n={t["n_mfw"]:<3} {t["did"]}')
        return

    lauf = {}
    warteschlange = list(offen)
    treffer = exakt = 0
    budget_bytes = int(a.gb * GB)
    belegt = 0

    def mitnehmen(signum, frame):
        # Without this, a shot-down parent leaves its children as
        # ORPHANS: measured 09.09.2026 eight of them with 2.4 GB each, which kept growing
        # while the parent was already gone. On this box that ends in a freeze,
        # not in the OOM kill (Disc #701).
        for kpid in list(lauf):
            try:
                os.kill(kpid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        print(f'  ABORT signal {signum}: {len(lauf)} children taken along', flush=True)
        os._exit(143)

    signal.signal(signal.SIGTERM, mitnehmen)
    signal.signal(signal.SIGINT, mitnehmen)
    while warteschlange or lauf:
        while warteschlange and len(lauf) < a.parallel:
            riegel, frist = budget(warteschlange[0]['c_ub'])
            if a.riegel_gb > 0:
                riegel = int(a.riegel_gb * GB)
            if lauf and belegt + riegel > budget_bytes:
                break              # wait until a child releases its latch
            t = warteschlange.pop(0)
            f = OUT / f'{t["did"].replace("/","_")}.json'
            pid = os.fork()
            if pid == 0:
                kind(t, riegel, f)
            belegt += riegel
            lauf[pid] = (t, time.time() + frist, f, riegel)
            print(f'  start {t["did"]} (c_ub={t["c_ub"]}, n={t["n_mfw"]}, '
                  f'{riegel // GB} GB, {frist // 3600} h) '
                  f'[{len(lauf)}/{a.parallel}, {belegt / GB:.0f}/{a.gb:.0f} GB]', flush=True)
        try:
            pid, status = os.waitpid(-1, os.WNOHANG)
        except ChildProcessError:
            break
        if pid and pid in lauf:
            t, _, f, riegel_frei = lauf.pop(pid)
            belegt -= riegel_frei
            # 'latch/signal' was a guess: two different causes in one message.
            # The exit status says which it was — so ask instead of assuming.
            if os.WIFSIGNALED(status):
                sig = os.WTERMSIG(status)
                name = signal.Signals(sig).name if sig in set(s.value for s in signal.Signals) else str(sig)
                v = (f'CHILD DEAD without artefact — signal {sig} ({name}); '
                     f'{"C++ abort on memory shortage, latch too small" if name == "SIGABRT" else "hard death"}')
            else:
                v = f'CHILD DEAD without artefact — exit {os.WEXITSTATUS(status)}'
            if f.exists():
                r = json.load(open(f))
                g = r.get('hfk', {}).get('seifert_genus')
                if g is not None:
                    s = bound(g, t['n_mfw'])
                    marke = ''
                    if s is not None and s == t['c_ub']:
                        marke = f'  ⟹ c={s} EXACT!'; exakt += 1
                    elif s is not None and s >= 21:
                        marke = f'  ⟹ c>={s}'; treffer += 1
                    v = f'g={g}{marke}'
                else:
                    v = f'FAIL: {str(r.get("fail"))[:80]}'
            print(f'  done {t["did"]}: {v}', flush=True)
            continue
        now = time.time()
        for pid in list(lauf):
            t, deadline, f, riegel_frei = lauf[pid]
            if now > deadline:
                os.kill(pid, signal.SIGKILL)
                os.waitpid(pid, 0)
                lauf.pop(pid)
                belegt -= riegel_frei
                json.dump({'did': t['did'], 'fail': 'DEADLINE (SIGKILL in parent)'},
                          open(f, 'w'), indent=1)
                print(f'  DEADLINE KILL {t["did"]}', flush=True)
        time.sleep(3)
    print(f'=== DONE: {exakt} exact crossing numbers, {treffer} further sharpenings ===', flush=True)
    print('Result file != corpus — the entry into the Pflicht records is a step of its own.', flush=True)


if __name__ == '__main__':
    main()
