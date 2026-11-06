#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Verifies promotion candidates of the Yamada route FRESH — nothing read back.

A candidate is an unresolved_c_le_20 representative for which c >= 2g + n - 1 >= 21
holds. That would be a class change with external effect (unresolved -> novel), so
NOTHING is taken over from corpus_db or the Pflicht record here except the PD code:
  * HFK recomputed (snappy)        -> g
  * HOMFLY recomputed (snappy)     -> MFW lower bound n
  * both held against the stored value; deviation = STOP, not an average
COUNTER-CASE: 4_1 (c=4, g=1, n=2) must BREAK the threshold. A checker that lets
everything through cannot be told apart from one that rejects everything.
"""
import site_config
import ast, json, os, sys, resource
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from yamada_c_bound import collect, je_cluster, bound, mfw_from_homfly
from span_lower_bound import laden

GB = 1024 ** 3
SCHWELLE = 21


def frisch(pd_str, riegel=6 * GB):
    """HFK genus and MFW fresh from the PD. Runs in the CHILD (latch), result via file."""
    import tempfile
    f = tempfile.NamedTemporaryFile(suffix='.json', delete=False).name
    pid = os.fork()
    if pid == 0:
        resource.setrlimit(resource.RLIMIT_AS, (riegel, riegel))
        out = {}
        try:
            import snappy
            L = snappy.Link(ast.literal_eval(pd_str))
            out['g'] = int(L.knot_floer_homology()['seifert_genus'])
        except Exception as e:
            out['g_fail'] = f'{type(e).__name__}: {e}'
        # snappy/spherogram can NOT do HOMFLY (measured 09.09.: no attribute).
        # Regina can, but sits in the base env -> two steps, no guessing.
        out['n_fail'] = 'in own step (regina/base)'
        json.dump(out, open(f, 'w'))
        os._exit(0)
    os.waitpid(pid, 0)
    try:
        return json.load(open(f))
    except Exception:
        return {'g_fail': 'child dead without artefact (latch/signal)'}


BASE_PY = site_config.PYTHON_BASE


def homfly_frisch(pd_str, riegel_s=900):
    """HOMFLY NEW with regina in the base env (snappy/spherogram cannot do it)."""
    import subprocess, textwrap
    code = textwrap.dedent(f"""
        import regina, json
        L = regina.Link.fromPD({pd_str!r})
        print(json.dumps({{'homfly': str(L.homfly())}}))
    """)
    try:
        r = subprocess.run([BASE_PY, '-c', code], capture_output=True, text=True, timeout=riegel_s)
        if r.returncode != 0:
            return None, f'regina rc={r.returncode}: {r.stderr.strip()[:90]}'
        return json.loads(r.stdout.strip().splitlines()[-1])['homfly'], None
    except Exception as e:
        return None, f'{type(e).__name__}: {e}'


def main():
    ks = laden()
    rows = collect(ks)
    tab, _ = je_cluster(rows)
    by_stem = {r.get('stem'): r for r in ks}

    kand = []
    for r in rows:
        if not (r['rep'] and r['rc'] == 'unresolved_c_le_20'):
            continue
        t = tab.get(r['did']) or {}
        s = bound(t.get('hfk_g'), t.get('bi_l'))
        if s is not None and s >= SCHWELLE:
            kand.append((r, s, t))

    print(f'COUNTER-CASE first — 4_1 (c=4, g=1, n=2) MUST break the threshold {SCHWELLE}')
    PD41 = '[(4,2,5,1),(8,6,1,5),(6,3,7,4),(2,7,3,8)]'
    g41 = frisch(PD41)
    h41, e41 = homfly_frisch(PD41)
    n41 = mfw_from_homfly(h41) if h41 else None
    # 4_1 has braid index 3 — the fresh MFW computation MUST hit that, otherwise
    # the HOMFLY route (convention/parser) is wrong and everything on it worthless.
    print(f'   {"OK " if n41 == 3 else "ERROR"} 4_1 MFW fresh = {n41} (braid index 3 known)'
          + (f'  [{e41}]' if e41 else ''))
    if n41 != 3:
        sys.exit('HOMFLY route not calibrated — no verification possible')
    s41 = bound(g41.get('g'), 2)
    ok41 = s41 is not None and s41 < SCHWELLE
    print(f'   {"OK " if ok41 else "ERROR"} 4_1: g={g41.get("g")} -> bound {s41} '
          f'({"below" if ok41 else "ABOVE"} the threshold)\n')
    if not ok41:
        sys.exit('Counter-case failed — checker unusable')

    print(f'{len(kand)} candidates, each recomputed FRESH:\n')
    bestaetigt, verworfen = [], []
    for r, s_alt, t in sorted(kand, key=lambda x: -x[1]):
        rec = by_stem.get(r['stem']) or {}
        pf = rec.get('pflicht_path')
        pd = None
        if pf and Path(pf).exists():
            pd = json.load(open(pf)).get('pd_list_str')
        print(f'  {str(r["did"])[:42]:<44} c_ub={r["c"]}  stored: g={t["hfk_g"]} n={t["bi_l"]} -> {s_alt}')
        if not pd:
            print('      ⛔ no PD — not verifiable, DISCARDED'); verworfen.append(r['did']); continue
        f = frisch(pd)
        h_neu, h_err = homfly_frisch(pd)
        n_neu = mfw_from_homfly(h_neu) if h_neu else None
        g_neu = f.get('g')
        print(f'      freshly computed: g={g_neu} n={n_neu}'
              + (f'   [g: {f.get("g_fail")}]' if f.get('g_fail') else '')
              + (f'   [homfly: {h_err}]' if h_err else ''))
        if n_neu is None:
            print('      ⛔ MFW not reproducible fresh — DISCARDED')
            verworfen.append(r['did']); continue
        if g_neu is None:
            print('      ⛔ HFK not reproducible fresh — DISCARDED'); verworfen.append(r['did']); continue
        if g_neu != t['hfk_g']:
            print(f'      🔴 GENUS DEVIATES ({t["hfk_g"]} stored vs {g_neu} fresh) — DISCARDED')
            verworfen.append(r['did']); continue
        n_use = n_neu if isinstance(n_neu, int) else t['bi_l']
        if isinstance(n_neu, int) and n_neu != t['bi_l']:
            print(f'      ⚠ MFW deviates ({t["bi_l"]} vs {n_neu}) — the SMALLER one applies (lower bound)')
            n_use = min(n_neu, t['bi_l'])
        s_neu = bound(g_neu, n_use)
        if s_neu is not None and s_neu >= SCHWELLE:
            exakt = ' == c_ub ⟹ CROSSING NUMBER EXACT' if s_neu == r['c'] else ''
            print(f'      ✅ CONFIRMED: c >= {s_neu}{exakt}')
            bestaetigt.append((r['did'], s_neu, r['c']))
        else:
            print(f'      ⛔ fresh only {s_neu} < {SCHWELLE} — DISCARDED'); verworfen.append(r['did'])
        print()
    print(f'=== {len(bestaetigt)} confirmed, {len(verworfen)} discarded ===')
    for d, s, c in bestaetigt:
        print(f'   {d}: c >= {s}' + ('  (exact)' if s == c else ''))
    print('\nA confirmed bound is NOT yet a promotion — the class change')
    print('is an external claim and remains Werner\'s decision.')


if __name__ == '__main__':
    main()
