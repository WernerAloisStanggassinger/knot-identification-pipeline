#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""M/E solver on ALL 102 distinct knots — also the 84 from k=4 (a=3, order 8).

JUSTIFICATION why this is allowed (Conant 2607.18655, read myself):
  Corollary 7.3 is FORMULATED for a=2, but the identity is the telescoping step
  itself and holds for 1 <= k <= a-1 (Thm 6.11 allows exactly this range):
        eps_l [l] nu_{k+1} = Q_k(u,2) Q_k(u,-2) = M_k^2 - 4 E_k^2
  At k = a-1 nu_a = nabla_K. So for a=3, k=2 -> the same equation on nabla_K.
  Integrality of M,E (M = sum_{j even} 2^j q_j, E = sum_{j odd} 2^{j-1} q_j)
  depends at no point on k=1. Likewise E(0) = (l^2-1)/8 mod 2: uses only M(0) odd
  and evaluation at w=0, where nu_k(0)=1.

The run is thus a TEST AGAINST US: Conant's theorem predicts solvability. If the
solver finds NO solution for a k=4 knot, the suspicion falls first on our data.

Fork per knot (footgun 13/17: isolation at process level, SymPy cache dies with the child),
parent with WNOHANG + deadline (footgun 18).
"""
import json, os, sys, time, importlib.util, collections

HIER = os.path.dirname(os.path.abspath(__file__))
sp = importlib.util.spec_from_file_location('me', os.path.join(HIER, 'me_rechnen.py'))
me = importlib.util.module_from_spec(sp); sp.loader.exec_module(me)   # __main__ guard: calibration
# does NOT run along -> explicitly forced below, otherwise nothing is computed.
spb = importlib.util.spec_from_file_location('b', os.path.join(HIER, 'sieb', 'breit.py'))
B = importlib.util.module_from_spec(spb); spb.loader.exec_module(B)

PAR = int(os.environ.get('PAR', '6'))
TMO = int(os.environ.get('TMO', '1800'))
OUT = os.path.join(HIER, 'me_ergebnis_a3.json')



def a_aus_k(k):
    """Tower height a from the curve parameter k: the rotary reflection has order 2k = 2^a.
    k=2 -> a=2 (order 4) · k=4 -> a=3 (order 8) · k=8 -> a=4 (order 16).
    BEFORE there stood here a table or a power-of-two expression that SILENTLY misrepresented k=8
    (me_lauf: a=2 instead of 4; bilanz: None -> the knot dropped out of every evaluation).
    Found 29.07.2026 while preparing the order-16 run, before data arose."""
    if k < 2 or (k & (k - 1)):
        raise ValueError(f'k must be a power of two >= 2, is {k}')
    return k.bit_length()

def eine_aufgabe(r):
    """computes one knot; return as dict (json-safe)."""
    nbz = B.nabla(r['alex'])
    if nbz is None:
        return dict(status='nabla_none')
    nab = me.nabla_w(str(nbz.as_expr()))
    ell = r['ell']
    vor = ((abs(ell) ** 2 - 1) // 8) % 2
    a = me.analyse(nab, ell)
    treffer = {f'{l_}|{e_}': v[0] for (l_, e_), v in a.items()
               if isinstance(v[0], int) and v[0] > 0}
    E0 = {f'{l_}|{e_}': v[1] for (l_, e_), v in a.items()
          if isinstance(v[0], int) and v[0] > 0}
    return dict(status='ok', ell=ell, deg_nabla=nbz.degree(), a_tower=a_aus_k(int(r['k'])),
                vorhersage_E0=vor, treffer=treffer, E0=E0,
                loesung_existiert=bool(treffer),
                E0_passt=all(v == [vor] for v in E0.values()) if E0 else None)


def kind(r, fd):
    try:
        res = eine_aufgabe(r)
    except Exception as ex:
        res = dict(status=f'FEHLER {type(ex).__name__}', detail=str(ex)[:120])
    res.update(k=r['k'], seed=r['seed'], nxy=r['nxy'], nz=r['nz'], dc=r['dc'],
               crot=r.get('crot', 1),
               vol=r.get('geo', {}).get('vol'))
    os.write(fd, json.dumps(res, default=str).encode())
    os.close(fd); os._exit(0)


if __name__ == '__main__':
    if not me.eichung():
        print('\n🔴 CALIBRATION FAILED — no computation on our knots.'); sys.exit(1)
    print('\n✓ Calibration passed (solver reproduces M,E of Conant for the EHW knot).\n')

    rows = json.load(open(os.path.join(HIER, 'sieb', 'distinkte_knoten.json')))
    c = collections.Counter(r['k'] for r in rows)
    print(f'=== M/E on {len(rows)} distinct knots, per k {dict(sorted(c.items()))} '
          f'(k=4 -> a=3), PAR={PAR}, deadline {TMO}s ===', flush=True)

    warte = list(rows); lebt = {}; fertig = []
    t0 = time.time()
    while warte or lebt:
        while warte and len(lebt) < PAR:
            r = warte.pop(0)
            rd, wr = os.pipe()
            pid = os.fork()
            if pid == 0:
                os.close(rd); kind(r, wr)
            os.close(wr)
            lebt[pid] = (r, rd, time.time())
        pid, st = os.waitpid(-1, os.WNOHANG)          # NEVER blocking (footgun 18)
        if pid and pid in lebt:
            r, rd, ts = lebt.pop(pid)
            try:
                roh = b''
                while True:
                    d = os.read(rd, 65536)
                    if not d: break
                    roh += d
                res = json.loads(roh) if roh else dict(status='leer', k=r['k'], seed=r['seed'])
            except Exception as ex:
                res = dict(status=f'lesefehler {type(ex).__name__}', k=r['k'], seed=r['seed'])
            os.close(rd)
            res['sek'] = round(time.time() - ts, 1)
            fertig.append(res)
            flag = ('✓' if res.get('loesung_existiert') else
                    ('—' if res.get('status') == 'ok' else '🔴'))
            print(f"  [{len(fertig):3d}/{len(rows)}] k={res.get('k')} seed={res.get('seed')} "
                  f"a={res.get('a_tower')} ell={res.get('ell')} deg={res.get('deg_nabla')} "
                  f"{flag} {res.get('status')} {res['sek']}s", flush=True)
            json.dump(fertig, open(OUT, 'w'), indent=1)
            continue
        for pid, (r, rd, ts) in list(lebt.items()):
            if time.time() - ts > TMO:
                print(f"  ⏱️  deadline: killing k={r['k']} seed={r['seed']} (pid {pid})", flush=True)
                try: os.kill(pid, 9)
                except ProcessLookupError: pass
        time.sleep(0.5)

    ok = [f for f in fertig if f.get('status') == 'ok']
    # Do NOT hard-code groups: on 29.07.2026 this balance printed only a=2 and a=3 and
    # concealed the seven order-16 knots, which very much stood in the output file
    # ("data complete, printout silent"). Now: every tower height that occurs gets a
    # row, and what would fall into NO group is counted loudly.
    NAME = {2: 'a=2 (k=2, Ordnung 4)', 3: 'a=3 (k=4, Ordnung 8)', 4: 'a=4 (k=8, Ordnung 16)'}  # Ordnung = order
    hoehen = sorted({f.get('a_tower') for f in ok if isinstance(f.get('a_tower'), int)})
    print(f'\n=== BALANCE after {time.time()-t0:.0f}s ===')
    for a in hoehen:
        grp = [f for f in ok if f.get('a_tower') == a]
        mit = [f for f in grp if f['loesung_existiert']]
        e0 = [f for f in grp if f.get('E0_passt') is True]
        print(f'  {NAME.get(a, f"a={a} (order {2**a})")}: {len(grp)} computed · '
              f'Solution exists {len(mit)}/{len(grp)} · E(0) bound satisfied {len(e0)}/{len(mit)}')
    ohne = [f for f in ok if not isinstance(f.get('a_tower'), int)]
    if ohne:
        print(f'  🔴 {len(ohne)} record(s) WITHOUT a_tower — fall into no group, '
              f'k values {sorted({f.get("k") for f in ohne})}. That is the error of 29.07.')
    if len(ok) != sum(1 for f in ok if isinstance(f.get('a_tower'), int)):
        print('  🔴 group sum != computed records.')
    schlecht = [f for f in fertig if f.get('status') != 'ok']
    if schlecht:
        print(f'  🔴 not computed: {len(schlecht)} '
              f'({collections.Counter(f.get("status") for f in schlecht)})')
    print(f'  -> {OUT}')
