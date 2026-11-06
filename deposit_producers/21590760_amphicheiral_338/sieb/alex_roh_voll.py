#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""ALEX-RAW-FULL — the FULL Alexander polynomial from the raw diagram, without simplify. Sieve v3.

BACKGROUND (01.08.2026, alex_roh_eichung.py): the 581-s median wall of the c=13 sieve
is 100 % `regina.simplify()`; the point identity carries 341/341 (median 1 s,
max 14 s at c_raw=3678). `fromSignedGauss` already reduces on embedding to
n_raw median 272 / max 1104. For UNDECIDED candidates the identity is
not enough, though — kawauchi_ok needs the COEFFICIENTS. This step here.

PROCEDURE (alex_punkt primitives, Disc #570):
  det A_raw(t) = ±t^j * Delta(t). Interpolation over a DEGREE LADDER: first D=512 (covers
  the measured stock: deg = j + span <= ~600), then 1024, then the safe
  bound n+1. Per (prime, stage): values at D+1 support points, Lagrange mod p,
  then 3 CHECK POINTS outside the support points — if the interpolated polynomial does not
  agree there, D was too small (or the computation is wrong): next stage. Only when
  ALL 4 primes are verified: CRT + symmetric lift.
  ⚠️ CRT cap 4 primes ~ 1e37: a coefficient overflow would NOT be seen by the per-prime
  verification. Stock: max ~5e11; the Delta(1) guard below catches wrap practically.

NORMALIZATION — READ OFF THE 341, NOT GUESSED (convention = calibration):
  regina: min degree 0, palindromic, leading coefficient POSITIVE (Delta(1) is ±1,
  in 8 of 341 truly −1 — the rule is "lead>0", NOT "Delta(1)=1").
  Two FREE GUARDS, because mathematically forced: (1) coefficients palindromic
  (Alexander symmetric), (2) Delta(1) = ±1 (knot, no link/split). Violated ->
  rej, never a record.

RECORD (--sieb): format like breit_c._ein_kandidat, but via='punkt_v3', c = n_roh
  (upper bound stays TRUE, only coarser than simplify-c — Disc #231 "<=N"),
  pd = raw PD. kaw via BC.nabla + BC.kawauchi_ok on the alex string (identical
  path as v1/v2). kaw=False -> confirmation NS=4320 as POINT IDENTITY against the
  own alex (16 points, seconds) instead of a full computation, field confirm_4320_punkt.
  results to breit_c_out/v3_*.jsonl — erledigt() sees them at the next start.

SELF-TEST (mandatory, --selbsttest, two-sided):
  (a) CARRY: trefoil+figure-eight exact · 3 smallest decided ones BYTE-EXACT
      against the stored regina string (format included).
  (b) BITE: forced D=64 on a case with deg>64 must at the check points
      fail (rej:verify) — otherwise the ladder checks nothing.

HARNESS: pattern ns_eichung/alex_roh_eichung (fork per candidate, WNOHANG+own clock,
drain pipes continuously, MemFloor without >4 bound, save atomically). v3 children are
LIGHT (~200-400 MB instead of 1-1.8 GB) — MemFloor loses its terror.
"""
import json, os, sys, time, glob, argparse

HIER = os.path.dirname(os.path.abspath(__file__))

ap = argparse.ArgumentParser()
ap.add_argument('--selbsttest', action='store_true')
ap.add_argument('--eiche', action='store_true', help='all decided ones: interpolate fully, compare byte-exactly')
ap.add_argument('--sieb', action='store_true', help='compute open candidates of the grid (v3 records)')
ap.add_argument('--seeds', type=int, default=100)
ap.add_argument('--kk', type=int, default=8)
ap.add_argument('--crot', type=int, default=13)
ap.add_argument('--n-probe', type=int, default=0)
ap.add_argument('--nw', type=int, default=6)
ap.add_argument('--tmo', type=int, default=1800)
ap.add_argument('--memfloor', type=int, default=3500)
AR = ap.parse_args()

import importlib.util                                                    # noqa: E402
spc = importlib.util.spec_from_file_location('bc', os.path.join(HIER, 'breit_c.py'))
BC = importlib.util.module_from_spec(spc); spc.loader.exec_module(BC)
spa = importlib.util.spec_from_file_location('apk', os.path.join(HIER, 'alex_punkt.py'))
AP = importlib.util.module_from_spec(spa); spa.loader.exec_module(AP)

import numpy as np                                                       # noqa: E402

NS = 1440
VARIANTE = 6
SCHLUESSEL = ('k', 'seed', 'nxy', 'nz', 'dc', 'crot')
LEITER = (384, 768, 1536)    # deg(det)=j+span: median ~200, k max 497 -> 768 covers the tail


def mem_mb():
    for ln in open('/proc/meminfo'):
        if ln.startswith('MemAvailable:'):
            return int(ln.split()[1]) // 1024
    return 1 << 30


def entschiedene(crot):
    out, gesehen = [], set()
    for f in glob.glob(os.path.join(HIER, 'breit_c_out', '*.jsonl')):
        for ln in open(f):
            try:
                r = json.loads(ln)
            except Exception:
                continue
            if r.get('crot') != crot or not r.get('alex'):
                continue
            key = tuple(r[s] for s in SCHLUESSEL)
            if key in gesehen:
                continue
            gesehen.add(key)
            out.append(r)
    return out


def roh_link(k, seed, nxy, nz, dc, crot, ns=NS):
    import regina
    P = BC.F.forge_c(k, seed, nxy, nz, dc, c=crot)
    idx = np.linspace(0, len(P), ns, endpoint=False).astype(int)
    sg = BC.A.signed_gauss(P[idx])
    if sg is None:
        return None, None, 'tangenz'
    L = regina.Link.fromSignedGauss(sg)
    pd = AP.parse_pd(str(L.pd()))
    return pd, len(sg) // 3, None


def volles_alex(pd, zwang_D=0):
    """Integer coefficients of det A(t), degree ladder + check points, CRT.

    Returns (co_canonical ascending, D_used) or raises ValueError('verify'/...).
    """
    n = len(pd)
    stufen = [zwang_D] if zwang_D else [d for d in LEITER if d < n + 1] + [n + 1]
    for D in stufen:
        xs = list(range(1, D + 2))                 # D+1 support points
        pruef = [D + 2, D + 3, D + 4]
        per_prime, ok = [], True
        for p in AP.PRIMES:
            ys = [AP.alexander_at(pd, x, VARIANTE, p) for x in xs]
            co_p = AP.lagrange_mod_p([x % p for x in xs], ys, p)
            for xv in pruef:                       # does the ladder bite?
                v = 0
                for c in reversed(co_p):
                    v = (v * xv + c) % p
                if v != AP.alexander_at(pd, xv, VARIANTE, p):
                    ok = False
                    break
            if not ok:
                break
            per_prime.append(co_p)
        if not ok:
            if zwang_D:
                raise ValueError('verify')         # (b) bite: forced D MUST land here
            continue
        co = [AP.crt_lift([pp[i] for pp in per_prime], AP.PRIMES)
              for i in range(len(per_prime[0]))]
        # canonical: strip t^j, leading coefficient positive (regina rule, read off)
        while co and co[0] == 0:
            co.pop(0)
        while co and co[-1] == 0:
            co.pop()
        if not co:
            raise ValueError('null')
        if co[-1] < 0:
            co = [-c for c in co]
        if co != co[::-1]:
            raise ValueError('nicht_palindrom')    # Alexander MUST be symmetric
        if sum(co) not in (1, -1):
            raise ValueError('delta1')             # a knot MUST have Delta(1)=+-1
        return co, D
    raise ValueError('leiter_erschoepft')


def format_regina(co):
    """ascending coefficients -> regina's print form ('x^84 - 89 x^83 + … - 41 x + 1')."""
    teile = []
    for e in range(len(co) - 1, -1, -1):
        c = co[e]
        if c == 0:
            continue
        mag = abs(c)
        if e == 0:
            t = str(mag)
        elif e == 1:
            t = 'x' if mag == 1 else f'{mag} x'
        else:
            t = f'x^{e}' if mag == 1 else f'{mag} x^{e}'
        if not teile:
            teile.append(t if c > 0 else f'- {t}')
        else:
            teile.append(('+ ' if c > 0 else '- ') + t)
    return ' '.join(teile)


def punkt_bestaetigung(pd_hi, alex_str):
    """NS=4320 confirmation as point identity (logic from alex_roh_eichung, candidate intersection)."""
    co = AP.parse_regina(alex_str)
    n = len(pd_hi)
    kand, geprueft = None, 0
    for p in AP.PRIMES[:2]:
        for x in (2, 3, 5, 7, 11, 13, 17, 19):
            v = 0
            for c in reversed(co):
                v = (v * x + c) % p
            if v == 0:
                continue
            d = AP.alexander_at(pd_hi, x, VARIANTE, p)
            if d == 0:
                return False
            q = d * pow(v, p - 2, p) % p
            ks, pw = set(), 1
            for k in range(n + len(co) + 8):
                if q == pw:
                    ks.add((+1, k))
                if q == p - pw:
                    ks.add((-1, k))
                pw = pw * x % p
            if not ks:
                return False
            kand = ks if kand is None else (kand & ks)
            if not kand:
                return False
            geprueft += 1
    return geprueft >= 4


def _kind_sieb(job, wr):
    k, seed, nxy, nz, dc, crot = job
    try:
        t0 = time.time()
        pd, c_raw, rej = roh_link(k, seed, nxy, nz, dc, crot)
        if rej:
            erg = dict(rej=rej)
        else:
            co, D = volles_alex(pd)
            alex = format_regina(co)
            pd_s = ' '.join(f'X[{a},{b},{c},{d}]' for a, b, c, d in pd)
            erg = dict(c=len(pd), alex=alex, pd=pd_s, c_raw=c_raw, n_roh=len(pd),
                       via='punkt_v3', D=D, t_alex=round(time.time() - t0, 1))
            try:
                nb = BC.nabla(alex)
                erg['kaw'] = None if nb is None else bool(BC.kawauchi_ok(nb))
                erg['nabla'] = None if nb is None else str(nb.as_expr())
            except Exception as e:
                erg['kaw'] = None
                erg['nabla_err'] = f'{type(e).__name__}'
            if erg.get('kaw') is False:
                pd_hi, _, rej_hi = roh_link(k, seed, nxy, nz, dc, crot, ns=4320)
                erg['confirm_4320_punkt'] = (None if rej_hi
                                             else punkt_bestaetigung(pd_hi, alex))
                erg['n_roh_4320'] = None if rej_hi else len(pd_hi)
    except ValueError as e:
        erg = dict(rej=f'alex_{e}')
    except Exception as e:
        erg = dict(rej=f'{type(e).__name__}')
    try:
        os.write(wr, json.dumps(erg).encode())
    finally:
        os.close(wr)
        os._exit(0)


def _kind_eiche(r, wr):
    try:
        t0 = time.time()
        pd, c_raw, rej = roh_link(*[r[s] for s in SCHLUESSEL])
        if rej:
            erg = dict(urteil=f'rej:{rej}')
        else:
            co, D = volles_alex(pd)
            s = format_regina(co)
            erg = dict(urteil='BYTE-GLEICH' if s == r['alex'] else 'WEICHT',
                       D=D, n_roh=len(pd), t=round(time.time() - t0, 1))
            if erg['urteil'] == 'WEICHT':
                erg['alex_neu'] = s[:120]
    except ValueError as e:
        erg = dict(urteil=f'rej:alex_{e}')
    except Exception as e:
        erg = dict(urteil=f'rej:{type(e).__name__}')
    try:
        os.write(wr, json.dumps(erg).encode())
    finally:
        os.close(wr)
        os._exit(0)


def fahre(aufgaben, kind, nw, tmo, memfloor, ziel=None, extra=None, jsonl=None):
    """Harness (pattern alex_roh_eichung). aufgaben=list; kind(job, wr).
    ziel: atomic overall JSON. jsonl: append file in breit_c row format (--sieb)."""
    erg, offen, lebt = [], list(aufgaben), {}
    t0, letzte = time.time(), 0.0
    fj = open(jsonl, 'a') if jsonl else None
    while offen or lebt:
        while offen and len(lebt) < nw:
            a = offen.pop(0)
            rd, wr = os.pipe()
            pid = os.fork()
            if pid == 0:
                os.close(rd)
                kind(a, wr)
            os.close(wr)
            os.set_blocking(rd, False)
            lebt[pid] = dict(a=a, rd=rd, start=time.time(), roh=b'')
        try:
            pid, _ = os.waitpid(-1, os.WNOHANG)
        except ChildProcessError:
            pid = 0
        for p, z in list(lebt.items()):
            try:
                d = os.read(z['rd'], 65536)
                if d:
                    z['roh'] += d
            except (BlockingIOError, OSError):
                pass
        if pid and pid in lebt:
            z = lebt.pop(pid)
            try:
                while True:
                    d = os.read(z['rd'], 65536)
                    if not d:
                        break
                    z['roh'] += d
            except Exception:
                pass
            os.close(z['rd'])
            try:
                neu = json.loads(z['roh']) if z['roh'] else dict(rej='leer')
            except Exception:
                neu = dict(rej='teilantwort')
            key = (dict(zip(SCHLUESSEL, z['a'])) if isinstance(z['a'], tuple)
                   else {s: z['a'][s] for s in SCHLUESSEL})
            neu = dict(**key, **neu, gesamt_s=round(time.time() - z['start'], 1))
            erg.append(neu)
            if fj:
                fj.write(json.dumps(neu) + '\n')
                fj.flush()
            if ziel:
                tmp = ziel + '.tmp'
                with open(tmp, 'w') as f:
                    json.dump(erg, f)
                os.replace(tmp, ziel)
            u = neu.get('urteil') or neu.get('rej') or (
                f"kaw={neu.get('kaw')} D={neu.get('D')}")
            print(f"  [{len(erg)}/{len(aufgaben)}] n_roh={neu.get('n_roh')} -> {u}"
                  f"  {neu['gesamt_s']:.0f}s", flush=True)
        jetzt = time.time()
        for p, z in list(lebt.items()):
            if jetzt - z['start'] > tmo:
                try:
                    os.kill(p, 9)
                except ProcessLookupError:
                    pass
                os.waitpid(p, 0)
                os.close(z['rd'])
                lebt.pop(p, None)
                key = (dict(zip(SCHLUESSEL, z['a'])) if isinstance(z['a'], tuple)
                       else {s: z['a'][s] for s in SCHLUESSEL})
                neu = dict(**key, rej='frist', frist_s=tmo)
                erg.append(neu)
                if fj:
                    fj.write(json.dumps(neu) + '\n')
                    fj.flush()
                print(f"  [{len(erg)}/{len(aufgaben)}] -> DEADLINE {tmo}s", flush=True)
        mv = mem_mb()
        if mv < memfloor and lebt:
            opfer = max(lebt, key=lambda p: lebt[p]['start'])
            z = lebt.pop(opfer)
            try:
                os.kill(opfer, 9)
            except ProcessLookupError:
                pass
            os.waitpid(opfer, 0)
            os.close(z['rd'])
            offen.append(z['a'])
            nw = max(1, nw - 1)
            print(f'🔴 MemFloor {mv} MB: back into the queue, workers -> {nw}', flush=True)
            time.sleep(5)
        if time.time() - letzte > 120:
            letzte = time.time()
            print(f'  … {len(erg)}/{len(aufgaben)}  {time.time()-t0:.0f}s'
                  f'  MemAvail {mv} MB  workers {len(lebt)}/{nw}', flush=True)
        time.sleep(0.05)
    if fj:
        fj.close()
    return erg


    # ⚠️ AP.FIG8 is deliberately NOT used here: the PD is INVALID as a PD
    # (edge 2 three times, edge 4 once; regina.fromPD: "strand appears more than
    # twice") — the docstring itself called it "candidate, not dogma", and the
    # doctrine says: calibrate against REAL data, not against what is remembered.
def selbsttest():
    print('(a1) trefoil exact:', flush=True)
    co, D = volles_alex(AP.TREFOIL, zwang_D=len(AP.TREFOIL) + 1)
    print(f'    trefoil: {tuple(co)} expected (1, -1, 1)', flush=True)
    ok_a1 = tuple(co) == (1, -1, 1)
    recs = sorted(entschiedene(AR.crot), key=lambda r: r.get('c_raw') or 0)
    print('\n(a2) 3 smallest decided ones + 1 counter-case BYTE-EXACT against the regina string:', flush=True)
    gg = next((r for r in recs if r.get('kaw') is False), None)
    probe = recs[:3] + ([gg] if gg else [])
    a = fahre(probe, _kind_eiche, nw=4, tmo=AR.tmo, memfloor=AR.memfloor)
    ok_a2 = all(e.get('urteil') == 'BYTE-GLEICH' for e in a)
    print('\n(b) BITE — forced D=64 must fail at the check points:', flush=True)
    pd, _, _ = roh_link(*[recs[0][s] for s in SCHLUESSEL])
    try:
        volles_alex(pd, zwang_D=64)
        ok_b = False
        print('    🔴 no error — the ladder checks NOTHING', flush=True)
    except ValueError as e:
        ok_b = str(e) == 'verify'
        print(f'    ValueError({e}) {"✓" if ok_b else "🔴 wrong reason"}', flush=True)
    print(f"\n(a1) known exact   : {'✅' if ok_a1 else '🔴'}")
    print(f"(a2) byte-equal     : {'✅' if ok_a2 else '🔴'}")
    print(f"(b)  ladder bites   : {'✅' if ok_b else '🔴'}")
    if ok_a1 and ok_a2 and ok_b:
        print('\n✅ Self-test passed.')
        return 0
    print('\n🔴 Self-test NOT passed — results worthless.')
    return 2


def main():
    if AR.selbsttest:
        return selbsttest()
    stamp = time.strftime('%Y%m%d_%H%M')
    if AR.eiche:
        recs = sorted(entschiedene(AR.crot), key=lambda r: r.get('c_raw') or 0)
        if AR.n_probe > 0:
            recs = recs[:AR.n_probe]
        ziel = os.path.join(HIER, f'alex_voll_eichung_{stamp}.json')
        print(f'=== FULL CALIBRATION: {len(recs)} decided ones byte-exact === -> {ziel}', flush=True)
        erg = fahre(recs, _kind_eiche, AR.nw, AR.tmo, AR.memfloor, ziel=ziel)
        gleich = sum(1 for e in erg if e.get('urteil') == 'BYTE-GLEICH')
        weicht = [e for e in erg if e.get('urteil') == 'WEICHT']
        rej = [e for e in erg if str(e.get('urteil', '')).startswith('rej')]
        sp = sorted(e['gesamt_s'] for e in erg if e.get('urteil') == 'BYTE-GLEICH')
        print(f'\n=== BALANCE === BYTE-GLEICH (byte-equal) {gleich}/{len(erg)} · WEICHT (differs) {len(weicht)}'
              f' · discarded {len(rej)}', flush=True)
        if sp:
            print(f'time: median {sp[len(sp)//2]:.0f}s · max {sp[-1]:.0f}s', flush=True)
        for e in weicht[:6]:
            print('  WEICHT (differs):', {s: e[s] for s in SCHLUESSEL}, e.get('alex_neu', '')[:80], flush=True)
        return 0 if (gleich == len(erg) and not weicht) else 1
    if AR.sieb:
        alle = [(AR.kk, s, nxy, nz, dc, AR.crot) for s in range(1, AR.seeds + 1)
                for nxy in (2, 3, 4) for nz in (2, 3) for dc in (0.45, 0.55, 0.7)]
        fertig = BC.erledigt()
        jobs = [j for j in alle if j not in fertig]
        if AR.n_probe > 0:
            jobs = jobs[:AR.n_probe]
        jl = os.path.join(HIER, 'breit_c_out', f'v3_{stamp}.jsonl')
        print(f'=== SIEVE v3 (point path): {len(jobs)} open of {len(alle)} === -> {jl}', flush=True)
        erg = fahre(jobs, _kind_sieb, AR.nw, AR.tmo, AR.memfloor, jsonl=jl)
        ent = [e for e in erg if e.get('alex')]
        gg = [e for e in ent if e.get('kaw') is False]
        rej = [e for e in erg if e.get('rej')]
        print(f'\n=== BALANCE === decided {len(ent)}/{len(erg)} · counter-cases {len(gg)}'
              f' · rej/frist {len(rej)}', flush=True)
        for e in gg:
            print('  COUNTER-CASE:', {s: e[s] for s in SCHLUESSEL},
                  'confirm_4320_punkt', e.get('confirm_4320_punkt'), flush=True)
        return 0
    print('nothing to do (--selbsttest | --eiche | --sieb)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
