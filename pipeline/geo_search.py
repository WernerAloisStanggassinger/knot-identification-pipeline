#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""geo_search.py — search for a GEOMETRIC triangulation instead of verifying expensively.

MUST run in an environment with snappy (gnottheorie_env or sage_stable; verify_hyperbolicity
needs sage ⟹ sage_stable, see --nur-suche).

WHY IT EXISTS (measured 2026-08-28, not assumed)
------------------------------------------------------
The escalation ladder (`zertifiziere_hyperbolizitaet.py`) calls `verify_hyperbolicity()` after EVERY randomize().
For the UEBER knots one round costs 240–800 s that way — and as measured almost all of it goes to the verification, not to the re-triangulation:

    UEBER_k6_s4 (931 crossings, 1003 tetrahedra):  randomize() = 2 s   solution_type() = 0 s
    the same round in the ladder on 27.08.      :              240 s

The verification CANNOT succeed there at all, though: `solution_type()` says in every
measured round *contains negatively oriented tetrahedra*, and the interval engine
(Krawczyk) only certifies solutions with consistently positively oriented tetrahedra.
⟹ 8 h of compute time on 27./28.08. went into 30 verifications that were
   hopeless before they started; in the same time ~14,000 re-triangulations are possible.

THE SEARCH HERE: many randomize() rounds, `solution_type()` as a filter (free),
verify_hyperbolicity() ONLY for a geometric solution — then with increasing precision.

⚠️ WHAT A NEGATIVE RESULT MEANS: "no geometric solution found in N rounds"
is a statement about the search, NOT about the knot. Not certified never means
not hyperbolic (footgun 14 / Disc #538).

--selbsttest runs counter-cases in BOTH directions (footgun 16): a searcher that only
confirms cannot be told apart from one that waves everything through.
"""
import site_config
import os
import json, os, sys, time, argparse, resource

GEOMETRISCH = ('all tetrahedra positively oriented',)


class FristAus(BaseException):
    """BaseException, not Exception — otherwise an `except Exception` swallows the deadline (Disc #733)."""


def _riegel(mb):
    resource.setrlimit(resource.RLIMIT_AS, (mb * 1024 * 1024, mb * 1024 * 1024))


def lade_pd(pfad):
    """Two origins, both allowed, the origin is reported:
       (a) output of pd_aus_kurve.py  (field `pd_vereinfacht`)
       (b) a pipeline record `*_pflicht_invarianten_complete.json` — then LITERALLY via
           `zertifiziere_hyperbolizitaet.pd_aus_record` (Werner's calibrated reader, Hard-Rule 2)."""
    d = json.load(open(pfad))
    if 'pd_vereinfacht' in d:
        return [tuple(c) for c in d['pd_vereinfacht']], d.get('stem', os.path.basename(pfad)), 'pd_aus_kurve'
    import importlib.util
    T = os.path.join(os.path.dirname(site_config.__file__), 'certify_hyperbolicity.py')   # released as pipeline/certify_hyperbolicity.py
    sp = importlib.util.spec_from_file_location('zh', T)
    zh = importlib.util.module_from_spec(sp); sp.loader.exec_module(zh)
    pd, quelle = zh.pd_from_record({'pflicht_path': pfad})
    if not pd:
        raise SystemExit(f'no usable PD in {pfad}')
    # The record often carries 'myknot' in the field `name` (pipeline default) — the FOLDER
    # carries the real stem. Hence folder first, `name` only as fallback.
    stem = os.path.basename(os.path.dirname(pfad)) or d.get('name')
    return [tuple(c) for c in pd], stem, f'Record ({quelle})' 


def suche(M, runden, log, verifizieren=True, bits=(None, 212, 500), tally=None):
    """randomize loop with solution_type as a filter. Returns (certified, info).

       `tally` is a dict that is written along DURING the run. Reason (28.08.):
       if the deadline is hit, the run leaves this function via the exception and the
       return value is lost — the truncated k6_s4 run reported
       `geometrisch_gefunden: None`, although it had seen FIVE geometric triangulations.
       A truncated run that conceals its findings understates."""
    import snappy  # noqa: F401  (only for the error message, in case the environment is wrong)
    gesehen = tally if tally is not None else {}
    t0 = time.time()
    bestes = None
    for i in range(0, runden + 1):
        if i:
            M.randomize()
        st = str(M.solution_type())
        gesehen[st] = gesehen.get(st, 0) + 1
        if st in GEOMETRISCH:
            bestes = i
            log(f"    [{time.time()-t0:6.0f}s] round {i:5d}: GEOMETRIC ({M.num_tetrahedra()} tetrahedra)")
            if not verifizieren:
                return False, dict(geometrisch_gefunden=True, runde=i, verifiziert=False,
                                   tetraeder=M.num_tetrahedra(), solution_types=gesehen)
            for b in bits:
                t = time.time()
                try:
                    ok, cert = M.verify_hyperbolicity() if b is None else M.verify_hyperbolicity(bits_prec=b)
                except Exception as e:
                    ok, cert = False, f'<{type(e).__name__}: {e}>'
                log(f"       verify bits_prec={b}: {ok}  ({time.time()-t:.0f}s)")
                if ok:
                    return True, dict(geometrisch_gefunden=True, runde=i, bits_prec=b,
                                      tetraeder=M.num_tetrahedra(), solution_types=gesehen)
        elif i % 50 == 0:
            log(f"    [{time.time()-t0:6.0f}s] round {i:5d}: {st}  ({M.num_tetrahedra()} tetrahedra)")
    return False, dict(geometrisch_gefunden=bestes is not None, runde=bestes,
                       runden_gefahren=runden, solution_types=gesehen,
                       sekunden=round(time.time() - t0, 1))


def pari_stack(gb):
    """Raise the PARI stack BEFORE snappy. Why this is here (measured 28.08.):
       `tetrahedra_shapes('rect', intervals=True)` throws at ~1000 tetrahedra
       `PariError: the PARI stack overflows (current size: 1073741824)`.
       `verify_hyperbolicity` catches that in `except (ValueError, RuntimeError)` and
       simply returns `False` — a MEMORY cap looks like a
       mathematical No there. sage's default: 8 MB stack, 1 GB cap."""
    try:
        from sage.all import pari
        vorher = pari.stacksize()
        pari.allocatemem(int(gb * 1024**3), int(4 * gb * 1024**3), silent=True)
        print(f"    PARI-Stack {vorher} -> {pari.stacksize()} (max {pari.stacksizemax()})", flush=True)
    except Exception as e:
        print(f"    PARI-Stack NOT raised: {type(e).__name__}: {e}", flush=True)


def lauf(pfad, runden, sek, hp, ziel, mb, gb_pari=4.0):
    import signal
    if gb_pari:
        pari_stack(gb_pari)
    import snappy
    pd, stem, quelle = lade_pd(pfad)
    print(f"[{stem}] PD {len(pd)} crossings from {quelle}", flush=True)
    signal.signal(signal.SIGALRM, lambda *a: (_ for _ in ()).throw(FristAus('Frist')))
    signal.alarm(sek)
    erg = dict(stem=stem, kreuzungen=len(pd), runden_soll=runden, hp=bool(hp))
    tally = {}                      # grows along, so that a deadline abort keeps the findings
    try:
        t = time.time()
        M = snappy.Link(pd).exterior()
        M.simplify()
        print(f"    exterior+simplify: {M.num_tetrahedra()} tetrahedra ({time.time()-t:.0f}s), "
              f"sol={M.solution_type()}", flush=True)
        erg['tetraeder_start'] = M.num_tetrahedra()
        try:
            erg['vol_kandidat'] = float(M.volume())
        except Exception:
            erg['vol_kandidat'] = None
        if hp:
            t = time.time()
            M = M.high_precision()
            print(f"    high_precision: sol={M.solution_type()} ({time.time()-t:.0f}s)", flush=True)
        ok, info = suche(M, runden, lambda s: print(s, flush=True), tally=tally)
        erg.update(zertifiziert=bool(ok), abgeschnitten=False, **info)
        if ok:
            try:
                erg['vol_rigoros'] = str(M.volume(verified=True))
            except Exception as e:
                erg['vol_fehler'] = f'{type(e).__name__}: {e}'
    except FristAus:
        print(f"    🔴 DEADLINE {sek}s hit — TRUNCATED. Statement about the time limit, "
              f"not about the knot.", flush=True)
        erg.update(zertifiziert=False, abgeschnitten=True, solution_types=dict(tally),
                   geometrisch_gefunden=bool(tally.get('all tetrahedra positively oriented')))
    finally:
        signal.alarm(0)
    print(f"[{stem}] certified={erg.get('zertifiziert')} "
          f"geometrisch_gefunden={erg.get('geometrisch_gefunden')}", flush=True)
    if ziel:
        json.dump(erg, open(ziel, 'w'), indent=1)
        print('written:', ziel, flush=True)
    return erg


def selbsttest():
    pari_stack(1.0)          # both directions must be right even with a raised stack
    import snappy
    ok = tot = 0

    def pruef(name, ist, soll):
        nonlocal ok, tot
        tot += 1
        gut = ist == soll
        ok += gut
        print(f"  {'✅' if gut else '🔴'} {name}: {ist!r} (expected {soll!r})")

    still = lambda s: None
    # 1+2 figure-eight knot: geometric AND certifiable — the yes case
    M = snappy.Manifold('4_1')
    pruef('4_1 solution_type geometric', str(M.solution_type()) in GEOMETRISCH, True)
    z, info = suche(M, 3, still)
    pruef('4_1 is certified', z, True)
    # 3 canary: the volume must be that of the figure-eight knot, not just any
    pruef('4_1 volume 2.0298832', round(float(M.volume()), 7), 2.0298832)
    # 4+5 trefoil: NOT hyperbolic — must NEVER be certified, not even after rounds
    T = snappy.Manifold('3_1')
    z, info = suche(T, 5, still)
    pruef('3_1 is NOT certified', z, False)
    pruef('3_1 finds no geometric solution', info.get('geometrisch_gefunden'), False)
    # 6 torus knot T(2,5) likewise
    z, _ = suche(snappy.Manifold('5_1'), 5, still)
    pruef('5_1 (torus) is NOT certified', z, False)
    # 7 counter-case against the filter itself: a non-geometric type must not count as geometric
    pruef('"negatively oriented" does not count as geometric',
          'contains negatively oriented tetrahedra' in GEOMETRISCH, False)
    print(f"\n  {ok}/{tot} probes")
    return ok == tot


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('pd_json', nargs='?')
    ap.add_argument('--runden', type=int, default=2000)
    ap.add_argument('--sek', type=int, default=7200)
    ap.add_argument('--mb', type=int, default=6000)
    ap.add_argument('--hp', action='store_true', help='high_precision instead of double precision')
    ap.add_argument('--pari-gb', type=float, default=4.0, help='PARI stack in GB (0 = do not touch)')
    ap.add_argument('--json', default=None)
    ap.add_argument('--selbsttest', action='store_true')
    a = ap.parse_args()
    _riegel(a.mb)
    if a.selbsttest:
        sys.exit(0 if selbsttest() else 1)
    if not a.pd_json:
        ap.error('pd_json missing')
    lauf(a.pd_json, a.runden, a.sek, a.hp, a.json, a.mb, a.pari_gb)
