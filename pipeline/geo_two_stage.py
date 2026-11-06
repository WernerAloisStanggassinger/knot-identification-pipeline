#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""geo_two_stage.py — the search cheap, the check expensive.

MEASURED 2026-08-29 on UEBER_k3_s2 (538 tetrahedra at start), 120 s per side:

    normal            353 rounds = 0.34 s/round · 327 geometric · smallest 500
    high_precision      3 rounds = 40.0 s/round ·   1 geometric · smallest 539

**Factor 118.** The run of the morning searched on `high_precision` and needed
76 minutes for 113 rounds (smallest 504) — the same work costs two minutes
in normal mode and finds more. This is the same thought as Disc #737
(`solution_type` as a free pre-filter), only one stage earlier:
**high precision belongs to the verification, not to the search.**

WHY SMALL MATTERS: `tetrahedra_shapes('rect', intervals=True)` blows up at
k3_s2 (528–537 tetrahedra) both a 4 GB and a 12 GB PARI stack, while the
only success so far (k3_s1p5) had **496** tetrahedra. Whether size
really is the lever, nobody knows — exactly that is what stage 2 tests, by trying
the smallest AND the largest collected triangulations.

🔴 THE COUNTER-CASE THAT CARRIES OR TOPPLES THE TWO-STAGE APPROACH:
a triangulation recognised as geometric at normal precision does NOT have to stay so
under `high_precision`. Hence every collected triangulation
is re-classified in stage 2, and the share that fails there is PRINTED.
If everything fails, the procedure is worthless — and that is then written there, instead of
an empty result looking like bad luck.

⚠️ A non-result is a statement about the search, never about the knot
(footgun 14 / Disc #538).
"""
import site_config
import argparse, importlib.util, json, os, resource, time

TOOLS = (site_config.TOOLS_DIR)
GEO = 'all tetrahedra positively oriented'


def lade_gs():
    sp = importlib.util.spec_from_file_location('gs', f'{TOOLS}/geo_search.py')
    gs = importlib.util.module_from_spec(sp); sp.loader.exec_module(gs)
    return gs


def stufe1(M0, sek, log, cache=None):
    """cheap: normal precision, many rounds, collect isosigs

    🔴 The result is written to disk IMMEDIATELY. On 29.08. exactly this search ran
    300 s, found 392 triangulations — and the 90-min `timeout` cut the run off later
    in the middle of stage 2. Because the JSON was only written at the end, EVERYTHING was gone:
    no result, no intermediate state, and the search had to start over. Disc #739 says it
    for the deadline abort, here it was the same family one stage higher.
    If the cache is already there, NO new search is made — stage 2 starts directly on it.
    """
    import snappy, json as _json, os as _os
    if cache and _os.path.exists(cache):
        gesehen = {k: v for k, v in _json.load(open(cache)).items()}
        log(f'  Stage 1 skipped — {len(gesehen)} triangulations from {cache}')
        return gesehen
    M = M0.copy()
    t0 = time.time(); runden = 0; gesehen = {}
    while time.time() - t0 < sek:
        M.randomize(); runden += 1
        if M.solution_type() == GEO:
            n = M.num_tetrahedra()
            iso = M.triangulation_isosig()
            if iso not in gesehen:
                gesehen[iso] = n
    log(f'  Stage 1 (normal): {runden} rounds in {time.time()-t0:.0f}s · '
        f'{len(gesehen)} distinct geometric triangulations · '
        f'sizes {min(gesehen.values()) if gesehen else "—"}..{max(gesehen.values()) if gesehen else "—"}')
    if cache and gesehen:
        import json as _json, os as _os
        _json.dump(gesehen, open(cache + '.tmp', 'w'))
        _os.replace(cache + '.tmp', cache)
        log(f'  Stage 1 saved: {cache}')
    return gesehen


def _sichern(ziel, ergebnis, hp_verloren):
    """Intermediate state to disk immediately. A run that dies at a deadline must be able to say
       what it had seen up to then (Disc #739)."""
    if not ziel:
        return
    import json as _json, os as _os
    _json.dump({'unvollstaendig': True, 'hp_verloren': hp_verloren, 'proben': ergebnis},
               open(ziel + '.tmp', 'w'), indent=1, ensure_ascii=False)
    _os.replace(ziel + '.tmp', ziel)


def stufe2(gesehen, n_klein, n_gross, pari_gb, log, teilziel=None):
    """expensive: high_precision + interval arithmetic, only on selected ones"""
    import snappy
    nach_groesse = sorted(gesehen.items(), key=lambda kv: kv[1])
    auswahl = [('klein', i, n) for i, n in nach_groesse[:n_klein]]
    auswahl += [('gross', i, n) for i, n in nach_groesse[-n_gross:]] if n_gross else []
    ergebnis = []; hp_verloren = 0
    for etikett, iso, n in auswahl:
        M = snappy.Manifold(iso).high_precision()
        sol = M.solution_type()
        if sol != GEO:
            hp_verloren += 1
            log(f'    [{etikett}] {n:4d} tetrahedra: under HP NO LONGER geometric ({sol[:34]})')
            ergebnis.append({'welche': etikett, 'tetraeder': n, 'hp_solution_type': sol,
                             'ausgang': 'HP_VERLOREN'})
            _sichern(teilziel, ergebnis, hp_verloren)
            continue
        t = time.time()
        eintrag = {'welche': etikett, 'tetraeder': n, 'isosig': iso, 'hp_solution_type': sol}
        try:
            sh = M.tetrahedra_shapes('rect', intervals=True)
            eintrag.update(ausgang='SHAPES_DURCH', n_shapes=len(sh))
            log(f'    [{etikett}] {n:4d} tetrahedra: shapes WENT THROUGH ({time.time()-t:.0f}s)')
            t2 = time.time()
            ok, cert = M.verify_hyperbolicity()
            eintrag['verify'] = bool(ok)
            log(f'         verify_hyperbolicity: {ok} ({time.time()-t2:.0f}s)')
            if ok:
                try: eintrag['vol_rigoros'] = str(M.volume(verified=True))
                except BaseException as e: eintrag['vol_err'] = str(e)
                ergebnis.append(eintrag)
                _sichern(teilziel, ergebnis, hp_verloren)
                return ergebnis, hp_verloren, True
        except BaseException as e:
            eintrag.update(ausgang='AUSNAHME', typ=type(e).__name__, text=str(e)[:160])
            log(f'    [{etikett}] {n:4d} tetrahedra: {type(e).__name__}: {str(e)[:80]} ({time.time()-t:.0f}s)')
        eintrag['sek'] = round(time.time() - t, 1)
        ergebnis.append(eintrag)
        _sichern(teilziel, ergebnis, hp_verloren)   # after EVERY probe, not at the end
    return ergebnis, hp_verloren, False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pd')
    ap.add_argument('--such-sek', type=int, default=300)
    ap.add_argument('--klein', type=int, default=6, help='how many of the SMALLEST to check')
    ap.add_argument('--gross', type=int, default=2, help='how many of the LARGEST (counter-check)')
    ap.add_argument('--mb', type=int, default=24000)
    ap.add_argument('--pari-gb', type=float, default=8.0)
    ap.add_argument('--json', default=None)
    a = ap.parse_args()

    resource.setrlimit(resource.RLIMIT_AS, (a.mb * 1024**2, a.mb * 1024**2))
    gs = lade_gs()
    gs.pari_stack(a.pari_gb)
    import snappy

    def log(s): print(s, flush=True)
    pd, stem, quelle = gs.lade_pd(a.pd)
    log(f'[{stem}] {len(pd)} crossings · search {a.such_sek}s normal, '
        f'then HP on {a.klein} smallest + {a.gross} largest · PARI {a.pari_gb} GB')
    M0 = snappy.Link(pd).exterior(); M0.simplify()
    log(f'  exterior+simplify: {M0.num_tetrahedra()} tetrahedra')

    cache = (a.json + '.stufe1.json') if a.json else None
    gesehen = stufe1(M0, a.such_sek, log, cache=cache)
    if not gesehen:
        log('  🔴 no geometric triangulation found — statement about the search, not the knot')
        return
    ergebnis, hp_verloren, treffer = stufe2(gesehen, a.klein, a.gross, a.pari_gb, log,
                                           teilziel=(a.json + '.teil.json') if a.json else None)

    log('\nVERDICT')
    log(f'  lost under HP: {hp_verloren} of {len(ergebnis)} checked')
    if hp_verloren == len(ergebnis):
        log('  ⚠️ ALL triangulations found at normal precision fail under HP —')
        log('     the two-stage search does NOT carry for this knot.')
    if treffer:
        log('  ✅ HYPERBOLICITY VERIFIED')
    else:
        durch = [e for e in ergebnis if e.get('ausgang') == 'SHAPES_DURCH']
        log(f'  no proof in this run · shapes came through for {len(durch)} of {len(ergebnis)}')
    if a.json:
        json.dump({'stem': stem, 'gesehen': len(gesehen), 'hp_verloren': hp_verloren,
                   'treffer': treffer, 'proben': ergebnis},
                  open(a.json, 'w'), indent=1, ensure_ascii=False)


if __name__ == '__main__':
    main()
