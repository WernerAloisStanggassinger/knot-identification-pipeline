#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""cusp_backfill.py — add the cusp sizes that wave 3 never reached.

WHY THIS EXISTS (measured 2026-08-29, not conjectured)
--------------------------------------------------------
448 of 794 records carry no `cusp_areas`, although Werner's
`enrich_optional_invariants.py --phase snappy` has long been able to. Two causes,
both measured — and neither of them is "the enricher cannot do it":

  (1) It dies BEFORE that, in wave 1 (trivial), on missing fields:
        out['nu_invariant'] = a.get('knot_floer_homology', {}).get('nu')
        AttributeError: 'NoneType' object has no attribute 'get'
      For large knots `knot_floer_homology` is missing, for others `conway`.
      Measured over 237 geometric records: **176 crashes, 61 through.**
      ⟹ Wave 3 is never even reached.
  (2) Where it gets through, snappy often reports
        'Manifold has non-geometric solution type contains negatively oriented tetrahedra'
      That is the SAME cause as with the OVER knots (Disc #737): badly
      triangulated, not un-hyperbolic. `randomize()` with `solution_type()` as a
      free filter opens them up (measured: fs_k5 after 2 rounds).

This tool does ONLY wave 3, from the isosig, with randomize beforehand — and leaves
Werner's enricher untouched (Hard Rule 2).

🔴 THE LATCH: the freshly computed volume must match the stored one
--------------------------------------------------------------------------
A triangulation from a foreign isosig looks exactly like the right one.
Therefore EVERY value is discarded whose `M.volume()` deviates relatively by more than
`--vol-toleranz` (default 1e-6) from the record volume. Without a stored
volume nothing is written. This is the same idea as the
Frankenstein latch in `certificate_record.py`.

WHAT IS WRITTEN: purely additive in `optional_invariants`, only keys that are not yet
there. An existing value is NEVER overwritten — if a freshly computed one
deviates from an existing one, that is REPORTED (`--laut`), not smoothed over.
Provenance is stated: `cusp_nachtrag_weg` = "isosig direkt" or "randomize#N".

`--selbsttest` runs four counter-cases (Hard Rule 3a):
  1 calibration — records that already have the value must reproduce it
  2 latch       — a foreign isosig must be REJECTED (otherwise the latch is decoration)
  3 additive    — outside optional_invariants not a byte may change
  4 canary      — a flat (non-hyperbolic) isosig must write NOTHING
"""
import site_config
import argparse, json, os, sys, time, hashlib

GEO = 'all tetrahedra positively oriented'
REPARIEREN = False   # set by --reparieren
CORPUS = (site_config.PIPELINE_DIR + '/corpus_db.json')
# 🔴 `cusp_areas_bits` / `cusp_translations_bits` MUST be in here. In the first
#    run on 29.08. they were missing: the precision escalation ran, entered its values —
#    and the note THAT increased precision was necessary fell by the wayside
#    (`cusp_areas_bits` stood at 0 over 794 records afterwards). A number without its
#    provenance can later no longer be told apart from a default-precision computation.
SCHLUESSEL = ('cusp_areas', 'cusp_area_total', 'cusp_areas_bits',
              'cusp_translations', 'cusp_translations_bits',
              'cusp_volume_max', 'cusp_volume_at_displacement_0', 'cusp_volume_max_via',
              'num_cusps', 'homology_str')


def pari_stack(gb):
    """A PariError instead of a result at higher precision is a stack problem, not a
       mathematical one (Disc #737). Without sage present: carry on silently."""
    try:
        from cypari2 import Pari
        Pari().allocatemem(int(gb * 1024**3), int(2 * gb * 1024**3), silent=True)
        return True
    except BaseException:
        return False


def rechne(isosig, runden):
    """isosig -> geometric triangulation -> cusp sizes. Errors are reported,
       not swallowed (the point where the enricher is blind)."""
    import snappy
    t0 = time.time()
    M = snappy.Manifold(isosig)
    start = M.solution_type()
    runde = 0
    while M.solution_type() != GEO and runde < runden:
        M.randomize(); runde += 1
    if M.solution_type() != GEO:
        return {'ok': False, 'grund': f'not geometric after {runde} rounds (start: {start})',
                'solution_type': M.solution_type(), 'sek': round(time.time() - t0, 2)}
    out = {'ok': True, 'solution_type': M.solution_type(), 'tetraeder': M.num_tetrahedra(),
           'cusp_nachtrag_weg': 'isosig direkt' if runde == 0 else f'randomize#{runde}'}
    # PRECISION ESCALATION (measured 29.08.2026): 75 of 249 geometric
    # triangulations throw at default precision
    #   InsufficientPrecisionError: Could neither verify that the two given tiles are ...
    # With bits_prec=212 they get through. Counter-check on records that already work
    # without: the same value to rel 1e-13 (fs_k5, gl_k4m0) — higher precision
    # does not change the result, it only makes it possible.
    for b in (None, 212, 500):
        try:
            ca = M.cusp_areas() if b is None else M.cusp_areas(bits_prec=b)
            out['cusp_areas'] = [float(x) for x in ca]
            out['cusp_area_total'] = float(sum(ca))
            if b: out['cusp_areas_bits'] = b
            out.pop('cusp_areas_err', None)
            break
        except BaseException as e:
            out['cusp_areas_err'] = f'{type(e).__name__}: {e}'
    for b in (None, 212, 500):
        try:
            tr = M.cusp_translations() if b is None else M.cusp_translations(bits_prec=b)
            out['cusp_translations'] = [{'meridian': str(t[0]), 'longitude': str(t[1])}
                                        for t in tr]
            if b: out['cusp_translations_bits'] = b
            out.pop('cusp_translations_err', None)
            break
        except BaseException as e:
            out['cusp_translations_err'] = f'{type(e).__name__}: {e}'
    # 15.10.2026: volume(0) right after cusp_neighborhood() is displacement 0 = the constant 0.3247595...
    # for EVERY one-cusp knot (Disc #745). This tool still did that after the enricher was fixed, and 39
    # records kept the constant. Same method + canaries as enrich_optional_invariants.py now, plus the
    # one-cusp identity max cusp volume = cusp area / 2.
    try:
        cn = M.cusp_neighborhood()
        v0 = float(cn.volume(0))
        cn.set_displacement(cn.stopping_displacement())
        vmax = float(cn.volume(0))
        vol = float(M.volume())
        if not (vmax == vmax) or vmax <= 0:
            out['cusp_volume_err'] = f'impossible maximum {vmax!r}'
        elif vmax < v0 - 1e-9:
            out['cusp_volume_err'] = f'maximum {vmax!r} below displacement-0 {v0!r}'
        elif vmax > vol + 1e-6:
            out['cusp_volume_err'] = f'cusp volume {vmax!r} above total volume {vol!r}'
        elif M.num_cusps() == 1 and 'cusp_area_total' in out and abs(vmax / out['cusp_area_total'] - 0.5) > 1e-6:
            out['cusp_volume_err'] = f'one cusp: max volume {vmax!r} != area/2 {out["cusp_area_total"] / 2!r}'
        else:
            out['cusp_volume_max'] = vmax
            out['cusp_volume_at_displacement_0'] = v0
            out['cusp_volume_max_via'] = 'set_displacement(stopping_displacement())'
    except BaseException as e:
        out['cusp_volume_err'] = f'{type(e).__name__}: {e}'
    try:
        out['num_cusps'] = int(M.num_cusps())
        out['homology_str'] = str(M.homology())
    except BaseException:
        pass
    try:
        out['volume_frisch'] = float(M.volume())
    except BaseException as e:
        out['volume_err'] = str(e)
    out['sek'] = round(time.time() - t0, 2)
    return out


def riegel(frisch, vol_record, tol):
    """TWO latches, both must hold.

    (1) Frankenstein: without a matching volume nothing is entered.
    (2) Plausibility: a cusp area <= 0 (or not finite) is impossible for a
        hyperbolic knot. 🔴 Found on 29.08. at
        `fleet_n8_2_3_3_5_5_5_7_7`: the same knot delivered on another
        geometric triangulation `cusp_areas = [0.0]` and translations of
        order `8.18 E-2172` — numerical collapse that looks like a result.
        It was caught ONLY by the conflict with an already existing value;
        with an empty field the zero would have moved in. The volume was fine,
        so latch (1) would have let it through.
    """
    import math
    ca = frisch.get('cusp_areas')
    if ca is not None:
        vals = [float(x) for x in ca]
        if not vals or any((v <= 0 or not math.isfinite(v)) for v in vals):
            return False, f'impossible cusp area {vals} — numerical collapse, nothing entered'
    if vol_record in (None, 0):
        return False, 'no stored volume — nothing entered'
    v = frisch.get('volume_frisch')
    if v is None:
        return False, 'no fresh volume computable'
    rel = abs(v - vol_record) / max(abs(vol_record), 1e-12)
    if rel > tol:
        return False, f'Volume deviates: fresh {v:.9f} against record {vol_record:.9f} (rel {rel:.2e})'
    return True, f'Volume matches (rel {rel:.1e})'


def eintragen(pfad, frisch, schreiben, laut):
    d = json.load(open(pfad))
    opt = d.setdefault('optional_invariants', {})
    neu, konflikt = {}, []
    for k in SCHLUESSEL + ('cusp_nachtrag_weg',):
        if k not in frisch:
            continue
        # --reparieren (15.10.2026): ONLY a cusp_volume_max that is exactly the displacement-0 constant
        # (bug signature, Disc #745) may be replaced; old value goes to the protocol before writing.
        if REPARIEREN and k == 'cusp_volume_max' and isinstance(opt.get(k), float) \
                and abs(opt[k] - 0.3247595264191645) < 1e-12:
            neu[k] = frisch[k]; continue
        if k in opt and opt[k] not in (None, [], '', {}):
            if opt[k] != frisch[k]:
                konflikt.append(k)
            continue
        neu[k] = frisch[k]
    if konflikt and laut:
        print(f'    ⚠️ existing values deviate, NOT overwritten: {konflikt}')
    if neu and schreiben:
        opt.update(neu)
        tmp = pfad + '.tmp_cusp'
        with open(tmp, 'w') as f:
            json.dump(d, f, indent=2, ensure_ascii=False)
        os.replace(tmp, pfad)
    return list(neu), konflikt


# ─────────────────────────── Counter-cases ───────────────────────────

def selbsttest():
    print('SELFTEST cusp_nachtrag — four counter-cases\n')
    db = json.load(open(CORPUS))
    vol_von = {r['stem']: r.get('vol') for r in db['knots']}
    mit, ohne_hyp = [], []
    for r in db['knots']:
        p = r.get('pflicht_path')
        if not p or not os.path.exists(p) or not r.get('isosig'):
            continue
        opt = (json.load(open(p)).get('optional_invariants') or {})
        if opt.get('cusp_areas') and len(mit) < 5:
            mit.append((r, opt['cusp_areas'], opt))
        if r.get('is_hyperbolic') is False and len(ohne_hyp) < 2:
            ohne_hyp.append(r)
        if len(mit) >= 5 and len(ohne_hyp) >= 2:
            break

    fehler = 0

    print('1) CALIBRATION — existing value must be reproduced')
    for r, soll, opt in mit:
        f = rechne(r['isosig'], 40)
        ist = f.get('cusp_areas')
        if not ist or len(ist) != len(soll):
            print(f'   🔴 {r["stem"]}: missing/malformed value'); fehler += 1; continue
        rel = max(abs(a - b) / max(abs(b), 1e-12) for a, b in zip(ist, soll))
        ok = rel < 1e-6
        print(f'   {"✅" if ok else "🔴"} {r["stem"][:40]:40s} rel {rel:.1e}')
        fehler += (not ok)
        # 15.10.2026: max cusp volume — never the displacement-0 constant, equal to a correct stored value
        vm = f.get('cusp_volume_max')
        if vm is None or abs(vm - 0.3247595264191645) < 1e-9:
            print(f'   🔴 {r["stem"][:40]}: cusp_volume_max {vm!r} (missing or displacement-0 constant)'); fehler += 1
        elif opt.get('cusp_volume_max_via', '').startswith('set_displacement'):
            relv = abs(vm - opt['cusp_volume_max']) / opt['cusp_volume_max']
            print(f'   {"✅" if relv < 1e-6 else "🔴"} {r["stem"][:40]:40s} cusp_volume_max rel {relv:.1e}')
            fehler += (relv >= 1e-6)

    print('\n2) LATCH — a FOREIGN isosig must be rejected')
    # 🔴 The partner is chosen by VOLUME, not by position in the list.
    #    The first draft took mit[1] — and that was `lj_4_9_11_ph2p356…` against
    #    `lj_4_9_11_ph0p785…`: the SAME knot in two parametrisations, same
    #    volume. The latch was reported green without ever having been tested
    #    (the checker measured itself).
    a = mit[0][0]
    va = vol_von.get(a['stem'])
    fremd = [r for r, _, _o in mit[1:]
             if vol_von.get(r['stem']) and va
             and abs(vol_von[r['stem']] - va) / va > 0.01]
    if not fremd:
        print('   ⚠️ no record with a clearly different volume in the sample — '
              'counter-case NOT run')
        fehler += 1
    else:
        b = fremd[0]
        f = rechne(b['isosig'], 40)                        # values of b ...
        ok, grund = riegel(f, va, 1e-6)                    # ... against the volume of a
        print(f'   {"🔴 WAVED THROUGH" if ok else "✅ rejected"}: {b["stem"][:30]} '
              f'(vol {vol_von[b["stem"]]:.4f}) against {a["stem"][:30]} (vol {va:.4f})')
        print(f'      {grund}')
        fehler += bool(ok)
    #  Opposite direction: the OWN isosig must get through, otherwise it latches everything out
    f2 = rechne(a['isosig'], 40)
    ok2, grund2 = riegel(f2, vol_von.get(a['stem']), 1e-6)
    print(f'   {"✅ own one gets through" if ok2 else "🔴 also latches out the OWN one"}: {grund2}')
    fehler += (not ok2)

    print('\n3) ADDITIVE — outside optional_invariants nothing may change')
    import shutil, tempfile
    quelle = mit[0][0]['pflicht_path']
    tmpd = tempfile.mkdtemp()
    kopie = os.path.join(tmpd, 'r.json'); shutil.copy(quelle, kopie)
    vor = json.load(open(kopie))
    f3 = rechne(mit[0][0]['isosig'], 40)
    eintragen(kopie, f3, True, False)
    nach = json.load(open(kopie))
    aussen_gleich = ({k: v for k, v in vor.items() if k != 'optional_invariants'}
                     == {k: v for k, v in nach.items() if k != 'optional_invariants'})
    print(f'   {"✅" if aussen_gleich else "🔴"} top level outside optional_invariants unchanged')
    fehler += (not aussen_gleich)
    behalten = all(nach['optional_invariants'].get(k) == v
                   for k, v in (vor.get('optional_invariants') or {}).items())
    print(f'   {"✅" if behalten else "🔴"} no existing optional value overwritten')
    fehler += (not behalten)
    shutil.rmtree(tmpd)

    print('\n4) CANARY — a non-hyperbolic knot must deliver NOTHING')
    if ohne_hyp:
        r = ohne_hyp[0]
        f4 = rechne(r['isosig'], 40)
        still = (not f4.get('ok')) or ('cusp_areas' not in f4)
        print(f'   {"✅" if still else "🔴 DELIVERS VALUES"} {r["stem"][:40]:40s} sol={f4.get("solution_type","?")[:34]}')
        fehler += (not still)
    else:
        print('   ⚠️ no non-hyperbolic record with isosig found — counter-case NOT run')
        fehler += 1

    print('\n5) COLLAPSE — cusp_areas = 0.0 must be rejected')
    #   The case from fleet_n8_2_3_3_5_5_5_7_7 (29.08.), re-enacted: volume matches,
    #   the area has collapsed. Latch (1) alone would wave it through.
    fake = {'cusp_areas': [0.0], 'volume_frisch': 33.532891941}
    ok5, grund5 = riegel(fake, 33.532891941, 1e-6)
    print(f'   {"🔴 WAVED THROUGH" if ok5 else "✅ rejected"}: {grund5}')
    fehler += bool(ok5)
    #   Opposite direction: a healthy area at the same volume must get through
    ok6, grund6 = riegel({'cusp_areas': [27.44], 'volume_frisch': 33.532891941}, 33.532891941, 1e-6)
    print(f'   {"✅ healthy area gets through" if ok6 else "🔴 also latches out the healthy"}: {grund6}')
    fehler += (not ok6)

    print(f'\n{"✅ ALL COUNTER-CASES PASSED" if fehler == 0 else f"🔴 {fehler} ERRORS"}')
    return 0 if fehler == 0 else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('liste', nargs='?', help='file with one record path per line')
    ap.add_argument('--runden', type=int, default=60)
    ap.add_argument('--vol-toleranz', type=float, default=1e-6)
    ap.add_argument('--schreiben', action='store_true', help='without this: dry run')
    ap.add_argument('--laut', action='store_true', default=True)
    ap.add_argument('--json', default=None)
    ap.add_argument('--selbsttest', action='store_true')
    ap.add_argument('--pari-gb', type=float, default=2.0)
    ap.add_argument('--reparieren', action='store_true', help='replace the displacement-0 constant in cusp_volume_max (only that)')
    a = ap.parse_args()
    global REPARIEREN; REPARIEREN = a.reparieren
    if not a.selbsttest:
        print('PARI stack raised:', pari_stack(a.pari_gb))
    if a.selbsttest:
        sys.exit(selbsttest())
    if not a.liste:
        ap.error('list missing (or --selbsttest)')

    db = json.load(open(CORPUS))
    vol_von = {r['pflicht_path']: r.get('vol') for r in db['knots'] if r.get('pflicht_path')}
    iso_von = {r['pflicht_path']: r.get('isosig') for r in db['knots'] if r.get('pflicht_path')}
    stem_von = {r['pflicht_path']: r.get('stem') for r in db['knots'] if r.get('pflicht_path')}

    pfade = [l.strip() for l in open(a.liste) if l.strip()]
    print(f'{len(pfade)} records · {"WRITING" if a.schreiben else "DRY"} · '
          f'vol tolerance {a.vol_toleranz} · max {a.runden} randomize rounds\n')
    bilanz = {'geschrieben': 0, 'abgelehnt_riegel': 0, 'nicht_geometrisch': 0,
              'kein_wert': 0, 'nichts_neu': 0, 'konflikte': 0}
    zeilen = []
    t0 = time.time()
    for i, p in enumerate(pfade):
        iso = iso_von.get(p)
        stem = stem_von.get(p, os.path.basename(p))[:44]
        if not iso:
            bilanz['kein_wert'] += 1
            print(f'  {stem:44s} — no isosig'); continue
        f = rechne(iso, a.runden)
        if not f.get('ok'):
            bilanz['nicht_geometrisch'] += 1
            print(f'  {stem:44s} — {f["grund"]}')
            zeilen.append({'pfad': p, 'stem': stem, **f}); continue
        ok, grund = riegel(f, vol_von.get(p), a.vol_toleranz)
        if not ok:
            bilanz['abgelehnt_riegel'] += 1
            print(f'  {stem:44s} 🔴 LATCH: {grund}')
            zeilen.append({'pfad': p, 'stem': stem, 'abgelehnt': grund, **f}); continue
        neu, konflikt = eintragen(p, f, a.schreiben, a.laut)
        bilanz['konflikte'] += bool(konflikt)
        if neu:
            bilanz['geschrieben'] += 1
            ca = f.get('cusp_areas')
            print(f'  {stem:44s} ✅ {len(neu)} fields · {f["cusp_nachtrag_weg"]:14s} '
                  f'· area {round(ca[0],4) if ca else "—"} · {f["sek"]:.1f}s')
        else:
            bilanz['nichts_neu'] += 1
        zeilen.append({'pfad': p, 'stem': stem, 'neu': neu, 'konflikt': konflikt, **f})

    print(f'\n=== BALANCE after {time.time()-t0:.0f}s ===')
    for k, v in bilanz.items():
        print(f'  {k:20s} {v:4d}')
    if a.json:
        json.dump({'bilanz': bilanz, 'zeilen': zeilen}, open(a.json, 'w'), indent=1, ensure_ascii=False)
        print('written:', a.json)


if __name__ == '__main__':
    main()
