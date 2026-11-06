#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""certificate_record.py — write a hyperbolicity certificate into the record.

Input is the output of `geo_search.py` (or `cert_failure.py`), target is the
pipeline record. The block is called `_HYPERBOLICITY_CERT_<date>` and carries the same
keys as the one from 27.08. (readers: `build_corpus_db._zertifikate`,
`build_web_export._zertifikat_bewiesen` — both search for the prefix and read the
verbatim verdict).

FOUR LATCHES, each checked with a counter-case (footgun 16). A tool that writes a
certificate without being able to refuse is a rubber-stamp machine:
  1. `zertifiziert` must be TRUE — no entry from a failed run.
  2. The rigorous volume must match the volume OF THE RECORD (default 1e-6).
     This is the Frankenstein latch: it catches the case that the certificate belongs to
     a DIFFERENT knot than the record (Disc #251/#252).
  3. An already existing `_HYPERBOLICITY_CERT_*` block is NEVER overwritten.
  4. The shape of the file stays: indentation preserved, `ensure_ascii=False` (Disc #736 (b)).
"""
import site_config
import json, os, sys, shutil, argparse, datetime

BACKUP = (site_config.BACKUP_DIR + '/20260828')


def _zahl(s):
    """'277.2961687388877473686184242?' -> float. Sage intervals end in '?'."""
    if s is None:
        return None
    t = str(s).strip().rstrip('?')
    try:
        return float(t)
    except ValueError:
        return None


def pruefe(erg, rec_vol, vorhandene_bloecke, toleranz=1e-6):
    """('schreiben', delta) or ('ablehnen', reason) — the verdict, without side effect."""
    if not erg.get('zertifiziert'):
        return 'ablehnen', 'run did NOT certify'
    if vorhandene_bloecke:
        return 'ablehnen', f'Record traegt schon {vorhandene_bloecke[0]} — nichts ueberschreiben'
    # First choice is the RIGOROUS volume. If it is missing, the numerical candidate may
    # carry the Frankenstein comparison — but ONLY the comparison; the block then says
    # explicitly that the volume is not certified (k3_s1p5, 28.08.:
    # `volume(verified=True)` threw `Can't take the argument of interval strictly
    # containing zero`, hyperbolicity had already been proven by then).
    v = _zahl(erg.get('vol_rigoros'))
    vol_zertifiziert = v is not None
    if v is None:
        v = _zahl(erg.get('vol_kandidat'))
    if v is None:
        return 'ablehnen', 'neither vol_rigoros nor vol_kandidat readable'
    if rec_vol is None:
        return 'ablehnen', 'Record carries no volume to compare'
    d = abs(v - float(rec_vol))
    if d > toleranz:
        return 'ablehnen', (f'volume deviates: certificate {v}, record {rec_vol}, |Δ|={d:.3e} '
                            f'> {toleranz:.0e} — the certificate possibly belongs to a '
                            f'DIFFERENT knot')
    return ('schreiben' if vol_zertifiziert else 'schreiben_ohne_volumen'), d


def eintragen(erg_json, record, schreiben=False, datum=None):
    erg = json.load(open(erg_json))
    if isinstance(erg, list):
        erg = erg[0]
    rohtext = open(record).read()
    d = json.loads(rohtext)
    H = (d.get('pflichtliste') or {}).get('HYPERBOLISCH', {}) or {}
    rec_vol = H.get('hyperbolic_volume')
    bloecke = [k for k in d if k.startswith('_HYPERBOLICITY_CERT')]
    art, info = pruefe(erg, rec_vol, bloecke)
    vol_zert = (art == 'schreiben')
    if art == 'ablehnen':
        print(f"🔴 REJECTED: {info}")
        return 1
    datum = datum or datetime.date.today().isoformat()
    # NOTE: the German text values stay German on purpose — build_corpus_db reads `verdikt` verbatim
    # and zertifikat_englisch.py matches these sentences literally (SAETZE).
    block = {
        'verdikt': 'HYPERBOLISCH BEWIESEN',
        'methode': 'snappy.verify_hyperbolicity (Intervallarithmetik, SnapPy in Sage)',
        'solution_type': 'all tetrahedra positively oriented',
        'vol_rigoros': str(erg.get('vol_rigoros')) if vol_zert else None,
        'vol_numerisch': erg.get('vol_kandidat'),
        'vol_valid': bool(vol_zert),
        'volume_certified': bool(vol_zert),
        'folgerung': (('Hyperbolicity rigorously proven; the stored volume is a '
                       'hyperbolisches Volumen.') if vol_zert else
                      ('Hyperbolizitaet rigoros bewiesen (verify_hyperbolicity = True). '
                       'The VOLUME is thus NOT certified: volume(verified=True) '
                       'failed at the interval arithmetic. The stored number is '
                       'a numerical candidate, not a proof. ⚠️ The corpus field is called '
                       '`vol_valid` and is derived from the verdict — here it carries '
                       'HYPERBOLICITY, not the volume.')),
        'when': datum,
        'blindtest': 'computed from the stored PD, no table lookup',
        'weg': (f"ManifoldHP+randomize#{erg.get('runde')}" if erg.get('hp')
                else f"Manifold+randomize#{erg.get('runde')}"),
        'suche': (f"{erg.get('tetraeder_start')} tetrahedra at the start, certified at "
                  f"{erg.get('tetraeder')}; solution_type-Zaehlung {erg.get('solution_types')}"),
        'vol_abweichung_zum_record': info,
        'werkzeug': 'snappy 3.3.2 in sage_stable, Suche via claude_setup_v2/tools/geo_suche.py',
        'beleg': erg_json,
    }
    schluessel = f"_HYPERBOLICITY_CERT_{datum.replace('-', '')}"
    print(f"✅ {os.path.basename(record)}: {schluessel} — |Δvol| = {info:.3e}")
    if not schreiben:
        print('   (dry run — write for real with --schreiben)')
        return 0
    os.makedirs(BACKUP, exist_ok=True)
    ziel = os.path.join(BACKUP, os.path.basename(record) + f'.bak_{datum.replace("-", "")}_hypcert')
    if not os.path.exists(ziel):
        shutil.copy2(record, ziel)
    d[schluessel] = block
    zeilen = rohtext.split('\n', 2)
    tiefe = len(zeilen[1]) - len(zeilen[1].lstrip(' ')) if len(zeilen) > 1 else 1
    tmp = record + '.tmp_cert'
    with open(tmp, 'w') as f:
        json.dump(d, f, indent=max(1, tiefe), ensure_ascii=False)
    os.replace(tmp, record)
    print(f"   written · backup {ziel}")
    return 0


def selbsttest():
    ok = tot = 0

    def pruef(name, ist, soll):
        nonlocal ok, tot
        tot += 1
        gut = ist == soll
        ok += gut
        print(f"  {'✅' if gut else '🔴'} {name}: {ist!r} (expected {soll!r})")

    gut = dict(zertifiziert=True, vol_rigoros='277.2961687388877473686184242?')
    pruef('clean run ⟹ schreiben', pruefe(gut, 277.29616873888773, [])[0], 'schreiben')
    # without rigorous volume: write, but explicitly WITHOUT a volume claim
    ohne = dict(zertifiziert=True, vol_kandidat=445.6796705444772)
    pruef('candidate only ⟹ schreiben_ohne_volumen',
          pruefe(ohne, 445.6796705444777, [])[0], 'schreiben_ohne_volumen')
    # 🔴 and the Frankenstein latch applies to the candidate TOO
    pruef('candidate of a foreign knot ⟹ ablehnen',
          pruefe(dict(zertifiziert=True, vol_kandidat=318.046), 445.6796705444777, [])[0],
          'ablehnen')
    pruef('not certified ⟹ ablehnen',
          pruefe(dict(gut, zertifiziert=False), 277.29616873888773, [])[0], 'ablehnen')
    pruef('block already there ⟹ ablehnen',
          pruefe(gut, 277.29616873888773, ['_HYPERBOLICITY_CERT_20260827'])[0], 'ablehnen')
    # Frankenstein canary: certificate of a DIFFERENT knot (k3_s3, vol 318.05)
    pruef('foreign volume ⟹ ablehnen',
          pruefe(dict(gut, vol_rigoros='318.0460626703397?'), 277.29616873888773, [])[0], 'ablehnen')
    pruef('no record volume ⟹ ablehnen', pruefe(gut, None, [])[0], 'ablehnen')
    pruef('unreadable vol_rigoros ⟹ ablehnen',
          pruefe(dict(gut, vol_rigoros='NaN-ish'), 277.29616873888773, [])[0], 'ablehnen')
    # Canary in the other direction: a deviation JUST below the bound must pass
    pruef('Δ=1e-9 ⟹ schreiben',
          pruefe(dict(gut, vol_rigoros='277.296168739887?'), 277.29616873888773, [])[0], 'schreiben')
    print(f"\n  {ok}/{tot} probes")
    return ok == tot


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('erg_json', nargs='?')
    ap.add_argument('record', nargs='?')
    ap.add_argument('--schreiben', action='store_true')
    ap.add_argument('--selbsttest', action='store_true')
    a = ap.parse_args()
    if a.selbsttest:
        sys.exit(0 if selbsttest() else 1)
    if not (a.erg_json and a.record):
        ap.error('erg_json and record required')
    sys.exit(eintragen(a.erg_json, a.record, a.schreiben))
