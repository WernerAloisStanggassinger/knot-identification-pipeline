#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
# Knot-Identification Pipeline — Werner Alois Stanggassinger & Claude Code (Anthropic).
"""
Certifies hyperbolicity of a corpus record with snappy verify_hyperbolicity()
(interval arithmetic = proof, not numerics).

MUST run in sage_stable — only there do snappy AND sage live together:
  PATH=/home/squeezy/miniforge3/envs/sage_stable/bin:$PATH \
  /home/squeezy/miniforge3/envs/sage_stable/bin/python zertifiziere_hyperbolizitaet.py <stem> [...]

Builds the manifold from the PD code OF THE RECORD (not from a table) —
blind-test doctrine. The stored isosig serves only as counter-check.

MANDATORY ESCALATION (HOT_FOOTGUNS footgun 14 / Disc #538):
  verify_hyperbolicity()[0] == False means "NOT CERTIFIED", NOT "not hyperbolic".
  Usually it is precision. Therefore before any non-hyperbolicity conclusion:
      ManifoldHP  ->  up to 30x randomize()  ->  canonical_retriangulation
  On 19.07. exactly this rescued 26 of 26.

--selbsttest runs known cases in BOTH directions (hyperbolic AND torus).
"""
import json, sys, os, argparse
import site_config

KORPUS = '@PIPELINE_DIR@/corpus_db.json'.replace('@PIPELINE_DIR@', site_config.PIPELINE_DIR)


def hole_record(stem):
    db = json.load(open(KORPUS))
    for eimer in ('knots', 'validation_knots'):
        for r in db.get(eimer, []):
            if r.get('stem') == stem:
                return r, eimer
    return None, None


def pd_aus_record(rec):
    """Pull the PD list from the record file. Several fields, because the records
       grew historically — the first usable one wins, the origin is reported."""
    d = json.load(open(rec['pflicht_path']))
    kandidaten = [
        ('pd_list_str', d.get('pd_list_str')),
        ('raw_subprocess_outputs.snappy.pd_input',
         (d.get('raw_subprocess_outputs') or {}).get('snappy', {}).get('pd_input')),
        ('pflichtliste.DIAGRAMM.pd_code',
         (d.get('pflichtliste') or {}).get('DIAGRAMM', {}).get('pd_code')),
    ]
    for name, v in kandidaten:
        if not v:
            continue
        s = str(v).strip()
        if s.startswith('PD['):                      # regina notation -> tuple
            s = s[3:-1].replace('X[', '(').replace(']', ')')
            s = '[' + s + ']'
        try:
            pd = eval(s, {'__builtins__': {}})
            if pd and all(len(x) == 4 for x in pd):
                return [tuple(int(y) for y in x) for x in pd], name
        except Exception:
            continue
    return None, None


def _versuch(M):
    """verify_hyperbolicity(), but an exception abort is a FAILED ATTEMPT,
       not a crash. For non-geometric solutions (e.g. connected sums) snappy
       throws ZeroDivisionError/LinAlgError from the Krawczyk engine — that is a
       normal result of the stage and must not end the escalation."""
    try:
        return M.verify_hyperbolicity()
    except Exception as e:
        return False, f'<{type(e).__name__}: {e}>'


def eskalation(M0, bauen, log):
    """Footgun 14: only after full escalation may 'not hyperbolic' be said."""
    ok, cert = _versuch(M0)
    log(f"    Stage 1  Manifold             -> {ok}")
    if ok:
        return True, cert, 'Manifold'

    try:
        MH = bauen(hoch=True)
        ok, cert = _versuch(MH)
        log(f"    Stage 2  ManifoldHP           -> {ok}")
        if ok:
            return True, cert, 'ManifoldHP'
    except Exception as e:
        log(f"    Stage 2  ManifoldHP           -> build error {type(e).__name__}: {e}")
        MH = None

    if MH is not None:
        for i in range(1, 31):
            try:
                MH.randomize()
            except Exception as e:
                log(f"    Stage 3  randomize() #{i:<2d}     -> error {type(e).__name__}")
                break
            ok, cert = _versuch(MH)
            if ok:
                log(f"    Stage 3  randomize() #{i:<2d}     -> True")
                return True, cert, f'ManifoldHP+randomize#{i}'
        else:
            log(f"    Stage 3  30x randomize()      -> False")

        try:
            MC = MH.canonical_retriangulation()
            ok, cert = _versuch(MC)
            log(f"    Stage 4  canonical_retri.   -> {ok}")
            if ok:
                return True, cert, 'canonical_retriangulation'
        except Exception as e:
            log(f"    Stage 4  canonical_retri.   -> error {type(e).__name__}: {e}")
    return False, None, 'voll eskaliert, nicht zertifiziert'


def eine(stem, log=print):
    import snappy
    rec, eimer = hole_record(stem)
    if rec is None:
        log(f"[{stem}] ❌ no record in corpus_db")
        return None
    pd, quelle = pd_aus_record(rec)
    if pd is None:
        log(f"[{stem}] ❌ no usable PD code in the record")
        return None

    log(f"[{stem}]  ({eimer})  PD from '{quelle}', {len(pd)} crossings")
    log(f"    Record says: is_hyperbolic={rec.get('is_hyperbolic')}  vol={rec.get('vol')}  "
        f"proof={rec.get('hyperbolicity_proof')}")

    def bauen(hoch=False):
        L = snappy.Link(pd)
        M = L.exterior()
        return M.high_precision() if hoch else M

    M = bauen()
    iso = M.triangulation_isosig()
    vol_roh = None
    try:
        vol_roh = float(M.volume())
    except Exception:
        pass
    log(f"    built from PD: {M.num_tetrahedra()} tetrahedra, vol(roh)={vol_roh}")

    ok, cert, weg = eskalation(M, bauen, log)

    erg = dict(stem=stem, eimer=eimer, pd_quelle=quelle, kreuzungen=len(pd),
               verify=bool(ok), weg=weg, vol_roh=vol_roh,
               vol_record=rec.get('vol'), isosig_neu=iso,
               isosig_record=rec.get('isosig'),
               proof_alt=rec.get('hyperbolicity_proof'),
               is_hyperbolic_alt=rec.get('is_hyperbolic'),
               zertifikat=(str(cert)[:400] if cert is not None else None),
               werkzeug=f"snappy {snappy.__version__} in sage_stable",
               aufruf="snappy.Link(<PD aus Record>).exterior().verify_hyperbolicity()")

    if ok and rec.get('vol') is not None and vol_roh:
        d = abs(float(rec['vol']) - vol_roh)
        erg['vol_abweichung'] = d
        log(f"    vol record vs. fresh: |Δ| = {d:.3e}  {'✅' if d < 1e-6 else '⚠️ DEVIATES'}")
    log(f"    ⟹ verify_hyperbolicity = {ok}  (path: {weg})")
    return erg


def selbsttest():
    """Does the checker bite? Known cases in BOTH directions."""
    import snappy
    print("=== SELFTEST (without corpus, snappy only) ===")
    alles_gut = True
    for name, soll in [('4_1', True), ('m004', True), ('6_2', True),
                       ('3_1', False), ('5_1', False), ('8_19', False)]:
        M = snappy.Manifold(name)
        ok, _ = M.verify_hyperbolicity()
        gut = (ok == soll)
        alles_gut &= gut
        print(f"  {name:6s} verify={str(ok):5s} expected={str(soll):5s} {'✅' if gut else '❌'}")
    print("SELFTEST", "GREEN" if alles_gut else "🔴 RED")
    return 0 if alles_gut else 1


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('stems', nargs='*')
    ap.add_argument('--selbsttest', action='store_true')
    ap.add_argument('--out', default=None, help='results as JSON to here')
    a = ap.parse_args()
    if a.selbsttest:
        sys.exit(selbsttest())
    erg = [e for e in (eine(s) for s in a.stems) if e]
    if a.out:
        json.dump(erg, open(a.out, 'w'), indent=2, default=str)
        print("written:", a.out)
