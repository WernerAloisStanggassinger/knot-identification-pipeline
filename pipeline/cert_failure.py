#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Hyperbolicity + rigorous volume for records that are NOT in corpus_db
(pipeline_failures). Feeder to Werner's tools/zertifiziere_hyperbolizitaet.py —
its functions are imported VERBATIM, nothing in it is touched (hard rule 2).
The only difference is the origin of the record: failure folder instead of corpus_db.

MUST run in sage_stable (snappy + sage together).

Occasion 27.08.2026: UEBER_k3_s2 carries is_hyperbolic=true / robust_hyperbolicity=true,
but vol=null — snappy ran into "timeout (300s)". In the same record, under
raw_subprocess_outputs.length_spectrum.volume, stands 462.1017807... — next to it, however,
"contains negatively oriented tetrahedra" and four failed Dirichlet constructions.
⚠️ This number is therefore a CANDIDATE, not a finding: it stems from a non-geometric
solution. Exactly for that reason the full escalation with verify_hyperbolicity
(interval arithmetic = proof) runs here, instead of taking the number over.
"""
import site_config
import os
import json, sys, os, importlib.util, resource, signal


class FristAus(BaseException):
    """Inherits from BaseException, NOT from Exception — otherwise the
       `except Exception` in the escalation loop catches it and the deadline vanishes silently."""
    pass

T = os.path.join(os.path.dirname(site_config.__file__), 'certify_hyperbolicity.py')   # released as pipeline/certify_hyperbolicity.py
sp = importlib.util.spec_from_file_location('zh', T)
zh = importlib.util.module_from_spec(sp); sp.loader.exec_module(zh)

def eskalation_laut(M0, bauen, log):
    """Like zh.escalation, but every stage AND every randomize round is logged.
       The order of the stages is taken over verbatim — only the logging is added."""
    import time
    t0 = time.time()
    def stempel():
        return f"[{time.strftime('%H:%M:%S')} +{time.time()-t0:6.0f}s]"

    ok, cert = zh._attempt(M0)
    log(f"    {stempel()} Stage 1  Manifold             -> {ok}")
    if ok:
        return True, cert, 'Manifold'
    try:
        MH = bauen(hoch=True)
        ok, cert = zh._attempt(MH)
        log(f"    {stempel()} Stage 2  ManifoldHP           -> {ok}")
        if ok:
            return True, cert, 'ManifoldHP'
    except Exception as e:
        log(f"    {stempel()} Stage 2  ManifoldHP  build error {type(e).__name__}: {e}")
        MH = None
    if MH is not None:
        for i in range(1, 31):
            tr = time.time()
            try:
                MH.randomize()
            except Exception as e:
                log(f"    {stempel()} Stage 3  randomize #{i:<2d} error {type(e).__name__}")
                break
            ok, cert = zh._attempt(MH)
            log(f"    {stempel()} Stage 3  randomize #{i:<2d} -> {ok}   ({time.time()-tr:.0f}s for the round)")
            if ok:
                return True, cert, f'ManifoldHP+randomize#{i}'
        else:
            log(f"    {stempel()} Stage 3  30x randomize      -> False")
        try:
            MC = MH.canonical_retriangulation()
            ok, cert = zh._attempt(MC)
            log(f"    {stempel()} Stage 4  canonical_retri.   -> {ok}")
            if ok:
                return True, cert, 'canonical_retriangulation'
        except Exception as e:
            log(f"    {stempel()} Stage 4  canonical_retri.  error {type(e).__name__}: {e}")
    return False, None, 'voll eskaliert, nicht zertifiziert'


def zert(pfad, stem, sek=3600):
    import snappy
    d = json.load(open(pfad))
    if 'pd_vereinfacht' in d:                    # output of tools/pd_aus_kurve.py
        pd = [tuple(c) for c in d['pd_vereinfacht']]
        quelle = 'pd_aus_kurve (regina-simplified, without scanner)'
        nh, ls = {}, {}
        stem = d.get('stem', stem)
    else:                                        # record from pipeline_failures
        pd, quelle = zh.pd_from_record({'pflicht_path': pfad})
        nh = d.get('non_hyperbolic_classification') or {}
        ls = (d.get('raw_subprocess_outputs') or {}).get('length_spectrum') or {}
    if pd is None or not pd:
        print(f"[{stem}] ❌ no usable PD"); return None
    print(f"[{stem}]  PD from '{quelle}', {len(pd)} crossings")
    print(f"    Record says: is_hyperbolic={nh.get('is_hyperbolic')} vol={nh.get('vol')} "
          f"decision={nh.get('decision')}")
    print(f"    Candidate from length_spectrum: vol={ls.get('volume')}  (NOT certified)")

    def bauen(hoch=False):
        M = snappy.Link(pd).exterior()
        return M.high_precision() if hoch else M

    signal.signal(signal.SIGALRM, lambda *a: (_ for _ in ()).throw(FristAus('Frist')))
    signal.alarm(sek)
    abgeschnitten = False
    try:
        M = bauen()
        ok, cert, weg = eskalation_laut(M, bauen, print)
    except FristAus:
        print(f"    🔴 DEADLINE {sek}s EXCEEDED — run TRUNCATED, not fully escalated.")
        print(f"       This is NO result about the knot, but one about the time limit.")
        ok, cert, weg, abgeschnitten = False, None, f'ABGESCHNITTEN nach Frist {sek}s', True
    finally:
        signal.alarm(0)

    erg = dict(stem=stem, zertifiziert=bool(ok), weg=weg, kreuzungen=len(pd),
               abgeschnitten=bool(abgeschnitten), vol_kandidat=ls.get('volume'))
    if ok:
        try:
            erg['solution_type'] = str(M.solution_type())
            erg['vol_rigoros'] = str(M.volume(verified=True))
            erg['vol_float'] = float(M.volume())
        except Exception as e:
            erg['vol_fehler'] = f'{type(e).__name__}: {e}'
        print(f"    ✅ CERTIFIED via {weg}")
        print(f"       vol_rigoros = {erg.get('vol_rigoros')}")
        print(f"       vol_float   = {erg.get('vol_float')}   candidate was {ls.get('volume')}")
    else:
        print(f"    🔴 NOT certified ({weg}) — that does NOT mean 'not hyperbolic'")
    return erg

if __name__ == '__main__':
    # Latch BEFORE the run (hard rule 3b): the box freezes on memory shortage instead of killing.
    mb = int(os.environ.get('MB', '6000'))
    resource.setrlimit(resource.RLIMIT_AS, (mb*1024*1024, mb*1024*1024))
    out = []
    for pfad in sys.argv[1:]:
        stem = os.path.basename(os.path.dirname(pfad))
        e = zert(pfad, stem, sek=int(os.environ.get('SEK', '3600')))
        if e: out.append(e)
    if out:
        z = os.environ.get('JSON')
        if z: json.dump(out, open(z, 'w'), indent=1); print('written:', z)
