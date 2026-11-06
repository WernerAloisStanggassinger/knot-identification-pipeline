#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Bring ALGEBRA.fibered up to date from the HFK dictionary lying in the same record.

FINDING (sector T, 19.09.2026): some records carry
    pflichtliste.ALGEBRA.knot_floer_homology.fibered = true
but NO flat key pflichtliste.ALGEBRA.fibered. The corpus_db builder reads
the flat key, therefore these knots stand at null there and count in no
Conj-19 balance. "A result in a result file is not in the corpus" — here even:
in the record, but not in the field.

THE RULE THAT IS CHECKED HERE BEFORE IT IS APPLIED:
  ALGEBRA.fibered == ALGEBRA.knot_floer_homology.fibered
If there is even ONE record where both exist and differ, the
derivation is wrong and NOTHING is written.

Call:  --pruefen (measure only, default)  |  --schreiben (with backup)
"""
import site_config
import json, os, sys, shutil, datetime

DB = (site_config.PIPELINE_DIR + '/corpus_db.json')

def hole(r):
    """(has_flat, value_flat, has_hfk, value_hfk)"""
    a = (r.get("pflichtliste") or {}).get("ALGEBRA") or {}
    kf = a.get("knot_floer_homology")
    hat_f = "fibered" in a
    hat_h = isinstance(kf, dict) and "fibered" in kf
    return hat_f, a.get("fibered"), hat_h, (kf.get("fibered") if hat_h else None)

def sammeln():
    db = json.load(open(DB))
    zeilen = []
    for k in db["knots"]:
        p = k.get("pflicht_path")
        if not p or not os.path.exists(p):
            zeilen.append((k.get("stem"), p, None, None, None, None, "DATEI_FEHLT"))
            continue
        try:
            r = json.load(open(p))
        except Exception as e:
            zeilen.append((k.get("stem"), p, None, None, None, None, f"LESEFEHLER {e}"))
            continue
        hf, wf, hh, wh = hole(r)
        if hf and hh:      lage = "BEIDE"
        elif hf and not hh: lage = "NUR_FLACH"
        elif hh and not hf: lage = "NUR_HFK"      # <- these are to be brought up to date
        else:               lage = "KEINES"
        zeilen.append((k.get("stem"), p, hf, wf, hh, wh, lage))
    return zeilen

def pruefen(zeilen, laut=True):
    from collections import Counter
    c = Counter(z[6] for z in zeilen)
    beide = [z for z in zeilen if z[6] == "BEIDE"]
    konflikt = [z for z in beide if z[3] != z[5]]
    if laut:
        print("=" * 78)
        print("STATE OF THE FIELD over all records")
        print("=" * 78)
        for lage, n in sorted(c.items()):
            print(f"  {lage:<12} {n:>5}")
        print()
        print(f"COUNTER-CASE (the decisive one): records in which BOTH exist: {len(beide)}")
        print(f"  of these CONTRADICTORY (flat != HFK): {len(konflikt)}")
        if konflikt:
            print("  🔴 The derivation rule is NOT covered. Examples:")
            for z in konflikt[:10]:
                print(f"     {z[0]}: flat={z[3]} hfk={z[5]}")
        else:
            print("  ✅ not a single contradiction — the rule is calibrated on all double cases.")
        nur = [z for z in zeilen if z[6] == "NUR_HFK"]
        print(f"\nTO BE BROUGHT UP TO DATE (HFK only, flat field missing): {len(nur)}")
        from collections import Counter as C2
        print("  Values in the HFK dict:", dict(C2(z[5] for z in nur)))
    return len(konflikt) == 0, [z for z in zeilen if z[6] == "NUR_HFK"], konflikt

def selbsttest():
    """Counter-cases for the measurer itself — on artificial records, without disk."""
    print("=" * 78); print("SELFTEST of the measurer"); print("=" * 78)
    faelle = [
        ({"pflichtliste":{"ALGEBRA":{"fibered":True,"knot_floer_homology":{"fibered":True}}}},  (True,True,True,True)),
        ({"pflichtliste":{"ALGEBRA":{"knot_floer_homology":{"fibered":True}}}},                 (False,None,True,True)),
        ({"pflichtliste":{"ALGEBRA":{"fibered":False}}},                                        (True,False,False,None)),
        ({"pflichtliste":{"ALGEBRA":{}}},                                                       (False,None,False,None)),
        # falsified probe: flat True, HFK False -> must become visible as a contradiction
        ({"pflichtliste":{"ALGEBRA":{"fibered":True,"knot_floer_homology":{"fibered":False}}}}, (True,True,True,False)),
    ]
    ok = True
    for i,(r,soll) in enumerate(faelle,1):
        ist = hole(r)
        gut = ist == soll
        ok &= gut
        print(f"[{'OK ' if gut else 'RED'}] Case {i}: {ist}  (expected {soll})")
    # and: does pruefen() recognise the contradiction?
    z = [("k1","/x",True,True,True,False,"BEIDE")]
    gedeckt, _, konf = pruefen(z, laut=False)
    gut = (not gedeckt) and len(konf) == 1
    ok &= gut
    print(f"[{'OK ' if gut else 'RED'}] falsified probe is reported as a contradiction -> covered={gedeckt}")
    print(f"Selftest: {'ALL GREEN' if ok else '🔴 RED'}\n")
    return ok

def schreiben(nur_hfk):
    stempel = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    sicher = f"{site_config.ROOT}/claude_setup_v2/tmp_durable/fibered_nachzug_{stempel}"
    os.makedirs(sicher, exist_ok=True)
    protokoll = []
    for stem, p, hf, wf, hh, wh, lage in nur_hfk:
        shutil.copy2(p, os.path.join(sicher, os.path.basename(p)))
        r = json.load(open(p))
        a = r["pflichtliste"]["ALGEBRA"]
        vorher = dict(a)
        a["fibered"] = wh
        a["fibered_herkunft"] = ("knot_floer_homology.fibered in the same record, brought up to date "
                                 "2026-09-19 (sector T / Disc #811 vicinity): the flat key "
                                 "was missing, therefore the knot stood at null in corpus_db")
        tmp = p + ".tmp"
        with open(tmp, "w") as f:
            json.dump(r, f, indent=2, ensure_ascii=False)
        os.replace(tmp, p)
        # check: only this one key has been added
        nach = json.load(open(p))["pflichtliste"]["ALGEBRA"]
        neu = set(nach) - set(vorher)
        geaendert = {k for k in vorher if k in nach and nach[k] != vorher[k]}
        protokoll.append({"stem": stem, "pfad": p, "wert": wh,
                          "neue_schluessel": sorted(neu), "geaenderte_schluessel": sorted(geaendert)})
    with open(os.path.join(sicher, "PROTOKOLL.json"), "w") as f:
        json.dump(protokoll, f, indent=2, ensure_ascii=False)
    return sicher, protokoll

if __name__ == "__main__":
    if not selbsttest():
        print("⛔ Selftest red — abort."); sys.exit(1)
    zeilen = sammeln()
    gedeckt, nur_hfk, konflikt = pruefen(zeilen)
    if "--schreiben" in sys.argv:
        if not gedeckt:
            print("\n⛔ There are contradictions — NOTHING will be written."); sys.exit(1)
        sicher, prot = schreiben(nur_hfk)
        schlecht = [p for p in prot if p["neue_schluessel"] != ["fibered","fibered_herkunft"] or p["geaenderte_schluessel"]]
        print(f"\nWRITTEN: {len(prot)} records · backup+log: {sicher}")
        print(f"Records in which more than the two new keys arose: {len(schlecht)}"
              + ("  🔴" if schlecht else "  ✅"))
    else:
        print("\n(checked only — with --schreiben it is brought up to date)")
