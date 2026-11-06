#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""
HARVEST: lower bounds for the unknotting number u(K) over the whole corpus.

Occasion: Myf Evans' *untangling number* for 3-periodic tangles (Potsdam GMI 2026).
The finite counterpart is the classical unknotting number u(K).
Work order: konversationen/potsdam_gmi2026/vom_laptop_20260810/MESSEN_ENTKNOTUNGSZAHL_20260809.md

WHAT HAPPENS HERE: only HARVESTING + COMPUTING FROM EXISTING FIELDS. No new invariant,
no Khovanov run, no writing into corpus_db (record = truth, DB would be generated).

The four bounds, all textbook:
    u(K) >= |sigma(K)| / 2            (Murasugi)
    u(K) >= |s(K)| / 2                (Rasmussen; via g_4 >= |s|/2 and u >= g_4)
    u(K) >= |tau(K)|                  (Ozsvath-Szabo)
    u(K) >= g_4(K)                    (smooth 4-genus, lower bound from the record)

Counter-checks that BITE (footgun 16 — a test without a counter-case is decoration):
    P1  thin  ==>  s == -sigma        (Lee/Rasmussen for thin knots)
    P2  u_lower <= u_upper            (diagram unknotting number must never undercut the algebra)
    P3  |s|/2 <= g_4_lower            (the g4 bound in the record must already contain s)
    P4  ground truth: trefoil u=1, figure-eight u=1
    --selbsttest deliberately flips the sign of s: P1 MUST then turn red.
"""
import site_config
import json
import os
import sys
from collections import Counter

WURZEL = site_config.ROOT
DB = os.path.join(WURZEL, "freeze_manifest/pipeline_scripts/corpus_db.json")
HIER = os.path.dirname(os.path.abspath(__file__))
ZIEL = os.path.join(HIER, "ernte_u_schranken.json")


def hole(d, *pfad, default=None):
    """Read nested — footgun 29: never guess on the wrong level."""
    cur = d
    for p in pfad:
        if not isinstance(cur, dict) or p not in cur:
            return default
        cur = cur[p]
    return cur


def ganzzahl(x):
    if x is None:
        return None
    try:
        return int(x)
    except (TypeError, ValueError):
        return None


def ernten(sabotage_s=False):
    db = json.load(open(DB))
    zeilen = []
    fehlend = Counter()
    for r in db["knots"]:
        pfad = r["pflicht_path"]
        if not os.path.exists(pfad):
            fehlend["record_datei_fehlt"] += 1
            continue
        d = json.load(open(pfad))
        A = hole(d, "pflichtliste", "ALGEBRA", default={})
        O = hole(d, "optional_invariants", default={})

        sig = ganzzahl(A.get("signature"))
        s = ganzzahl(A.get("rasmussen_s"))
        if s is None:
            s = ganzzahl(O.get("rasmussen_s"))
        if sabotage_s and s is not None:
            s = -s  # --selbsttest: P1 must react to this
        tau = ganzzahl(A.get("tau_invariant"))
        g4l = ganzzahl(A.get("smooth_4_genus_lower_bound"))
        thin = O.get("khovanov_thin")
        if thin is None:
            thin = A.get("khovanov_thin")
        width = ganzzahl(O.get("khovanov_width"))

        u_diag = ganzzahl(A.get("unknotting_number"))
        u_sl = ganzzahl(A.get("unknotting_number_strict_lower"))
        u_su = ganzzahl(A.get("unknotting_number_strict_upper"))
        u_ub = ganzzahl(A.get("unknotting_number_upper_bound"))
        witness = A.get("unknotting_witness_subset")
        # 🔴 ONLY witness-backed upper bounds count. `unknotting_number` alone is NOT an
        # upper bound: for `knots-of-love` it says 6 there, while the record explicitly
        # "u > 4 on this diagram" and carries neither upper_bound nor witness — the 6
        # is there an assignment, not a proof. (Exactly one record in the corpus, checked 11.08.)
        u_obere = u_ub if u_ub is not None else u_su
        if u_obere is not None and not witness:
            u_obere = None
        u_obere_quelle = ("upper_bound+witness" if u_obere is not None else
                          ("VERWORFEN: unknotting_number ohne Zeugen" if u_diag is not None else None))

        kand = {}
        if sig is not None:
            kand["sigma"] = abs(sig) // 2
        if s is not None:
            kand["s"] = abs(s) // 2
        if tau is not None:
            kand["tau"] = abs(tau)
        if g4l is not None:
            kand["g4"] = g4l
        u_lower = max(kand.values()) if kand else None
        sieger = sorted([k for k, v in kand.items() if v == u_lower]) if kand else []

        if sig is None:
            fehlend["sigma"] += 1
        if s is None:
            fehlend["rasmussen_s"] += 1
        if tau is None:
            fehlend["tau"] += 1
        if thin is None:
            fehlend["khovanov_thin"] += 1

        zeilen.append(dict(
            stem=r["stem"], folder=r.get("folder"), category=r.get("category"),
            c=r.get("c"), det=r.get("det"),
            distinct_knot_id=r.get("distinct_knot_id"),
            ist_vertreter=bool(r.get("is_distinct_representative")),
            n_realisierungen=r.get("n_realisierungen"),
            knotinfo_id=r.get("knotinfo_id"),
            sig_record=sig, sig_db=ganzzahl(r.get("sig")),
            s=s, tau=tau, g4_lower=g4l,
            thin=thin, kh_width=width,
            b_sigma=kand.get("sigma"), b_s=kand.get("s"),
            b_tau=kand.get("tau"), b_g4=kand.get("g4"),
            u_lower=u_lower, sieger=sieger,
            u_obere=u_obere, u_obere_quelle=u_obere_quelle, u_diag=u_diag,
            u_strict_lower_record=u_sl,
            pflicht_path=pfad,
        ))
    return zeilen, fehlend


def pruefen(zeilen):
    """The four counter-checks. Returns (report, n_red)."""
    b = []
    rot = 0

    # P1: thin ==> s == -sigma
    p1 = [z for z in zeilen if z["thin"] is True and z["s"] is not None and z["sig_record"] is not None]
    p1_bad = [z for z in p1 if z["s"] != -z["sig_record"]]
    b.append(f"P1  thin ==> s == -sigma : {len(p1)-len(p1_bad)}/{len(p1)} ok, {len(p1_bad)} Verletzung(en)")
    rot += len(p1_bad)
    for z in p1_bad[:5]:
        b.append(f"    ROT {z['stem']}: sigma={z['sig_record']} s={z['s']}")

    # P2: u_lower <= u_obere
    p2 = [z for z in zeilen if z["u_lower"] is not None and z["u_obere"] is not None]
    p2_bad = [z for z in p2 if z["u_lower"] > z["u_obere"]]
    b.append(f"P2  u_lower <= u_obere   : {len(p2)-len(p2_bad)}/{len(p2)} ok, {len(p2_bad)} Verletzung(en)")
    rot += len(p2_bad)
    for z in p2_bad[:8]:
        b.append(f"    ROT {z['stem']}: u_lower={z['u_lower']} > u_obere={z['u_obere']} (c={z['c']})")

    # P3: |s|/2 <= g4_lower
    p3 = [z for z in zeilen if z["b_s"] is not None and z["g4_lower"] is not None]
    p3_bad = [z for z in p3 if z["b_s"] > z["g4_lower"]]
    b.append(f"P3  |s|/2 <= g4_lower    : {len(p3)-len(p3_bad)}/{len(p3)} ok, {len(p3_bad)} Verletzung(en)")
    rot += len(p3_bad)
    for z in p3_bad[:5]:
        b.append(f"    ROT {z['stem']}: |s|/2={z['b_s']} > g4_lower={z['g4_lower']}")

    # P3b: sigma from record == sigma from DB (consistency record<->DB)
    p3b = [z for z in zeilen if z["sig_record"] is not None and z["sig_db"] is not None]
    p3b_bad = [z for z in p3b if z["sig_record"] != z["sig_db"]]
    b.append(f"P3b sigma Record == DB   : {len(p3b)-len(p3b_bad)}/{len(p3b)} ok, {len(p3b_bad)} Abweichung(en)")
    rot += len(p3b_bad)
    for z in p3b_bad[:5]:
        b.append(f"    ROT {z['stem']}: Record {z['sig_record']} != DB {z['sig_db']}")

    # P4: ground truth. Torus knots: u(T(p,q)) = (p-1)(q-1)/2 (Kronheimer-Mrowka, THEOREM —
    # computed here, not copied from a table). Plus the figure-eight (u=1).
    torus = {"K3_1": (2, 3), "K5_1": (2, 5), "K7_1": (2, 7), "K9_1": (2, 9),
             "K8_19": (3, 4), "K10_124": (3, 5)}
    grund = {name: (p - 1) * (q - 1) // 2 for name, (p, q) in torus.items()}
    grund["K4_1"] = 1
    for name in list(grund):           # mirrors carry the same unknotting number
        grund[name + "m"] = grund[name]
    treffer, p4_bad = 0, []
    for z in zeilen:
        ki = (z["knotinfo_id"] or "").strip()
        if ki in grund and z["u_lower"] is not None:
            treffer += 1
            if z["u_lower"] > grund[ki]:
                p4_bad.append((z["stem"], ki, z["u_lower"], grund[ki]))
    scharf = sum(1 for z in zeilen
                 if (z["knotinfo_id"] or "").strip() in grund and z["u_lower"] is not None
                 and z["u_lower"] == grund[(z["knotinfo_id"] or "").strip()])
    b.append(f"P4  Grundwahrheit u_lower <= u_wahr : {treffer-len(p4_bad)}/{treffer} ok, "
             f"{len(p4_bad)} Verletzung(en); davon SCHARF (u_lower == u_wahr): {scharf}")
    rot += len(p4_bad)
    for t in p4_bad[:5]:
        b.append(f"    ROT {t[0]} ({t[1]}): u_lower={t[2]} > u_wahr={t[3]}")

    return b, rot


def main():
    selbsttest = "--selbsttest" in sys.argv
    zeilen, fehlend = ernten()
    bericht, rot = pruefen(zeilen)

    print(f"Records harvested: {len(zeilen)}")
    print("missing fields:", dict(fehlend))
    print()
    for z in bericht:
        print(z)
    print()

    if selbsttest:
        z2, _ = ernten(sabotage_s=True)
        b2, rot2 = pruefen(z2)
        print("--- SELF-TEST: sign of s deliberately flipped ---")
        for z in b2[:4]:
            print(z)
        beisst = rot2 > rot
        print(f"--- Self-test {'GREEN (the test bites)' if beisst else 'RED (the test is decoration!)'}: "
              f"red {rot} -> {rot2} ---")
        if not beisst:
            sys.exit(2)

    json.dump({"gebaut_am": "2026-08-11", "n_records": len(zeilen),
               "fehlende_felder": dict(fehlend), "proben": bericht,
               "n_rot": rot, "zeilen": zeilen},
              open(ZIEL + ".tmp", "w"), indent=1)
    os.replace(ZIEL + ".tmp", ZIEL)  # footgun 28: atomic
    print(f"\nwritten: {ZIEL}  ({os.path.getsize(ZIEL)} bytes)")
    print(f"RED total: {rot}")


if __name__ == "__main__":
    main()
