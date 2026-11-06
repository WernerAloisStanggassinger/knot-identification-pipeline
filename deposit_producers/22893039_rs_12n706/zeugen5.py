#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""WITNESS SEARCH at c=5 for the twelve open cases (GEGENWIND §4).

§4: "prove r_s = 5: one witness in R^s_5 suffices (targeted search)."
A witness needs NO complete enumeration — it is an object one presents.
So these twelve are the EASIER case, not the harder one.

Procedure as for 12n706, only with c=5 diagrams:
  Delta sites (from R^s_5) -> real labels -> volume filter -> isometry_signature(verified).
The identity is checked VERIFIED (interval arithmetic), not just numerically (§2).
"""
import site_config
import itertools, json, os, sys, time
HIER = (site_config.ROOT + '/maschinen_ideen/RS4_12N706_20260919')
sys.path.insert(0, HIER)

ARBEITER = 6
FRIST = 180
GRENZE = int(sys.argv[2]) if len(sys.argv) > 2 else 3
VORHER = int(sys.argv[3]) if len(sys.argv) > 3 else 0   # bound already searched
ZIEL = sys.argv[1]
AUS = os.path.join(HIER, f"zeuge_{ZIEL}.json")

OFFEN = {
    "11a58":[2,-9,18,-23,18,-9,2], "11a164":[1,-7,20,-35,43,-35,20,-7,1],
    "11a326":[1,-6,19,-36,45,-36,19,-6,1], "12n49":[4,-20,33,-20,4],
    "12n504":[1,-6,15,-24,29,-24,15,-6,1], "12a3":[2,-14,40,-57,40,-14,2],
    "12a173":[1,-7,20,-35,43,-35,20,-7,1], "12a279":[2,-15,40,-55,40,-15,2],
    "12a646":[2,-9,21,-33,39,-33,21,-9,2], "12a715":[4,-18,38,-49,38,-18,4],
    "12a786":[2,-15,40,-55,40,-15,2], "12a1083":[2,-9,21,-33,39,-33,21,-9,2],
    # CALIBRATION CASE: r_s = 5 is PROVEN here (Table 3, Thm. 1.5). If the machinery finds there
    # no witness, it measures nothing — then all negative findings are worthless.
    "12n4":[1,-7,19,-27,19,-7,1],
}


def arbeite(auf):
    from snappy import Link as _L
    from spherogram import Link
    from pruef_aglet_flaeche import karten
    from band_rand import band_rand_sequenz
    from dt_bau import link_aus_sequenz
    wort, ueber, labels, kidx = auf
    w = [int(c) for c in wort]
    ueber = {int(k): v for k, v in ueber.items()}
    labels = {int(k): v for k, v in labels.items()}
    seq = band_rand_sequenz(w, ueber, karten(w)[kidx], labels, True, True)
    L, dt = link_aus_sequenz(seq)
    if L is None or len(L.link_components) != 1:
        return {"stufe": "kein_knoten"}
    L.simplify('global')
    ziel = Link('K' + ZIEL)
    ze = ziel.exterior()
    zv = float(ze.volume())
    try:
        E = L.exterior(); v = float(E.volume())
    except Exception:
        return {"stufe": "nicht_hyperbolisch"}
    if abs(v - zv) > 1e-6:
        return {"stufe": "volumen_verfehlt"}
    # VERIFIED (GEGENWIND §2): isometry_signature with interval arithmetic
    try:
        s1 = E.isometry_signature(verified=True); s2 = ze.isometry_signature(verified=True)
        echt = bool(s1 and s1 == s2)
    except Exception:
        echt = False
    return {"stufe": "geprueft", "wort": wort, "labels": labels, "c": len(L.crossings),
            "dt": list(dt), "iso_verified": echt}


def in_kind(auf):
    r, w = os.pipe()
    pid = os.fork()
    if pid == 0:
        try:
            os.close(r); os.write(w, json.dumps(arbeite(auf)).encode()); os.close(w)
        except BaseException:
            pass
        finally:
            os._exit(0)
    os.close(w)
    return pid, r


if __name__ == "__main__":
    from eiche_union import koeffs, koeffs_von_liste
    from enum_bogendiagramme import enumerieren
    from band_rand import alexander_von, strands
    from abschluss import aglet_ueber

    k = koeffs_von_liste(OFFEN[ZIEL])
    stellen = []
    for t in sorted(enumerieren(5, laut=False)):
        wort = "".join(map(str, t)); w = list(t)
        fest, vork = aglet_ueber(w)
        innere = sorted(c for c in vork if c not in fest)
        for wahl in itertools.product(*[[vork[c][0], vork[c][1]] for c in innere]):
            ueber = dict(fest); ueber.update(dict(zip(innere, wahl)))
            _, n_str = strands(w, ueber)
            frei = [j for j in range(n_str) if j not in (0, n_str - 1)]
            for bits in itertools.product((0, 1), repeat=len(frei)):
                lab = dict(zip(frei, bits))
                for s in alexander_von(w, ueber, lab, True, True):
                    if koeffs(s) == k:
                        stellen.append((wort, ueber, lab)); break
    print(f"TARGET {ZIEL}: {len(stellen)} Delta sites in R^s_5, |a| <= {GRENZE}"
          + (f" (only NEW compared to |a| <= {VORHER})" if VORHER else ""), flush=True)
    if not stellen:
        print("  no Delta site — this case is not reachable via R^s_5 at all"); sys.exit(0)

    ung = [a for a in range(-GRENZE, GRENZE + 1) if a % 2]
    ger = [a for a in range(-GRENZE, GRENZE + 1) if a % 2 == 0]
    aufgaben = []
    for wort, ueber, lab0 in stellen:
        frei = sorted(lab0)
        for kidx in (0, 1):
            for werte in itertools.product(*[(ung if lab0[j] % 2 else ger) for j in frei]):
                if VORHER and all(abs(v) <= VORHER for v in werte):
                    continue          # already checked at |a| <= VORHER, do not compute twice
                aufgaben.append((wort, {str(a): b for a, b in ueber.items()},
                                 {str(a): b for a, b in zip(frei, werte)}, kidx))
    print(f"{len(aufgaben)} tasks", flush=True)

    zustand = {"ziel": ZIEL, "grenze": GRENZE, "aufgaben": len(aufgaben),
               "zeugen": [], "stufen": {}, "erledigt": 0, "tot": 0}
    laufend, i = {}, 0
    while i < len(aufgaben) or laufend:
        while len(laufend) < ARBEITER and i < len(aufgaben):
            pid, r = in_kind(aufgaben[i]); laufend[pid] = (r, time.time()); i += 1
        for pid in list(laufend):
            r, ts = laufend[pid]
            fertig, _ = os.waitpid(pid, os.WNOHANG)
            if fertig == 0:
                if time.time() - ts > FRIST:
                    os.kill(pid, 9); os.waitpid(pid, 0); os.close(r)
                    del laufend[pid]; zustand["erledigt"] += 1; zustand["tot"] += 1
                continue
            daten = b""
            while True:
                st = os.read(r, 65536)
                if not st: break
                daten += st
            os.close(r); del laufend[pid]; zustand["erledigt"] += 1
            try:
                res = json.loads(daten.decode())
                zustand["stufen"][res.get("stufe","?")] = zustand["stufen"].get(res.get("stufe","?"),0)+1
                if res.get("iso_verified"):
                    zustand["zeugen"].append(res)
                    print(f"  WITNESS {res['wort']} labels={res['labels']} c={res['c']} "
                          f"(isometry_signature verified)", flush=True)
                    json.dump(zustand, open(AUS,"w"), indent=1)
                    if len(zustand["zeugen"]) >= 2:      # two suffice; one already proves
                        i = len(aufgaben)                 # start no new tasks
            except Exception:
                zustand["tot"] += 1
            if zustand["erledigt"] % 300 == 0:
                print(f"  [{zustand['erledigt']}/{len(aufgaben)}] witnesses {len(zustand['zeugen'])} "
                      f"· stages {zustand['stufen']}", flush=True)
        time.sleep(0.01)
    json.dump(zustand, open(AUS,"w"), indent=1)
    print(f"\nDONE {ZIEL}: {zustand['erledigt']} tasks · witnesses {len(zustand['zeugen'])} "
          f"· dead {zustand['tot']} · stages {zustand['stufen']}")
    if zustand["zeugen"]:
        z = zustand["zeugen"][0]
        print(f"  ==> r_s({ZIEL}) <= 5, proven by: {z['wort']} labels={z['labels']}")
        print(f"      DT {tuple(z['dt'])}")
    else:
        print(f"  ==> no witness in |a| <= {GRENZE}. NO proof of r_s >= 6.")
