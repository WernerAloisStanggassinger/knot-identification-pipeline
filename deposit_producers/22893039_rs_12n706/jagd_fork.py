#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Decide 12n23 — with ISOLATION AT PROCESS LEVEL.

WHY FORK: the first attempt (offene_restliche.py, |a|<=5) vanished without traceback in the middle of the
run. No OOM (12 GB free, journal empty). Cause is the SnapPea C code: a hard crash
comes as a SIGNAL, not as an exception, and is caught by no try/except (Footgun 13).
Cure: every task in a THROWAWAY CHILD (os.fork), result via pipe, deadline in the parent via
an OWN CLOCK with os.waitpid(..., WNOHANG) -- never blocking (Footgun 18/24) and
explicitly WITHOUT signal.alarm (Footgun 20: SIGALRM without a handler kills the watchdog itself).

Result is saved atomically after EVERY hit (tmp + os.replace, Footgun 28).
"""
import site_config
import itertools, json, os, sys, time

HIER = (site_config.ROOT + '/maschinen_ideen/RS4_12N706_20260919')
sys.path.insert(0, HIER)

ZIEL_NAME = sys.argv[1] if len(sys.argv) > 1 else "12n23"
GRENZE = int(sys.argv[2]) if len(sys.argv) > 2 else 3
ARBEITER = 6                      # Conj-1.6 runs alongside (Footgun 14): not all 24
FRIST = 120                       # seconds per task
MAX_C = 40                        # now only an emergency brake; the sharp filter is the VOLUME
AUS = os.path.join(HIER, f"ergebnis_{ZIEL_NAME}.json")

DELTA = {"12n23": [2, -5, 2],
         "12n706": [1, -4, 6, -8, 11, -8, 6, -4, 1],   # calibration case: the hit is known
         "12n49": [4, -20, 33, -20, 4],
         "12n504": [1, -6, 15, -24, 29, -24, 15, -6, 1]}


def arbeite(aufgabe):
    """Runs IN THE CHILD. Builds everything itself so that the parent stays untouched."""
    from snappy import Link as _L            # noqa: F401  (makes exterior() usable)
    from spherogram import Link
    from pruef_aglet_flaeche import karten
    from band_rand import band_rand_sequenz
    from dt_bau import link_aus_sequenz
    wort, ueber, labels, kidx = aufgabe
    w = [int(c) for c in wort]
    ueber = {int(k): v for k, v in ueber.items()}
    labels = {int(k): v for k, v in labels.items()}
    seq = band_rand_sequenz(w, ueber, karten(w)[kidx], labels, True, True)
    L, dt = link_aus_sequenz(seq)
    if L is None:
        return {"stufe": "kein_dt"}
    if len(L.link_components) != 1:
        return {"stufe": "mehrkomponentig"}
    L.simplify('global')
    nc = len(L.crossings)
    if nc > MAX_C:
        return {"stufe": "zu_gross", "c": nc}      # only an emergency brake against outliers
    ziel = Link('K' + ZIEL_NAME)
    # VOLUME FIRST: exact invariant, much sharper and cheaper than the old
    # size filter "c <= 16 after simplify" (that was heuristic, because simplify does not
    # guarantee to find the minimum, and discarded ~80% unchecked).
    ze = ziel.exterior()
    zv = float(ze.volume())
    try:
        v = float(L.exterior().volume())
    except Exception:
        return {"stufe": "nicht_hyperbolisch", "c": nc}
    if abs(v - zv) > 1e-6:
        return {"stufe": "volumen_verfehlt", "c": nc}
    if L.jones_polynomial() not in (ziel.jones_polynomial(), ziel.mirror().jones_polynomial()):
        return {"stufe": "jones_verfehlt", "c": nc}
    iso = bool(L.exterior().is_isometric_to(ze))
    return {"stufe": "geprueft", "wort": wort, "labels": labels,
            "c": nc, "dt": list(dt), "iso": iso}


def in_kind(aufgabe):
    r, w = os.pipe()
    pid = os.fork()
    if pid == 0:
        try:
            os.close(r)
            res = arbeite(aufgabe)
            os.write(w, json.dumps(res).encode())
            os.close(w)
        except BaseException:
            pass
        finally:
            os._exit(0)
    os.close(w)
    return pid, r


def ernte(pid, r, start):
    """Poll non-blocking; terminate by SIGKILL after the DEADLINE."""
    while True:
        fertig, status = os.waitpid(pid, os.WNOHANG)
        if fertig == pid:
            break
        if time.time() - start > FRIST:
            try:
                os.kill(pid, 9)
            except ProcessLookupError:
                pass
            os.waitpid(pid, 0)
            os.close(r)
            return "FRIST"
        time.sleep(0.02)
    daten = b""
    while True:
        st = os.read(r, 65536)
        if not st:
            break
        daten += st
    os.close(r)
    if os.WIFSIGNALED(status):
        return f"SIGNAL{os.WTERMSIG(status)}"
    try:
        return json.loads(daten.decode())
    except Exception:
        return "LEER"


def sichern(zustand):
    tmp = AUS + ".tmp"
    with open(tmp, "w") as f:
        json.dump(zustand, f, indent=1)
    os.replace(tmp, AUS)


if __name__ == "__main__":
    from eiche_union import koeffs, koeffs_von_liste, faelle
    from band_rand import alexander_von
    from pruef_aglet_flaeche import WOERTER

    k = koeffs_von_liste(DELTA[ZIEL_NAME])
    stellen = []
    for wort in WOERTER[4]:
        for w, ueber, lab, innere, n_str in faelle(wort):
            for s in alexander_von(w, ueber, lab, True, True):
                if koeffs(s) == k:
                    stellen.append((wort, ueber, lab))
    print(f"TARGET {ZIEL_NAME}: {len(stellen)} Delta sites, |a| <= {GRENZE}, "
          f"{ARBEITER} workers, deadline {FRIST}s", flush=True)

    ung = [a for a in range(-GRENZE, GRENZE + 1) if a % 2]
    ger = [a for a in range(-GRENZE, GRENZE + 1) if a % 2 == 0]
    aufgaben = []
    for wort, ueber, lab0 in stellen:
        frei = sorted(lab0)
        for kidx in (0, 1):
            for werte in itertools.product(*[(ung if lab0[j] % 2 else ger) for j in frei]):
                aufgaben.append((wort, {str(a): b for a, b in ueber.items()},
                                 {str(a): b for a, b in zip(frei, werte)}, kidx))
    print(f"{len(aufgaben)} tasks", flush=True)

    zustand = {"ziel": ZIEL_NAME, "grenze": GRENZE, "aufgaben": len(aufgaben),
               "treffer": [], "signal_tot": 0, "frist_tot": 0, "erledigt": 0,
               "stufen": {}}      # <- WHERE the tasks get stuck (Footgun 42)
    laufend, i = {}, 0
    while i < len(aufgaben) or laufend:
        while len(laufend) < ARBEITER and i < len(aufgaben):
            pid, r = in_kind(aufgaben[i])
            laufend[pid] = (r, time.time(), aufgaben[i]); i += 1
        for pid in list(laufend):
            r, t0, auf = laufend[pid]
            fertig, _ = os.waitpid(pid, os.WNOHANG)
            if fertig == 0 and time.time() - t0 <= FRIST:
                continue
            res = ernte(pid, r, t0) if fertig == 0 else None
            if fertig == pid:
                daten = b""
                while True:
                    st = os.read(r, 65536)
                    if not st: break
                    daten += st
                os.close(r)
                try: res = json.loads(daten.decode())
                except Exception: res = "SIGNAL/LEER"
            del laufend[pid]
            zustand["erledigt"] += 1
            if isinstance(res, str):
                zustand["signal_tot" if res.startswith(("SIGNAL", "LEER")) else "frist_tot"] += 1
            elif isinstance(res, dict):
                zustand["stufen"][res.get("stufe", "?")] = \
                    zustand["stufen"].get(res.get("stufe", "?"), 0) + 1
            if isinstance(res, dict) and res.get("iso"):
                zustand["treffer"].append(res)
                print(f"  ✅ HIT {res['wort']} labels={res['labels']} c={res['c']}", flush=True)
                sichern(zustand)
            if zustand["erledigt"] % 200 == 0:
                print(f"  [{zustand['erledigt']}/{len(aufgaben)}] "
                      f"hits {len(zustand['treffer'])} · hard {zustand['signal_tot']} · "
                      f"deadline {zustand['frist_tot']} · stages {zustand['stufen']}", flush=True)
                sichern(zustand)
        time.sleep(0.02)
    sichern(zustand)
    print(f"\nDONE: {zustand['erledigt']} tasks · hits {len(zustand['treffer'])} · "
          f"died hard {zustand['signal_tot']} · deadline {zustand['frist_tot']}")
    print(f"WHERE the tasks ended up: {zustand['stufen']}")
    g = (zustand['stufen'].get('geprueft', 0) + zustand['stufen'].get('jones_verfehlt', 0)
         + zustand['stufen'].get('volumen_verfehlt', 0))
    print(f"  ==> {g} tasks got as far as the volume/Jones comparison. If that is 0, the")
    print(f"      search checked NOTHING and 'no hit' means nothing (Footgun 42).")
    for t in zustand["treffer"][:5]:
        print(f"   {t['wort']} labels={t['labels']} c={t['c']}  DT {tuple(t['dt'])}")
    if not zustand["treffer"]:
        print(f"   NO hit. That is NO proof of r_s >= 5 — only: in this search space")
        print(f"   (|a| <= {GRENZE}, c <= {MAX_C}) not present.")
