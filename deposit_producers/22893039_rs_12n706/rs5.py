#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""COMPUTE R^s_5 — the set that the paper explicitly leaves open.

  Akash et al., arXiv:2606.02390v1, S. 3: "If r_s(K) = 5, the set R^s_5 is too large to
  compute via our methods."

Our machinery can do it: the enumerator delivers 1/2/8/26 reduced arc diagrams for
c=2/3/4/5 (the first three calibrated against Lemma 3.1, 3.2, Prop. 3.3), and the union constructor
reproduces R^s_2, R^s_3 and R^s_4 EXACTLY (27 of 27, Table 1).

At c=5: 2 aglet crossings fixed + 3 inner free (2^3) x 4 free labels mod 2 (2^4) = 128
per diagram, so 26*128 = 3328 computations.

TWO INDEPENDENT CONTROLS, both from the source:
 (1) POSITIVE: every knot from Table 3 with PROVEN r_s = 5 must have its Delta in R^s_5.
     If one is missing, the enumeration is incomplete.
 (2) NEGATIVE: Theorem 1.6 says r_s(K)=5 ==> det(K) <= 169. So NO polynomial in R^s_5 may
     carry a determinant > 169. That is a prediction of the paper about our set.
"""
import site_config
import itertools, json, os, sys, time
HIER = (site_config.ROOT + '/maschinen_ideen/RS4_12N706_20260919')
sys.path.insert(0, HIER)

ARBEITER = 6
FRIST = 90
AUS = os.path.join(HIER, "rs5_ergebnis.json")

# Table 3: knots with PROVEN r_s = 5 (column "lower" = Thm. 1.5)
RS5_BEWIESEN = {
    "12n4":[1,-7,19,-27,19,-7,1], "12n43":[1,-5,10,-15,19,-15,10,-5,1],
    "12n106":[1,-4,10,-16,19,-16,10,-4,1], "12n170":[6,-20,29,-20,6],
    "12n380":[2,-9,18,-23,18,-9,2], "12n399":[1,-7,19,-27,19,-7,1],
    "12n420":[1,-7,19,-27,19,-7,1], "12n440":[2,-9,18,-23,18,-9,2],
    "12n553":[4,-20,33,-20,4], "12n556":[4,-20,33,-20,4],
    "12n636":[1,-7,19,-27,19,-7,1], "12n657":[1,-4,9,-16,21,-16,9,-4,1],
    "12n676":[2,-2,-4,9,-4,-2,2], "12n702":[2,-12,28,-37,28,-12,2],
    "12n782":[2,-8,18,-25,18,-8,2], "12n802":[1,-5,14,-25,31,-25,14,-5,1],
    "12n876":[2,-8,18,-25,18,-8,2],
}
# 🔴 CORRECTION 20.09.: R^s_5 decides ONLY IN ONE DIRECTION.
#   Delta NOT IN R^s_5  ==>  r_s != 5      (valid: contraposition of the definition)
#   Delta     IN R^s_5  ==>  NOTHING AT ALL    (Delta is not a complete invariant; there is
#                                            just SOME knot with r_s=5 and this Delta)
# The first version of this script printed "YES ==> r_s = 5" — a fallacy.
# A case can only be decided positively by an EXPLICIT disc (as with 12n706).
# The open cases with det <= 169, to which the obstruction is applicable at all:
OFFEN = {
    "11a58":[2,-9,18,-23,18,-9,2], "11a164":[1,-7,20,-35,43,-35,20,-7,1],
    "11a326":[1,-6,19,-36,45,-36,19,-6,1], "12n49":[4,-20,33,-20,4],
    "12n504":[1,-6,15,-24,29,-24,15,-6,1], "12a3":[2,-14,40,-57,40,-14,2],
    "12a173":[1,-7,20,-35,43,-35,20,-7,1], "12a279":[2,-15,40,-55,40,-15,2],
    "12a646":[2,-9,21,-33,39,-33,21,-9,2], "12a715":[4,-18,38,-49,38,-18,4],
    "12a786":[2,-15,40,-55,40,-15,2], "12a1083":[2,-9,21,-33,39,-33,21,-9,2],
    "12a1202":[9,-42,67,-42,9],
}


def arbeite(auf):
    """IN THE CHILD: one (diagram, crossing info, labels) -> Alexander coefficients."""
    from eiche_union import koeffs
    from band_rand import alexander_von
    wort, ueber, lab = auf
    w = [int(c) for c in wort]
    ueber = {int(k): v for k, v in ueber.items()}
    lab = {int(k): v for k, v in lab.items()}
    aus = []
    for s in alexander_von(w, ueber, lab, True, True):
        k = koeffs(s)
        aus.append(list(k) if k else None)
    return aus


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
    from enum_bogendiagramme import enumerieren
    from abschluss import aglet_ueber
    from band_rand import strands

    diagramme = sorted(enumerieren(5, laut=False))
    print(f"c=5: {len(diagramme)} reduced arc diagrams", flush=True)
    aufgaben = []
    for t in diagramme:
        wort = "".join(map(str, t))
        w = list(t)
        fest, vork = aglet_ueber(w)
        innere = sorted(c for c in vork if c not in fest)
        for wahl in itertools.product(*[[vork[c][0], vork[c][1]] for c in innere]):
            ueber = dict(fest); ueber.update(dict(zip(innere, wahl)))
            _, n_str = strands(w, ueber)
            frei = [j for j in range(n_str) if j not in (0, n_str - 1)]
            for bits in itertools.product((0, 1), repeat=len(frei)):
                aufgaben.append((wort, {str(a): b for a, b in ueber.items()},
                                 {str(a): b for a, b in zip(frei, bits)}))
    print(f"{len(aufgaben)} tasks (expected 26 x 8 x 16 = 3328)", flush=True)

    menge, laufend, i, erledigt, tot = set(), {}, 0, 0, 0
    t0 = time.time()
    while i < len(aufgaben) or laufend:
        while len(laufend) < ARBEITER and i < len(aufgaben):
            pid, r = in_kind(aufgaben[i]); laufend[pid] = (r, time.time()); i += 1
        for pid in list(laufend):
            r, ts = laufend[pid]
            fertig, _ = os.waitpid(pid, os.WNOHANG)
            if fertig == 0:
                if time.time() - ts > FRIST:
                    os.kill(pid, 9); os.waitpid(pid, 0); os.close(r)
                    del laufend[pid]; erledigt += 1; tot += 1
                continue
            daten = b""
            while True:
                st = os.read(r, 65536)
                if not st: break
                daten += st
            os.close(r); del laufend[pid]; erledigt += 1
            try:
                for k in json.loads(daten.decode()):
                    if k: menge.add(tuple(k))
            except Exception:
                tot += 1
            if erledigt % 400 == 0:
                print(f"  [{erledigt}/{len(aufgaben)}] polynomials {len(menge)} · dead {tot} "
                      f"· {time.time()-t0:.0f}s", flush=True)
        time.sleep(0.01)

    det = lambda v: abs(sum(c * (-1)**j for j, c in enumerate(v)))
    print(f"\nDONE: {erledigt} tasks, {tot} without result, {len(menge)} polynomials in {time.time()-t0:.0f}s")
    json.dump({"polynome": sorted(map(list, menge)), "aufgaben": erledigt, "tot": tot},
              open(AUS, "w"), indent=1)
    dets = {}
    for v in menge: dets[det(list(v))] = dets.get(det(list(v)), 0) + 1
    print(f"Determinants in R^s_5: {dict(sorted(dets.items()))}")
    print()
    def norm(v):
        k = [v, v[::-1], [-x for x in v], [-x for x in v][::-1]]
        return min(tuple(x) for x in k)
    print("CONTROL 1 (positive): Table-3 knots with PROVEN r_s = 5 must be in it")
    fehlt = [n for n, d in RS5_BEWIESEN.items() if norm(d) not in menge]
    print(f"  {len(RS5_BEWIESEN)-len(fehlt)}/{len(RS5_BEWIESEN)} found"
          + (f"   🔴 MISSING: {fehlt}" if fehlt else "   ALL IN"))
    print()
    print("CONTROL 2 (negative): Theorem 1.6 — no polynomial in R^s_5 may have det > 169")
    zuviel = [d for d in dets if d > 169]
    print(f"  largest determinant {max(dets)}   " +
          (f"🔴 ABOVE 169: {sorted(zuviel)}" if zuviel else "<= 169 ✅ Theorem 1.6 confirmed"))
    print()
    print("THE OPEN CASES — R^s_5 acts ONLY as an obstruction (one direction!):")
    entschieden = []
    for n, d in sorted(OFFEN.items()):
        drin = norm(d) in menge
        if drin:
            print(f"  {n:9s} det {det(d):4d}  Delta in R^s_5  ==> NO statement "
                  f"(Delta is not a complete invariant) — stays open")
        else:
            entschieden.append(n)
            print(f"  {n:9s} det {det(d):4d}  Delta NOT in R^s_5  ==> r_s != 5  "
                  f"==> DECIDED with the upper bound from the table")
    print()
    print(f"DECIDED by the obstruction: {entschieden if entschieden else '(none)'}")
    print("All others need an EXPLICIT disc, as 12n706 got one.")
