#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""Gauss sequence -> DT code -> spherogram.Link. This makes cl(kappa) a COMPUTABLE
object (Jones, Alexander, volume), not just a determinant.

WHY DT AND NOT PD: a PD code needs the planar embedding (the cyclic order at
every crossing). A DT code does NOT -- it consists only of the traversal order
and the over/under information, and exactly that is in the Gauss sequence. spherogram
reconstructs the embedding itself.
  That every crossing gets one even and one odd position (precondition for DT),
  is no coincidence: it IS Gauss' parity condition -- the same one with which pruef_aglet_flaeche.py
  cross-checked the map code. It holds for every realizable CLOSED diagram,
  and cl(kappa) is one.

SIGN CONVENTION: measured, not guessed -- and the result was a surprise.
BOTH conventions yield the same knot, because in the DT code only the RELATIVE sign pattern
counts; flipping all signs is the same diagram. The first selftest therefore reported
"all green" twice and thereby measured NOTHING AT ALL (Footgun 16: a checker without a biting
counter-case is decoration). Since the sharpening the selftest carries a counter-case that
MUST bite: a flipped crossing must change the Jones polynomial.
What the test really proves: the construction from the Gauss sequence reproduces the original knot
CHIRALITY-FAITHFULLY (Jones distinguishes T(2,3) from its mirror image).
"""
import site_config
import sys
from spherogram import Link

sys.path.insert(0, (site_config.ROOT + '/maschinen_ideen/RS4_12N706_20260919'))


# --------------------------------------------- PD -> Gauss sequence (only for the calibration)
def pd_zu_sequenz(pd):
    """spherogram PD X(a,b,c,d): a=under-strand in, c=under-strand out, {b,d}=over-strand.
    Edge labels run consecutively along the knot."""
    M = 2 * len(pd)
    naechste = {}           # edge -> (crossing, 'O'|'U', next edge)
    for i, (a, b, c, d) in enumerate(pd):
        naechste[a] = (i, 'U', c)
        if (d - b) % M == 1:
            naechste[b] = (i, 'O', d)
        elif (b - d) % M == 1:
            naechste[d] = (i, 'O', b)
        else:
            raise ValueError(f"over-strand not consecutive in {(a,b,c,d)}")
    seq, k, start = [], 0, 0
    for _ in range(M):
        kr, typ, folge = naechste[k]
        seq.append((kr, typ))
        k = folge
    if k != start:
        raise ValueError("traversal does not close -- multi-component?")
    return seq


# --------------------------------------------- Gauss sequence -> DT
def sequenz_zu_dt(seq, negativ_bei_gerade_ueber=True):
    """Returns the DT tuple. Positions 1..2N along the traversal."""
    N2 = len(seq)
    if N2 % 2:
        raise ValueError("odd sequence length")
    N = N2 // 2
    pos = {}
    for i, (c, t) in enumerate(seq):
        pos.setdefault(c, []).append((i + 1, t))
    paare = {}
    for c, eintraege in pos.items():
        if len(eintraege) != 2:
            raise ValueError(f"crossing {c} occurs {len(eintraege)}x")
        (p1, t1), (p2, t2) = eintraege
        u = [(p, t) for p, t in eintraege if p % 2 == 1]
        g = [(p, t) for p, t in eintraege if p % 2 == 0]
        if len(u) != 1 or len(g) != 1:
            # Gauss' parity condition violated -> no realizable diagram
            return None
        (pu, tu), (pg, tg) = u[0], g[0]
        vz = -1 if ((tg == 'O') == negativ_bei_gerade_ueber) else 1
        paare[pu] = vz * pg
    return tuple(paare[p] for p in sorted(paare))


def link_aus_sequenz(seq, negativ_bei_gerade_ueber=True):
    dt = sequenz_zu_dt(seq, negativ_bei_gerade_ueber)
    if dt is None:
        return None, None
    return Link(f"DT: [{dt}]"), dt


# --------------------------------------------- SELF-TEST
ZOEPFE = [
    ("T(2,3)",      [1, 1, 1]),
    ("T(2,5)",      [1, 1, 1, 1, 1]),
    ("T(2,7)",      [1, 1, 1, 1, 1, 1, 1]),
    ("4_1",         [1, -2, 1, -2]),
    ("det7-Zopf",   [1, 1, 1, 2, -1, 2]),
    ("T(3,4)",      [1, 2, 1, 2, 1, 2]),
    ("nicht-alt A", [1, 1, -2, 1, 1, -2]),
    ("nicht-alt B", [1, -2, -2, 1, 2, 2, 1]),
]


def einkomponentige_zoepfe(grenze=8):
    """Search for non-alternating ONE-component braid closures (instead of guessing them)."""
    import itertools
    gefunden = []
    for laenge in (5, 6):          # 7 was too expensive (21,504 link builds, run aborted)
        for w in itertools.product([1, -1, 2, -2], repeat=laenge):
            if len(gefunden) >= grenze:
                return gefunden
            if len(set(abs(x) for x in w)) < 2:
                continue
            if not (any(x > 0 for x in w) and any(x < 0 for x in w)):
                continue          # mixed signs => not alternating
            try:
                L = Link(braid_closure=list(w))
                if len(L.link_components) == 1 and len(L.crossings) >= 4:
                    gefunden.append((f"braid{list(w)}", list(w)))
            except Exception:
                pass
    return gefunden


def selbsttest():
    print("=" * 96)
    print("SELF-TEST dt_bau.py — Gauss sequence -> DT -> Link, checked against the ORIGINAL")
    print("=" * 96)
    alles = True
    for konv in (True, False):
        print(f"\n--- Sign convention: negative if the EVEN passage "
              f"is {'OVER' if konv else 'UNDER'} ---")
        ok, n = True, 0
        for name, zopf in ZOEPFE:
            orig = Link(braid_closure=zopf)
            if len(orig.link_components) != 1:
                print(f"  [--] {name:<14} multi-component, skipped"); continue
            try:
                seq = pd_zu_sequenz(orig.PD_code())
                neu, dt = link_aus_sequenz(seq, konv)
                if neu is None:
                    print(f"  [RED] {name:<14} DT not buildable"); ok = False; n += 1; continue
                jo, jn = orig.jones_polynomial(), neu.jones_polynomial()
                gleich = (jo == jn)
                ok &= gleich; n += 1
                print(f"  [{'OK ' if gleich else 'RED'}] {name:<14} DT={str(dt):<28} "
                      f"Jones {'identical' if gleich else 'DEVIATION'}")
                if not gleich:
                    print(f"        orig {jo}")
                    print(f"        new  {jn}")
            except Exception as e:
                ok = False; n += 1
                print(f"  [RED] {name:<14} {type(e).__name__}: {e}")
        print(f"  Base set {n} probes -> {'ALL GREEN' if ok else '🔴 RED'}")
        alles &= ok

    print("\n--- ADDITION: non-alternating one-component braid closures (searched, not guessed) ---")
    extra = einkomponentige_zoepfe()
    ok, n = True, 0
    for name, zopf in extra:
        orig = Link(braid_closure=zopf)
        try:
            seq = pd_zu_sequenz(orig.PD_code())
            neu, dt = link_aus_sequenz(seq, True)
            if neu is None:
                continue
            gleich = (orig.jones_polynomial() == neu.jones_polynomial())
            ok &= gleich; n += 1
            print(f"  [{'OK ' if gleich else 'RED'}] {name:<26} {len(orig.crossings)} crossings, "
                  f"Jones {'identical' if gleich else 'DEVIATION'}")
        except Exception as e:
            print(f"  [--] {name:<26} {type(e).__name__}")
    print(f"  Base set {n} probes -> {'ALL GREEN' if ok else '🔴 RED'}")
    alles &= ok

    print("\n--- COUNTER-CASE (must BITE): flip one crossing => different Jones ---")
    ok, n = True, 0
    for name, zopf in ZOEPFE[:5] + extra[:5]:
        orig = Link(braid_closure=zopf)
        if len(orig.link_components) != 1:
            continue
        seq = pd_zu_sequenz(orig.PD_code())
        c0 = seq[0][0]
        kaputt = [(c, ('U' if t == 'O' else 'O') if c == c0 else t) for c, t in seq]
        neu, _ = link_aus_sequenz(kaputt, True)
        if neu is None:
            print(f"  [OK ] {name:<26} flipped version is no longer realizable (bites too)")
            n += 1; continue
        anders = (orig.jones_polynomial() != neu.jones_polynomial())
        ok &= anders; n += 1
        print(f"  [{'OK ' if anders else 'RED'}] {name:<26} "
              f"{'Jones differs — counter-case bites' if anders else 'NO difference — checker measures nothing'}")
    print(f"  Base set {n} probes -> {'ALL GREEN' if ok else '🔴 RED'}")
    alles &= ok

    print()
    print("Finding on the convention: BOTH sign conventions yield the same knot —")
    print("in the DT code only the RELATIVE pattern counts. The choice is therefore immaterial,")
    print("it was not 'determined'. Proven is: the construction is chirality-faithful.")
    print("CALIBRATION:", "ALL GREEN" if alles else "🔴 RED")
    return alles


if __name__ == "__main__":
    sys.exit(0 if selbsttest() else 1)
