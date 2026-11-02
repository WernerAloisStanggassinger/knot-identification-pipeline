# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Flype-orbit chirality test for prime alternating knots (01.11.2026).
# Menasco–Thistlethwaite 1993: any two reduced alternating diagrams of the same (prime alternating) knot on S² are related by flypes.
# ⟹ K amphichiral ⟺ the reflected diagram mD lies in the flype orbit of D (up to isotopies of S² and the 3-D rotation, Regina sig
#   with allowReflection=False, allowReversal=True, allowRotation=True). mD not in the COMPLETE orbit ⟹ K chiral (combinatorial proof,
#   independent of SnapPy's isometry signature).
# Flypes as rewrites of the signed Gauss word (our A343358 deposit, DOI 10.5281/zenodo.23855791, DEFINITION_CHECK_20261030.md §5, cases A/B/C): T = crossings whose letters fill two
#   blocks U, V; c adjacent. A: c1 U X c2 V Y -> U c1 X V c2 Y · B: c1 U X V c2 Y -> U c1 X c2 V Y · C: c1 U c2 X V Y -> U X c1 V c2 Y.
#   T is rotated by π: letters of T keep their sign and swap O/U; c keeps its sign, its O/U is fixed by alternation. Applied to the word and
#   to its reverse (covers the inverse moves). OVER-GENERATION IS SAFE (extra members can only make the test inconclusive), MISSING flypes
#   are not — hence the counter-cases: every amphichiral alternating knot among them must find its mirror in the orbit.
#   Kept only if Regina says classical (planar on S²), alternating, same crossing number, same Jones polynomial.
import sys, re, json, bz2, csv
import regina, snappy
dt = lambda s: [(ord(ch.lower()) - 96) * 2 * (1 if ch.islower() else -1) for ch in s]
tok = re.compile(r'([OU])(\d+)([+-])')
def word(L): return [(m.group(1), int(m.group(2)), m.group(3)) for m in tok.finditer(L.signedGauss())]
def link(w): return regina.Link.fromSignedGauss(''.join(f'{o}{l}{s}' for o, l, s in w))
def flipT(blk): return [('U' if o == 'O' else 'O', l, s) for o, l, s in blk]
def closed(U, V):
    from collections import Counter
    c = Counter(l for _, l, _ in U + V); return all(v == 2 for v in c.values())
def moves(w):
    n2 = len(w)
    for ww in (w, w[::-1]):
        labels = {l for _, l, _ in ww}
        for cl in labels:
            occ = [i for i, t in enumerate(ww) if t[1] == cl]
            for start in occ:
                r = ww[start:] + ww[:start]; c1 = r[0]; p = next(i for i in range(1, n2) if r[i][1] == cl); c2 = r[p]
                # A: c1 U X c2 V Y
                for i in range(2, p + 1):
                    U = r[1:i]; X = r[i:p]
                    for j in range(p + 2, n2 + 1):
                        V = r[p + 1:j]; Y = r[j:]
                        if closed(U, V) and not ({l for _, l, _ in U + V} & {l for _, l, _ in X + Y}):
                            yield flipT(U) + [c1] + X + flipT(V) + [c2] + Y
                # B: c1 U X V c2 Y
                for i in range(2, p):
                    U = r[1:i]
                    for k in range(i, p):
                        X = r[i:k]; V = r[k:p]; Y = r[p + 1:]
                        if V and closed(U, V) and not ({l for _, l, _ in U + V} & {l for _, l, _ in X + Y}):
                            yield flipT(U) + [c1] + X + [c2] + flipT(V) + Y
                # C: c1 U c2 X V Y
                U = r[1:p]
                if U:
                    for k in range(p + 1, n2):
                        for m in range(k + 1, n2 + 1):
                            X = r[p + 1:k]; V = r[k:m]; Y = r[m:]
                            if closed(U, V) and not ({l for _, l, _ in U + V} & {l for _, l, _ in X + Y}):
                                yield flipT(U) + X + [c1] + flipT(V) + [c2] + Y
def fix_alt(w):
    # set O/U of every letter by alternation from position 0's neighbour pattern; return None if no alternating assignment matches the
    # O/U of the T-letters and outside letters (only c's letters are free) — simplest: require the word to alternate after re-setting c.
    for s0 in ('O', 'U'):
        ok = all(t[0] == (s0 if i % 2 == 0 else ('U' if s0 == 'O' else 'O')) for i, t in enumerate(w) if True)
        if ok: return w
    return None
def realign_c(w, cl):
    # c's two letters get the O/U that makes the word alternate (one must be O, the other U)
    for a in ('O', 'U'):
        b = 'U' if a == 'O' else 'O'; seen = 0; ww = []
        for o, l, s in w:
            if l == cl: ww.append((a if seen == 0 else b, l, s)); seen += 1
            else: ww.append((o, l, s))
        if fix_alt(ww): return ww
    return None
def orbit(L0, cap=200000):
    n = L0.size(); J = L0.jones(); seen = {L0.sig(False, True, True): word(L0)}; todo = [word(L0)]; rejected = 0
    while todo:
        w = todo.pop()
        for raw in moves(w):
            # the moved crossing c is the one whose relative position changed; realign each label's O/U by alternation trying every label
            cand = fix_alt(raw)
            if cand is None:
                cand = next((x for x in (realign_c(raw, l) for l in {t[1] for t in raw}) if x), None)
            if cand is None: rejected += 1; continue
            try: M = link(cand)
            except Exception: rejected += 1; continue
            if not (M.isClassical() and M.isAlternating() and M.size() == n): rejected += 1; continue
            sg = M.sig(False, True, True)
            if sg in seen: continue
            if M.jones() != J: rejected += 1; continue
            seen[sg] = cand; todo.append(cand)
            if len(seen) > cap: return seen, rejected, False
    return seen, rejected, True
def test(L0):
    S, rej, complete = orbit(L0); Mr = regina.Link(L0); Mr.reflect()
    return dict(orbit_size=len(S), rejected=rej, complete=complete, mirror_in_orbit=Mr.sig(False, True, True) in S)
def from_dt(code): return regina.Link.fromPD([[x + 1 for x in c] for c in snappy.Link('DT:' + str([tuple(dt(code))])).PD_code()])
def from_name(nm): return regina.Link.fromPD([[x + 1 for x in c] for c in snappy.Link(nm).PD_code()])

if __name__ == '__main__':
    import os as _os
    sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
    from census_codes import main
    AMPH = ('4_1', '6_3', '8_3', '8_9', '8_12', '8_17', '8_18', '10_17', '10_33', '10_37', '10_43', '10_45', '10_79', '10_81', '10_88',
            '10_99', '10_109', '10_115', '10_118', '10_123')
    counter = [(nm, (lambda nm=nm: from_name(nm)), (lambda r, e=e: r['complete'] and r['mirror_in_orbit'] is e)) for nm, e in
               [(k, True) for k in AMPH] + [(k, False) for k in ('3_1', '5_2', '8_19', '9_42', '10_1', '10_48', '10_71', '10_104')]
               if from_name(nm).isAlternating()]
    from census_codes import is_reduced
    def witness(L0):                   # verdict wrapper (driver): needs a REDUCED alternating diagram of a prime knot
        r = test(L0); r['alternating_diagram'] = L0.isAlternating(); r['reduced_diagram'] = is_reduced(L0.pdData())
        usable = r['alternating_diagram'] and r['reduced_diagram'] and r['complete']
        r['verdict'] = ('CHIRAL (mirror diagram not in the complete generated flype orbit)' if usable and not r['mirror_in_orbit'] else
                        ('mirror diagram in the generated orbit (no proof of amphichirality)' if usable else 'not usable (orbit incomplete or diagram not reduced alternating)'))
        return r
    main(witness, counter, from_dt)
