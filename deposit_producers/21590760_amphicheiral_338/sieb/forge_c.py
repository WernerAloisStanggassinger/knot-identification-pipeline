#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""MACHINE 115-c — generalisation to Conant's rotation parameter c.

WHY. Conant's standard model (2607.18655, §1.2) is
      rho(z1,z2) = (e^{2*pi*i*c/2^a} z1, conj(z2)),  a >= 2, c ODD,
and the constant congruence from Corollary 7.3, E(0) = (l^2-1)/8 (mod 2), is only
sharp if l = +-3 (mod 8) — then E(0) must be ODD.

FINDING 28.07.: machine 115 is the case c = 1. Its xy modes are ns = [1] + [1 +- 2kj],
i.e. all n = 1 (mod 2k), and the winding number l is one of the occurring modes.
Also l = 1 (mod 2^a) BY CONSTRUCTION. Measured: a=2 -> l in {1,-3,5,-7,9}, all = 1 mod 4;
a=3 -> l in {1,-7,9,17}, all = 1 mod 8. That is why the sharp congruence ran empty at a=3:
not because it does not hold, but because the machine never enters its branch.

CURE. sigma_c = R_z(2*pi*c/2^a) o (z -> -z) is to act on the curve as a parameter shift
by 2*pi/2^a. That requires:
    xy modes:  n * (2*pi/2^a) = 2*pi*c/2^a  (mod 2*pi)   =>  n = c   (mod 2^a)
    z modes:   n * (2*pi/2^a) = pi          (mod 2*pi)   =>  n = 2^{a-1} (mod 2^a)
The z condition is unchanged (k = 2^{a-1}); only the xy basis moves from 1 to c.
For c = 1 forge_c falls back exactly to forge — that is checked below, bit-exactly.

⚠️ ATTRIBUTION + TRAP FOR a>=4 (added 28.07., prior-art hit against ourselves):
The rule is NOT ours. Conant Lemma 4.2 reads sharper and the other way round

        c * l = 1  (mod 2^a) ,

from which l odd follows. For a <= 3 that is the same as l = c, because in (Z/4)* and (Z/8)*
every element is its own inverse (3*3=9=1, 5*5=25=1, 7*7=49=1 mod 8). FROM a=4 NO LONGER:
mod 16, c=3 has l = 11, not 3. Whoever uses this file for k=8 (a=4, order 16),
must compute l = c^{-1} (mod 2^a).

🔴 CORRECTION 31.07.: here stood "Our measurements below CONFIRM Conant's Lemma". That is
WRONG and always was. forge_c sets the xy modes to n = c (mod 2^a) (lines 53-55), and
l IS one of these modes ==> l = c (mod 2^a) is TRUE BY CONSTRUCTION, not measured.
The column c*l = 1 is thus equivalent to c^2 = 1 (mod 2^a) — for a <= 3 true for EVERY
odd c (pure arithmetic). A test that cannot fail proves nothing.
This file can neither confirm nor refute Conant's Lemma 4.2. What is measured below is
ONLY the sigma_c error; everything else is the machine's consistency with itself.
The same kind of sentence flew out of the letter on 29.07. and out of bilanz.py on 30.07.

For a=3 it means: c=3 -> l = 3 (mod 8) -> E(0) MUST be odd.
"""
import sys, os
import numpy as np

HIER = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HIER)
import importlib.util
sp = importlib.util.spec_from_file_location('a', os.path.join(HIER, 'stufe_a_pd.py'))
A = importlib.util.module_from_spec(sp); sp.loader.exec_module(A)

N = A.N


def forge_c(k, seed, n_xy=3, n_z=2, decay=0.55, c=1):
    """Machine 115 with rotation parameter c (odd). c=1 == forge() bit-exactly."""
    if c % 2 == 0:
        raise ValueError('c must be odd (Conant Lemma 4.2: l is odd)')
    rng = np.random.default_rng(seed)
    t = np.linspace(0, 2 * np.pi, N, endpoint=False)
    ns_xy = [c]
    for j in range(1, n_xy + 1):
        ns_xy += [c + 2 * k * j, c - 2 * k * j]
    cc = np.zeros(N, complex)
    for n in ns_xy:
        amp = (1.0 if n == c else decay ** abs((n - c) // (2 * k)))
        a = amp * (rng.normal() + 1j * rng.normal()) * (0.45 if n != c else 1.2)
        cc += a * np.exp(1j * n * t)
    z = np.zeros(N)
    for j in range(n_z):
        n = k + 2 * k * j
        z += (decay ** j) * 0.7 * (rng.normal() * np.cos(n * t) + rng.normal() * np.sin(n * t))
    return np.column_stack([cc.real, cc.imag, z])


def sigma_c_ok(P, k, c, tol=1e-12):
    """Certificate: sigma_c(P(t)) == P(t + pi/k) bit-exactly, with sigma_c = R_z(c*pi/k) o (z->-z)."""
    ang = c * np.pi / k
    R = np.array([[np.cos(ang), -np.sin(ang), 0],
                  [np.sin(ang), np.cos(ang), 0],
                  [0, 0, -1.0]])
    return float(np.abs(P @ R.T - np.roll(P, -(N // (2 * k)), axis=0)).max())


def windung(P):
    z = P[:, 0] + 1j * P[:, 1]
    return int(np.round((np.unwrap(np.angle(np.append(z, z[0])))[-1] - np.angle(z[0])) / (2 * np.pi)))


if __name__ == '__main__':
    print('=== CALIBRATION 1: c=1 must reproduce forge() BIT-EXACTLY ===')
    schlecht = 0
    for k in (2, 4):
        for seed in (7, 29, 715, 1347):
            for nxy in (2, 3, 4):
                P0 = A.forge(k, seed, nxy, 2, 0.55)
                P1 = forge_c(k, seed, nxy, 2, 0.55, c=1)
                d = float(np.abs(P0 - P1).max())
                if d != 0.0:
                    schlecht += 1
                    print(f'  🔴 k={k} seed={seed} nxy={nxy}: max|diff| = {d:.3e}')
    print(f'  {"✓ bit-exactly identical" if not schlecht else f"🔴 {schlecht} deviations"}')
    if schlecht:
        sys.exit(1)

    print('\n=== CALIBRATION 2: sigma_c is a symmetry of the curve (construction certificate) ===')
    print('⚠️  This table checks the MACHINE, not the Lemma 4.2 of Conant.')
    print('    forge_c sets the xy modes to n = c (mod 2^a); l is one of these modes')
    print('    ==> column "l=c?" is TRUE BY CONSTRUCTION, not a measurement result.')
    print('    ==> column "c*l=1?" is thus the same as c^2 = 1 (mod 2^a) — pure arithmetic.')
    print('    ONLY the sigma_c error is MEASURED here.')
    print(f'{"k":>2} {"a":>2} {"c":>2} {"seed":>5} {"sigma error":>14} {"l":>4} {"l mod 2^a":>10} '
          f'{"l=c?":>6} {"c*l mod 2^a":>12} {"c*l=1?":>7} {"E(0) must":>10}')
    fehler = 0            # real violation: sigma_c or construction consistency
    n_zeilen = 0
    conant_ok = 0         # only COUNTED, never counted as receipt
    conant_bruch = {}     # a -> set of the c for which c*l != 1
    for k, a in ((2, 2), (4, 3), (8, 4)):      # a=4 = the COUNTER-CASE (footgun 16)
        M = 2 ** a
        for c in range(1, M, 2):
            for seed in (7, 29, 715):
                P = forge_c(k, seed, 3, 2, 0.55, c=c)
                s = sigma_c_ok(P, k, c)
                l = windung(P)
                ok_konstr = (l - c) % M == 0          # machine against itself
                ok_conant = (c * l - 1) % M == 0      # Conant's condition, only observed
                e0 = 'ungerade' if ((abs(l) ** 2 - 1) // 8) % 2 else 'gerade'
                n_zeilen += 1
                if s > 1e-11 or not ok_konstr:
                    fehler += 1
                if ok_conant:
                    conant_ok += 1
                else:
                    conant_bruch.setdefault(a, set()).add(c)
                print(f'{k:>2} {a:>2} {c:>2} {seed:>5} {s:>14.2e} {l:>4} {l % M:>10} '
                      f'{"✓" if ok_konstr else "🔴":>6} {(c * l) % M:>12} '
                      f'{"✓" if ok_conant else "🔴":>7} {e0:>10}')

    print()
    if fehler:
        print(f'  🔴 {fehler} violations (sigma_c error or l != c mod 2^a) — machine broken')
        sys.exit(1)
    print(f'  ✓ MEASURED: sigma_c certificate holds in {n_zeilen}/{n_zeilen} rows (< 1e-11).')
    print(f'  ✓ CONSTRUCTION: l = c (mod 2^a) in {n_zeilen}/{n_zeilen} — self-consistency, NO receipt.')
    print(f'  ℹ  OBSERVED (no receipt): c*l = 1 (mod 2^a) holds in {conant_ok}/{n_zeilen} rows.')
    for a in sorted(conant_bruch):
        cs = ', '.join(str(c) for c in sorted(conant_bruch[a]))
        print(f'     a={a}: violated at c = {cs}  (c^2 != 1 mod {2 ** a})')
    if conant_bruch:
        print('     EXPECTED, and precisely that is the counter-case: the machine delivers l = c,')
        print('     Conant requires l = c^{-1}. From a=4 on the two fall apart. Whoever wants Conant\'s')
        print('     whoever wants the sharp branch at order 16 must aim at l = c^{-1} (mod 16).')
    else:
        print('     🔴 NO counter-case occurred — then this table checks nothing (footgun 16).')
        sys.exit(1)
