#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""pd_to_theta_input.py

PD code -> (Cs, phi) converter for van der Veen's Theta.sage
(Bar-Natan/van der Veen, "Perturbed Gaussian generating series").

Background (Discovery #133):
 - Theta.sage expects input `[Cs, phi]` where
     * Cs is an n x 3 crossing matrix with rows (sign, i, j),
       sign in {+1, -1}, i = half-strand index of the incoming over-strand,
       j = half-strand index of the incoming under-strand;
       half-strand indices 0..2n-1 (with outgoing = incoming+1).
     * phi is a rotation vector of length 2n+1.
 - This file reproduces exactly the Rot[] function from
   external_tools/bar_natan_theta/Theta.m (lines 22-40).

Input convention: standard PD list of 4-tuples, each
(a, b, c, d) -- counter-clockwise starting at the incoming under-strand
(KnotTheory` convention). PD edge labels are 1..2n.

Direct call:
    python pd_to_theta_input.py            # self-test on the trefoil

Library:
    from pd_to_theta_input import pd_to_cs_phi
    cs_rows, phi = pd_to_cs_phi([(1,4,2,5),(3,6,4,1),(5,2,6,3)])
"""
from __future__ import annotations
import json
import sys
from typing import List, Sequence, Tuple


def _is_positive_crossing(pd: Sequence[int], two_n: int) -> bool:
    """PositiveQ in the KnotTheory` sense: True iff b - d == 1 (mod 2n).

    For X[a,b,c,d] with 1-based PD edge labels.
    """
    a, b, c, d = pd
    return ((b - d) % two_n) == 1


def pd_to_cs_phi(pd_list: Sequence[Sequence[int]]) -> Tuple[List[Tuple[int, int, int]], List[int]]:
    """Convert a standard PD code to (Cs rows, phi).

    Args:
        pd_list: list of 4-tuples (a,b,c,d) with 1-based PD edge labels.

    Returns:
        (cs_rows, phi):
          cs_rows -- list of (sign, i, j) triples in the order of the
                     given PD crossings, 0-based half-strand indices
                     (0..2n-1). i = incoming over-strand, j = incoming under-strand.
          phi     -- list of length 2n+1, rotation numbers per
                     half-strand index (position 0..2n).

    Note: 1:1 translation of the Rot[] function from Theta.m.
    Half-strand indices are computed 1-based internally here (as
    in the Mathematica source) and shifted by -1 at the end, so that
    they match the Sage theta convention (0-based, compatible with
    knots.sage entries).
    """
    pd_list = [tuple(c) for c in pd_list]
    n = len(pd_list)
    two_n = 2 * n
    if two_n == 0:
        return [], [0]

    # Step 1: build xs = list of Xp/Xm markers (1-based internally).
    # Xp[d,a]  for a positive crossing
    # Xm[b,a]  for a negative crossing
    # We store (sign, i, j) with i=over-in, j=under-in.
    xs: List[List[int]] = []  # mutable list per crossing
    for c in pd_list:
        a, b, _c, d = c
        if _is_positive_crossing(c, two_n):
            xs.append([+1, d, a])
        else:
            xs.append([-1, b, a])

    # Step 2: front propagation with rotation accumulation.
    # rots is 1-based as in Mathematica: rots[1..2n+1] -- we
    # allocate length 2n+2 and ignore index 0.
    rots = [0] * (two_n + 2)
    front: List[int] = [1]

    def _match_and_replace(token: int) -> List[int]:
        """Find the crossing in xs that has `token` as i (over-in)
        or j (under-in), and build the replacement list.
        Mathematica pattern:
          Xp[k,l_] | Xm[l_,k]   -> {l+1, k+1, -l}
          Xp[l_,k] | Xm[k,l_]   -> ++rots[[l]]; {-l, k+1, l+1}
        token = current k.
        """
        for sign, i, j in xs:
            # Case A: Xp[k,l] corresponds to (sign=+1, i=k, j=l);
            #         Xm[l,k] corresponds to (sign=-1, i=l, j=k).
            if sign == +1 and i == token:
                l = j
                return [l + 1, token + 1, -l]
            if sign == -1 and j == token:
                l = i
                return [l + 1, token + 1, -l]
            # Case B: Xp[l,k] corresponds to (sign=+1, i=l, j=k);
            #         Xm[k,l] corresponds to (sign=-1, i=k, j=l).
            if sign == +1 and j == token:
                l = i
                rots[l] += 1  # rots[[l]] in Mathematica
                return [-l, token + 1, l + 1]
            if sign == -1 and i == token:
                l = j
                rots[l] += 1
                return [-l, token + 1, l + 1]
        return []  # token not found (should not happen)

    for k in range(1, two_n + 1):
        if -k not in front:
            new_front: List[int] = []
            replaced = False
            for elt in front:
                if elt == k and not replaced:
                    new_front.extend(_match_and_replace(k))
                    replaced = True
                else:
                    new_front.append(elt)
            front = new_front
        else:
            # `-k` is in front. Mathematica equivalent:
            #   Cases[front, k|-k] /. {k,-k} :> --rots[[k]]
            # `Cases[front, k|-k]` extracts all occurrences of k or
            # -k in the order in which they appear in front.
            # The list pattern `{k,-k}` then matches exactly the
            # 2-element sequence [k, -k] (k before -k). Exactly then
            # rots[[k]] is decremented by 1.
            cases = [x for x in front if x == k or x == -k]
            if cases == [k, -k]:
                rots[k] -= 1

    # Step 3: rots is 1-based (index 1..2n+1). Convert to a
    # 0-based phi vector of length 2n+1 by taking rots[1..2n+1].
    # The knots.sage entry also uses length 2n+1.
    phi = rots[1:two_n + 2]  # length 2n+1

    # Step 4: make the Cs indices 0-based (Sage convention).
    cs_rows: List[Tuple[int, int, int]] = [
        (int(sign), int(i) - 1, int(j) - 1) for sign, i, j in xs
    ]
    return cs_rows, phi


def cs_phi_to_sage_literal(cs_rows: Sequence[Sequence[int]], phi: Sequence[int]) -> str:
    """Format (Cs, phi) as a Sage literal that can be passed to Theta(...)."""
    mat_rows = ", ".join(
        "[{}, {}, {}]".format(s, i, j) for s, i, j in cs_rows
    )
    phi_str = ", ".join(str(int(x)) for x in phi)
    return "[matrix([{}]), [{}]]".format(mat_rows, phi_str)


# --------------------------------------------------------------------
# Self-test
# --------------------------------------------------------------------
def _selftest() -> int:
    # Standard trefoil PD from KnotTheory`:
    trefoil_pd = [(1, 4, 2, 5), (3, 6, 4, 1), (5, 2, 6, 3)]
    cs, phi = pd_to_cs_phi(trefoil_pd)
    print("=== 3_1 (Trefoil) ===")
    print("PD  =", trefoil_pd)
    print("Cs  =", cs)
    print("phi =", phi)
    print("Sage literal:", cs_phi_to_sage_literal(cs, phi))

    # Expected from knots.sage:
    expected_cs = [(-1, 3, 0), (-1, 5, 2), (-1, 1, 4)]
    expected_phi = [0, 0, 0, -1, 0, 0, 0]
    ok_cs = cs == expected_cs
    ok_phi = phi == expected_phi
    print("Cs  match knots.sage?", ok_cs)
    print("phi match knots.sage?", ok_phi)
    overall = ok_cs and ok_phi

    # Unknot with a single twist (Reidemeister I) -- a constructive
    # test: the PD of a twisted unknot over two crossings (R1 + R1^-1) is
    # not trivial, so we leave it out here and rely on an
    # additional end-to-end test against knots.sage in a Sage script.

    return 0 if overall else 1


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--selftest":
        sys.exit(_selftest())
    if len(sys.argv) > 1:
        # Accept the PD as JSON on stdin or as a path to a JSON file.
        with open(sys.argv[1], "r") as fh:
            pd = json.load(fh)
        cs, phi = pd_to_cs_phi(pd)
        print(json.dumps({"cs": cs, "phi": phi}))
    else:
        sys.exit(_selftest())
