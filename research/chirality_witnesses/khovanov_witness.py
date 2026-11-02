# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Usage (environment with khoca): KH_PD=pds.json [KH_ALG=sl2Z|sl3Z|BNQ|BNF2] python -I khovanov_witness.py NAME   (pds.json from make_kh_pd.py)
# KH_ALG=sl3Z -> khoca sl(3) homology over Z (Frobenius algebra Z[X]/X^3); default sl2Z.
# Khovanov chirality witness, one knot: khoca (Lewark) over Z (InteractiveCalculator(0)), reduced + unreduced homology of K and of its
# mirror (SnapPy PDs from make_kh_pd.py). Amphichiral ⟹ K and mK are isotopic ⟹ equal Khovanov homology; we compare the homology of the
# two diagrams directly; counter-cases (Z coefficients, in the PD file from make_kh_pd.py counter): 3_1, 9_42, 10_125 differ; 4_1 equal.
# differ ⟹ K not isotopic to mK ⟹ chiral (independent of SnapPy's isometry signature). equal ⟹ no information.
import json, sys, re, time, collections, os
from khoca import InteractiveCalculator
term = re.compile(r'(\d*)\s*\*?\s*t\^(-?\d+)\s*q\^(-?\d+)(?:\[(\d+)\])?')
def parse(s): 
    c = collections.Counter()
    for m in term.finditer(s): c[(int(m.group(2)), int(m.group(3)), int(m.group(4) or 0))] += int(m.group(1) or 1)
    return c
name = sys.argv[1]; rec = next(x for x in json.load(open(os.environ.get('KH_PD','pd_kh_targets.json'))) if x[0] == name)
ALG = os.environ.get("KH_ALG", "sl2Z")
CALC = {"sl2Z": (0,), "sl3Z": (0, (0, 0, 0), 0), "BNQ": (1, (0, 1), 0), "BNF2": (2, (0, 1), 0)}
# BNQ / BNF2: Bar-Natan deformation Q[X]/(X^2+X) resp. F2[X]/(X^2+X) — its homology carries the Rasmussen invariant s
# (over Q) resp. the F2 version; s(mK) = -s(K), so a difference K vs mK = chirality witness; output groups compared as a whole.
KH = InteractiveCalculator(*CALC[ALG]); t0 = time.time(); H = []
for pd in rec[2:4]:
    _, m = KH([list(c) for c in pd], print_messages=True); H.append([parse(x) for x in m if 't^' in x])
eq = H[0] == H[1]
print(json.dumps(dict(name=name, kind=rec[1], kh_equal_mirror=eq, n_groups=[len(h) for h in H], total_rank_K=[sum(v for k, v in h.items() if k[2] == 0) for h in H[0]],
                      alg=ALG, verdict=f'CHIRAL ({ALG} homology differs from mirror)' if not eq else f'{ALG} blind (equal to mirror)', sec=round(time.time() - t0, 1))), flush=True)
