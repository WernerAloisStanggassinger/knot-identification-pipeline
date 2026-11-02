# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Cross-check, not a witness: SnapPy's NUMERICAL symmetry group (not certified) says amphicheiral or not, for census names.
# Usage: python -I snappy_symmetry_check.py names.txt   (BURTON_CENSUS_DIR as for the witnesses)
import sys, os, snappy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from census_codes import codes, DT
names = [l.split()[0] for l in open(sys.argv[1]) if l.strip()]
for n, v in codes(names).items():
    print(n, 'SnapPy symmetry_group().is_amphicheiral() (numerical):',
          snappy.Link('DT:' + str([tuple(DT(v))])).exterior().symmetry_group().is_amphicheiral())
