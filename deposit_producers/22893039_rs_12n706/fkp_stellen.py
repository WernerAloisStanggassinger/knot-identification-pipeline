#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""FKP threshold for EVERY Delta site — the step that `STAND.md` lists as missing.

STAND.md, step 5: „the threshold is computed at *one* place; for the **24** Delta sites
of 12n23 it must be redetermined for each". This script does exactly that.

Per site (arc diagram + crossing info + label parities) and per map:
  1. Build M = knot with TWIST CIRCLES instead of twist crossings (`komplement.band_mit_kreisen`).
  2. Verify hyperbolicity and compute vol(M) with interval arithmetic.
  3. Solve for the FKP threshold. Futer-Kalfagianni-Purcell:  vol(M(a)) >= (1 - (2pi/l)^2)^{3/2} vol(M).
     Sought is the smallest l from which the right-hand side lies above vol(target knot):
         l  >=  2pi / sqrt(1 - (vol_ziel / vol_M)^{2/3})
     Slopes that are longer are thereby PROVEN excluded.
  4. `short_slopes(length=l)` per cusp = the finite remainder list, which then has to be checked individually.

⚠️ TWO HONESTIES that the result does NOT exceed:
  · The twist circle only generates EVEN twists (measured, STAND.md). The odd family is
    still missing — a complete remainder list here therefore does not yet prove r_s(12n23) >= 5.
  · The threshold excludes LONG slopes. What remains is a finite list, not a proof;
    the list must afterwards be worked through piece by piece.

CALIBRATION BUILT IN: the same run first drives the sites of **12n706**, where the machinery on
20.09. yielded vol(M) = 30.507105 and l > 9.930. If that deviates, the remainder is worthless.

Usage:  python3 fkp_stellen.py [target ...]      (default: 12n706 12n23)
Each site runs in its own process (SnapPea dies by signal, not by exception).
"""
import json
import math
import os
import subprocess
import sys
import time

HIER = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HIER)

# target volumes: verified, measured on 20.09. (STAND.md, step 5)
VOL_ZIEL = {'12n23': 14.166394, '12n706': 13.417374}
DELTA = {'12n23': [2, -5, 2],
         '12n706': [1, -4, 6, -8, 11, -8, 6, -4, 1]}
FRIST = 900


def stellen(ziel):
    """All (word, crossing info, label parities) with Delta = Delta(target)."""
    from eiche_union import koeffs, koeffs_von_liste, faelle
    from band_rand import alexander_von
    from pruef_aglet_flaeche import WOERTER
    k = koeffs_von_liste(DELTA[ziel])
    aus = []
    for wort in WOERTER[4]:
        for w, ueber, lab, innere, n_str in faelle(wort):
            for s in alexander_von(w, ueber, lab, True, True):
                if koeffs(s) == k:
                    aus.append((wort, ueber, lab))
    return aus


KIND = r'''
import json, math, sys
sys.path.insert(0, %(hier)r)
wort, ueber, lab, kidx, vol_ziel = json.loads(%(arg)r)
ueber = {int(a): b for a, b in ueber.items()}
lab = {int(a): b for a, b in lab.items()}
w = [int(c) for c in wort]
from pruef_aglet_flaeche import karten
from komplement import band_mit_kreisen, dt_mit_drehung
from snappy import Link as _SL      # erst damit funktioniert spherograms exterior()
from spherogram import Link

rec = {'wort': wort, 'kidx': kidx, 'kreise': sorted(lab)}
komp = band_mit_kreisen(w, ueber, karten(w)[kidx], sorted(lab), True, True)
dt, off = dt_mit_drehung(komp)
if dt is None:
    rec['status'] = 'kein_dt'; print('REC' + json.dumps(rec)); raise SystemExit
rec['dt'] = [list(t) for t in dt]
dtstr = 'DT: [' + ','.join('(' + ','.join(str(v) for v in t) + ')' for t in dt) + ']'
rec['dt_string'] = dtstr
L = Link(dtstr)
M = L.exterior()
rec['cusps'] = M.num_cusps()
M.simplify()
ok, iv = M.verify_hyperbolicity()
rec['hyperbolisch'] = bool(ok)
if not ok:
    rec['status'] = 'nicht_hyperbolisch'; print('REC' + json.dumps(rec)); raise SystemExit
vol = M.volume(verified=True)
rec['vol_M'] = float(vol.center()) if hasattr(vol, 'center') else float(vol)
rec['vol_M_intervall'] = str(vol)
q = (vol_ziel / rec['vol_M']) ** (2.0 / 3.0)
if q >= 1.0:
    rec['status'] = 'schwelle_unmoeglich'   # vol(M) too small: FKP does not apply here
    print('REC' + json.dumps(rec)); raise SystemExit
rec['l_schwelle'] = 2 * math.pi / math.sqrt(1.0 - q)
sl = M.short_slopes(length=rec['l_schwelle'])
rec['short_slopes'] = [[list(map(int, s)) for s in cusp] for cusp in sl]
rec['restliste'] = [len(c) for c in rec['short_slopes']]
rec['status'] = 'ok'
print('REC' + json.dumps(rec))
'''


def eine(wort, ueber, lab, kidx, vol_ziel):
    arg = json.dumps([wort, {str(a): b for a, b in ueber.items()},
                      {str(a): b for a, b in lab.items()}, kidx, vol_ziel])
    code = KIND % {'hier': HIER, 'arg': arg}
    t0 = time.time()
    try:
        pr = subprocess.run([sys.executable, '-c', code], capture_output=True,
                            text=True, timeout=FRIST)
        for z in pr.stdout.splitlines():
            if z.startswith('REC'):
                r = json.loads(z[3:])
                r['sek'] = round(time.time() - t0, 1)
                return r
        return {'wort': wort, 'kidx': kidx, 'status': 'abbruch', 'rc': pr.returncode,
                'fehler': (pr.stderr or '')[-300:], 'sek': round(time.time() - t0, 1)}
    except subprocess.TimeoutExpired:
        return {'wort': wort, 'kidx': kidx, 'status': 'frist', 'sek': round(time.time() - t0, 1)}


def main():
    ziele = sys.argv[1:] or ['12n706', '12n23']
    for ziel in ziele:
        st = stellen(ziel)
        aus = os.path.join(HIER, f'fkp_stellen_{ziel}.jsonl')
        print(f'=== {ziel}: {len(st)} Delta sites x 2 maps = {2*len(st)} tasks '
              f'· vol(target) = {VOL_ZIEL[ziel]} · output {os.path.basename(aus)}', flush=True)
        with open(aus, 'a') as f:
            for i, (wort, ueber, lab) in enumerate(st, 1):
                for kidx in (0, 1):
                    r = eine(wort, ueber, lab, kidx, VOL_ZIEL[ziel])
                    r['ziel'] = ziel
                    f.write(json.dumps(r) + '\n')
                    f.flush()
                    kurz = (f"vol {r.get('vol_M'):.6f} · l {r.get('l_schwelle'):.3f} "
                            f"· Restliste {r.get('restliste')}") if r.get('status') == 'ok' else r.get('status')
                    print(f"  [{i}/{len(st)}] {wort} map {kidx}: {kurz} ({r['sek']}s)", flush=True)
        print(f'--- {ziel} done {time.strftime("%F %T")}', flush=True)


if __name__ == '__main__':
    main()
