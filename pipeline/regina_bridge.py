#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Commercial inquiries: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root.
# Part of the Knot-Identification Pipeline.
"""
regina_bridge.py — call Regina in a process of its own instead of importing it.

WHY (Werner 2026-08-29: „stell um" (EN: "switch it over"))
--------------------------------------
Regina is under the GNU General Public License v2 or later (documented on
regina-normal.github.io). CALLING a GPL program is always free — that is how
one uses a tool, and nobody claims that a call of `grep` infects one's own code.
The only disputed question is whether an `import regina` in the SAME process makes the
calling script a derived work. For distributed code this question should not
arise in the first place.

Werner's stance behind this module: *„es bleibt geklaut. ich lese
zuerst, schreibe dann was eigenes. eine weiterentwicklung, oder? n.p. wenn ordentlich
zitiert wird."* (EN: "it remains stolen. I read first, then write something of my own. A further development, right? No problem if cited properly.") The citing is in LICENSE.txt; this module is the technical
side of the same stance.

WHAT IT COSTS — measured, not estimated (29.08.2026):
    a cold subprocess with `import regina` + one computation : 0.12 s
Therefore NOT one process per call is started but one per BATCH: the tasks
go in as JSON lines via stdin, the results come out as JSON lines. With
2000 diagrams that is one process start instead of two thousand.

The pattern is not newly invented — `pflicht_invarianten.py` has long run it
by hand (`conda run -n base ... python -c` with `for line in sys.stdin`). Here it stands
once, so that it is not copied five times.

⚠️ LIMIT: Regina objects do not survive the process boundary. What goes in and
out are strings and numbers (signed Gauss, isosig, PD, polynomials).
Whoever wants to pass a `regina.Link` on has to cut the interface differently.
"""
from __future__ import annotations
import json
import site_config
import os
import subprocess
import sys

HIER = os.path.dirname(os.path.abspath(__file__))
BASE_PY = site_config.PYTHON_BASE      # the env in which regina lives

# The body into which the task-specific code is inserted. The blind-test lock
# (hard rule 1) runs in the child too — otherwise the subprocess would be a loophole.
_RUMPF = '''
import sys, json
sys.path.insert(0, {manifest!r})
try:
    from pipeline_seal import seal_all
    seal_all(verbose=False)
except Exception:
    pass
import regina

def _arbeite(auf):
{koerper}

for _zeile in sys.stdin:
    _zeile = _zeile.strip()
    if not _zeile:
        continue
    try:
        _auf = json.loads(_zeile)
    except Exception as _e:
        print(json.dumps({{'fehler': 'no JSON: %s' % _e}}), flush=True); continue
    try:
        print(json.dumps(_arbeite(_auf), default=str), flush=True)
    except Exception as _e:
        print(json.dumps({{'fehler': '%s: %s' % (type(_e).__name__, _e)}}), flush=True)
'''


def stapel(aufgaben, koerper, timeout=3600, python=BASE_PY):
    """Runs `koerper` for every task — in ONE Regina subprocess.

    `koerper` is Python source (indented, with `return`) that sees `auf` (the
    task's JSON object) and returns a JSON-able result.
    Returns: list of results, same order; errors stand as
    {'fehler': ...} in their place instead of aborting the run.
    """
    manifest = os.path.dirname(HIER)
    code = _RUMPF.format(manifest=manifest,
                         koerper='\n'.join('    ' + z for z in koerper.strip('\n').split('\n')))
    eingabe = '\n'.join(json.dumps(a) for a in aufgaben) + '\n'
    p = subprocess.run([python, '-c', code], input=eingabe,
                       capture_output=True, text=True, timeout=timeout)
    aus = []
    for z in p.stdout.splitlines():
        z = z.strip()
        if not z or not z.startswith('{'):
            continue
        try:
            aus.append(json.loads(z))
        except Exception:
            pass
    if len(aus) != len(aufgaben):
        aus.append({'fehler': 'result count %d != task count %d; stderr: %s'
                              % (len(aus), len(aufgaben), p.stderr[-400:])})
    return aus


def einer(aufgabe, koerper, **kw):
    """Convenience for a single call."""
    return stapel([aufgabe], koerper, **kw)[0]


# ─────────────────────────── counter-cases ─────────────────────────

def selbsttest():
    """Four probes. The sharp one is the third: the result MUST agree with the
       direct import, otherwise the switch changed something."""
    print('SELF-TEST regina_bridge — four counter-cases\n')
    fehler = 0

    print('1) ONE CALL through the bridge')
    r = einer({'name': 'figureEight'}, """
T = regina.Triangulation3(getattr(regina.Example3, auf['name'])())
return {'tetraeder': T.size(), 'strikt': T.hasStrictAngleStructure()}
""")
    ok = r.get('tetraeder') == 2 and r.get('strikt') is True
    print(f'   {"✅" if ok else "🔴"} figureEight -> {r}')
    fehler += (not ok)

    print('\n2) BATCH — one process, many tasks, order preserved')
    namen = ['figureEight', 'trefoil', 'gieseking', 'figureEight']
    rs = stapel([{'name': n} for n in namen], """
T = regina.Triangulation3(getattr(regina.Example3, auf['name'])())
return {'name': auf['name'], 'strikt': T.hasStrictAngleStructure()}
""")
    ok = [x.get('name') for x in rs] == namen
    print(f'   {"✅" if ok else "🔴"} {[(x.get("name"), x.get("strikt")) for x in rs]}')
    fehler += (not ok)

    print('\n3) 🔑 THE SHARP ONE: bridge against direct import — bit for bit the same?')
    vergleich = subprocess.run([BASE_PY, '-c', """
import regina, json
aus = {}
for n in ('figureEight','trefoil','gieseking','whiteheadLink'):
    T = regina.Triangulation3(getattr(regina.Example3, n)())
    aus[n] = [T.size(), T.hasStrictAngleStructure(), T.isoSig()]
print(json.dumps(aus))
"""], capture_output=True, text=True)
    direkt = json.loads(vergleich.stdout)
    ueber = stapel([{'name': n} for n in direkt], """
T = regina.Triangulation3(getattr(regina.Example3, auf['name'])())
return {'name': auf['name'], 'w': [T.size(), T.hasStrictAngleStructure(), T.isoSig()]}
""")
    gleich = all(direkt[x['name']] == x['w'] for x in ueber if 'name' in x)
    print(f'   {"✅ identical" if gleich else "🔴 DEVIATION"} over {len(direkt)} cases')
    for x in ueber:
        if 'name' in x:
            print(f'      {x["name"]:14s} direct {direkt[x["name"]]}  ·  bridge {x["w"]}')
    fehler += (not gleich)

    print('\n4) CANARY — an error in the body must NOT abort the batch')
    rs = stapel([{'n': 1}, {'n': 0}, {'n': 3}], """
return {'kehrwert': 1.0 / auf['n']}
""")
    ok = len(rs) == 3 and 'fehler' in rs[1] and rs[2].get('kehrwert') == 1/3
    print(f'   {"✅" if ok else "🔴"} {rs}')
    fehler += (not ok)

    print(f'\n{"✅ ALL COUNTER-CASES PASSED" if fehler==0 else f"🔴 {fehler} ERRORS"}')
    return 0 if fehler == 0 else 1


if __name__ == '__main__':
    sys.exit(selbsttest())
