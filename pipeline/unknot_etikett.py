#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""unknot_etikett.py — writes the unknot verdict (label) to where it will be found.

A Stage-2 abort with the message "no active cells" is TWO different things: for
a curve with a crossing-free projection it is a CORRECT rejection
(unknot, nothing to identify), for a twisted curve it is a gap.
Without a label both look the same and get recomputed next time.

Written only where TWO independent lines agree:
  1. kinematic    — scan_triage.py: is there a frame without a crossing?
  2. topological  — unknot_zweitmeinung.py: regina minimal crossing number == 0
On a contradiction NOTHING is written; it is reported instead.
"""
import sys, os, json, glob, argparse, datetime

def laden(p):
    return json.load(open(p))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--triage', required=True, help='scan_triage.py --json')
    ap.add_argument('--regina', required=True, help='unknot_zweitmeinung.py --json')
    ap.add_argument('--stand', required=True, help='root of the candidate folders')
    ap.add_argument('--registry', help='central registry file')
    ap.add_argument('--datum', required=True)
    ap.add_argument('--schreiben', action='store_true', help='write; without it = dry run')
    a = ap.parse_args()

    tri = {r['stamm']: r for r in laden(a.triage) if r.get('sorte') == 'UNKNOTEN'}
    reg = {}
    for r in laden(a.regina):
        reg[os.path.basename(os.path.dirname(r['path']))] = r

    agree, contradiction, without = [], [], []
    for stem, t in sorted(tri.items()):
        r = reg.get(stem)
        if r is None:
            # duplicate: the same curve under a different folder name
            hits = [v for k, v in reg.items() if v['stem'] == t['stamm'].split('__')[-1]]
            r = hits[0] if hits else None
        if r is None or 'fehler' in r:
            without.append(stem); continue
        if not r['unknot']:
            contradiction.append((stem, r['reduced'])); continue
        agree.append((stem, t, r))

    print(f"unknot candidates: {len(tri)}")
    print(f"  ✅ both lines agree   : {len(agree)}")
    print(f"  🔴 contradiction      : {len(contradiction)}  {contradiction}")
    print(f"  ❔ no second opinion  : {len(without)}  {without}")
    if contradiction:
        print("\n🔴 contradiction ⟹ NOTHING is labelled. Resolve it first.")
        return 1

    entries = []
    for stem, t, r in agree:
        v = {
            "stem": stem,
            "verdict": "unknot",
            "verdict_date": a.datum,
            "meaning": ("The Stage-2 abort for this curve is a CORRECT rejection, "
                        "not a pipeline failure: there is nothing to identify. "
                        "Do not re-queue."),
            "evidence": {
                "kinematic_scan": {
                    "tool": "claude_setup_v2/tools/scan_triage.py",
                    "n_crossings_total": t['n'],
                    "n_frames": t['frames'],
                    "min_crossings_per_frame": t['minc'],
                    "median_crossings_per_frame": t['medc'],
                    "fraction_of_crossing_free_frames": round(t['nullanteil'], 6),
                    "argument": ("A closed curve with a crossing-free projection is "
                                 "embedded in the plane, hence unknotted.")},
                "topological_gate": {
                    "tool": "claude_setup_v2/tools/unknot_zweitmeinung.py",
                    "method": "curve -> pyknotid planar diagram -> regina minimal crossing number",
                    "reduced_crossing_number": r['reduced'],
                    "calibration": "circle->0, trefoil->3, figure-eight->4 (self-test)",
                    "source": "maschinen_ideen/_topo_gate_audit/audit_topology.py (Disc #264)"}},
            "independence": ("The two lines share no code, no camera model and no "
                             "(ti,tj) grid; they cannot fail the same way.")}
        entries.append(v)
        target = os.path.join(a.stand, stem, f"{stem}_unknot_verdict.json")
        if a.schreiben:
            tmp = target + '.tmp'
            json.dump(v, open(tmp, 'w'), indent=1)
            os.replace(tmp, target)     # atomic, never a half-written file
    print(f"\n{'written' if a.schreiben else 'WOULD write'}: "
          f"{len(entries)} × <stamm>_unknot_verdict.json under {a.stand}")

    if a.registry and a.schreiben:
        os.makedirs(os.path.dirname(a.registry), exist_ok=True)
        alt = laden(a.registry) if os.path.exists(a.registry) else {"entries": []}
        existing = {e['stem'] for e in alt['entries']}
        new_entries = [e for e in entries if e['stem'] not in existing]
        alt['entries'] = sorted(alt['entries'] + new_entries, key=lambda e: e['stem'])
        alt['updated'] = a.datum
        alt['note'] = ("Curves that the identification pipeline correctly refused to "
                       "process because they are unknotted. Two independent lines of "
                       "evidence per entry. These are NOT pipeline failures.")
        tmp = a.registry + '.tmp'
        json.dump(alt, open(tmp, 'w'), indent=1)
        os.replace(tmp, a.registry)
        print(f"Registry: +{len(new_entries)} new, now {len(alt['entries'])} entries → {a.registry}")
    return 0

if __name__ == '__main__':
    sys.exit(main())
