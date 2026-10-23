#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""sort_and_render.py — stage 9 SORT+RENDER for pipeline_run.sh.

Reads a knot folder (with pflicht_*_complete.json + crosscheck.json),
decides the sort target according to [[feedback_render_and_sort_on_new_knots]]
and [[project_naming_doctrine]], renames/moves it into the right folder and
triggers a JPG render if none exists yet and there is no KnotInfo hit.

Decision matrix:
  1. Failure (rc≠0 or missing mandatory core fields) → erzeugte_knoten/pipeline_failures/
  2. KnotInfo hit (crosscheck.knotinfo_match.status startswith "matched") →
     erzeugte_knoten/bekannt/<klassisch>[__<original>]/
  3. Cluster match (fingerprint c<c>_d<det>_<sym>_<chir> exists) →
     erzeugte_knoten/<fingerprint>/
  4. Otherwise: stays a top-level singleton (or erzeugte_knoten/singletons/ depending on convention)

Classical notation (Discovery #175):
  c≤10: K prefix stays (K0_1, K3_1, K8_19, K10_43)
  c≥11: no K prefix, with underscore (12n_837, 11n_133)
  Mirror: suffix `m`
  Multi-hit: <klassisch>__<original>

own_knots/ sculptures are NOT renamed (artist names are kept).

Usage:
  sort_and_render.py <knot_dir>                      # process one knot
  sort_and_render.py --batch [<root>]                # all knots under <root> (default: erzeugte_knoten/)
  sort_and_render.py <knot_dir> --no-render          # sort only, no render
  sort_and_render.py <knot_dir> --dry-run            # execute nothing, only log
"""
from __future__ import annotations
import argparse
import site_config
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

# === Paths ===
ROOT = Path(site_config.ROOT)
ERZ = Path(site_config.GENERATED_KNOTS_DIR)
BEKANNT = ERZ / 'bekannt'
FAILURES = ERZ / 'pipeline_failures'
SINGLETONS = ERZ / 'new_found_knots' / 'singletons_novelty'
FREEZE = Path(site_config.FREEZE_DIR)
RENDER_TEMPLATE = Path(site_config.RENDER_TEMPLATE)
RENDER_SCRIPT = Path(site_config.PIPELINE_DIR) / 'render_knot_blender.py'
BLENDER = site_config.BLENDER

# === Classical-name detection ===
# c≤10 folders have a K prefix; c≥11 do not.
KNOTINFO_PREFIX_RE = re.compile(r'^(K\d+_\d+m?|\d+[an]_\d+m?)')

# Mandatory core fields that must not all be None (failure type b detection)
CORE_FIELDS = ('signature', 'symmetry_group', 'hyperbolic_volume', 'amphicheiral')


def is_knotinfo_folder(name: str) -> bool:
    """True if folder name matches klassisch KnotInfo prefix."""
    return bool(KNOTINFO_PREFIX_RE.match(name))


def find_pflicht_complete(d: Path) -> Path | None:
    """Find <stem>_pflicht_invarianten*_complete.json in folder."""
    candidates = sorted(d.glob('*_pflicht_invarianten*_complete.json'))
    return candidates[0] if candidates else None


def load_pflicht_data(pfl: Path) -> dict | None:
    try:
        return json.load(open(pfl))
    except Exception:
        return None


def load_crosscheck_data(d: Path, stem: str) -> dict | None:
    cc = d / f'{stem}_crosscheck.json'
    if cc.exists():
        try: return json.load(open(cc))
        except: return None
    return None


def is_failure(d: Path, pfl: Path | None, pdata: dict | None) -> tuple[bool, str]:
    """Failure detection.
    Type a: no _complete.json
    Type b: _complete.json exists but the mandatory core fields are all None.
    """
    if pfl is None or pdata is None:
        return True, 'type_a (no _complete.json)'
    pl = pdata.get('pflichtliste', {})
    alg = pl.get('ALGEBRA', {}) or {}
    hyp = pl.get('HYPERBOLISCH', {}) or {}
    vals = [
        alg.get('signature'),
        hyp.get('symmetry_group'),
        hyp.get('hyperbolic_volume'),
        hyp.get('is_amphicheiral'),
    ]
    if all(v is None for v in vals):
        return True, 'type_b (mandatory core fields all None)'
    return False, ''


def get_knotinfo_match(cdata: dict | None) -> str | None:
    """Return the classical name from crosscheck.knotinfo_match.best_match if high-confidence."""
    if not cdata:
        return None
    km = cdata.get('knotinfo_match', {}) or {}
    status = km.get('status', '')
    best = km.get('best_match')
    if status.startswith('matched') and best:
        return best
    # Snappy fallback (also classical after the pipeline_crosscheck.py patch)
    sm = cdata.get('snappy_census_match', {}) or {}
    if sm.get('status', '').startswith('matched') and sm.get('best_match'):
        return sm['best_match']
    return None


def fingerprint_from_pflicht(pdata: dict) -> tuple[int|None, int|None, str|None, str|None]:
    """Compute (c, det, sym, chir) for cluster-fingerprint."""
    pl = pdata.get('pflichtliste', {})
    diag = pl.get('DIAGRAMM', {}) or {}
    alg = pl.get('ALGEBRA', {}) or {}
    hyp = pl.get('HYPERBOLISCH', {}) or {}
    c = diag.get('crossing_number_min')
    det = alg.get('determinant')
    sym_full = hyp.get('symmetry_group') or ''
    # Normalise sym: D6, D3, Z/2, Z/3, etc.
    sym = re.sub(r'\s+', '', sym_full).replace('+', '_').replace('/', '')
    sym = sym or None
    amphi = hyp.get('is_amphicheiral')
    chir = 'amphi' if amphi else ('chiral' if amphi is not None else None)
    return c, det, sym, chir


def cluster_folder_name(c, det, sym, chir) -> str | None:
    """Stanggassinger class folder name (Discovery #177).

    Format: S_<c>_d<det>_<sym>_<chir>
    Optional suffixes (manual/later):
      _m      mirror of an existing S_ class
      _v<n>   volume discriminator if (c,det,sym,chir) has several Theta classes
    """
    if c is None or det is None or not sym or not chir:
        return None
    return f'S_{c}_d{det}_{sym}_{chir}'


def find_existing_cluster(c, det, sym, chir) -> Path | None:
    """Look for existing S_<c>_d<det>_<sym>_<chir> folder in erzeugte_knoten/.

    _v<n> suffixes and _m are NOT recognised automatically either — the cluster
    match branch looks only for the primary S_ folder. The Theta check decides
    later whether it can be promoted.
    """
    fname = cluster_folder_name(c, det, sym, chir)
    if not fname:
        return None
    p = ERZ / fname
    return p if p.exists() and p.is_dir() else None


def load_theta_numeric(folder: Path, stem: str) -> str | None:
    """Load theta_at_22_7__21_13 as a string (for an exact-equality comparison)."""
    tn = folder / f'{stem}_theta_numeric.json'
    if not tn.exists():
        return None
    try:
        d = json.load(open(tn))
        # Field name varies: theta_total or theta
        return str(d.get('theta_total') or d.get('theta') or '')
    except Exception:
        return None


def theta_match(folder_a: Path, folder_b: Path) -> bool:
    """True if both folders have theta values AND they are identical."""
    stem_a = re.sub(r'_pflicht_invarianten.*_complete\.json$', '',
                    find_pflicht_complete(folder_a).name) if find_pflicht_complete(folder_a) else folder_a.name
    stem_b = re.sub(r'_pflicht_invarianten.*_complete\.json$', '',
                    find_pflicht_complete(folder_b).name) if find_pflicht_complete(folder_b) else folder_b.name
    ta = load_theta_numeric(folder_a, stem_a)
    tb = load_theta_numeric(folder_b, stem_b)
    if ta is None or tb is None:
        return False  # missing theta → no match
    return ta == tb


def find_singleton_match(c, det, sym, chir, exclude: Path) -> Path | None:
    """Look for an existing singleton (in singletons/ or erzeugte_knoten/ Toplevel)
    with the same fingerprint. If found, this is a 2nd-singleton and should be
    promoted to a Family-Folder.
    """
    if c is None or det is None:
        return None
    candidates = []
    # Search singletons/
    if SINGLETONS.exists():
        candidates.extend([p for p in SINGLETONS.iterdir() if p.is_dir()])
    # Search Toplevel of erzeugte_knoten/ (excluding bekannt/, pipeline_failures/, singletons/, family-folders)
    for p in ERZ.iterdir():
        if not p.is_dir():
            continue
        if p.name in ('bekannt', 'pipeline_failures', 'singletons', 'new_found_knots', 'singletons_unknots'):
            continue
        if p.name.startswith('c') and re.match(r'^c\d+_d\d+_', p.name):
            continue  # is a family folder
        if p.name.startswith('.'):
            continue  # hidden
        candidates.append(p)
    for cand in candidates:
        if cand == exclude:
            continue
        # Load pflicht of candidate to compare fingerprint
        pfl = find_pflicht_complete(cand)
        if not pfl:
            continue
        pdata = load_pflicht_data(pfl)
        if not pdata:
            continue
        cc, dd, ss, ch = fingerprint_from_pflicht(pdata)
        if cc == c and dd == det and ss == sym and ch == chir:
            return cand
    return None


# --- folder/file rename helper -----------------------------------------

def rename_folder_and_files(old_dir: Path, new_dir: Path, dry_run: bool = False) -> Path:
    """Rename folder old_dir -> new_dir, also rename all interior files to use new prefix."""
    old_name = old_dir.name
    new_name = new_dir.name
    log_prefix = '[DRY] ' if dry_run else ''
    if old_dir == new_dir:
        return new_dir
    print(f'{log_prefix}MOVE {old_dir} -> {new_dir}')
    if dry_run:
        return new_dir
    new_dir.parent.mkdir(parents=True, exist_ok=True)
    if new_dir.exists():
        raise FileExistsError(f'Target {new_dir} exists')
    # First rename interior files (while old_dir still exists)
    for f in list(old_dir.iterdir()):
        if not f.is_file():
            continue
        if f.name.startswith(old_name):
            nf = old_dir / (new_name + f.name[len(old_name):])
            f.rename(nf)
    # Then rename folder
    old_dir.rename(new_dir)
    return new_dir


def patch_content_refs(target_dir: Path, old_name: str, new_name: str, dry_run: bool = False) -> int:
    """Replace all occurrences of old_name with new_name in JSON/MD/log/txt content
    inside target_dir.
    """
    if old_name == new_name:
        return 0
    n = 0
    for f in target_dir.rglob('*'):
        if not f.is_file():
            continue
        if f.suffix.lower() not in ('.json', '.md', '.log', '.txt'):
            continue
        if '.bak' in f.name:
            continue
        try:
            content = f.read_text(encoding='utf-8', errors='ignore')
        except Exception:
            continue
        if old_name not in content:
            continue
        new_content = content.replace(old_name, new_name)
        if new_content == content:
            continue
        if dry_run:
            print(f'[DRY] PATCH {f}')
        else:
            bak = f.with_suffix(f.suffix + '.bak-pre-sort')
            if not bak.exists():
                shutil.copy2(f, bak)
            f.write_text(new_content, encoding='utf-8')
        n += 1
    return n


def klassisch_target_name(klassisch_id: str, original_stem: str) -> str:
    """Compose a multi-hit folder name or a singleton name.

    Rules:
      - klassisch_id is Pipeline-Crosscheck-Output (already klassisch).
      - If original_stem starts with klassisch_id (= already correctly named):
          single name.
      - Else: <klassisch_id>__<original_stem>.
    """
    if original_stem == klassisch_id:
        return klassisch_id
    if original_stem.startswith(klassisch_id + '__'):
        return original_stem  # already correctly Multi-Hit
    return f'{klassisch_id}__{original_stem}'


def determine_target(folder: Path, dry_run: bool = False) -> tuple[str, Path, str]:
    """Determine target action.
    Returns (kind, target_path, reason). Kind ∈ {'failure','knotinfo','cluster','cluster_promote','singleton','noop','own_knots_skip'}.
    """
    # own_knots/ sculptures are off-limits
    rel = folder.relative_to(ROOT) if folder.is_relative_to(ROOT) else folder
    if any(part == 'own_knots' for part in rel.parts):
        return 'own_knots_skip', folder, 'own_knots: artist names stay (doctrine)'

    pfl = find_pflicht_complete(folder)
    pdata = load_pflicht_data(pfl) if pfl else None

    # Failure check
    fail, reason = is_failure(folder, pfl, pdata)
    if fail:
        target = FAILURES / folder.name
        return 'failure', target, reason

    # KnotInfo check
    stem_match = re.match(r'(.+?)_pflicht_invarianten.*_complete\.json$', pfl.name)
    stem = stem_match.group(1) if stem_match else folder.name
    cdata = load_crosscheck_data(folder, stem)
    ki_id = get_knotinfo_match(cdata)
    if ki_id:
        new_name = klassisch_target_name(ki_id, folder.name)
        target = BEKANNT / new_name
        return 'knotinfo', target, f'KnotInfo hit: {ki_id}'

    # Cluster fingerprint + Theta check (Discovery #177)
    c, det, sym, chir = fingerprint_from_pflicht(pdata)
    if c is not None and det is not None and sym and chir:
        existing = find_existing_cluster(c, det, sym, chir)
        if existing:
            # Theta match is mandatory: check against an existing member
            members = [m for m in existing.iterdir() if m.is_dir()]
            if members:
                ref_member = members[0]
                if theta_match(folder, ref_member):
                    target = existing / klassisch_target_name(existing.name, folder.name)
                    return 'cluster', target, f'cluster match (theta-verified): {existing.name}'
                else:
                    # Mandatory fingerprint matched but Theta did not — separate class
                    return 'singleton', SINGLETONS / folder.name, f'mandatory-list match {existing.name}, but theta differs → own class, stays a singleton'
            else:
                # Cluster folder empty → simply move in
                target = existing / klassisch_target_name(existing.name, folder.name)
                return 'cluster', target, f'cluster match (empty cluster): {existing.name}'
        # Check for a 2nd singleton with Theta verification
        sib = find_singleton_match(c, det, sym, chir, exclude=folder)
        if sib:
            if theta_match(folder, sib):
                new_cluster = ERZ / cluster_folder_name(c, det, sym, chir)
                return 'cluster_promote', new_cluster, f'2nd singleton (theta-verified), promoted to an S_ cluster with sibling: {sib.name}'
            else:
                return 'singleton', SINGLETONS / folder.name, f'mandatory-list match {sib.name}, but theta differs → no promotion, both stay singletons'

    # Default: Singleton in singletons/
    target = SINGLETONS / folder.name
    return 'singleton', target, 'no DB hit + no cluster match'


def do_render(knot_dir: Path, dry_run: bool = False) -> bool:
    """Run render_knot_blender.py on each .fseries in folder."""
    fseries_files = list(knot_dir.glob('*.fseries'))
    if not fseries_files:
        return False
    ok = True
    for fs in fseries_files:
        stem = fs.stem
        jpg = knot_dir / f'{stem}_render.jpg'
        if jpg.exists():
            print(f'  RENDER skip (exists): {stem}')
            continue
        if dry_run:
            print(f'  [DRY] RENDER: {fs} -> {jpg}')
            continue
        cmd = [
            BLENDER, '-b', str(RENDER_TEMPLATE),
            '--python', str(RENDER_SCRIPT), '--',
            '--fseries', str(fs),
            '--jpg', str(jpg),
        ]
        print(f'  RENDER: {stem}')
        rc = subprocess.run(cmd, capture_output=True).returncode
        if rc != 0 or not jpg.exists():
            ok = False
            print(f'    FAIL rc={rc}')
    return ok


def process_folder(folder: Path, render: bool = True, dry_run: bool = False) -> dict:
    folder = folder.resolve()
    print(f'\n=== {folder} ===')
    kind, target, reason = determine_target(folder, dry_run=dry_run)
    print(f'  Decision: {kind}  ({reason})')
    print(f'  Target:   {target}')

    result = {'folder': str(folder), 'kind': kind, 'target': str(target), 'reason': reason}

    if kind == 'own_knots_skip':
        # Just render if needed
        if render:
            do_render(folder, dry_run=dry_run)
        return result

    # No-op if already in target location
    if target.resolve() == folder.resolve():
        print(f'  → already at target, no move')
        if render and kind not in ('knotinfo',):
            do_render(folder, dry_run=dry_run)
        elif render and kind == 'knotinfo':
            print(f'  RENDER skip: KnotInfo hit ({folder.name})')
        return result

    if kind == 'noop':
        return result

    # Move folder
    old_name = folder.name
    if kind == 'cluster_promote':
        # Need to:
        #  1. Create new family folder = target
        #  2. Move BOTH the current folder + sibling singleton into it
        sib = None
        c, det, sym, chir = fingerprint_from_pflicht(load_pflicht_data(find_pflicht_complete(folder)))
        sib = find_singleton_match(c, det, sym, chir, exclude=folder)
        if not sib:
            print('  ERROR: sibling not found')
            return result
        if not dry_run:
            target.mkdir(parents=True, exist_ok=True)
        # Move folder
        new_path = target / folder.name
        if not dry_run:
            shutil.move(str(folder), str(new_path))
        print(f'  -> {new_path}')
        # Move sibling
        sib_new = target / sib.name
        if not dry_run:
            shutil.move(str(sib), str(sib_new))
        print(f'  -> {sib_new}')
        result['promoted_with'] = str(sib)
        result['target_folder'] = str(target)
        if render:
            do_render(new_path, dry_run=dry_run)
            do_render(sib_new, dry_run=dry_run)
        return result

    # Standard rename + move
    if kind == 'knotinfo':
        # bekannt/<klassisch>[__<original>] with interior file rename
        if not dry_run:
            target.parent.mkdir(parents=True, exist_ok=True)
        # First rename interior files within OLD location
        old_name_inside = folder.name
        new_name_inside = target.name
        if old_name_inside != new_name_inside:
            if not dry_run:
                for f in list(folder.iterdir()):
                    if f.is_file() and f.name.startswith(old_name_inside):
                        nf = folder / (new_name_inside + f.name[len(old_name_inside):])
                        f.rename(nf)
                # Then move folder
                shutil.move(str(folder), str(target))
                # Patch content refs
                patch_content_refs(target, old_name_inside, new_name_inside, dry_run=False)
            else:
                print(f'  [DRY] interior rename + move')
        else:
            if not dry_run:
                shutil.move(str(folder), str(target))
            else:
                print(f'  [DRY] move')
        # Skip render for KnotInfo hits (doctrine)
        if render:
            print(f'  RENDER skip: KnotInfo hit ({target.name})')
    elif kind == 'failure':
        if not dry_run:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(folder), str(target))
        if render:
            do_render(target if not dry_run else folder, dry_run=dry_run)
    elif kind == 'cluster':
        if not dry_run:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(folder), str(target))
        if render:
            do_render(target if not dry_run else folder, dry_run=dry_run)
    elif kind == 'singleton':
        if not dry_run:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(folder), str(target))
        if render:
            do_render(target if not dry_run else folder, dry_run=dry_run)

    return result


CLUSTER_RE = re.compile(r'^(S_\d+_d\d+_|c\d+_d\d+_)')  # S_ classes + legacy

def batch_process(root: Path, render: bool = True, dry_run: bool = False) -> list[dict]:
    """Iterate all knot-folders under root that need processing.
    Skips folders already in their correct location (bekannt/, failure/, family/).
    """
    candidates = []
    for p in sorted(root.rglob('*.fseries')):
        if not p.is_file():
            continue
        if any(part.startswith('.') for part in p.relative_to(root).parts):
            continue  # hidden
        if 'pipeline_failures' in p.parts:
            continue  # already in failures
        if 'bekannt' in p.parts:
            continue  # already KnotInfo-classified
        if 'singletons' in p.parts:
            continue  # already singleton
        # Skip if any parent dir is a family-cluster (c<c>_d<det>_<sym>_<chir>)
        kdir = p.parent
        if any(CLUSTER_RE.match(part) for part in p.relative_to(root).parts):
            continue  # already in family folder
        if kdir == root:
            continue
        candidates.append(kdir)

    candidates = list(dict.fromkeys(candidates))
    print(f'{len(candidates)} candidates to process.')
    results = []
    for kd in candidates:
        try:
            r = process_folder(kd, render=render, dry_run=dry_run)
            results.append(r)
        except Exception as e:
            print(f'  ERROR: {e}')
            results.append({'folder': str(kd), 'error': str(e)})
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('folder', type=Path, nargs='?', help='Knot folder (or omit it and use --batch)')
    ap.add_argument('--batch', action='store_true', help='Iterate over all candidates')
    ap.add_argument('--root', type=Path, default=ERZ, help='Root for --batch (default: erzeugte_knoten/)')
    ap.add_argument('--no-render', action='store_true', help='Skip render')
    ap.add_argument('--dry-run', action='store_true', help='Execute nothing, only log')
    args = ap.parse_args()

    render = not args.no_render

    if args.batch:
        batch_process(args.root, render=render, dry_run=args.dry_run)
    elif args.folder:
        process_folder(args.folder, render=render, dry_run=args.dry_run)
    else:
        ap.error('Need either <folder> or --batch')


if __name__ == '__main__':
    main()
