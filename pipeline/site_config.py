# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
"""Where the pipeline finds its own modules, external programs and data.

Every location the pipeline uses is a name in this module; no code outside it contains an absolute path
(comments and docstrings may still name the authors' folders, as documentation).

Two profiles:
- portable (default): own modules from this folder, programs from PATH, data below the repository root.
- authors (KNOT_PIPELINE_PROFILE=authors): exactly the locations used in the authors' runs of 2024-2026, listed in AUTHORS.
  Every published result was computed with these; the release gate (tools/path_equivalence.py) proves that, with the authors'
  profile, each released script is the same program as the one that produced those results.

Any single name can be overridden by the environment variable KNOT_PIPELINE_<NAME>, in either profile.
"""
import os, shutil, sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

_K = '/media/squeezy/kauffman_storage'
_CONDA = '/home/squeezy/miniforge3'
AUTHORS = {
    'ROOT': _K,
    'FREEZE_DIR': _K + '/freeze_manifest',
    'PIPELINE_DIR': _K + '/freeze_manifest/pipeline_scripts',
    'TOOLS_DIR': _K + '/claude_setup_v2/tools',
    'TOPO_GATE_DIR': _K + '/maschinen_ideen/_topo_gate_audit',
    'HECKE_DIR': _K + '/maschinen_ideen/_mutation_separation_build_2026_06_16',
    'KAWAUCHI_SIEVE_DIR': _K + '/konversationen/kawauchi_20260725/sieb',
    'MACHINE_YARD': _K + '/maschinen_ideen',
    'GENERATED_KNOTS_DIR': _K + '/erzeugte_knoten',
    'WEB_EXPORT_DIR': _K + '/web_export',
    'THETA_DIR': _K + '/external_tools/bar_natan_theta',
    'KNOTINFO_INDEX': _K + '/external_tools/knotinfo_index/index_c_det_sig.json',
    'KNOTINFO_RECORDS': _K + '/external_tools/knotinfo_index/records.json',
    'FREMLIN_CACHE': _K + '/external_tools/fremlin_cache/fremlin_cache.json',
    'CAMERA': _K + '/freeze_manifest/camera_131072.npz',
    'RENDER_TEMPLATE': _K + '/freeze_manifest/render_template.blend',
    'KNOTJOB_JAR': '/home/squeezy/tools/knotjob/KnotJob/KnotJob.jar',
    'JAVA_BIN': _CONDA + '/envs/sage_stable/bin/java',
    'SAGE_BIN': _CONDA + '/envs/sage_stable/bin/sage',
    'PYTHON_MAIN': _CONDA + '/envs/gnottheorie_env/bin/python3',
    'PYTHON_BASE': _CONDA + '/bin/python',
    'CONDA_BIN': _CONDA + '/bin/conda',
    'BLENDER': '/home/squeezy/.local/bin/blender',
}

def _which(prog):
    return shutil.which(prog) or prog

PORTABLE = {
    'ROOT': REPO,
    'FREEZE_DIR': HERE,
    'PIPELINE_DIR': HERE,
    'TOOLS_DIR': HERE,
    'TOPO_GATE_DIR': HERE,
    'HECKE_DIR': HERE,
    'KAWAUCHI_SIEVE_DIR': HERE,
    'MACHINE_YARD': os.path.join(REPO, 'machines'),
    'GENERATED_KNOTS_DIR': os.path.join(REPO, 'output'),
    'WEB_EXPORT_DIR': os.path.join(REPO, 'web_export'),
    'THETA_DIR': os.path.join(REPO, 'third_party', 'bar_natan_theta'),
    'KNOTINFO_INDEX': os.path.join(REPO, 'data', 'knotinfo_index', 'index_c_det_sig.json'),
    'KNOTINFO_RECORDS': os.path.join(REPO, 'data', 'knotinfo_index', 'records.json'),
    'FREMLIN_CACHE': os.path.join(REPO, 'data', 'fremlin_cache', 'fremlin_cache.json'),
    'CAMERA': os.path.join(REPO, 'data', 'camera_131072.npz'),
    'RENDER_TEMPLATE': os.path.join(REPO, 'third_party', 'render_template.blend'),
    'KNOTJOB_JAR': os.path.join(REPO, 'third_party', 'KnotJob', 'KnotJob.jar'),
    'JAVA_BIN': _which('java'),
    'SAGE_BIN': _which('sage'),
    'PYTHON_MAIN': sys.executable,
    'PYTHON_BASE': sys.executable,
    'CONDA_BIN': _which('conda'),
    'BLENDER': _which('blender'),
}
assert set(AUTHORS) == set(PORTABLE)

PROFILE = os.environ.get('KNOT_PIPELINE_PROFILE', 'portable')
if PROFILE not in ('portable', 'authors'):
    raise ValueError(f"KNOT_PIPELINE_PROFILE must be 'portable' or 'authors', not {PROFILE!r}")
_BASE = AUTHORS if PROFILE == 'authors' else PORTABLE

for _name, _value in _BASE.items():
    globals()[_name] = os.environ.get('KNOT_PIPELINE_' + _name, _value)

def describe():
    """One line per name: value and where it came from (for logs)."""
    out = []
    for n in AUTHORS:
        src = 'env' if ('KNOT_PIPELINE_' + n) in os.environ else PROFILE
        out.append(f'{n:20s} {globals()[n]}  [{src}]')
    return '\n'.join(out)

if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--get':   # for shell scripts: PIPE=$(python3 site_config.py --get PIPELINE_DIR)
        print(globals()[sys.argv[2]] if sys.argv[2] in AUTHORS else sys.exit(f'unknown name {sys.argv[2]}'))
    else:
        print(describe())
