# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
"""
pipeline_seal.py
================
Disables every knot-table / census access in Sage / snappy / spherogram,
so that the pipeline works under blind-test conditions.

Pipeline doctrine: only our own computations and pipeline outputs, no
table lookups. This is central for pipeline work at K>16, where tables do not
exist anyway: the pipeline must not lean on tables in a hidden way, and
an 8_11 test must give the same result as later K=18 tests.

Activation at the start of every pipeline script:

    import sys
    sys.path.insert(0, '/media/squeezy/kauffman_storage/freeze_manifest')
    from pipeline_seal import seal_all
    seal_all()

What is blocked:
  1. Sage `sage.knots.knotinfo.KnotInfo` (demo stub enum with K0_1..K7_2 + links)
  2. Sage `sage.databases.knotinfo_db.KnotInfoDataBase` (`create_filecache`,
     `read`, `read_column_dict`, `read_row_dict`, `knot_list`, `link_list`)
  3. snappy `Manifold(name)` for census / knot-table names (e.g. 'm003',
     '8_11', 'K8a9', 'o9_42258', '8a9'). Triangulation isosigs (e.g.
     'kLLPwAPkccdefghjijjhggaaaaaaof') stay allowed; so do PD code lists,
     DT code strings and file paths.
  4. snappy `Manifold.identify()` (based on the census DB)
  5. snappy census classes (`OrientableCuspedCensus`, `LinkExteriors`,
     `CensusKnots` etc.), replaced by sealed proxies
  6. spherogram `Link(name)` for census names; PD code / DT code stay allowed

What is NOT blocked:
  - our own pipeline computations (volume, polynomials from PD/braid, HFK, Khovanov)
  - triangulation isosigs (snappy's internal geometric format, not a
    knot-table lookup)
  - PD code, DT code, Gauss code, oriented Gauss code, braid word
  - Marin data (`database_cubic_hecke`) for the CHA Markov trace: algebraic
    structure constants, not a knot table
  - KnotJob algorithms (Khovanov, s-invariant etc.); the knot tables
    (Rolfsen.kjb) are already sealed via `knotjob_SEQUESTERED_BLINDTEST/`

As of: 2026-05-08
"""
import re


class PipelineSealViolation(RuntimeError):
    pass


# Table / census name patterns. Strings that match one of these patterns
# are treated as lookup input and blocked.
CENSUS_NAME_PATTERNS = [
    re.compile(r'^[mstv]\d+$'),         # m003, s116, v123, t456
    re.compile(r'^o\d+_\d+$'),          # o9_42258 (OrientableCuspedCensus)
    re.compile(r'^K\d+[an]\d+$'),       # K8a9 (KnotInfo knot ID)
    re.compile(r'^L\d+[an]\d+(_\d+)?$'),# L4a1, L4a1_0 (LinkInfo)
    re.compile(r'^\d+_\d+$'),           # 8_11 (Rolfsen)
    re.compile(r'^\d+[an]\d+$'),        # 8a9 (HTW-Style)
    re.compile(r'^[Tt]refoil', re.IGNORECASE),
    re.compile(r'^[Ff]igure[- ]?[Ee]ight'),
    re.compile(r'^[Uu]nknot$'),
]


def looks_like_census_name(spec):
    if not isinstance(spec, str):
        return False
    return any(p.match(spec.strip()) for p in CENSUS_NAME_PATTERNS)


def _blocked_factory(message):
    def _blocked(*args, **kwargs):
        raise PipelineSealViolation(message)
    return _blocked


def seal_sage_knotinfo():
    """Block Sage's KnotInfo demo stub and KnotInfoDataBase."""
    try:
        import sage.knots.knotinfo as ki
        import sage.databases.knotinfo_db as kdb
    except ImportError:
        return False

    blocked = _blocked_factory(
        "Pipeline seal: Sage KnotInfo table access blocked. "
        "Use your own computation."
    )

    class _SealedKnotInfo:
        def __getattr__(self, name):
            blocked()
        def __getitem__(self, key):
            blocked()
        def __iter__(self):
            blocked()
        def __call__(self, *a, **k):
            blocked()
        def __repr__(self):
            return '<KnotInfo: SEALED by pipeline_seal>'

    ki.KnotInfo = _SealedKnotInfo()

    for method in ('create_filecache', 'read', 'read_column_dict',
                   'read_row_dict', 'read_num_knots', 'knot_list', 'link_list'):
        if hasattr(kdb.KnotInfoDataBase, method):
            setattr(kdb.KnotInfoDataBase, method, blocked)
    return True


def seal_snappy():
    """Block snappy census lookups and Manifold.identify()."""
    try:
        import snappy
    except ImportError:
        return False

    original_init = snappy.Manifold.__init__

    def sealed_init(self, spec=None, *args, **kwargs):
        if looks_like_census_name(spec):
            raise PipelineSealViolation(
                f"Pipeline-Seal: snappy.Manifold('{spec}') census lookup blocked. "
                "Use a PD code, DT code or triangulation isosig."
            )
        original_init(self, spec, *args, **kwargs)

    snappy.Manifold.__init__ = sealed_init

    snappy.Manifold.identify = _blocked_factory(
        "Pipeline-Seal: snappy.Manifold.identify() blockiert (Census-DB-Lookup)."
    )

    for census_attr in ('OrientableClosedCensus', 'NonorientableClosedCensus',
                        'OrientableCuspedCensus', 'NonorientableCuspedCensus',
                        'LinkExteriors', 'CensusKnots',
                        'AlternatingKnotExteriors', 'NonalternatingKnotExteriors',
                        'HTLinkExteriors'):
        if hasattr(snappy, census_attr):
            class _SealedCensus:
                _name = census_attr
                def __init__(self_inner, name):
                    self_inner._name = name
                def _err(self_inner):
                    raise PipelineSealViolation(
                        f"Pipeline-Seal: snappy.{self_inner._name} census DB blocked."
                    )
                def __getitem__(self_inner, k): self_inner._err()
                def __iter__(self_inner):       self_inner._err()
                def __call__(self_inner, *a, **k): self_inner._err()
                def __len__(self_inner):        self_inner._err()
                def __contains__(self_inner, x): self_inner._err()
            setattr(snappy, census_attr, _SealedCensus(census_attr))
    return True


def seal_spherogram():
    """Block spherogram.Link(name) for census names."""
    try:
        import spherogram
    except ImportError:
        return False

    original_link_init = spherogram.Link.__init__

    def sealed_link_init(self, *args, **kwargs):
        spec = args[0] if args else kwargs.get('crossings', kwargs.get('name'))
        if looks_like_census_name(spec):
            raise PipelineSealViolation(
                f"Pipeline-Seal: spherogram.Link('{spec}') census lookup blocked. "
                "Use a PD code list or a DT code."
            )
        original_link_init(self, *args, **kwargs)

    spherogram.Link.__init__ = sealed_link_init
    return True


def seal_all(verbose=True):
    """Activate all three seals. Envs that cannot be imported are skipped."""
    results = {
        'sage_knotinfo': seal_sage_knotinfo(),
        'snappy':        seal_snappy(),
        'spherogram':    seal_spherogram(),
    }
    if verbose:
        for name, ok in results.items():
            status = 'sealed' if ok else 'env not loaded (skipped)'
            print(f'pipeline_seal: {name}: {status}')
    return results


if __name__ == '__main__':
    seal_all()
