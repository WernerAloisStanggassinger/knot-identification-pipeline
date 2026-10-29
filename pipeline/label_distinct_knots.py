# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
"""label_distinct_knots.py — distinct-knot bookkeeping in corpus_db.json.

Werner directive 2026-05-30: „nicht löschen, sondern die Distinkt-Knoten-Identität
explizit machen ... Buchhaltung muß passen." (EN: do not delete, but make the distinct-knot identity
explicit ... the bookkeeping must fit.) The many machine-parameter records
(different sh/tw/helix) are NOT duplicates to be deleted — they are the C_n/S_n
family multiplicity (anchor-2 evidence). They ALL remain; every record
gets:
  distinct_knot_id          — canonical name of the realised distinct knot
  n_realisierungen          — how many records realise the same distinct knot
  is_distinct_representative — exactly one record per distinct knot = True
  realizes_class            — known | novel_hyperbolic | novel_nonhyp

Identity = manifold level (det,sig,alex,theta@22/7,vol@6), NOT isosig.
Idempotent. Backup before writing. meta bookkeeping must balance (assert).

Call:  python3 label_distinct_knots.py [--apply]
Without --apply only a dry-run report.
"""
import sys, json, shutil, re
import site_config
from pathlib import Path
from collections import defaultdict, Counter

from datetime import date

sys.path.insert(0, site_config.TOOLS_DIR)
from span_lower_bound import (span_in_t as _span_in_t,   # calibrated, not copied
                                  alternating_status as _alt_status,
                                  beyond_every_enumeration as _jenseits,
                                  ENUM_ALTERNATING)
# Mirror detection: ONE place, one truth (Disc #808). The theta route for sigma = 0
# is hardened against three counter-cases (conway_cw/conway_kt are NOT a mirror pair) and
# the parser calibrated against sympy — maschinen_ideen/KNOWN_DOPPELZAEHLUNG_20260919/.
from theta_mirror import are_mirrors
# Alexander normal form (Disc #809): take the ±-convention out of the identity key.
# Counter-case: maschinen_ideen/KNOWN_DOPPELZAEHLUNG_20260919/alex_normalform_gegenfall.py
# — merges EXACTLY the one documented pair, no other (509 -> 508 groups).
from alexander_normal_form import canonical as alex_kanonisch
# One place, one truth (doctrine l.175): the genus and braid-index ingredients come
# from the calibrated tool, they are NOT rebuilt here.
from yamada_c_bound import (genus_lower_from_alexander as _g_alex,
                               mfw_from_homfly as _mfw)

ROOT = Path(site_config.ROOT)
DB = ROOT / "freeze_manifest/pipeline_scripts/corpus_db.json"
# Maguire-Extended-Retractions (c14-17, Disc #512/#516): distinct_knot_ids whose
# novelty is refuted against Maguire's c<=17 hyp-Jones census (volume-confirmed,
# is_isometric_to). SEPARATE from crosscheck_knotinfo (KnotInfo c<=13). Field:
# realizes_tabulated_extended. Durable + reproducible (no hand fix that rots).
RETRACT_FILE = ROOT / "freeze_manifest/pipeline_scripts/maguire_retractions.json"
# SECOND source (2026-09-04): Burton's complete census. Own file and own
# record field, so that the provenance of a retraction stays readable — whoever asks later
# WHY a knot is no longer novel must see the source, not only the result.
# Occasion: the Burton run of 10.07.2026 found 17 hits, 14 of them were still standing seven weeks
# later as novel_hyperbolic in corpus_db — among them `sundown`(=trainstation,
# built in steel) with c=18 and a NOVELTY label on the public site.
BURTON_FILE = ROOT / "freeze_manifest/pipeline_scripts/burton_retractions.json"
# FREEZE REGISTRY (2026-09-05, Werner: "isosig-registry bauen" (EN: build the isosig registry) — Disc #755/#758):
# distinct_knot_id hung on the cluster FOLDER and broke at every re-sorting
# (S_18_d27_D3_chiral_sig6 -> _m, 04.09.). The stored isosig field is NOT suitable as a
# substitute (measured 05.09.: 66 distinct knots carry SEVERAL isosigs — triangulation,
# not canonical signature; 68 records without; a canonical signature would be mirror-blind).
# Hence: the FIRST assignment is frozen here, key = ident() tuple (the
# true identity). canon_id() now only coins a name for NEW identities.
REGISTRY_FILE = ROOT / "freeze_manifest/pipeline_scripts/distinct_id_registry.json"
# PROVEN identities (2026-08-09): for non-hyperbolic knots the key has no
# geometric discriminator, and canon_id() falls back to the alphabetically first STEM —
# i.e. to a machine parameter name instead of an identity. This file carries per group
# the POSITIVE proof (braid conjugation), never an inference from equality of invariants.
IDENT_FILE = ROOT / "freeze_manifest/pipeline_scripts/proven_identities.json"


def load_burton_retractions():
    if BURTON_FILE.exists():
        return json.load(open(BURTON_FILE)).get("retractions", {})
    return {}


# Census ABSENCE (c19 run 06.09.2026): target det computed completely in the c<=19 census,
# no volume hit => knot NOT in the census => c>=20 PROVEN (census complete,
# knot hyperbolic prime). Raises ONLY the lower bound to 20 — the novelty bar
# stays 21 (c=20 gap, Disc #753), classes NEVER change through this.
ABSENCE_FILE = ROOT / "freeze_manifest/pipeline_scripts/burton_absence_c19.json"


def load_burton_absence():
    if ABSENCE_FILE.exists():
        d = json.load(open(ABSENCE_FILE))
# Union of stage-B absences (9), tighten_65 absences (31) and
        # stubborn_crusher absences (28, 06.09. evening: the last non-geometric
        # census blockers computed rigorously, 0 hits). Absence targets may meanwhile be
        # novel (genus promotion) — then the absence changes nothing about cl,
        # but counts as arrived.
        return (set(d.get("absent_c_le_19_rigoros_UNRESOLVED", []))
                | set(d.get("absent_c_le_19_rigoros_TIGHTEN_20260906", []))
                | set(d.get("absent_c_le_19_rigoros_CRUSHER_20260906", []))
# det hunt 07.-08.09.2026: target det does not occur in the COMPLETE c<=19 census
                # (matches=0, all target numbers exact, detfails=0) — stronger
                # than the vol route, no volume window needed. Details in _meta.
                | set(d.get("absent_c_le_19_rigoros_DETJAGD_20260908", [])))
    return set()


def load_retractions():
    try:
        return json.load(open(RETRACT_FILE)).get("retractions", {})
    except FileNotFoundError:
        return {}


DISTINCT_FILE = ROOT / "freeze_manifest/pipeline_scripts/proven_distinctions.json"


def load_distinctions():
    """stem -> entry from proven_distinctions.json (PROVEN distinctness).

    The distinct key is blind to mutation: a mutant pair agrees in
    det, sig, Alexander, theta and vol and lands in ONE group. Where the
    complements are demonstrably NOT isometric (Gordon-Luecke), this
    entry splits the stem off into its own distinct_knot_id."""
    try:
        return json.load(open(DISTINCT_FILE)).get("distinctions", {})
    except FileNotFoundError:
        return {}


# stem -> case name; loaded once, because ident() runs per record
PROVEN_DISTINCT = {v["abgespalten"]: k for k, v in load_distinctions().items()
                   if v.get("abgespalten")}


def ident(m):
    vol = m.get('vol')
    vol_key = round(vol, 6) if isinstance(vol, (int, float)) else None
    # Disc #538: for NON-hyperbolic knots `vol` is not an invariant but snappy garbage
    # (sum of the JSJ pieces for connected sums, companion geometry for satellites). As a
    # key component this garbage then SEPARATES identical knots: 5_2#m5_2 lay 3x as its
    # own "distinct knot" in the corpus, because the garbage volume
    # differed in the 7th decimal (5.656247650 / 5.656244176 / 5.656244177), 4_1#4_1 2x. So: non-hyperbolic
    # => vol out of the key, what remains is det+sig+Alexander+theta (for non-hyperbolic
    # knots the load-bearing invariants).
    if m.get('is_hyperbolic') is False:
        vol_key = None
    # 2026-08-09: the condition above alone is NOT enough. A record with is_hyperbolic=None and
    # vol=0.0 kept volk=0.0, its identical twin with is_hyperbolic=False got volk=None
    # -> (0.0 != None) separated two realisations of the SAME knot. Measured on
    # st_dn5_s11_x5z2d55 / st_dn6_s20_x4z2d45: both proven left-handed T(5,6) (braid conjugation),
    # bit-identical in det/sig/alex/theta, torn apart solely by the unverified flag.
    # Rule: vol belongs in the key only if it is a REAL hyperbolic volume —
    # same threshold as is_hyp() (>0.5), otherwise it is snappy garbage (Disc #538) or simply 0.
    # A missing/unverified flag must never act as a distinguishing feature.
    if isinstance(vol, (int, float)) and vol <= 0.5:
        vol_key = None
    # Disc #809: NOT str(alex) — Δ is only defined up to ±t^k, and the corpus
    # carried the same Δ once as "-x^18 + 13x^17 - ..." and once as
    # "x^18 - 13x^17 + ...". That separated poly_s08 from kf332_polys08_sh+0.5_tw-1, although
    # Jones and HOMFLY are bit-identical, theta = 0 on both sides and the
    # complements are ISOMETRIC. The normal form only flips the overall sign and
    # sorts — it DELIBERATELY does NOT normalise the t^k shift (not a
    # single such case in the corpus; normalise only what is measured).
    basis = (str(m.get('det')), str(m.get('sig')), alex_kanonisch(m.get('alex')),
             str(m.get('theta_at_22_7__21_13')), vol_key)
    # PROVEN distinctness (proven_distinctions.json): the key is
    # mutant-blind, the isometry of the complements is not. A stem listed here
    # gets its own key element and thereby its own ID.
    fall = PROVEN_DISTINCT.get(m.get('stem'))
    if fall:
        return basis + ('proven-distinct:' + fall,)
    return basis


def load_identities():
    """stem -> proven identity name (see IDENT_FILE)."""
    try:
        raw = json.load(open(IDENT_FILE)).get("identities", {})
    except FileNotFoundError:
        return {}
    out = {}
    for name, e in raw.items():
        for s in e.get("stems", []):
            out[s] = name
    return out


PROVEN_IDENT = load_identities()


def is_composite_id(m):
    """Is the identity a COMPOSITE/satellite description instead of a
    table entry? (Disc #538) `5_2#m5_2`, `4_1#4_1`, `Whitehead-Double(3_1)` are
    known objects, but appear in NO prime-knot table -> they belong in the bucket
    'known non-hyp (torus/sum)', not in 'known (table)'.
    Otherwise a connected sum would count as a table hit and inflate distinct_known."""
    s = str(m.get('crosscheck_knotinfo') or '')
    # '#' = connected sum · 'Double'/'satellite' = satellite · 'C_{p,q}(K)' = cabling
    return '#' in s or 'Double' in s or 'satellite' in s.lower() or s.startswith('C_{')


def is_known(m):
    # crosscheck_knotinfo = KnotInfo c<=13; realizes_tabulated_extended = Maguire c14-17
    # (Disc #512/#516). Both sources kept separate, but both = tabulated (not novel).
    if is_composite_id(m):
        return False
    return bool(m.get('knotinfo_id') or m.get('crosscheck_knotinfo')
                or m.get('crosscheck_snappy') or m.get('realizes_tabulated_extended')
                or m.get('realizes_tabulated_burton'))


def c_lower_bound(m):
    """Proven LOWER bound for the minimal crossing number — or None.

    Werner 2026-09-04: "untere schranke = k 20. no novelty." (EN: lower bound = c 20. no novelty.) The prime knots are enumerated up to
    AND INCLUDING c=20 (Thistlethwaite, AGT 25 (2025) 329-344). A knot
    for which we cannot PROVE c>=21 is no novelty candidate — it is a
    knot whose table entry we merely do not have (c=20 not public) or have not
    yet searched (c=19 run open). That is a state of our receipts, not a
    property of the knot.

    Basis: span V(K) <= c(K) for EVERY knot (Kauffman-Murasugi-Thistlethwaite),
    with equality exactly for reduced alternating diagrams.
    ⚠️ One-sided: a SMALL bound proves nothing — span=12 does not mean "c is small",
    but "this route is not enough here".

    🔴 NO COPY OF THE FORMULA HERE. The first draft on 04.09. copied the computation
    and forgot the branch for the CONSTANT term — 25 knots
    got a bound too small by 1, and because "too small" here means "more cautious", it
    would have looked like diligence. Caught only because the copy was checked against the
    original. One place, one truth: span_lower_bound.py is calibrated against 250
    alternating records (span == c, zero deviations), this function calls it.
    """
    j = m.get('jones')
    if not isinstance(j, str) or not j.strip():
        return None
    return _span_in_t(j)


def hfk_exact_genus(m):
    """Seifert genus ONLY if HFK-exactly documented — otherwise None.

    BEFUND_13 lesson (Galway return 07.08.): the corpus field `genus` can be a
    bound-as-value. It counts here only if the pflicht record contains
    knot_floer_homology.seifert_genus AND it agrees with the corpus field
    (Ozsvath-Szabo: HFK detects the Seifert genus — that is a theorem, not an
    algorithm artefact). Receipt of 05.09.: kf373 carries genus=13 WITHOUT an HFK receipt
    and rightly fails here.
    """
    g = m.get('genus')
    if not isinstance(g, int):
        return None
    pf = m.get('pflicht_path')
    if not pf or not Path(pf).exists():
        return None
    try:
        d = json.load(open(pf))
        sg = (d.get('pflichtliste', {}).get('ALGEBRA', {})
              .get('knot_floer_homology') or {}).get('seifert_genus')
    except Exception:
        return None
    return g if sg == g else None


def _nontrivial_receipted(m):
    """The braid-index floor s>=2 holds only for a PROVABLY non-trivial knot
    (b(K)=1 <=> unknot). Non-triviality is demanded as DOCUMENTED, never assumed:
    det != 1 (Alexander(-1) != +-1) or verified hyperbolicity."""
    d = m.get('det')
    try:
        if d is not None and abs(int(d)) != 1:
            return True
    except (TypeError, ValueError):
        pass
    return m.get('is_hyperbolic') is True or bool(m.get('hyperbolicity_proof'))


TURAEV_FILE = ROOT / "freeze_manifest/pipeline_scripts/turaev_bounds.json"
_TURAEV_APPLIED = set()


def _turaev_bounds():
    """Route-3 receipts {stem: {span_V, w_HF, w_KH_Q, receipt}} — empty if the file does not exist."""
    if not TURAEV_FILE.exists():
        return {}
    return json.load(open(TURAEV_FILE)).get('bounds', {})


def c_lower_bound_with_via(m):
    """Best proven lower bound + provenance, as (value, via) — or (None, None).

    TWO independent routes (AUFTRAG_KREUZUNGSSCHRANKEN_20260903 from the Potsdam laptop,
    genus route built in on Werner's go 05.09.):
      (1) span V(K) <= c(K)  — Kauffman-Murasugi-Thistlethwaite, calibrated tool.
      (2) c >= 2g+1          — Seifert algorithm on a minimal diagram:
          g <= (c - s + 1)/2 with s Seifert circles, and s >= 2 for every non-
          trivial knot (Yamada 1987: min. Seifert circles = braid index >= 2)
          ==> c >= 2g + s - 1 >= 2g + 1. Holds only with HFK-EXACT g (hfk_exaktes_genus).
    The larger bound wins; via names the winner honestly.
    """
    best, via = None, None
    s = c_lower_bound(m)
    if s is not None:
        best, via = s, 'span V(K) <= c (Kauffman-Murasugi-Thistlethwaite)'

    # --- Route 2, since 09.09.2026 in THREE stages instead of one --------------
    # c >= 2g + s - 1 is monotone in BOTH arguments, so a LOWER bound may stand
    # anywhere and the result remains a lower bound.
    #   g:  HFK-exact (theorem)  >  Alexander degree (g >= span(Delta)/2)
    #   s:  measured braid index (MFW from HOMFLY)  >  floor s>=2 (Yamada: s = b(K),
    #       and b(K)=1 <=> unknot — a theorem, not an estimate)
    # Before, only (HFK-exact, s=2) stood here. That left ingredients unused that
    # had long been in the record: 7 of the 54 'unresolved_c_le_20' leave the
    # c=20 gap with them, six of them without a single new computation (Disc #773).
    g_hfk = hfk_exact_genus(m)
    g_alx = _g_alex(m.get('alex'))
    if isinstance(g_hfk, int):
        g, g_via = g_hfk, 'g HFK-exact, Ozsvath-Szabo'
    elif isinstance(g_alx, int):
        g, g_via = g_alx, 'g >= span(Alexander)/2 (Seifert inequality)'
    else:
        g, g_via = None, None

    n_mfw = _mfw(m.get('homfly'))
    if isinstance(n_mfw, int) and n_mfw >= 2:
        n, n_via = n_mfw, 'braid index >= %d (MFW from HOMFLY)' % n_mfw
    elif g is not None and _nontrivial_receipted(m):
        n, n_via = 2, 'braid index >= 2 (b(K)=1 <=> unknot; non-triviality documented)'
    else:
        n, n_via = None, None

    if g is not None and n is not None:
        candidate = 2 * g + n - 1
        if best is None or candidate > best:
            best = candidate
            via = 'c >= 2g+n-1 (Seifert algorithm/Yamada; %s; %s)' % (g_via, n_via)
    # --- Route 3 (since 28.10.2026, Disc #943): Turaev genus ---------------------
    #   c(K) >= span V(K) + g_T(K)   (Champanerkar-Kofman arXiv:1406.1945v1 eq. (1))
    #   g_T(K) >= w_HF(K) - 1        (Lowrance arXiv:0709.0720v1 Thm 1.1)
    #   g_T(K) >= w_KH(K) - 1        (reduced Khovanov width, [CK] Cor. 4.5)
    # Widths are measured, not derivable from the record: they come from the receipt file
    # turaev_bounds.json (tools/turaev_bound.py + witnesses). Latch: the span stored there
    # must be the span of THIS record's Jones, else the entry belongs to another knot.
    tb = _turaev_bounds().get(str(m.get('stem') or ''))
    if tb:
        sv = c_lower_bound(m)
        if sv != tb['span_V']:
            raise AssertionError(
                f"turaev_bounds.json: span V {tb['span_V']} != record span {sv} at {m.get('stem')} — "
                f"the receipt does not belong to this record.")
        w = max(tb['w_HF'], tb.get('w_KH_Q', 1))
        candidate = int(tb['span_V'] + w - 1)
        _TURAEV_APPLIED.add(str(m.get('stem')))
        if best is None or candidate > best:
            best = candidate
            via = ('c >= span V + g_T, g_T >= w - 1 (Turaev, Champanerkar-Kofman arXiv:1406.1945v1 eq.(1); '
                   'w_HF=%d Lowrance arXiv:0709.0720v1 Thm 1.1; w_KH(Q)=%s [CK] Cor. 4.5; receipt %s)'
                   % (tb['w_HF'], tb.get('w_KH_Q', '-'), tb.get('receipt', 'turaev_bounds.json')))
    # Consistency latch: a lower bound above the reported UPPER bound c
    # would be a contradiction in the record — fail loudly, never smooth over.
    c_ub = m.get('c')
    if best is not None and isinstance(c_ub, int) and best > c_ub:
        raise AssertionError(
            f"c_lower_bound {best} ({via}) > c upper bound {c_ub} at {m.get('stem')} — "
            f"contradiction: genus or c of this record is wrong.")
    return best, via


# From which proven crossing number on a knot CAN be a candidate at all.
# 21, because up to 20 is enumerated completely. A number, no feeling.
# 🔴 Disc #814 (20.09.2026): this is only the FIRST of two limits. The ALTERNATING
# prime knots are published-enumerated up to c=23 (Rankin/Flint/Schermann, JKTR 13 (2004))
# and computed up to c=24 (Fontaine 2007) — verbatim from Thistlethwaite himself, AGT 25
# (2025) p.2: "This has indeed been accomplished quite dramatically for n <= 23
# [Rankin et al. 2004]". ENUM_ALTERNIEREND = 24 comes from span_untere_schranke.
ENUMERATION_LIMIT = 20


def is_hyp(m):
    # Disc #323: the connected-sum gate (enrich_non_hyperbolic) delivers is_hyperbolic=False
    # AUTHORITATIVELY for connected sums — vol alone lies there (snappy sums the
    # JSJ pieces to one volume). is_hyperbolic=False therefore overrules vol>0.5.
    # For True/None the old criterion stays: real volume>0.5 (an old,
    # unverified is_hyperbolic=True with vol=None must NOT wrongly promote a
    # non-hyperbolic knot to novel — regression protection).
    if m.get('is_hyperbolic') is False:
        return False
    v = m.get('vol')
    return isinstance(v, (int, float)) and v > 0.5


def canon_id(ms):
    # PROVEN identity beats every fallback (2026-08-09). If the file contradicts itself
    # (two names in ONE group), abort loudly instead of silently choosing one — a gate
    # that picks its result is not a gate.
    named = sorted({PROVEN_IDENT[s] for m in ms
                    if (s := m.get('stem')) in PROVEN_IDENT})
    if len(named) > 1:
        raise SystemExit(f"🔴 proven_identities.json contradictory: one group carries {named} "
                         f"(stems: {[m.get('stem') for m in ms]})")
    if named:
        return named[0]
    # true S-cluster (S_<c>_d...) as ID; NOT the bucket folders
    # 'bekannt'/'singletons_novelty'/'singletons_unknots' (those are not knot names)
    cfs = sorted(set(m.get('cluster_folder') for m in ms
                     if str(m.get('cluster_folder') or '').startswith('S_')))
    if cfs:
        return cfs[0]
    kis = sorted(set(str(m.get('knotinfo_id') or m.get('crosscheck_knotinfo'))
                     for m in ms if (m.get('knotinfo_id') or m.get('crosscheck_knotinfo'))))
    if kis:
        return kis[0]
    return sorted((m.get('stem') or '?') for m in ms)[0]


def registry_key(key):
    """ident() tuple -> deterministic JSON string (volk is rounded to 6 places,
    json.dumps of the same float is deterministic in Python)."""
    return json.dumps(list(key), ensure_ascii=False)


def load_registry():
    try:
        return json.load(open(REGISTRY_FILE))
    except FileNotFoundError:
        return {"_doc": ("Freeze registry for distinct_knot_id (Disc #755/#758). "
                         "Key = ident() tuple as JSON list [det, sig, alex, "
                         "theta@22/7;21/13, vol@6]. An ID once assigned stays — "
                         "re-sorting the cluster folders never changes it again. "
                         "NEVER edit by hand; the labeler registers new "
                         "identities itself. Conflicts break LOUDLY."),
                "frozen": {}}


def assign_ids(grp, frozen):
    """Assign IDs: frozen wins, new identities are coined + registered.

    Returns: list of the NEWLY registered (rk, entry). Mutates g['id'] and frozen.
    Conflicts (duplicate registry IDs, proven name contradicts frozen ID) -> SystemExit.
    """
    ids = [e['id'] for e in frozen.values()]
    if len(ids) != len(set(ids)):
        dup = sorted({i for i in ids if ids.count(i) > 1})
        raise SystemExit(f"🔴 distinct_id_registry corrupt: duplicate IDs {dup}")
    used = set(ids)  # ALL frozen IDs are taken — also those currently not in the corpus
    new = []
    open_ = []
    for g in grp:
        e = frozen.get(registry_key(g['key']))
        if e:
            # a proven name (proven_identities) must never contradict a frozen ID
            hat_proven = any((m.get('stem') in PROVEN_IDENT) for m in g['ms'])
            if hat_proven and g['base'] != e['id']:
                raise SystemExit(
                    f"🔴 Registry conflict: identity {g['key']} is frozen as '{e['id']}' "
                    f"but carries the proven name '{g['base']}' "
                    f"(proven_identities). Migrate the registry by hand, do not guess.")
            g['id'] = e['id']
        else:
            open_.append(g)
    base_counts = defaultdict(int)
    for g in open_:
        base_counts[g['base']] += 1
    for g in open_:
        if base_counts[g['base']] > 1:
            sig = g['rep'].get('sig')
            base = f"{g['base']}_sig{sig}"
        else:
            base = g['base']
        cand = base; k = 1
        while cand in used:
            k += 1; cand = f"{base}#{k}"
        used.add(cand)
        g['id'] = cand
        entry = {'id': cand, 'first_seen': date.today().isoformat(),
                   'witness_stem': g['rep'].get('stem')}
        frozen[registry_key(g['key'])] = entry
        new.append((registry_key(g['key']), entry))
    return new


def mirror_collapsed_count(reps, class_name='novel_hyperbolic'):
    """Distinct knots of a class after mirror merging.

    Disc #808 — THE GAP THAT WAS HERE: the condition used to read inline
    `si == -ei and si != 0`. That excluded sigma = 0 BY CONSTRUCTION, and theta
    was never read. Exactly the case from Disc #807 (kf350_spiral_amp0.3_freq4 <->
    S_18_d1_D3_chiral: sigma = 0 on both sides, complements isometric, only
    theta = -theta' shows the mirror) was the only one the function COULD
    not catch. The detection now lives in theta_spiegel.sind_spiegel() — ONE place, one
    truth, hardened against three counter-cases and calibrated against sympy.

    `klasse` selects the class (Disc #808: for 'known' this number did not exist,
    although KnotInfo and Burton list ONE representative per chirality pair — the
    asymmetry was itself the source of error, cf. Disc #779/#804).
    """
    group = [r for r in reps if r['_class'] == class_name]
    by_dv = defaultdict(list)
    for r in group:
        m = r['_rep']
        v = m.get('vol')
        by_dv[(str(m.get('det')), round(v, 4) if v else None)].append(m)
    drop = 0
    routes = Counter()
    for k, ms in by_dv.items():
        seen = set()
        for d in ms:
            for e in ms:
                if d is e or id(d) in seen or id(e) in seen:
                    continue
                weg = are_mirrors(d, e)
                if weg:
                    seen.add(id(d)); seen.add(id(e)); drop += 1
                    routes[weg] += 1
                    break
    return len(group) - drop, dict(routes)


def main(apply=False):
    db = json.load(open(DB))
    knots = db['knots']
    retractions = load_retractions()
    groups = defaultdict(list)
    for m in knots:
        groups[ident(m)].append(m)

    # 1. Base ID per group, stably sorted (for deterministic disambiguation)
    grp = []
    for key, ms in groups.items():
        rep = ms[0]
        grp.append({'key': key, 'ms': ms, 'base': canon_id(ms), 'rep': rep})
    grp.sort(key=lambda g: (g['base'], str(g['rep'].get('det')),
                            str(g['rep'].get('sig')), str(g['rep'].get('vol'))))
    # 2. Assign IDs: FREEZE REGISTRY wins (Disc #755/#758), only new
    #    identities are coined (suffix logic as before: _sig<sig>, then #k;
    #    fix 2026-06-19 against the mutated used counter lives on in assign_ids).
    registry = load_registry()
    frozen = registry.setdefault('frozen', {})
    newly_registered = assign_ids(grp, frozen)

    reps = []
    n_retracted_applied = 0
    n_burton_applied = 0
    n_absence_applied = 0
    burton = load_burton_retractions()
    absence = load_burton_absence()
    # det-hunt subset: different proof route (det absence instead of vol window) => own
    # via text; the evidence line must not cite the 06.09. vol runs (false bit).
    absence_det_hunt = (set(json.load(open(ABSENCE_FILE)).get(
        "absent_c_le_19_rigoros_DETJAGD_20260908", []))
        if ABSENCE_FILE.exists() else set())
    for g in grp:
        ms, did, rep = g['ms'], g['id'], g['rep']
        bt = burton.get(did)
        if bt:
            n_burton_applied += 1
            for m in ms:
                m['realizes_tabulated_burton'] = bt
        # Maguire-Extended-Retraction: distinct_knot_id tabulated against the c<=17 census.
        # Stamp the field on ALL records of the group, BEFORE is_known -> counts as known.
        rt = retractions.get(did)
        if rt:
            n_retracted_applied += 1
            for m in ms:
                m['realizes_tabulated_extended'] = rt
        if is_known(rep):
            cls = 'known'
        elif is_hyp(rep):
            # Werner 2026-09-04: "untere schranke = k 20. no novelty." (EN: lower bound = c 20. no novelty.)
            # A hyperbolic knot without a table hit is a novelty candidate ONLY if
            # c >= 21 is PROVEN. Otherwise it lies in the enumerated range, and
            # that we do not find it says something about our tables, not
            # about it.
            cl, cl_via = c_lower_bound_with_via(rep)
            if did in absence:
                n_absence_applied += 1
                if cl is None or cl < 20:
                    if did in absence_det_hunt:
                        cl, cl_via = 20, ('target det does not occur in the COMPLETE c<=19 census '
                                          '(det exactly via regina |Alexander(-1)|, all target numbers '
                                          'exact, matches=0/detfails=0; run work_detjagd_20260907, '
                                          '2026-09-07/08; det with 2 witnesses each bit-exact, Disc #762)')
                    else:
                        cl, cl_via = 20, ('not in the c<=19 census (Burton complete: 3-16 AND '
                                          '17-19 layer streamed, 352,151,858 hyp knots, '
                                          'detfails=0/vol-fails=0; runs work_c19_20260904 + '
                                          'low_stage_a/b, 2026-09-06)')
                    c_ub = rep.get('c')
                    assert not (isinstance(c_ub, int) and cl > c_ub), (
                        f"Absence bound {cl} > c upper bound {c_ub} at {rep.get('stem')} — "
                        f"contradiction: the knot would be in the census, where it is proven absent.")
                for m in ms:
                    m['burton_absent_c_le_19'] = True
            for m in ms:
                m['c_lower_bound_proven'] = cl
                m['c_lower_bound_via'] = cl_via
            # Disc #814: TWO limits. A novelty candidate is only one who exceeds both —
            # either c >= 25, or c >= 21 AND proven non-alternating. Alternating ones
            # with c = 21..24 lie in Rankin/Flint/Schermann's enumeration; that we do not have
            # their table is the same situation as with c = 20.
            _alt_state = _alt_status(rep, cl)
            for m in ms:
                m['alternating_proven'] = _alt_state
            if cl is None or cl <= ENUMERATION_LIMIT:
                cls = 'unresolved_c_le_20'
            elif _jenseits(rep, cl):
                cls = 'novel_hyperbolic'
            else:
                cls = 'unresolved_alternating_c_le_24'
        else:
            # Non-hyperbolic + no prime-table match ⇒ torus / satellite / connected-sum:
            # a KNOWN family, NEVER a novelty (doctrine: novelty = hyperbolic-prim + no-match).
            # The census (snappy/Fremlin) is hyperbolic-focused, so torus knots (T(3,7), T(6,7))
            # and connected sums (5_2#m5_2, 6×trefoil) miss it although they are classical.
            # Werner 2026-06-12 ("technisch wahr und trotzdem falsch" (EN: technically true and still wrong)). NOT novel_nonhyp.
            cls = 'known_nonhyp'
        for i, m in enumerate(ms):
            m['distinct_knot_id'] = did
            m['n_realisierungen'] = len(ms)
            m['is_distinct_representative'] = (i == 0)
            m['realizes_class'] = cls
        reps.append({'id': did, 'n': len(ms), '_class': cls, '_rep': rep})

    n_dist = len(reps)
    n_known = sum(1 for r in reps if r['_class'] == 'known')
    n_nh = sum(1 for r in reps if r['_class'] == 'novel_hyperbolic')
    n_kn = sum(1 for r in reps if r['_class'] == 'known_nonhyp')  # known non-hyp (torus/sum)
    # Werner 2026-09-04: hyperbolic, no table hit, but c>=21 NOT proven.
    # Neither known nor novelty candidate — open, and for a reason that lies with
    # US: up to c=20 is enumerated, we merely do not have the c=20 table.
    n_un = sum(1 for r in reps if r['_class'] == 'unresolved_c_le_20')
    # Disc #814: hyperbolic, no hit, c >= 21 proven — but alternating (or
    # alternation undecided) at c <= 24, i.e. in the range that Rankin/Flint/Schermann
    # enumerated in 2004. No novelty candidate, for the same reason as
    # unresolved_c_le_20: the table exists, we do not have it.
    n_ua = sum(1 for r in reps if r['_class'] == 'unresolved_alternating_c_le_24')
    n_nn = 0  # novel_nonhyp abolished (non-hyperbolic is never novelty); field stays =0 for backward compat
    n_mc, mc_routes = mirror_collapsed_count(reps)
    # Disc #808: the same number for 'known'. KnotInfo and Burton list ONE representative
    # per chirality pair — 128 counts chirality, the tables do not. Both numbers are
    # true, for different questions. The asymmetry (novel had the number, known did not)
    # was itself the source of error, cf. Disc #779/#804.
    n_known_mc, known_mc_routes = mirror_collapsed_count(reps, 'known')
    n_kn_mc, _ = mirror_collapsed_count(reps, 'known_nonhyp')
    n_un_mc, _ = mirror_collapsed_count(reps, 'unresolved_c_le_20')
    n_ua_mc, _ = mirror_collapsed_count(reps, 'unresolved_alternating_c_le_24')
    sum_real = sum(r['n'] for r in reps)

    # On-demand dial split (Werner 2026-07-01: „eigene on-demand-familie, nicht in die
    # outreach-headline" (EN: own on-demand family, not in the outreach headline)). Dialled knots (origin=on_demand_dial, set in build_corpus_db
    # from the family marker) are real distinct corpus knots (stay in n_dist/n_nh/n_mc),
    # but do NOT count into the outreach sweep headline. Both numbers remain derivable.
    def _is_dial(r):
        return r['_rep'].get('origin') == 'on_demand_dial'
    reps_sweep = [r for r in reps if not _is_dial(r)]
    n_nh_sweep = sum(1 for r in reps_sweep if r['_class'] == 'novel_hyperbolic')
    n_mc_sweep, _ = mirror_collapsed_count(reps_sweep)
    dial_fams = Counter(r['_rep'].get('origin_family') or '?'
                        for r in reps if _is_dial(r) and r['_class'] == 'novel_hyperbolic')
    n_dial_nh = sum(dial_fams.values())

    # BOOKKEEPING MUST FIT
    assert sum_real == len(knots), f"Σ n_realisierungen {sum_real} != records {len(knots)}"
    assert n_known + n_nh + n_kn + n_un + n_ua == n_dist, "Class sum != distinct"
    ids = [r['id'] for r in reps]
    assert len(set(ids)) == n_dist, (
        f"ID collision: {len(ids)} groups but only {len(set(ids))} unique IDs "
        f"-> {[i for i in ids if ids.count(i) > 1][:5]}")
    # Distinctions gate: every proven distinctness MUST have taken effect. A silent
    # no-op would look like success here (the number would simply stay) — exactly the trap
    # from Disc #266. What is checked is not that the file was read, but that the
    # partners end up carrying DIFFERENT distinct_knot_id.
    _dist = load_distinctions()
    if _dist:
        by_stem = {}
        for m in knots:
            if m.get('stem'):
                by_stem[m['stem']] = m.get('distinct_knot_id')
        for fall, e in _dist.items():
            a, b = e.get('abgespalten'), e.get('verbleibt')
            ida, idb = by_stem.get(a), by_stem.get(b)
            assert ida is not None, f"proven_distinctions '{fall}': stem {a} not in the corpus"
            assert idb is not None, f"proven_distinctions '{fall}': stem {b} not in the corpus"
            assert ida != idb, (
                f"proven_distinctions '{fall}' INEFFECTIVE: {a} and {b} still carry "
                f"the same distinct_knot_id {ida!r} — split did not arrive")
            for stem, expected_field in ((a, 'erwartete_id_abgespalten'), (b, 'erwartete_id_verbleibt')):
                expected_ = e.get(expected_field)
                actual = by_stem.get(stem)
                assert expected_ is None or actual == expected_, (
                    f"proven_distinctions '{fall}': {stem} got id {actual!r}, expected {expected_!r}")
        print(f"  ✓ proven_distinctions: {len(_dist)} split(s) arrived, "
              f"partners carry different IDs")

    # Turaev gate (Disc #943): every route-3 receipt must have reached a representative.
    # A stem that drifted (other representative chosen) would silently drop the bound.
    _tb = _turaev_bounds()
    if _tb:
        lost = sorted(set(_tb) - _TURAEV_APPLIED)
        assert not lost, f"turaev_bounds.json: {len(lost)} receipt(s) did not arrive: {lost}"
        print(f"  ✓ turaev_bounds: {len(_tb)}/{len(_tb)} receipts arrived")

    # Retraction gate: EVERY retraction must land on a distinct_knot_id. A silent
    # no-op (did drift) would NOT lower the headline and would look like success (Disc #266).
    if retractions:
        idset = set(ids)
        missing = [k for k in retractions if k not in idset]
        assert not missing, (
            f"Maguire retraction does NOT land: {len(missing)} did(s) without match "
            f"-> {missing[:8]} (label/adjudication out of step)")
        if burton:
            missing_b = [k for k in burton if k not in idset]
            assert not missing_b, f"Burton retraction IDs not in the corpus: {missing_b}"
            assert n_burton_applied == len(burton), (
                f"Burton retraction count: applied {n_burton_applied} != file {len(burton)}")
        assert n_retracted_applied == len(retractions), (
            f"Retraction count: applied {n_retracted_applied} != file {len(retractions)}")
    # Absence gate: the same did-drift trap as with retractions (Disc #266) — every
    # absence ID must arrive, a silent no-op would look like success.
    if absence:
        idset = set(ids)
        missing_a = [k for k in absence if k not in idset]
        assert not missing_a, f"Burton absence IDs not in the corpus: {missing_a}"
        assert n_absence_applied == len(absence), (
            f"Absence count: applied {n_absence_applied} != file {len(absence)}")

    print("== DISTINCT-KNOT BOOKKEEPING ==")
    print(f"  corpus_db Records:                 {len(knots)}")
    print(f"  Σ n_realisierungen (= records):    {sum_real}   {'OK' if sum_real==len(knots) else 'ERROR'}")
    print(f"  distinct knots (manifold):         {n_dist}")
    print(f"    known (table):                   {n_known}")
    print(f"    novel hyperbolic (chiral):       {n_nh}")
    print(f"       of which mirror-merged:       {n_mc}  (true total of novel-hyp)")
    print(f"         → sweep-only (OUTREACH HEADLINE): {n_mc_sweep}")
    print(f"         → on-demand dial (own families, NOT headline): {n_dial_nh}  {dict(dial_fams) if dial_fams else ''}")
    print(f"    known non-hyp (torus/sum):       {n_kn}  (novel_nonhyp abolished = 0)")
    print(f"    OPEN, c>=21 not proven:          {n_un}  ← up to c=20 is enumerated; these")
    print(f"                                         lie in the enumerated range, we")
    print(f"                                         merely lack the table (c=20 not")
    print(f"                                         public, c=19 run open)")
    print(f"    OPEN, alternating c<=24:         {n_ua}  ← Disc #814: the alternating")
    print(f"                                         prime knots are published up to c=23")
    print(f"                                         enumerated (Rankin/Flint/Schermann")
    print(f"                                         2004), computed up to c=24 (Fontaine")
    print(f"                                         2007). Same situation as c=20.")
    print(f"  Class sum {n_known}+{n_nh}+{n_kn}+{n_un}+{n_ua} = {n_known+n_nh+n_kn+n_un+n_ua}  {'OK' if n_known+n_nh+n_kn+n_un+n_ua==n_dist else 'ERROR'}")
    # Disc #808: mirror-collapsed per class. ADDITIONAL KEY FIGURES, they replace nothing —
    # the left column counts distinct KNOTS, the right the TABLE ENTRIES they
    # realise (KnotInfo/Burton list one representative per chirality pair).
    # ⛔ Not to the outside without Werner's gate.
    print(f"  ── Mirror-collapsed per class (Disc #808, ADDITIONAL key figure) ──")
    for _lab, _raw, _mc in (("known", n_known, n_known_mc),
                            ("novel hyperbolic", n_nh, n_mc),
                            ("known non-hyp", n_kn, n_kn_mc),
                            ("open c<=20", n_un, n_un_mc),
                            ("open altern. c<=24", n_ua, n_ua_mc)):
        print(f"    {_lab:22s} {_raw:4d} distinct  ->  {_mc:4d} mirror-collapsed "
              f"({_raw - _mc} mirror pair(s))")
    _w = {**{f"novel/{k}": v for k, v in mc_routes.items()},
          **{f"known/{k}": v for k, v in known_mc_routes.items()}}
    print(f"    detected via: {_w if _w else 'nothing collapsed'}")
    if known_mc_routes.get('theta') or mc_routes.get('theta'):
        print(f"    ↑ the theta route took effect — sigma = 0 on both sides, only theta separates "
              f"(the gap from Disc #807/#808)")
    if retractions:
        if burton:
            print(f"  Burton census retractions applied: {n_burton_applied}/{len(burton)}  "
                  f"(c<=19 complete; c=20 exists, data not public)")
        print(f"  Maguire-Extended retractions applied: {n_retracted_applied}/{len(retractions)}  "
              f"(c14-17 tabulated, included in 'known'; separate from KnotInfo c<=13)")
    if absence:
        print(f"  Burton census ABSENCE applied: {n_absence_applied}/{len(absence)}  "
              f"(c>=20 proven via absence in the c<=19 census; classes unchanged, bar 21)")

    # Top multiplicities (family evidence)
    print("\n  largest realisation multiplicities (machine families):")
    for r in sorted(reps, key=lambda x: -x['n'])[:8]:
        print(f"    {r['id']:34} {r['n']:>3} realisations  [{r['_class']}]")

    if apply:
        # Backup centrally in _backups/ (Werner 2026-05-30: „ein ordner für die backups" (EN: one folder for the backups)),
        # rolling (no timestamp pile-up).
        bak = ROOT / "_backups/files/freeze_manifest/pipeline_scripts/corpus_db.json.bak-pre-label"
        bak.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(DB, bak)
        # Pull the base bookkeeping along (in case records were removed directly, without a rebuild)
        from collections import Counter as _C
        db['meta']['total_knots'] = len(knots)
        db['meta']['categories'] = dict(_C(m.get('category') for m in knots if m.get('category')))
        db['meta']['distinct_knots'] = n_dist
        db['meta']['distinct_known'] = n_known
        db['meta']['distinct_novel_hyperbolic'] = n_nh
        db['meta']['distinct_novel_hyperbolic_mirror_collapsed'] = n_mc  # true total (incl. dial)
        # Disc #808: mirror-collapsed numbers for ALL classes, not only novel_hyperbolic.
        # ADDITIONAL KEY FIGURES — they replace nothing. distinct_known remains the number of
        # distinct knots; *_mirror_collapsed is the number of TABLE ENTRIES that they
        # realise (KnotInfo/Burton list one representative per chirality pair).
        # ⛔ NOT for the outside presentation without Werner's gate — that would be a new claim.
        db['meta']['mirror_collapsed_by_class'] = {
            'known': n_known_mc,
            'novel_hyperbolic': n_mc,
            'known_nonhyp': n_kn_mc,
            'unresolved_c_le_20': n_un_mc,
            'unresolved_alternating_c_le_24': n_ua_mc,
            'wege_novel_hyperbolic': mc_routes,
            'wege_known': known_mc_routes,
            'note': ('Disc #808. Mirror pairs are detected via TWO routes: sig = -sig\' != 0, '
                     'or — when sig is ZERO on both sides — theta = -theta\' with theta != 0, '
                     'filtered at the point and confirmed on the full polynomial. The second route was missing '
                     'until 19.09.2026 (the condition read `si == -ei and si != 0`, sigma = 0 was '
                     'excluded by construction). Counter-cases in '
                     'maschinen_ideen/KNOWN_DOPPELZAEHLUNG_20260919/gegenfall_spiegel.py: '
                     'conway_cw/conway_kt (mutants, same det/sig/vol) must NOT collapse. '
                     'ADDITIONAL KEY FIGURE, replaces no existing number.'),
        }
        # Outreach sweep headline = novel-hyp-mirror-collapsed WITHOUT on-demand-dialled families
        # (Werner 2026-07-01). THIS is the number for mails/balances/web.
        db['meta']['outreach_headline_sweep'] = n_mc_sweep
        db['meta']['distinct_novel_hyperbolic_sweep'] = n_nh_sweep
        db['meta']['on_demand_dial'] = {
            'total_members': n_dial_nh,
            'families': dict(dial_fams),
            'note': ('On-demand dialled families (novelty_dial, Disc #384/#387). Real distinct '
                     'corpus knots (included in total/distinct/novel_hyp), but NOT in the '
                     'outreach sweep headline (Werner 2026-07-01). Headline = outreach_headline_sweep.'),
        }
        db['meta']['distinct_known_nonhyp'] = n_kn       # known non-hyp (torus/connected sum)
        db['meta']['distinct_unresolved_c_le_20'] = {
            'count': n_un,
            'rule': 'hyperbolic + no table hit, but c>=21 NOT proven',
            'warum': ('The prime knots are enumerated up to AND INCLUDING c=20 '
                      '(Thistlethwaite, Algebraic & Geometric Topology 25 (2025) 329-344; '
                      '1,847,319,428 items). These knots therefore lie in a COMPLETELY '
                      'tabulated range — that we do not find them is a statement '
                      'about our tables, not about the knots. The c=20 data are '
                      'not public; Burton\'s census up to c=19 is here, the comparison '
                      'against c=19 was still outstanding on 04.09.2026.'),
            'schranke': ('TWO routes, the larger wins (field c_lower_bound_proven, '
                         'provenance in c_lower_bound_via): span V(K) <= c '
                         '(Kauffman-Murasugi-Thistlethwaite) AND c >= 2g+1 '
                         '(Seifert algorithm, only with HFK-exact genus; '
                         'genus route built in 2026-09-05 on Werner\'s go)'),
            'werner': 'untere schranke = k 20. no novelty. (2026-09-04) · genus. (2026-09-05)',
        }
        db['meta']['distinct_unresolved_alternating_c_le_24'] = {
            'count': n_ua,
            'rule': ('hyperbolic + no table hit + c >= 21 PROVEN, but alternating '
                     '(or alternation undecided) at c <= 24'),
            'warum': ('Disc #814, 20.09.2026. There are TWO enumeration limits, not one. All '
                      'prime knots are enumerated up to c = 20 (Thistlethwaite 2025); the ALTERNATING ones, '
                      'however, are published up to c = 23 (Rankin, Flint & Schermann, J. Knot Theory Ramif. 13 '
                      '(2004) 57-150; 25,182,878,921 items at c=23, finished 14.03.2004) and computed up to c = 24 '
                      '(Fontaine, Rankin & Flint 2007; OEIS A002864 a(24) = 128,564,665,125). '
                      'Verbatim from Thistlethwaite himself, AGT 25 (2025) p. 2: "This has indeed been '
                      'accomplished quite dramatically for n <= 23 [Rankin et al. 2004]." An '
                      'alternating knot with c = 21..24 thus lies in an existing '
                      'enumeration — that we do not find it is a statement about our tables '
                      '(Hoste, Handbook: "not yet publicly available"), not about the knot. '
                      'Same situation as unresolved_c_le_20.'),
            'status': ('Field alternating_proven per record: True = alternating PROVEN (flag from '
                       'Regina OR span V == c_ub by Kauffman-Murasugi-Thistlethwaite), None = '
                       'undecided (span V == lower bound < c_ub). False never occurs here: '
                       'proven non-alternating ones with c >= 21 are novel_hyperbolic.'),
            'beleg': 'konversationen/alternierende_enumeration_20260920/',
        }
        db['meta']['distinct_novel_nonhyp'] = n_nn       # abolished = 0 (non-hyp is never novelty)
        db['meta']['burton_absence_c19'] = {
            'applied': n_absence_applied,
            'quelle': 'burton_absence_c19.json (run work_c19_20260904, 2026-09-06)',
            'bedeutung': ('c>=20 PROVEN via absence in the complete c<=19 census '
                          '(350,449,955 knots, all categories exact, detfails=0). '
                          'Raises only c_lower_bound_proven; classes unchanged (bar 21).'),
        }
        db['meta']['distinct_retracted_burton_census'] = {
            'applied': n_burton_applied,
            'source': 'burton_census_c_le_19',
            'url': 'https://regina-normal.github.io/data.html',
            'boundary': 'c<=19 checked completely',
            'note': ('distinct_knot_ids that appear in Burton\'s complete census up to c=19 '
                     '(is_isometric_to-confirmed, census name in the receipt). Field '
                     'realizes_tabulated_burton. ⚠️ OPEN LIMIT: the prime knots are enumerated up to '
                     'AND INCLUDING c=20 (Thistlethwaite, AGT 25 (2025) 329-344, '
                     '1,847,319,428 items), but the c=20 data are not public. Whoever has '
                     'span V(K)>=21 is proven beyond independently of that; the rest '
                     'carry the c=20 gap.'),
        }
        db['meta']['distinct_retracted_maguire_extended'] = {
            'applied': n_retracted_applied,
            'source': 'maguire_hyp_jones_c_le_17',
            'disc': [512, 516],
            'boundary': 'beyond-17',
            'note': ('distinct_knot_ids tabulated against Maguire\'s c<=17 hyp-Jones census '
                     '(volume-confirmed, is_isometric_to). Field realizes_tabulated_extended, '
                     'SEPARATE from crosscheck_knotinfo (KnotInfo c<=13). Soft floor: 90 '
                     'unresolved no-Jones + 7 collisions remain provisionally novel; c18/19 open.'),
        }
        db['meta']['labeled_at'] = date.today().isoformat()
        json.dump(db, open(DB, 'w'), indent=2, ensure_ascii=False)
        print(f"\n== WRITTEN: {DB}  (backup {bak.name}) ==")
        if newly_registered:
            if REGISTRY_FILE.exists():
                shutil.copy2(REGISTRY_FILE, str(REGISTRY_FILE) + ".bak-pre-label")
            json.dump(registry, open(REGISTRY_FILE, 'w'), indent=1, ensure_ascii=False)
            print(f"== REGISTRY: {len(newly_registered)} new identity(ies) frozen "
                  f"({REGISTRY_FILE.name}, total {len(frozen)}) ==")
        else:
            print(f"== REGISTRY: nothing new (total {len(frozen)} frozen) ==")
    else:
        print("\n== DRY RUN (no --apply) ==")
        print(f"== REGISTRY: {len(newly_registered)} identity(ies) WOULD be frozen "
              f"(total afterwards {len(frozen)}) ==")
    print("LABEL_DISTINCT_KNOTS_DONE")


def self_check():
    """Four probes, with counter-case in BOTH directions (footgun 16):
    (a) re-sorting WITH registry: ID stays.
    (b) the same re-sorting WITHOUT registry: ID flips — proves that (a) checks something.
    (c) new identity: gets a fresh ID and is registered.
    (d) corrupt registry (duplicate ID): aborts loudly.
    Everything in memory, NOTHING is written."""
    db = json.load(open(DB))
    knots = json.loads(json.dumps(db['knots']))  # deep copy
    groups = defaultdict(list)
    for m in knots:
        groups[ident(m)].append(m)
    grp = [{'key': k, 'ms': ms, 'base': canon_id(ms), 'rep': ms[0]} for k, ms in groups.items()]
    grp.sort(key=lambda g: (g['base'], str(g['rep'].get('det')),
                            str(g['rep'].get('sig')), str(g['rep'].get('vol'))))
    # seed the registry in memory
    frozen = {}
    assign_ids(grp, frozen)
    before = {g['ms'][0]['stem']: g['id'] for g in grp}
    # victim: a group whose ID comes from the S_ folder (no proven name)
    victim = next(g for g in grp if str(g['id']).startswith('S_')
                 and not any(m.get('stem') in PROVEN_IDENT for m in g['ms']))
    stem = victim['ms'][0]['stem']
    # (a) simulate re-sorting: change the folder name, reassign IDs — WITH registry
    for m in victim['ms']:
        m['cluster_folder'] = str(m.get('cluster_folder')) + '_umsortiert'
    grp2 = [{'key': k, 'ms': ms, 'base': canon_id(ms), 'rep': ms[0]} for k, ms in groups.items()]
    grp2.sort(key=lambda g: (g['base'], str(g['rep'].get('det')),
                             str(g['rep'].get('sig')), str(g['rep'].get('vol'))))
    assign_ids(grp2, frozen)
    after = {g['ms'][0]['stem']: g['id'] for g in grp2}
    a_ok = after[stem] == before[stem]
    # (b) WITHOUT registry the same re-sorting flips the ID — otherwise (a) checks nothing
    grp3 = [{'key': k, 'ms': ms, 'base': canon_id(ms), 'rep': ms[0]} for k, ms in groups.items()]
    grp3.sort(key=lambda g: (g['base'], str(g['rep'].get('det')),
                             str(g['rep'].get('sig')), str(g['rep'].get('vol'))))
    assign_ids(grp3, {})
    without = {g['ms'][0]['stem']: g['id'] for g in grp3}
    b_ok = without[stem] != before[stem]
    # (c) new identity
    fake = {'stem': 'selbsttest_fake', 'det': '999999999999991', 'sig': '0',
            'alex': 't-selbsttest', 'theta_at_22_7__21_13': 'x', 'vol': 12.345678,
            'cluster_folder': 'S_99_d999_TEST', 'c': 99}
    gf = {'key': ident(fake), 'ms': [fake], 'base': canon_id([fake]), 'rep': fake}
    n_new = assign_ids([gf], frozen)
    c_ok = len(n_new) == 1 and gf['id'] == 'S_99_d999_TEST'
    # (d) corrupt registry
    try:
        assign_ids([], {'k1': {'id': 'DOPPELT'}, 'k2': {'id': 'DOPPELT'}})
        d_ok = False
    except SystemExit:
        d_ok = True
    for name, ok, txt in [('a', a_ok, f'Re-sorting with registry: {stem} stays {before[stem]}'),
                          ('b', b_ok, f'without registry it flips ({without[stem]} != {before[stem]}) — test bites'),
                          ('c', c_ok, 'new identity: fresh ID + registered'),
                          ('d', d_ok, 'duplicate registry ID aborts loudly')]:
        print(f"  ({name}) {'✅' if ok else '🔴'} {txt}")
    if not (a_ok and b_ok and c_ok and d_ok):
        print("🔴 SELFTEST RED"); return 1
    print("✅ SELFTEST GREEN (4 probes, nothing written)"); return 0


if __name__ == "__main__":
    if '--self-check' in sys.argv:
        sys.exit(self_check())
    main(apply='--apply' in sys.argv)
