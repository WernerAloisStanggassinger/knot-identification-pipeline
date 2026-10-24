#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
# pipeline_run.sh — end-to-end pipeline B for one knot.
#
# Call from the project folder that contains the single .fseries:
#     cd /path/to/myknot_project
#     bash /media/squeezy/kauffman_storage/freeze_manifest/pipeline_run.sh
#
# Optional arguments:
#     bash pipeline_run.sh [<fseries-file>] [--skip-scanner] [--skip-tracker]
#
# Detects automatically:
#   - the single .fseries in the current directory (preflight.py)
#   - pi antisymmetry / pi periodicity (Discovery #70)
#   - target_count + frame_cap from the n_cross distribution (detect_target.py)
#   - aggressive_simplify n_min/n_max
#
# Pipeline stages:
#   1. knot_scanner.py                   → <knot>_crossings.parquet
#   2. knot_tracker_dbscan.py            → <knot>_tracks_dbscan.parquet
#   3. kinematic_resolve.py              → <knot>_kinematic_result.json
#   4. aggressive_simplify.py            → <knot>_aggressive_simplify.json
#   5. mandatory_invariants.py            → <knot>_mandatory_invariants.json
#   6. mandatory_invariants_extra.py      → <knot>_mandatory_invariants_complete.json
#   7. 7x enrich_*.py + repair          (in place on complete.json,
#                                        incl. enrich_non_hyperbolic for
#                                        torus + connected-sum classification,
#                                        Discovery #144)
#   8. summary_writer.py                 → <knot>_summary.md
#
# At every step: abort if rc != 0; error message with the stage index.
set -euo pipefail

# REGINA LATCH (Werner's go 2026-08-22 "set the latch to 120").
# WHY 120 s: regina's Alexander doubles every ~16 crossings — measured on our own
# corpus c=96 17.3s · 114 36.5s · 128 65.3s · 144 124.4s · 192 >1200s.
# 120 s thus covers everything up to c≈143 (largest c in the corpus: 128); above that
# regina is >20x slower than the fallback (point Alexander, variant 6,
# 557/557 corpus determinants hit: c=128 in 4.4s, c=192 in 11s).
# The earlier 1200 s latch cost three runs on 22.08. 20 minutes of waiting each
# and delivered NOTHING in the end.
# ${VAR:-} still lets a value set from outside win.
export REGINA_TIMEOUT="${REGINA_TIMEOUT:-120}"

# All locations come from site_config.py (portable by default; KNOT_PIPELINE_PROFILE=authors = the authors' runs).
SC_PY="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/site_config.py"
PIPE=$(python3 -I "$SC_PY" --get PIPELINE_DIR)
CAMERA=$(python3 -I "$SC_PY" --get CAMERA)
PROJ_DIR="$(pwd)"

log()   { echo -e "\033[1;34m[PIPELINE]\033[0m $*"; }
err()   { echo -e "\033[1;31m[PIPELINE ERROR]\033[0m $*" >&2; }
ok()    { echo -e "\033[1;32m[PIPELINE OK]\033[0m $*"; }
stage() { echo -e "\n\033[1;33m=== STAGE $1: $2 ===\033[0m"; }

# timing
T0=$(date +%s)

# ─────────────────────────────────────────────
# Stage 9 as a function + EXIT trap.
# Stage 9 (SORT+RENDER) MUST run even when earlier stages aborted with
# exit 3 (e.g. length_spectrum timeout) — otherwise broken knots stay
# behind as top-level folders and block the next
# wave. sort_and_render.py detects type a (no _complete.json) and
# type b (mandatory core fields all None) automatically and sorts them into
# erzeugte_knoten/pipeline_failures/.
# On by default, can be switched off via SKIP_SORT_RENDER=1.
# ─────────────────────────────────────────────
STAGE9_RAN=0

run_stage9() {
    if [ "${SKIP_SORT_RENDER:-0}" = "1" ]; then
        log "Stage 9 (Sort+Render) skipped — SKIP_SORT_RENDER=1"
        STAGE9_RAN=1
        return 0
    fi
    [ "$STAGE9_RAN" = "1" ] && return 0
    STAGE9_RAN=1

    stage "9" "SORT + RENDER"
    local T9
    T9=$(date +%s)
    local SORT_ARGS=()
    [ "${SORT_NO_RENDER:-0}" = "1" ] && SORT_ARGS+=(--no-render)
    [ "${SORT_DRY_RUN:-0}" = "1" ] && SORT_ARGS+=(--dry-run)
    conda run -n gnottheorie_env python "$PIPE/sort_and_render.py" "$PROJ_DIR" "${SORT_ARGS[@]}" \
        || err "Stage 9 (Sort+Render) failed — the mandatory data are saved, though"
    ok "Sort+Render done in $(($(date +%s)-T9))s"

    # If stage 9 moved the folder, find the new PFLC/SUMMARY
    if [ -n "${PFLC:-}" ] && [ ! -f "$PFLC" ] && [ -n "${KNOT:-}" ]; then
        local NEW_PFLC
        NEW_PFLC=$(find "$(python3 -I "$SC_PY" --get GENERATED_KNOTS_DIR)" -name "*${KNOT}*_mandatory_invariants*_complete.json" 2>/dev/null | head -1)
        if [ -n "$NEW_PFLC" ] && [ -f "$NEW_PFLC" ]; then
            log "Folder moved → new PFLC: $NEW_PFLC"
            PFLC=$NEW_PFLC
            PROJ_DIR=$(dirname "$NEW_PFLC")
            SUMMARY=$PROJ_DIR/${KNOT}_summary.md
        fi
    fi
}

on_exit() {
    local rc=$?
    if [ $rc -ne 0 ] && [ "$STAGE9_RAN" = "0" ]; then
        echo >&2
        err "Pipeline aborted (rc=$rc) — stage 9 cleanup runs anyway..."
        run_stage9 || true
    fi
    exit $rc
}
trap on_exit EXIT

# Signal traps for a watchdog kill (Disc #217 + kf373 salvage 2026-05-19).
# mem_watchdog.sh sends SIGTERM, waits 3s, then SIGKILL (uncatchable).
# Within that 3s window stage 9 must get the folder out of the top level,
# otherwise the knot stays behind as a wave blocker (kf373 → pipeline_failures/
# wave32_watchdog_killed/ pattern).
# SORT_NO_RENDER=1: Blender alone would already exceed the 3s grace period.
on_signal() {
    local signame=$1
    local exitcode=$2
    if [ "$STAGE9_RAN" = "0" ]; then
        echo >&2
        err "Signal $signame received — stage 9 cleanup (--no-render) running..."
        SORT_NO_RENDER=1 run_stage9 || true
    fi
    exit $exitcode
}
trap 'on_signal SIGTERM 143' TERM
trap 'on_signal SIGINT  130' INT
trap 'on_signal SIGHUP  129' HUP

# ─────────────────────────────────────────────
# 0. find the .fseries
# ─────────────────────────────────────────────
if [[ $# -ge 1 && -f "$1" ]]; then
    FSERIES="$1"
    shift
else
    mapfile -t FS_FILES < <(find "$PROJ_DIR" -maxdepth 1 -name '*.fseries' -type f | sort)
    if [[ ${#FS_FILES[@]} -eq 0 ]]; then
        err "No .fseries found in the project folder ($PROJ_DIR)"
        exit 2
    elif [[ ${#FS_FILES[@]} -gt 1 ]]; then
        err "Several .fseries found — please call them one at a time:"
        for f in "${FS_FILES[@]}"; do err "  $f"; done
        exit 2
    fi
    FSERIES="${FS_FILES[0]}"
fi

KNOT="$(basename "$FSERIES" .fseries)"
log "Knot: $KNOT"
log "Source: $FSERIES"
log "Pipeline scripts: $PIPE"
log "Standard camera:  $CAMERA"

# output paths
CROSS=$PROJ_DIR/${KNOT}_crossings.parquet
TRACKS=$PROJ_DIR/${KNOT}_tracks_dbscan.parquet
KIN=$PROJ_DIR/${KNOT}_kinematic_result.json
AGG=$PROJ_DIR/${KNOT}_aggressive_simplify.json
PFL=$PROJ_DIR/${KNOT}_mandatory_invariants.json
PFLC=$PROJ_DIR/${KNOT}_mandatory_invariants_complete.json
SUMMARY=$PROJ_DIR/${KNOT}_summary.md
PREFLIGHT_JSON=$PROJ_DIR/${KNOT}_preflight.json
DETECT_JSON=$PROJ_DIR/${KNOT}_detect_target.json

# ─────────────────────────────────────────────
# 0a. Preflight: pi symmetry + recommendations
# ─────────────────────────────────────────────
stage "0a" "PREFLIGHT — symmetry detection"
conda run -n gnottheorie_env python "$PIPE/preflight.py" "$FSERIES" --json > "$PREFLIGHT_JSON" \
    || { err "preflight.py failed"; exit 3; }

# parse flags from the JSON
FOLD_PI=$(conda run -n gnottheorie_env python -c \
    "import json;d=json.load(open('$PREFLIGHT_JSON'));print('1' if d['tracker_fold_pi_antisym_flag'] else '0')")
PERIOD_PI=$(conda run -n gnottheorie_env python -c \
    "import json;d=json.load(open('$PREFLIGHT_JSON'));print('1' if d['scanner_period_pi_flag'] else '0')")
SYM_CLASS=$(conda run -n gnottheorie_env python -c \
    "import json;d=json.load(open('$PREFLIGHT_JSON'));print(d['symmetry_class'])")

log "Symmetry class: $SYM_CLASS"
log "  → Scanner --period-pi:        $([ "$PERIOD_PI" = "1" ] && echo 'YES' || echo 'NO')"
log "  → Tracker --fold-pi-antisym:  $([ "$FOLD_PI" = "1" ] && echo 'YES' || echo 'NO')"

SCANNER_PERIOD_FLAG=""
TRACKER_FOLD_FLAG=""
RESOLVE_FOLD_FLAG=""
[ "$PERIOD_PI" = "1" ] && SCANNER_PERIOD_FLAG="--period-pi"
[ "$FOLD_PI"   = "1" ] && TRACKER_FOLD_FLAG="--fold-pi-antisym"
[ "$FOLD_PI"   = "1" ] && RESOLVE_FOLD_FLAG="--fold-pi-antisym"

# ─────────────────────────────────────────────
# RECIPE #3 (error log 2026-06-04 ~22:40, Werner's go 2026-06-05): POLYGON_PD=1
# replaces stages 1–4a completely with polygon_resolve.py (stage P) — polygon
# diagram directly from the .fseries (own projection, Nyquist gate) instead of
# the raw scanner diagram. For polyc07-class curves whose scanner path dies on RAM
# (tr14/tr15e: stage 4 balloon >28G or stage 5 regina rc=137), while
# the exact layer runs light on RAM on the same curve. Runs in the base env
# (regina+numpy). Stage 5+ consumes kinematic_result.json unchanged.
# Consequence as with SKIP_AGGRESSIVE_SIMPLIFY: no multi-frame sweep, but
# the diagram is regina-simplified + Nyquist-confirmed.
# ─────────────────────────────────────────────
if [ "${POLYGON_PD:-0}" = "1" ]; then

stage "P" "POLYGON-RESOLVE — recipe #3, replaces stages 1–4a (POLYGON_PD=1)"
TP=$(date +%s)
conda run -n base python "$PIPE/polygon_resolve.py" "$FSERIES" \
    --name "$KNOT" \
    --output "$KIN" \
    ${POLYGON_AXIS:+--axis "$POLYGON_AXIS"} \
    || { err "Stage P (Polygon-Resolve) failed"; exit 3; }
ok "Polygon-Resolve done in $(($(date +%s)-TP))s"

else

# ─────────────────────────────────────────────
# 1. Scanner
# ─────────────────────────────────────────────
stage "1" "SCANNER — CUDA, 131072 cameras"
T1=$(date +%s)
conda run -n gnottheorie_env python "$PIPE/knot_scanner.py" "$FSERIES" \
    --camera "$CAMERA" \
    --output "$CROSS" \
    --N 4096 \
    --max-cross-per-frame 2048 \
    $SCANNER_PERIOD_FLAG \
    || { err "Stage 1 (Scanner) failed"; exit 3; }
ok "Scanner done in $(($(date +%s)-T1))s"

# ─────────────────────────────────────────────
# 2. DBSCAN-Tracker
# ─────────────────────────────────────────────
stage "2" "DBSCAN-TRACKER on (ti, tj)"
T2=$(date +%s)
# ── LADDER INSTEAD OF ONE ATTEMPT (Disc #670, 2026-08-23) ──────────────────────
# --min-cell-count 2000 is NOT a density threshold. Measured over the re-run of 22.08.
# the separation is overlap-free: green from 1,967,552 crossings,
# red up to 1,775,870. The threshold is in effect a hidden LOWER BOUND ON
# THE CROSSING COUNT and thus on the winding density of the curve — roomy,
# simple knots failed, FIVE TREFOILS among them (wp_trefoil_mneg4,
# 2,032,127 crossings, is in the corpus as K3_1; its sibling
# wp_trefoil_m0, 533,109 crossings, failed — same machine).
# Hence on "0 active cells" the threshold is lowered step by step.
#
# GATE BEFORE IT, otherwise it would be Frankenstein all over again: for a curve with a
# crossing-free projection the ladder would lower until SOMETHING
# passes. Such a curve is embedded in the plane ⟹ UNKNOT ⟹ the
# rejection is CORRECT and stands (60 such cases in the re-run).
# MIN_CELL from outside wins and switches the ladder off (exactly one attempt).
# TRACKER_EPS from outside (2026-08-24, Disc #703): eps is the GRID SPACING, not
# the threshold. For very dense diagrams (UEBER_k3_s1p5: min 958 crossings per
# frame) 48 % of the cells are active at eps=0.01, they percolate, and everything
# falls into ONE cluster. Measured on this knot:
#   eps=0.01  ->    1 cluster (100 % of the points in it)  unusable
#   eps=0.005 ->   56 clusters (largest 99.7 %)            still one lump
#   eps=0.002 -> 1059 clusters (largest  1.1 %)            percolation broken
#   eps=0.001 ->  809 clusters (largest  0.5 %)
# Default stays 0.01 — that is right for the whole corpus so far.
TRIAGE="$(python3 -I "$SC_PY" --get TOOLS_DIR)/scan_triage.py"
TRACKER_ERR="$PROJ_DIR/${KNOT}_tracker_stderr.log"
: > "$TRACKER_ERR"
tracker_versuch() {
    # The child's stderr is RECORDED. Before, it vanished, and
    # "See above for error" pointed at nothing (Disc #702).
    conda run -n gnottheorie_env python "$PIPE/knot_tracker_dbscan.py" "$CROSS" \
        --output "$TRACKS" \
        --eps "${TRACKER_EPS:-0.01}" \
        --min-cell-count "$1" \
        $TRACKER_FOLD_FLAG 2>>"$TRACKER_ERR"
    local rc=$?
    if [ $rc -ne 0 ] && [ -s "$TRACKER_ERR" ]; then
        sed 's/^/[TRACKER-STDERR] /' "$TRACKER_ERR" | tail -20 >&2
    fi
    return $rc
}
MIN_CELL_USED=""
for MC in ${MIN_CELL:-2000 1000 500 300 200 100}; do
    # NOT a bare call: under set -e a bare rc=1 killed the script before the ladder could lower (broken 24.08.-23.10.2026).
    TRACKER_RC=0
    tracker_versuch "$MC" || TRACKER_RC=$?
    if [ $TRACKER_RC -eq 0 ]; then
        MIN_CELL_USED="$MC"
        break
    fi
    # rc=1 is the ONLY answer that means "0 active cells"
    # (knot_tracker_dbscan.py:249, sys.exit(1) when n_active == 0).
    # Everything else is a crash and must NOT let the ladder continue —
    # otherwise a kill gets printed six times as a mathematical verdict (Disc #702).
    if [ $TRACKER_RC -ne 1 ]; then
        err "Stage 2: tracker ABORTED with rc=$TRACKER_RC at --min-cell-count $MC."
        err "         This is NOT an empty result. The ladder does NOT go lower here,"
        err "         because nothing was measured about the curve."
        if [ $TRACKER_RC -ge 128 ]; then
            err "         rc>=128 ⟹ killed from outside, signal $((TRACKER_RC-128))."
        fi
        err "         Suspect memory. Check:  journalctl --no-pager -u earlyoom -n 20"
        err "         Crossing table: $(du -h "$CROSS" 2>/dev/null | cut -f1)"
        err "         Child's stderr: $TRACKER_ERR"
        exit 3
    fi
    if [ -n "${MIN_CELL:-}" ]; then
        break                       # explicit request from outside: no follow-up
    fi
    if conda run -n gnottheorie_env python "$TRIAGE" --gate "$CROSS"; then
        err "Stage 2: UNKNOT — the threshold is not too high, the curve is unknotted."
        err "         No retry. Write the label with"
        err "         pipeline/unknot_label.py (Disc #670)."
        exit 3
    fi
    log "Tracker without active cells at --min-cell-count $MC — the curve IS knotted, lowering further."
done
if [ -z "$MIN_CELL_USED" ]; then
    err "Stage 2 (Tracker) failed — even the lowest threshold gave 0 active cells"
    err "         (rc=1 on every rung, so the tracker DID run through every time)."
    exit 3
fi
if [ "$MIN_CELL_USED" != "2000" ]; then
    log "⚠️ Tracker needed --min-cell-count $MIN_CELL_USED instead of 2000 (Disc #670)."
fi
printf '{"min_cell_count_used": %s, "eps": 0.01, "provenance": "tracker ladder, Disc #670"}\n' \
    "$MIN_CELL_USED" > "$PROJ_DIR/${KNOT}_tracker_params.json"
ok "Tracker done in $(($(date +%s)-T2))s (--min-cell-count $MIN_CELL_USED)"

# ─────────────────────────────────────────────
# 2a. target count + agg range detection
# ─────────────────────────────────────────────
stage "2a" "TARGET-DETECTION (n_cross distribution)"
conda run -n gnottheorie_env python "$PIPE/detect_target.py" "$CROSS" --json > "$DETECT_JSON" \
    || { err "detect_target.py failed"; exit 3; }

TARGET=$(conda run -n gnottheorie_env python -c \
    "import json;d=json.load(open('$DETECT_JSON'));print(d['target_count'])")
FRAME_CAP=$(conda run -n gnottheorie_env python -c \
    "import json;d=json.load(open('$DETECT_JSON'));print(d['frame_cap'])")
AGG_MIN=$(conda run -n gnottheorie_env python -c \
    "import json;d=json.load(open('$DETECT_JSON'));print(d['agg_n_min'])")
AGG_MAX=$(conda run -n gnottheorie_env python -c \
    "import json;d=json.load(open('$DETECT_JSON'));print(d['agg_n_max'])")
ALL_EVEN=$(conda run -n gnottheorie_env python -c \
    "import json;d=json.load(open('$DETECT_JSON'));print('1' if d['all_even_n_cross'] else '0')")

log "Target-Count: $TARGET   Frame-Cap: $FRAME_CAP"
log "Aggressive-Simplify-Range: [$AGG_MIN, $AGG_MAX]"
if [ "$ALL_EVEN" = "1" ]; then
    log "All frame n_cross even → pi antisymmetry signature confirmed."
fi

# ─────────────────────────────────────────────
# 3. Kinematic-Resolve
# ─────────────────────────────────────────────
stage "3" "KINEMATIC-RESOLVE — frame choice + signed Gauss + Regina"
T3=$(date +%s)
conda run -n gnottheorie_env python "$PIPE/kinematic_resolve.py" "$CROSS" \
    --tracks "$TRACKS" \
    --fseries "$FSERIES" \
    --camera "$CAMERA" \
    --N 4096 \
    --target-count "$TARGET" \
    --frame-cap "$FRAME_CAP" \
    --match-tolerance 1.0 \
    $RESOLVE_FOLD_FLAG \
    --output "$KIN" \
    || { err "Stage 3 (Resolve) failed"; exit 3; }
ok "Resolve done in $(($(date +%s)-T3))s"

# ─────────────────────────────────────────────
# 4. Aggressive-Simplify (mandatory per Discovery #91, optionally skippable)
# ─────────────────────────────────────────────
# ENV override SKIP_AGGRESSIVE_SIMPLIFY=1 skips stage 4 + 4a
# completely — the kinematic_resolve diagram goes straight to stage 5+.
# Useful for the Sundowner class (c_median > 100), where stage 4 runs into
# Regina simplify timeouts even with AGG_MAX_FRAMES + AGG_MAX_C_MULT.
# Consequence: crossing_number_min is not minimal (single-frame diagram),
# but the algebra (det, Alex, Jones, HOMFLY, vol, sym) is correct and
# suffices for cluster building. Discovery #147.
#
# STAGE 4 HIGH-c CAP (Werner's go 2026-06-17, analogous to STAGE6_CMAX): at high c
# (c≈54+) stage 4 (Regina aggressive-simplify multi-frame sweep) becomes a
# global_timeout swamp (~880s) and then delivers NO usable best — merge_aggressive
# falls back to the kinematic diagram anyway since the root fix (2026-06-17).
# Instead of paying for the swamp first: for c > STAGE4_CMAX skip stages 4+4a right away
# (= automatic SKIP_AGGRESSIVE_SIMPLIFY for exactly the knots that need it).
# c = kinematic crossing count. FIX 2026-07-04 (Disc #437b): in POLYGON_PD mode NO 'regina'
# field is written into $KIN -> regina.simplified_size=None -> C4=0 -> the cap NEVER fired (verified:
# 5/5 wave knots regina.simplified_size=None, real c in chosen.n_cross_simplified={43,127,29,...}).
# Now chosen.n_cross_simplified first, fallback regina.simplified_size (non-polygon path). Default 0.
C4=$(python3 -c "import json;d=json.load(open('$KIN'));c=d.get('chosen',{}).get('n_cross_simplified');c=c if c is not None else d.get('regina',{}).get('simplified_size');print(int(c) if c is not None else 0)" 2>/dev/null || echo 0)
case "$C4" in (*[!0-9]*) C4=0;; esac
SKIP_S4=0; SKIP_S4_REASON=""
if [ "${SKIP_AGGRESSIVE_SIMPLIFY:-0}" = "1" ]; then
    SKIP_S4=1; SKIP_S4_REASON="SKIP_AGGRESSIVE_SIMPLIFY=1 (Discovery #147)"
elif [ "${STAGE4_CMAX:-0}" -gt 0 ] && [ "$C4" -gt "${STAGE4_CMAX}" ]; then
    SKIP_S4=1; SKIP_S4_REASON="c=$C4 > STAGE4_CMAX=${STAGE4_CMAX} (Hoch-c-Deckel 2026-06-17)"
fi
if [ "$SKIP_S4" = "1" ]; then
    log "Stage 4 + 4a SKIPPED — $SKIP_S4_REASON"
    log "  kinematic_resolve diagram goes straight to stage 5"
    log "  WARNING: crossing_number_min not minimal (upper bound, Disc #231)"
else
    stage "4" "AGGRESSIVE-SIMPLIFY — multi-frame sweep"
    T4=$(date +%s)
    # OOM brake stage 4 (hole A2, OOM 2026-07-01, Hard Rule 6 exception):
    # aggressive_simplify reads crossings.parquet (here 641 MB) into the
    # main process and forks it via mp.Pool(--processes) → every worker
    # COW-balloons to several GB on dense diagrams. At c=42 (raw
    # diagrams median 182 crossings) = OOM/reboot. Throttle pool + frame cap
    # by raw diagram density (TARGET=min(n_cross) as a proxy). Env AGG_*
    # overrides (Werner/r3 keep full control: := assignment).
    if [[ "${TARGET:-}" =~ ^[0-9]+$ ]] && [ "$TARGET" -gt 40 ]; then
        if [ "$TARGET" -gt 80 ]; then
            : "${AGG_PROCESSES:=3}"; : "${AGG_MAX_WORKERS:=6}"
            : "${AGG_MAX_FRAMES:=500}"; : "${AGG_MAX_C_MULT:=1.2}"
        else
            : "${AGG_PROCESSES:=6}"; : "${AGG_MAX_WORKERS:=8}"
            : "${AGG_MAX_FRAMES:=2000}"
        fi
        export AGG_PROCESSES AGG_MAX_WORKERS AGG_MAX_FRAMES
        [ -n "${AGG_MAX_C_MULT:-}" ] && export AGG_MAX_C_MULT
        # #385 GAP (2026-07-01): the aggressive_simplify/merge_aggressive worker calls
        # regina.homfly() — on dense high-c diagrams that explodes exponentially
        # (15GB/50min/frame → OOM/orphan via the timeout path, which #385's signal killer
        # does NOT catch). homfly/jones are ONLY extras there (best_homfly/best_jones), the
        # frame choice runs via size_simp_regina → skip = identical choice, no hog.
        # Only active in the density-throttle regime (TARGET>40); env overrides (:= wins).
        : "${R3_SKIP_HOMFLY:=1}"; export R3_SKIP_HOMFLY
        log "Stage 4 OOM throttle (dense diagram, raw c≈$TARGET): AGG_PROCESSES=$AGG_PROCESSES AGG_MAX_WORKERS=$AGG_MAX_WORKERS AGG_MAX_FRAMES=$AGG_MAX_FRAMES AGG_MAX_C_MULT=${AGG_MAX_C_MULT:-<unset>} R3_SKIP_HOMFLY=$R3_SKIP_HOMFLY"
    fi
    # ENV overrides for dense diagrams (sundowner/poly_c09 class,
    # Discovery #138). At c_ref>100 stage 4 otherwise becomes a 3.5h swamp.
    # Recommendation: AGG_MAX_FRAMES=500 AGG_MAX_C_MULT=1.2 for extremely dense
    # diagrams (c_median>50). If even that is not enough: SKIP_AGGRESSIVE_SIMPLIFY=1.
    AGG_MAX_FRAMES_ARG="${AGG_MAX_FRAMES:-5000}"
    AGG_EXTRA_ARGS=""
    if [ -n "${AGG_MAX_C_MULT:-}" ]; then
        AGG_EXTRA_ARGS="--max-c-multiplier ${AGG_MAX_C_MULT}"
    fi
    log "Stage 4 args: --max-frames $AGG_MAX_FRAMES_ARG $AGG_EXTRA_ARGS"
    conda run -n gnottheorie_env python "$PIPE/aggressive_simplify.py" "$CROSS" \
        --fseries "$FSERIES" \
        --N 4096 \
        --n-min "$AGG_MIN" \
        --n-max "$AGG_MAX" \
        --max-frames "$AGG_MAX_FRAMES_ARG" \
        --processes "${AGG_PROCESSES:-12}" \
        --max-workers "${AGG_MAX_WORKERS:-12}" \
        --output "$AGG" \
        $AGG_EXTRA_ARGS \
        || { err "Stage 4 (Aggressive-Simplify) failed"; exit 3; }
    ok "Aggressive-Simplify done in $(($(date +%s)-T4))s"

    # Discovery #91 / 2026-05-10 Saturn finding: aggressive_simplify often finds
    # a smaller diagram than kinematic_resolve (single-frame
    # stagnation). Patch kinematic_result.json with the aggressive best, if
    # better, BEFORE mandatory_invariants.py (stage 5) starts. Otherwise
    # stage 5 computes on the non-minimal diagram and reports a wrong
    # crossing_number_min.
    stage "4a" "MERGE-AGGRESSIVE — mandatory selection Discovery #91"
    conda run -n gnottheorie_env python "$PIPE/merge_aggressive.py" \
        "$KIN" "$AGG" --in-place \
        || { err "Stage 4a (Merge-Aggressive) failed"; exit 3; }
fi

fi  # Ende POLYGON_PD-Verzweigung (Rezept #3)

# ─────────────────────────────────────────────
# 5. Mandatory invariants (mandatory_invariants; level 4 in the doctrine)
# ─────────────────────────────────────────────
# STAGE 5 POINT-EVAL CAP (Werner's go 2026-06-17, Disc #248 doctrine "value instead of symbol"):
# regina homfly/jones/bracket + snappy HFK + sage symbolic homfly/jones + Khovanov are
# state sums 2^c -> at c≈54-104 timeout (1200s) / OOM (rc=137); EXACTLY the wall at which the
# 17.06 night run lost all 6 stage 5 knots (0 harvest). For c > STAGE5_CMAX stage 5 runs
# in point_eval mode: polynomial-time det (Alexander Δ(-1) + Goeritz + Seifert) / signature / vol /
# sym / Alexander fingerprint — the exponential forms left out. det/θ suffice for
# high-c novelty (type a, Disc #231/#147). C4 = kinematic simplified_size (line 300, in scope).
# Default 0 = OFF (backward compatible). Recommended together with STAGE6_CMAX: STAGE5_CMAX=24.
# Self-supplying (reads c fresh from $KIN, exists in BOTH paths — also POLYGON_PD,
# where C4 from line 300 is not in scope). FIX 2026-07-04 (Disc #437b): chosen.n_cross_simplified
# first (POLYGON_PD has no regina field -> simplified_size=None -> the cap NEVER fired), fallback regina.
C5=$(python3 -c "import json;d=json.load(open('$KIN'));c=d.get('chosen',{}).get('n_cross_simplified');c=c if c is not None else d.get('regina',{}).get('simplified_size');print(int(c) if c is not None else 0)" 2>/dev/null || echo 0)
case "$C5" in (*[!0-9]*) C5=0;; esac
MANDATORY_PE_ARG=""
if [ "${STAGE5_CMAX:-0}" -gt 0 ] && [ "$C5" -gt "${STAGE5_CMAX}" ]; then
    log "Stage 5 POINT-EVAL mode — c=$C5 > STAGE5_CMAX=${STAGE5_CMAX} (Disc #248: homfly/jones/bracket/HFK/khovanov skipped; det via Alexander+Goeritz+Seifert)"
    MANDATORY_PE_ARG="--point-eval"
fi

stage "5" "PFLICHT-INVARIANTEN (mandatory invariants) — multi-env orchestrator"
T5=$(date +%s)
conda run -n gnottheorie_env python "$PIPE/mandatory_invariants.py" \
    --kinematic-result "$KIN" \
    --name "$KNOT" \
    --output "$PFL" \
    ${SKIP_KNOTJOB:+--skip-knotjob} \
    ${MANDATORY_PE_ARG} \
    || { err "Stage 5 (Pflicht) failed"; exit 3; }
ok "Mandatory invariants done in $(($(date +%s)-T5))s"

# STAGE 6 HIGH-c CAP (Werner's go 2026-06-15, monster-timeout lesson): stage 6 (HFK/khovanov/
# length_spectrum/unknotting) is exponential. Skip for c > STAGE6_CMAX — det/θ from stage 5
# suffice for novelty; the knot becomes type a (no _complete.json → sort_and_render line 56), θ via salvage.
# Default 0 = OFF (backward compatible); only active when STAGE6_CMAX is set.
C6=$(python3 -c "import json,ast;d=json.load(open('$PFL'));pd=d.get('pd_list_str');print(len(ast.literal_eval(pd)) if pd else 0)" 2>/dev/null || echo 0)
if [ "${STAGE6_CMAX:-0}" -gt 0 ] && [ "$C6" -gt "${STAGE6_CMAX}" ]; then
    log "Stage 6+ SKIPPED — c=$C6 > STAGE6_CMAX=${STAGE6_CMAX} (high-c cap, det/θ from stage 5). Type a."
    exit 0
fi

# R3 DET ADJUDICATOR: only up to stage 5 (det quartet regina/Jones@i/HFK is in $PFL) —
# stage 6+ (homfly/khovanov/length_spectrum) are exponential + unnecessary for the det
# adjudication. STOP_AFTER_PFLICHT=1 -> clean exit 0 (backward compatible: unset -> continue).
if [ "${STOP_AFTER_PFLICHT:-0}" = "1" ]; then
    log "Stage 6+ skipped — STOP_AFTER_PFLICHT=1 (R3 det adjudicator, det in $PFL)"
    exit 0
fi

# ─────────────────────────────────────────────
# 6. Mandatory invariants, extra (mandatory_invariants_extra)
# ─────────────────────────────────────────────
stage "6" "PFLICHT-EXTRA — khovanov, length_spectrum, s, unknotting, kauffman"
T6=$(date +%s)
# OOM brake (hole A, OOM 2026-07-01, Hard Rule 6 exception = documented path):
# the unknotting pool in mandatory_invariants_extra grows on densified
# high-c knots to 5-9 GB/worker → 12 workers are guaranteed to blow 32 Gi. Throttle
# the worker count by crossing number (TARGET≈c). Env MANDATORY_EXTRA_PROCESSES
# (e.g. via r3_reresolve) overrides and also wins inside the script itself.
if [ -n "${MANDATORY_EXTRA_PROCESSES:-}" ]; then
    PE_PROC="$MANDATORY_EXTRA_PROCESSES"          # only for the log; the script reads the env itself
elif ! [[ "${TARGET:-}" =~ ^[0-9]+$ ]]; then
    PE_PROC=12                                   # TARGET unbekannt → alter Default
elif [ "$TARGET" -gt 36 ]; then
    PE_PROC=2
elif [ "$TARGET" -gt 30 ]; then
    PE_PROC=4
elif [ "$TARGET" -gt 24 ]; then
    PE_PROC=6
else
    PE_PROC=12
fi
log "Pflicht-Extra unknotting pool: $PE_PROC workers (TARGET-c=${TARGET:-?})"
conda run -n gnottheorie_env python "$PIPE/mandatory_invariants_extra.py" "$PFL" \
    --max-unknotting-subset "${R3_MAX_UNKNOT:-4}" \
    --processes "$PE_PROC" \
    || { err "Stage 6 (Pflicht-Extra) failed"; exit 3; }
ok "Pflicht-Extra done in $(($(date +%s)-T6))s"

# ─────────────────────────────────────────────
# 7. Enrichments (in place on _complete.json)
# ─────────────────────────────────────────────
stage "7" "ENRICHMENTS — repair_khovanov + 6x enrich_*"
T7=$(date +%s)
for SCRIPT in repair_khovanov_parsing.py \
              enrich_braid_index_mfw.py \
              enrich_unknotting_signed_gauss.py \
              enrich_bridge_number_bounds.py \
              enrich_length_spectrum.py \
              enrich_optional_invariants.py \
              enrich_non_hyperbolic.py; do
    log "  $SCRIPT …"
    # Disc #289: enrichments are nice-to-have. A failure (monster timeout in
    # length_spectrum/khovanov at c>~40) must NOT discard the knot — the
    # core invariants (stage 5+6) are valid. Before: || exit 3 -> a whole
    # night's harvest of monster knots was lost. Now: log + continue.
    conda run -n gnottheorie_env python "$PIPE/$SCRIPT" "$PFLC" --in-place > /dev/null \
        || err "Stage 7 ($SCRIPT) failed — enrichment, NOT fatal (Disc #289), continuing"
done
ok "Enrichments done in $(($(date +%s)-T7))s"

# ─────────────────────────────────────────────
# 8. Summary
# ─────────────────────────────────────────────
stage "8" "SUMMARY — Markdown overview"
conda run -n gnottheorie_env python "$PIPE/summary_writer.py" "$PFLC" \
    --output "$SUMMARY" \
    || { err "Stage 8 (Summary) failed"; exit 3; }

# ─────────────────────────────────────────────
# 8a. Bar-Natan / van der Veen Theta — Discovery #133, Roadmap 4
# ─────────────────────────────────────────────
# Second, independent identification path (PD -> Cs/phi -> Theta(t1,t2)).
# Default since Discovery #148: numeric path (specialisation t1=22/7, t2=21/13).
# On by default since the "fetch everything that can be fetched: it is mandatory" doctrine 2026-05-14.
# Stage 9 sort_and_render needs theta data for the cluster match.
# Symbolic (full polynomial) explodes in RAM at c>~20 — so symbolic can ONLY
# be selected via THETA_SYMBOLIC=1 (Discovery #157).
# Can be switched off via SKIP_THETA=1.
if [ "${SKIP_THETA:-0}" != "1" ]; then
    stage "8a" "THETA-NUMERIC — Bar-Natan/vdV (sage_stable)"
    THETA_JSON=$PROJ_DIR/${KNOT}_theta_numeric.json
    T8A=$(date +%s)
    if [ "${THETA_SYMBOLIC:-0}" = "1" ]; then
        THETA_SCRIPT="theta_crosscheck.py"
        THETA_JSON=$PROJ_DIR/${KNOT}_theta.json
    else
        THETA_SCRIPT="theta_numeric_crosscheck.py"
    fi
    conda run -n gnottheorie_env python "$PIPE/$THETA_SCRIPT" "$PFLC" \
        --output "$THETA_JSON" \
        || { err "Stage 8a (Theta) failed — mandatory data saved"; }
    ok "Theta crosscheck done in $(($(date +%s)-T8A))s   →  $THETA_JSON   ($THETA_SCRIPT)"
else
    log "Stage 8a (Theta) skipped — SKIP_THETA=1"
fi

# ─────────────────────────────────────────────
# 8b. KnotInfo / snappy / Fremlin Crosscheck — Discoveries #133/#134/#168/#175
# ─────────────────────────────────────────────
# Produces <knot>_crosscheck.json with knotinfo_match, snappy_census_match,
# fremlin_match. Stage 9 sort_and_render needs best_match for
# sorting into bekannt/ (known). On by default since "fetching everything is mandatory".
# Can be switched off via SKIP_CROSSCHECK=1, the snappy subprocess via CROSSCHECK_SKIP_SNAPPY=1.
if [ "${SKIP_CROSSCHECK:-0}" != "1" ]; then
    stage "8b" "CROSSCHECK — KnotInfo + snappy-Census + Fremlin"
    CC_JSON=$PROJ_DIR/${KNOT}_crosscheck.json
    T8B=$(date +%s)
    CC_ARGS=()
    [ "${CROSSCHECK_SKIP_SNAPPY:-0}" = "1" ] && CC_ARGS+=(--skip-snappy)
    conda run -n gnottheorie_env python "$PIPE/pipeline_crosscheck.py" "$PFLC" \
        --output "$CC_JSON" "${CC_ARGS[@]}" \
        || { err "Stage 8b (Crosscheck) failed — mandatory data saved"; }
    ok "Crosscheck done in $(($(date +%s)-T8B))s   →  $CC_JSON"
else
    log "Stage 8b (Crosscheck) skipped — SKIP_CROSSCHECK=1"
fi

# ─────────────────────────────────────────────
# Stage 9: SORT + RENDER (function run_stage9, defined above)
# Sorts the folder according to the mandatory/crosscheck result into:
#   pipeline_failures/  (rc!=0 or empty mandatory core fields)
#   bekannt/<klassisch>__<original>  (KnotInfo hit)
#   S_<c>_d<det>_<sym>_<chir>/       (Stanggassinger cluster, theta-verified)
#   singletons/                      (otherwise)
# Render PNG + GLB if a fseries exists and it is not a KnotInfo knot.
# The EXIT trap ensures that stage 9 also runs after a premature abort.
# Discovery #175/#177 classical + S_ prefix; [[feedback_render_and_sort_on_new_knots]].
# ─────────────────────────────────────────────
run_stage9

# ─────────────────────────────────────────────
# final summary
# ─────────────────────────────────────────────
TOTAL=$(($(date +%s)-T0))
echo
echo -e "\033[1;32m============================================\033[0m"
echo -e "\033[1;32m  PIPELINE DONE  (Total: ${TOTAL}s)\033[0m"
echo -e "\033[1;32m============================================\033[0m"
log "Knot: $KNOT"
log "Outputs:"
log "  $PFLC"
log "  $SUMMARY"
echo
conda run -n gnottheorie_env python -c "
import json
with open('$PFLC') as f: d=json.load(f)
def cl(o):
    n=p=0
    def rec(x):
        nonlocal n,p
        if isinstance(x,dict):
            for v in x.values(): rec(v)
        else:
            n+=1
            if x is not None and x!='' and x!=[] and x!={}: p+=1
    rec(o); return p,n
pp,pn = cl(d['mandatory_list'])
op,on = cl(d['optional_invariants'])
print(f'  Pflichtliste:        {pp}/{pn} populated  ({pp/pn*100:.1f}%)')
print(f'  Optional:            {op}/{on} populated  ({op/on*100:.1f}%)')
print(f'  Konsistenz:          {\"OK\" if d.get(\"consistency_ok\") else \"FAIL\"}')
print()
alg = d['mandatory_list']['ALGEBRA']; diag = d['mandatory_list']['DIAGRAM']; hyp = d['mandatory_list']['HYPERBOLIC']
print(f'  c={diag.get(\"crossing_number_min\")}, det={alg.get(\"determinant\")}, '
      f'vol={hyp.get(\"hyperbolic_volume\"):.6f}, sym={hyp.get(\"symmetry_group\")}, '
      f'amphi={hyp.get(\"is_amphicheiral\")}, fibered={alg.get(\"fibered\")}')
"
