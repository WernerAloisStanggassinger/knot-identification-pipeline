#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# batch_r3.sh — det=19 (14) + det=15 (13) + det=11 (12) = 39 knots, on-demand-dial.
# Per tuple: .fseries (Fidelity-Loop) → point-eval pipeline_run.sh (blind det/vol/sig) → standalone theta →
# Salvage record (cp complete + provenance, family from spec) → singletons_novelty. Afterwards clear parquet.
# Fix vs #389: salvage inline via `find` (no grep on path), folder without 'complete' in the name.
set -u
D=${KNOT_PIPELINE_ROOT:?set KNOT_PIPELINE_ROOT to the project root}/maschinen_ideen/191_twisted_torus_hybrid/ondemand_e2e_2026-06-29
GNOTPY=${KNOT_PIPELINE_PYTHON_MAIN:-python3}
PIPE=${KNOT_PIPELINE_ROOT:?set KNOT_PIPELINE_ROOT to the project root}/freeze_manifest/pipeline_scripts
PR=${KNOT_PIPELINE_ROOT:?set KNOT_PIPELINE_ROOT to the project root}/freeze_manifest/pipeline_run.sh
SN=${KNOT_PIPELINE_ROOT:?set KNOT_PIPELINE_ROOT to the project root}/erzeugte_knoten/new_found_knots/singletons_novelty
B="$D/batch_r3"; mkdir -p "$B"
V="$B/verdicts.jsonl"; : > "$V"
TS=$(date -Iseconds)

TUPLES=(
 "5 3 2 10 twisted_torus_det19" "7 3 2 9 twisted_torus_det19" "11 3 2 10 twisted_torus_det19" "13 3 2 9 twisted_torus_det19" "7 5 2 9 twisted_torus_det19" "13 5 2 10 twisted_torus_det19" "9 7 2 10 twisted_torus_det19" "9 5 2 10 twisted_torus_det19" "11 5 2 9 twisted_torus_det19" "11 9 2 9 twisted_torus_det19" "11 7 2 10 twisted_torus_det19" "13 11 2 10 twisted_torus_det19" "13 7 2 10 twisted_torus_det19" "13 9 2 9 twisted_torus_det19"
 "5 3 2 8 twisted_torus_det15" "11 3 2 8 twisted_torus_det15" "13 3 2 7 twisted_torus_det15" "7 5 2 7 twisted_torus_det15" "13 5 2 8 twisted_torus_det15" "9 7 2 8 twisted_torus_det15" "9 5 2 8 twisted_torus_det15" "11 5 2 7 twisted_torus_det15" "11 9 2 7 twisted_torus_det15" "11 7 2 8 twisted_torus_det15" "13 11 2 8 twisted_torus_det15" "13 9 2 7 twisted_torus_det15" "13 7 2 8 twisted_torus_det15"
 "11 3 2 6 twisted_torus_det11" "7 5 2 5 twisted_torus_det11" "13 5 2 6 twisted_torus_det11" "9 7 2 6 twisted_torus_det11" "13 8 2 4 twisted_torus_det11" "11 5 2 5 twisted_torus_det11" "9 5 2 6 twisted_torus_det11" "11 9 2 5 twisted_torus_det11" "11 7 2 6 twisted_torus_det11" "13 11 2 6 twisted_torus_det11" "13 9 2 5 twisted_torus_det11" "13 7 2 6 twisted_torus_det11"
)

build_fseries() {  # p q r s outpath -> echo "n err"
  local p=$1 q=$2 r=$3 s=$4 out=$5
  local n=$(( 96 > 4*((p-1)*q+(r-1)*r*s)/p ? 96 : 4*((p-1)*q+(r-1)*r*s)/p ))
  local o err
  while [ "$n" -le 2048 ]; do
    o=$(cd "$D" && "$GNOTPY" braid_to_fseries.py "$p" "$q" "$r" "$s" "$n" "$out" 2>/dev/null)
    err=$(echo "$o" | grep -oE 'roundtrip_relerr=[0-9.]+' | cut -d= -f2)
    awk "BEGIN{exit !($err<0.03)}" && break
    n=$((n*2))
  done
  echo "$n $err"
}

echo "BATCH_R3 START $(date '+%F %T')  (${#TUPLES[@]} Tupel)"
for spec in "${TUPLES[@]}"; do
  read p q r s fam <<< "$spec"
  stem="T_${p}_${q}_${r}_${s}"
  wd="$B/$stem"; mkdir -p "$wd"
  echo "==== $stem ($fam) START $(date '+%T') ===="
  fs="$wd/${stem}.fseries"
  read n err < <(build_fseries "$p" "$q" "$r" "$s" "$fs")
  ( cd "$wd" && STAGE5_CMAX=24 STAGE6_CMAX=24 SKIP_SORT_RENDER=1 bash "$PR" > "$wd/pipeline.log" 2>&1 )
  rc=$?
  pj=$(find "$wd" -maxdepth 1 -name "${stem}_pflicht_invarianten.json" | head -1)
  if [ -n "$pj" ]; then
    "$GNOTPY" "$PIPE/theta_numeric_crosscheck.py" "$pj" --output "$wd/${stem}_theta_numeric.json" > "$wd/theta.log" 2>&1
    cp "$pj" "$wd/${stem}_pflicht_invarianten_complete.json"
    cat > "$wd/salvaged_from_failure.json" <<EOF
{"salvaged_from_stage5": true, "salvage_date": "$TS", "mode": "point-eval + standalone theta", "reason": "on-demand-dial novel-hyp; Stage-6/HOMFLY ausgelassen, det/vol/sig blind + theta kanonisch (Type-a Disc #231/#147).", "family": "$fam"}
EOF
    dst="$SN/$stem"; mkdir -p "$dst"
    cp "$fs" "$wd/${stem}_pflicht_invarianten.json" "$wd/${stem}_pflicht_invarianten_complete.json" "$wd/${stem}_theta_numeric.json" "$wd/salvaged_from_failure.json" "$dst/" 2>/dev/null
    det=$("$GNOTPY" -c "import json;d=json.load(open('$pj'))['pflichtliste']['ALGEBRA'];print(d.get('determinant_alexander') or d.get('determinant'))" 2>/dev/null)
    vol=$("$GNOTPY" -c "import json;print(json.load(open('$pj'))['pflichtliste']['HYPERBOLISCH'].get('hyperbolic_volume'))" 2>/dev/null)
    sig=$("$GNOTPY" -c "import json;print(json.load(open('$pj'))['pflichtliste']['ALGEBRA'].get('signature'))" 2>/dev/null)
    echo "{\"stem\":\"$stem\",\"fam\":\"$fam\",\"det\":\"$det\",\"vol\":\"$vol\",\"sig\":\"$sig\",\"n\":$n,\"rt\":\"$err\"}" >> "$V"
    echo "==== $stem END det=$det vol=$vol sig=$sig $(date '+%T') ===="
  else
    echo "{\"stem\":\"$stem\",\"fam\":\"$fam\",\"error\":\"no pflicht\",\"rc\":$rc}" >> "$V"
    echo "==== $stem END rc=$rc NO-PFLICHT $(date '+%T') ===="
  fi
  rm -f "$wd"/*.parquet
done
echo "BATCH_R3 DONE $(date '+%F %T')  verdicts: $(wc -l < "$V")/${#TUPLES[@]}"
