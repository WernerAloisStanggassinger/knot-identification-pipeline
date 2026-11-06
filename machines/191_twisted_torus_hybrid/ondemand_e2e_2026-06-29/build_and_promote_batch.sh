#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# build_and_promote_batch.sh — det=13 completion (5) + det=17 family (16).
# Per tuple: build .fseries (fidelity loop) → point-eval pipeline_run.sh (blind det/vol/sig) →
# standalone theta → salvage record (cp complete + provenance) → into singletons_novelty.
# Rebuild+label+naming afterwards separately (foreground). Werner 2026-07-01: „vervollständige det=13 auf 12, dann det=17" (EN: complete det=13 to 12, then det=17).
set -u
D=${KNOT_PIPELINE_ROOT:?set KNOT_PIPELINE_ROOT to the project root}/maschinen_ideen/191_twisted_torus_hybrid/ondemand_e2e_2026-06-29
GNOTPY=${KNOT_PIPELINE_PYTHON_MAIN:-python3}
PIPE=${KNOT_PIPELINE_ROOT:?set KNOT_PIPELINE_ROOT to the project root}/freeze_manifest/pipeline_scripts
PR=${KNOT_PIPELINE_ROOT:?set KNOT_PIPELINE_ROOT to the project root}/freeze_manifest/pipeline_run.sh
SN=${KNOT_PIPELINE_ROOT:?set KNOT_PIPELINE_ROOT to the project root}/erzeugte_knoten/new_found_knots/singletons_novelty
B="$D/batch_complete"; mkdir -p "$B"
V="$B/verdicts.jsonl"; : > "$V"
TS=$(date -Iseconds)

# "p q r s family"
TUPLES=(
 "13 3 2 6 twisted_torus_det13" "13 5 2 7 twisted_torus_det13" "13 7 2 7 twisted_torus_det13"
 "13 9 2 6 twisted_torus_det13" "13 11 2 7 twisted_torus_det13"
 "5 3 2 9 twisted_torus_det17" "7 3 2 8 twisted_torus_det17" "11 3 2 9 twisted_torus_det17"
 "13 3 2 8 twisted_torus_det17" "7 5 2 8 twisted_torus_det17" "13 5 2 9 twisted_torus_det17"
 "13 10 2 3 twisted_torus_det17" "9 7 2 9 twisted_torus_det17" "9 5 2 9 twisted_torus_det17"
 "13 8 2 5 twisted_torus_det17" "11 5 2 8 twisted_torus_det17" "11 9 2 8 twisted_torus_det17"
 "11 7 2 9 twisted_torus_det17" "13 11 2 9 twisted_torus_det17" "13 7 2 9 twisted_torus_det17"
 "13 9 2 8 twisted_torus_det17"
)

build_fseries() {  # p q r s outpath -> echo "n err c"
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

echo "BUILD+PROMOTE START $(date '+%F %T')  (${#TUPLES[@]} tuples)"
for spec in "${TUPLES[@]}"; do
  read p q r s fam <<< "$spec"
  stem="T_${p}_${q}_${r}_${s}"
  wd="$B/$stem"; mkdir -p "$wd"
  echo "==== $stem (fam=$fam) START $(date '+%T') ===="
  # 1) build .fseries
  fs="$wd/${stem}.fseries"
  read n err < <(build_fseries "$p" "$q" "$r" "$s" "$fs")
  echo "  fseries n=$n roundtrip=$err"
  # 2) point-eval pipeline (blind), skip sort/render (we place it ourselves)
  ( cd "$wd" && STAGE5_CMAX=24 STAGE6_CMAX=24 SKIP_SORT_RENDER=1 bash "$PR" > "$wd/pipeline.log" 2>&1 )
  rc=$?
  pj=$(ls "$wd"/*_pflicht_invarianten.json 2>/dev/null | grep -v complete | head -1)
  # 3) standalone theta
  if [ -n "$pj" ]; then
    "$GNOTPY" "$PIPE/theta_numeric_crosscheck.py" "$pj" --output "$wd/${stem}_theta_numeric.json" > "$wd/theta.log" 2>&1
    # 4) salvage complete + provenance
    cp "$pj" "$wd/${stem}_pflicht_invarianten_complete.json"
    cat > "$wd/salvaged_from_failure.json" <<EOF
{"salvaged_from_stage5": true, "salvage_date": "$TS", "mode": "point-eval + standalone theta", "reason": "high-c novel-hyp on-demand-dial; exponentielle Stage-6-Felder ausgelassen, det/vol/sig blind + theta kanonisch. Type-a Disc #231/#147.", "family": "$fam"}
EOF
    # 5) into singletons_novelty
    dst="$SN/$stem"; mkdir -p "$dst"
    cp "$fs" "$dst/" 2>/dev/null
    cp "$wd/${stem}_pflicht_invarianten.json" "$wd/${stem}_pflicht_invarianten_complete.json" "$wd/${stem}_theta_numeric.json" "$wd/salvaged_from_failure.json" "$dst/" 2>/dev/null
    det=$("$GNOTPY" -c "import json;d=json.load(open('$pj'));print(d['pflichtliste']['ALGEBRA'].get('determinant_alexander') or d['pflichtliste']['ALGEBRA'].get('determinant'))" 2>/dev/null)
    vol=$("$GNOTPY" -c "import json;d=json.load(open('$pj'));print(d['pflichtliste']['HYPERBOLISCH'].get('hyperbolic_volume'))" 2>/dev/null)
    sig=$("$GNOTPY" -c "import json;d=json.load(open('$pj'));print(d['pflichtliste']['ALGEBRA'].get('signature'))" 2>/dev/null)
    echo "{\"stem\":\"$stem\",\"fam\":\"$fam\",\"rc\":$rc,\"det\":\"$det\",\"vol\":\"$vol\",\"sig\":\"$sig\",\"n_harm\":$n,\"roundtrip\":\"$err\"}" >> "$V"
    echo "==== $stem END rc=$rc det=$det vol=$vol sig=$sig $(date '+%T') ===="
  else
    echo "{\"stem\":\"$stem\",\"fam\":\"$fam\",\"rc\":$rc,\"error\":\"no pflicht.json\"}" >> "$V"
    echo "==== $stem END rc=$rc NO-PFLICHT $(date '+%T') ===="
  fi
done
echo "BUILD+PROMOTE DONE $(date '+%F %T')"
