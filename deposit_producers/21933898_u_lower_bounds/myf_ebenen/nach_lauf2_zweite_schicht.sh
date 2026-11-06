#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# When lauf2 is through, its 12 cores do NOT go back to theta (that stands on nice 19),
# but to the second layer of lauf3 — the same heap, unrolled from the back.
# Absolute paths, no CWD footgun (rule 5).
set -u

HIER=${KNOT_PIPELINE_ROOT:-/media/squeezy/kauffman_storage}/konversationen/potsdam_gmi2026/myf_ebenen
PY=${KNOT_PIPELINE_SAGE_ENV_BIN:-/home/squeezy/miniforge3/envs/sage_stable/bin}/python
ARBEITER=${1:-6}
BUDGET=${2:-900}
LOG="$HIER/lauf3/zweite_schicht.log"

mkdir -p "$HIER/lauf3"
echo "[$(date '+%F %T')] waiting for the end of lauf2 (neue_obere_schranken2.py)" >> "$LOG"
while pgrep -f "neue_obere_schranken2.py" > /dev/null; do
    sleep 60
done
echo "[$(date '+%F %T')] run 2 is done — starting ${ARBEITER} workers backwards, budget ${BUDGET}s" >> "$LOG"

for w in $(seq 0 $((ARBEITER-1))); do
    PATH=${KNOT_PIPELINE_SAGE_ENV_BIN:-/home/squeezy/miniforge3/envs/sage_stable/bin}:$PATH \
    nohup "$PY" "$HIER/zeugen_andere_diagramme.py" "$w" "$ARBEITER" "$BUDGET" rueckwaerts \
        > "$HIER/lauf3/log_r${w}.txt" 2>&1 &
    echo "  r$w PID $!" >> "$LOG"
done
echo "[$(date '+%F %T')] zweite Schicht steht" >> "$LOG"
wait
echo "[$(date '+%F %T')] ZWEITE SCHICHT FERTIG" >> "$LOG"
