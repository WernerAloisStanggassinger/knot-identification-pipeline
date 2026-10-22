#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# Contact: stanggassinger@stainlesssteel4u.de
# Full license: LICENSE.txt at repository root (GNU General Public License v3 or later).
# Part of the Knot-Identification Pipeline.
# mem_watchdog.sh — memory watchdog for pipeline sweeps.
#
# Background: Discovery #214 (2026-05-17). systemd-oomd reacts only at the
# OOM condition and takes the whole vte terminal scope with it. We want to
# terminate early, before that happens: at MemAvailable < THRESH
# the whole process group is killed cleanly + a watchdog log is
# written.
#
# Two modes:
#
# (A) Standalone CLI:
#     bash mem_watchdog.sh <target_pid> [threshold_gib] [interval_sec]
#       - default threshold: 3 GiB MemAvailable
#       - default interval: 5 sec
#       - log: /tmp/mem_watchdog_<target_pid>.log
#
# (B) Library (in sweep_wave*.sh):
#     source /media/squeezy/kauffman_storage/freeze_manifest/pipeline_scripts/mem_watchdog.sh
#     start_watchdog 3   # threshold 3 GiB
#     # ... do work ...
#     stop_watchdog
#
# Doctrine: the memory watchdog is *early* (3 GiB margin), not reactive.
# systemd-oomd only comes in at <100 MiB MemAvailable, by then it is too late.

set -uo pipefail

# --- Library mode: only export functions ---

start_watchdog() {
    # Args: [threshold_gib=3] [interval_sec=5]
    local thresh_gib=${1:-3}
    local interval=${2:-5}
    local target_pgid=$$
    local logfile="/tmp/mem_watchdog_${target_pgid}.log"
    export MEM_WATCHDOG_LOG="$logfile"
    echo "=== mem_watchdog START $(date -Is) target_pgid=$target_pgid thresh=${thresh_gib}G interval=${interval}s ===" > "$logfile"
    (
        local thresh_kib=$((thresh_gib * 1024 * 1024))
        local soft_kib=$((thresh_gib * 1024 * 1024 * 5 / 3))  # warn at 5/3*thresh
        while sleep "$interval"; do
            if ! kill -0 "$target_pgid" 2>/dev/null; then
                echo "[$(date -Is)] target $target_pgid gone, watchdog exiting" >> "$logfile"
                exit 0
            fi
            local avail_kib
            avail_kib=$(awk '/^MemAvailable:/{print $2; exit}' /proc/meminfo)
            if [[ -z "$avail_kib" ]]; then
                echo "[$(date -Is)] WARN: could not read MemAvailable" >> "$logfile"
                continue
            fi
            local avail_g=$(awk -v k="$avail_kib" 'BEGIN{printf "%.2f", k/1024/1024}')
            if (( avail_kib < thresh_kib )); then
                echo "[$(date -Is)] CRITICAL: MemAvailable=${avail_g}G < ${thresh_gib}G — killing process group $target_pgid" >> "$logfile"
                # Send TERM to the entire process group (negative PID)
                kill -TERM -- -"$target_pgid" 2>/dev/null || true
                sleep 3
                kill -KILL -- -"$target_pgid" 2>/dev/null || true
                echo "[$(date -Is)] kill complete, watchdog exiting" >> "$logfile"
                exit 1
            elif (( avail_kib < soft_kib )); then
                echo "[$(date -Is)] WARN: MemAvailable=${avail_g}G nearing threshold ${thresh_gib}G" >> "$logfile"
            fi
        done
    ) &
    export MEM_WATCHDOG_PID=$!
    echo "[mem_watchdog] started PID=$MEM_WATCHDOG_PID log=$logfile thresh=${thresh_gib}G" >&2
}

stop_watchdog() {
    if [[ -n "${MEM_WATCHDOG_PID:-}" ]] && kill -0 "$MEM_WATCHDOG_PID" 2>/dev/null; then
        kill -TERM "$MEM_WATCHDOG_PID" 2>/dev/null || true
        wait "$MEM_WATCHDOG_PID" 2>/dev/null || true
        echo "[$(date -Is)] watchdog stopped by parent" >> "${MEM_WATCHDOG_LOG:-/dev/null}"
        echo "[mem_watchdog] stopped" >&2
    fi
    unset MEM_WATCHDOG_PID MEM_WATCHDOG_LOG
}

# --- Standalone CLI mode ---

# If sourced (BASH_SOURCE[0] != $0), just expose functions and return.
# If executed directly, run CLI mode.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    if [[ $# -lt 1 ]]; then
        echo "Usage: $0 <target_pid> [threshold_gib=3] [interval_sec=5]" >&2
        echo "Library mode: source $0; start_watchdog [thresh] [interval]; ... ; stop_watchdog" >&2
        exit 2
    fi
    target_pid="$1"
    thresh_gib="${2:-3}"
    interval="${3:-5}"
    if ! kill -0 "$target_pid" 2>/dev/null; then
        echo "Target PID $target_pid not running" >&2
        exit 1
    fi
    logfile="/tmp/mem_watchdog_${target_pid}.log"
    thresh_kib=$((thresh_gib * 1024 * 1024))
    soft_kib=$((thresh_gib * 1024 * 1024 * 5 / 3))
    echo "=== mem_watchdog START $(date -Is) target_pid=$target_pid thresh=${thresh_gib}G interval=${interval}s ===" | tee "$logfile"
    while sleep "$interval"; do
        if ! kill -0 "$target_pid" 2>/dev/null; then
            echo "[$(date -Is)] target $target_pid gone, watchdog exiting" | tee -a "$logfile"
            exit 0
        fi
        avail_kib=$(awk '/^MemAvailable:/{print $2; exit}' /proc/meminfo)
        avail_g=$(awk -v k="$avail_kib" 'BEGIN{printf "%.2f", k/1024/1024}')
        if (( avail_kib < thresh_kib )); then
            echo "[$(date -Is)] CRITICAL: MemAvailable=${avail_g}G < ${thresh_gib}G — killing $target_pid + group" | tee -a "$logfile"
            # In CLI mode we don't know the group; try killing PID first, then any descendants via pkill
            kill -TERM "$target_pid" 2>/dev/null || true
            pkill -TERM -P "$target_pid" 2>/dev/null || true
            sleep 3
            kill -KILL "$target_pid" 2>/dev/null || true
            pkill -KILL -P "$target_pid" 2>/dev/null || true
            echo "[$(date -Is)] kill complete, watchdog exiting" | tee -a "$logfile"
            exit 1
        elif (( avail_kib < soft_kib )); then
            echo "[$(date -Is)] WARN: MemAvailable=${avail_g}G nearing threshold ${thresh_gib}G" | tee -a "$logfile"
        fi
    done
fi
