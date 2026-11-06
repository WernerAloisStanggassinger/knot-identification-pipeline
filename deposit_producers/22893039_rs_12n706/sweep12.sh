#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# The twelve open cases one after another. One witness per case suffices (GEGENWIND §4).
cd ${KNOT_PIPELINE_ROOT:-/media/squeezy/kauffman_storage}/maschinen_ideen/RS4_12N706_20260919
P=${KNOT_PIPELINE_SAGE_ENV_BIN:-/home/squeezy/miniforge3/envs/sage_stable/bin}
for Z in 11a58 11a164 11a326 12n49 12n504 12a3 12a173 12a279 12a646 12a715 12a786 12a1083; do
  echo "=========== $Z  $(date +%H:%M:%S) ==========="
  PATH=$P:$PATH PYTHONUNBUFFERED=1 $P/python zeugen5.py "$Z" 3 2>&1 | grep -E "ZIEL|Aufgaben$|ZEUGE|FERTIG|==>"
done
echo "SWEEP FERTIG $(date +%H:%M:%S)"
