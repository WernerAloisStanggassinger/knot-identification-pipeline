#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2024-2026 Werner Alois Stanggassinger.
# Authored in collaboration with Claude Code (Anthropic).
# The eleven still open cases with |a| <= 5. 12n504 is decided and drops out.
# VORHER=3: the combinations already checked at |a| <= 3 are skipped.
cd ${KNOT_PIPELINE_ROOT:-/media/squeezy/kauffman_storage}/maschinen_ideen/RS4_12N706_20260919
P=${KNOT_PIPELINE_SAGE_ENV_BIN:-/home/squeezy/miniforge3/envs/sage_stable/bin}
for Z in 12a279 12a786 12a646 12a1083 11a326 12n49 12a3 12a715 11a164 12a173 11a58; do
  echo "=========== $Z  $(date +%H:%M:%S) ==========="
  PATH=$P:$PATH PYTHONUNBUFFERED=1 $P/python zeugen5.py "$Z" 5 3 2>&1 \
    | stdbuf -oL grep --line-buffered -E "ZIEL|Aufgaben$|ZEUGE|FERTIG|==>|\[[0-9]+/"
done
echo "SWEEP-G5 FERTIG $(date +%H:%M:%S)"
