#!/usr/bin/env bash
# A mini database of close relatives: species 0.5% from their genus (about 1% apart), strains 0.2%.
set -euo pipefail
B=${B:-$HOME/protal-lr-build}  # holds src/ (a checkout), build/protal, mini_db/
cd "$B/src"
PROTAL=$B/build/protal bash scripts/mini_db/build_mini_db.sh "$B/mini_db_close" --species_divergence 0.005 --strain_divergence 0.002 > "$B/minidb_close.log" 2>&1 || { tail -20 "$B/minidb_close.log"; exit 1; }
tail -2 "$B/minidb_close.log"
