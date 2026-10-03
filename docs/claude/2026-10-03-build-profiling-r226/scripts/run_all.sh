#!/bin/bash
# Runs the measurements of this report in WSL, niced: sim_cost.sh, then long_cost.py on its first 1000-pair
# community. usage: run_all.sh OUT
set -uo pipefail
out=${1:-$HOME/bprof26}
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$out"
uptime > "$out/load_before.txt"
nice -n 10 bash "$here/sim_cost.sh" "$out/sim" > "$out/sim_cost.out" 2>&1
nice -n 10 python3 "$here/long_cost.py" "$out/sim/p1000/manifest.tsv" "$out/long" --bases 30000000 > "$out/long_cost.out" 2>&1
uptime > "$out/load_after.txt"
echo ALL_DONE >> "$out/long_cost.out"
