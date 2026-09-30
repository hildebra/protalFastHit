#!/usr/bin/env bash
# Simulate the reads of runs A, B, C (short and long fragments) and run protal (with qcmsa) on each.
# Usage: run_sim.sh [run ...]   (default: A B Cs Cl)
set -euo pipefail
ACC=$HOME/audit5/accuracy
BIN=$HOME/audit5/bin
DB=$HOME/audit5/world/protal_db
export PROTAL_QCMSA_SCRIPT=$HOME/audit5/src/scripts/qcmsa.py
T=2
runs=("$@"); [ ${#runs[@]} -eq 0 ] && runs=(A B Cs Cl)
for r in "${runs[@]}"; do
  case $r in
    A)  man=A; frag="--fragment_mean 350 --fragment_stdev 50";;
    B)  man=B; frag="--fragment_mean 350 --fragment_stdev 50";;
    Cs) man=C; frag="--fragment_mean 200 --fragment_stdev 20";;
    Cl) man=C; frag="--fragment_mean 600 --fragment_stdev 50";;
    *) echo "unknown run $r"; exit 1;;
  esac
  sim=$ACC/sim_$r; out=$ACC/prot_$r
  if [ ! -f "$sim/protal.meta" ]; then
    rm -rf "$sim"
    /usr/bin/time -f "sim $r: %e s" $BIN/simulate_metagenomes --from_manifest $ACC/design/$man.manifest.tsv \
      -o "$sim" --protal_metafile "$out" -t $T $frag > $ACC/logs/sim_$r.log 2>&1
  fi
  if [ ! -d "$out/strains" ]; then
    /usr/bin/time -f "protal $r: %e s" $BIN/protal --db $DB --map "$sim/protal.meta" -t $T > $ACC/logs/protal_$r.log 2>&1 || echo "protal $r exit $?"
  fi
done
