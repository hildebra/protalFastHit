#!/usr/bin/env bash
# Profiling the 5M-pair SAM: f5ee645 (par2) against the parsing + decompression build (parse) with 0-3 decompression
# threads (PROTAL_EXPERIMENT_DECOMPRESS_THREADS; 0: the reading thread decompresses, as before), alternated, niced.
set -uo pipefail
cp $HOME/mt-work/parse/src/build/protal $HOME/mt-work/bin/protal-parse
W=$HOME/mt-work/speed2; mkdir -p $W
SAM=$HOME/mt-audit/runs/b6/s1.sam.zst
DB=$HOME/bench071/V071/protal_db
[ -f $W/runs.tsv ] || printf "binary\tthreads\tdecompress\trep\tprofiling_s\tload\n" > $W/runs.tsv
seconds() { awk '{ s = 0; for (i = 1; i <= NF; i++) { v = $i; if (v ~ /ms$/) { sub(/ms/, "", v); s += v / 1000 } else if (v ~ /m$/) { sub(/m/, "", v); s += v * 60 } else if (v ~ /s$/) { sub(/s/, "", v); s += v } } printf "%.3f", s }'; }
run() {  # run BIN THREADS DECOMPRESS REP
  local bin=$1 t=$2 d=$3 rep=$4 o=$W/out
  local load=$(cut -d' ' -f1 /proc/loadavg)
  rm -rf $o
  if [ "$d" = "-" ]; then nice -n 5 $HOME/mt-work/bin/protal-$bin --db $DB --profile_only $SAM --prefix pe5M -o $o -t $t --no_qcmsa > $W/log 2>&1
  else PROTAL_EXPERIMENT_DECOMPRESS_THREADS=$d nice -n 5 $HOME/mt-work/bin/protal-$bin --db $DB --profile_only $SAM --prefix pe5M -o $o -t $t --no_qcmsa > $W/log 2>&1; fi
  local p=$(grep 'Profiling took' $W/log | sed 's/.*took //' | seconds)
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" $bin $t $d $rep $p $load | tee -a $W/runs.tsv
}
for rep in 1 2 3; do
  run par2 1 - $rep; run parse 1 0 $rep
  run par2 3 - $rep; run parse 3 0 $rep; run parse 3 1 $rep
  run par2 6 - $rep; run parse 6 0 $rep; run parse 6 1 $rep; run parse 6 2 $rep; run parse 6 3 $rep
done
