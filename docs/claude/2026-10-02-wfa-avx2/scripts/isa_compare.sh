#!/usr/bin/env bash
# The baseline build (protal: x86-64 plus its run-time AVX2 kernels) against protal_avx2 (main.cpp for x86-64-v3),
# both of f07acb3's code (the samroomhead tree: 5c7d64c + 8798733's change; WFA2's portable kernels in both),
# alternated, niced, after the full check of the WFA2 change:
#  - aligning the 500k-pair sample on one thread (--no_profile), five times each;
#  - aligning the 5M-pair sample and Nanopore 90 Mb on six threads (--no_profile), five times each;
#  - profiling the 5M-pair SAM on six threads (--profile_only), seven times each.
# Rows into ~/mt-work/isacmp/runs.tsv.
set -uo pipefail
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad
until grep -q "WFACHECK DONE" $SP/wfa_check.out 2>/dev/null; do sleep 30; done
W=$HOME/mt-work/isacmp; rm -rf $W; mkdir -p $W
(cd $HOME/mt-work/samroomhead/src && nice cmake --build build --target protal protal_avx2 -j 4 > $W/build.log 2>&1) || { echo "FAIL build"; exit 1; }
declare -A BIN=( [baseline]=$HOME/mt-work/samroomhead/src/build/protal [avx2]=$HOME/mt-work/samroomhead/src/build/protal_avx2 )
DB=$HOME/bench071/V071/protal_db; P=$HOME/bench071/samples/points; D=$HOME/bench071/samples_deep/points
declare -A ARGS=(
  [pe500k]="-1 $P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz -2 $P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz --no_profile"
  [pe5M]="-1 $D/rl150_p5000000/sim/reads/rl150_p5000000_s_1_R1.fq.gz -2 $D/rl150_p5000000/sim/reads/rl150_p5000000_s_1_R2.fq.gz --no_profile"
  [ont90M]="-1 $P/ont_b90000000/sim/reads/ont_b90000000_s_1.fq.gz --read_type ont --no_profile"
  [prof5M]="--profile_only $HOME/mt-audit/runs/b6/s1.sam.zst" )
seconds() { awk '{ s = 0; for (i = 1; i <= NF; i++) { v = $i; if (v ~ /ms$/) { sub(/ms/, "", v); s += v / 1000 } else if (v ~ /m$/) { sub(/m/, "", v); s += v * 60 } else if (v ~ /s$/) { sub(/s/, "", v); s += v } } printf "%.3f", s }'; }
printf "case\tthreads\tbinary\trep\tstage_s\twall_s\tuser_s\tload\n" > $W/runs.tsv
run() {  # run CASE THREADS BINARY REP
  local c=$1 t=$2 b=$3 rep=$4 d=$W/run; rm -rf $d
  local load=$(cut -d' ' -f1 /proc/loadavg)
  nice -n 5 /usr/bin/time -f "%e %U" -o $W/time ${BIN[$b]} --db $DB ${ARGS[$c]} --prefix s -o $d -t $t --no_qcmsa > $W/log 2>&1 || echo "FAIL $c $b"
  local stage=$(grep -E 'Aligning reads took|Profiling took' $W/log | head -1 | sed 's/.*took //' | seconds)
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $c $t $b $rep $stage "$(tail -1 $W/time | tr ' ' '\t')" $load | tee -a $W/runs.tsv
}
for rep in 1 2 3 4 5; do
  for b in baseline avx2; do run pe500k 1 $b $rep; done
  for c in pe5M ont90M; do for b in baseline avx2; do run $c 6 $b $rep; done; done
done
for rep in 1 2 3 4 5 6 7; do for b in baseline avx2; do run prof5M 6 $b $rep; done; done
echo ISACMP DONE
