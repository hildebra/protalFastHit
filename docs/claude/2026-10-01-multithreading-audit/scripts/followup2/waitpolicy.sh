#!/usr/bin/env bash
# Input-path audit: do threads that spin on critical(reader) starve the inflating threads? protal-inexp (a212559,
# batch 32) at 6 threads on the 6 vCPUs, libgomp's default wait policy (spin up to GOMP_SPINCOUNT 300k before
# sleeping, as OpenMP counts 6 threads on 6 CPUs and not the 2 inflating std::threads) against
# OMP_WAIT_POLICY=passive (sleep at once), alternated; mini and full database, BGZF / gzip / uncompressed.
# Rows into ~/mt-work/input/waitpolicy.tsv.
set -uo pipefail
W=$HOME/mt-work/input
B=$HOME/bench071
declare -A R1=( [bgzf]=$B/samples_deep/points/rl150_p5000000/sim/reads/rl150_p5000000_s_1_R1.fq.gz
                [gzip]=$HOME/mt-audit/plain/R1.std.gz [plain]=$HOME/mt-audit/plain/R1.fq )
declare -A R2=( [bgzf]=$B/samples_deep/points/rl150_p5000000/sim/reads/rl150_p5000000_s_1_R2.fq.gz
                [gzip]=$HOME/mt-audit/plain/R2.std.gz [plain]=$HOME/mt-audit/plain/R2.fq )
declare -A DB=( [mini]=$HOME/mt-work/rbfull/mini_db/protal_db [full]=$B/V071/protal_db )
[ -f $W/waitpolicy.tsv ] || printf "db\tinput\tthreads\tpolicy\trep\taligning_s\tpairs_per_s\twall_s\tuser_s\tsys_s\tload\treader_s_per_thread\n" > $W/waitpolicy.tsv
seconds() { awk '{ s = 0; for (i = 1; i <= NF; i++) { v = $i; if (v ~ /ms$/) { sub(/ms/, "", v); s += v / 1000 } else if (v ~ /m$/) { sub(/m/, "", v); s += v * 60 } else if (v ~ /s$/) { sub(/s/, "", v); s += v } } printf "%.3f", s }'; }
quiet() { local i; for i in $(seq 6); do awk -v l=$(cut -d" " -f1 /proc/loadavg) "BEGIN { exit !(l < 6) }" && return; sleep 10; done; }
run() {  # run DB INPUT THREADS POLICY REP
  quiet
  local db=$1 in=$2 t=$3 policy=$4 rep=$5 d=$W/run
  rm -rf $d; mkdir -p $d
  local load=$(cut -d' ' -f1 /proc/loadavg)
  local envs=(PROTAL_EXPERIMENT_BATCH=32)
  [ $policy = passive ] && envs+=(OMP_WAIT_POLICY=passive)
  env "${envs[@]}" nice -n 5 /usr/bin/time -f "%e %U %S" -o $d.time $HOME/mt-work/bin/protal-inexp --db ${DB[$db]} \
      -1 ${R1[$in]} -2 ${R2[$in]} --prefix s -o $d -t $t --no_profile > $d.log 2>&1
  local a=$(grep 'Aligning reads took' $d.log | sed 's/.*took //' | seconds)
  local reader=$(awk -F'\t' '$1 == "Sequence reader" { print $4 }' $d/misc/s_runtime.tsv)
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $db $in $t $policy $rep $a $(awk -v a=$a 'BEGIN { printf "%.0f", 4995812 / a }') \
      $(tail -1 $d.time) $load $reader | tee -a $W/waitpolicy.tsv
}
for rep in 1 2 3; do
  cat ${R1[@]} ${R2[@]} > /dev/null
  for in in bgzf gzip plain; do for p in default passive; do run mini $in 6 $p $rep; done; done
  for in in bgzf gzip plain; do for p in default passive; do run mini $in 4 $p $rep; done; done
  for in in bgzf gzip plain; do for p in default passive; do run full $in 6 $p $rep; done; done
done
echo WAITPOLICY DONE
