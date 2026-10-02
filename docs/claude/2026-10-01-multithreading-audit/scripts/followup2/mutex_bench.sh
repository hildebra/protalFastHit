#!/usr/bin/env bash
# The reader lock as a std::mutex (mutex.patch on a212559) against a212559 itself (protal-inexp at batch 32, the
# same code), at libgomp's default wait policy, alternated, niced: mini database at 4 and 6 threads, full database
# at 6 threads, BGZF / gzip / uncompressed. Rows into ~/mt-work/input/mutex.tsv.
set -uo pipefail
SP=$(cd "$(dirname "$0")" && pwd)  # this folder
W=$HOME/mt-work/input
bash $SP/build_wt.sh mutex $SP/mutex.patch a212559 protal > $W/mutex_build.out 2>&1 || { cat $W/mutex_build.out; echo "FAIL build"; exit 1; }
cp $HOME/mt-work/mutex/src/build/protal $HOME/mt-work/bin/protal-mutex
B=$HOME/bench071
declare -A R1=( [bgzf]=$B/samples_deep/points/rl150_p5000000/sim/reads/rl150_p5000000_s_1_R1.fq.gz
                [gzip]=$HOME/mt-audit/plain/R1.std.gz [plain]=$HOME/mt-audit/plain/R1.fq )
declare -A R2=( [bgzf]=$B/samples_deep/points/rl150_p5000000/sim/reads/rl150_p5000000_s_1_R2.fq.gz
                [gzip]=$HOME/mt-audit/plain/R2.std.gz [plain]=$HOME/mt-audit/plain/R2.fq )
declare -A DB=( [mini]=$HOME/mt-work/rbfull/mini_db/protal_db [full]=$B/V071/protal_db )
declare -A BIN=( [omp]=$HOME/mt-work/bin/protal-inexp [mutex]=$HOME/mt-work/bin/protal-mutex )
[ -f $W/mutex.tsv ] || printf "db\tinput\tthreads\tlock\trep\taligning_s\tpairs_per_s\twall_s\tuser_s\tsys_s\tload\treader_s_per_thread\n" > $W/mutex.tsv
seconds() { awk '{ s = 0; for (i = 1; i <= NF; i++) { v = $i; if (v ~ /ms$/) { sub(/ms/, "", v); s += v / 1000 } else if (v ~ /m$/) { sub(/m/, "", v); s += v * 60 } else if (v ~ /s$/) { sub(/s/, "", v); s += v } } printf "%.3f", s }'; }
quiet() { local i; for i in $(seq 6); do awk -v l=$(cut -d" " -f1 /proc/loadavg) "BEGIN { exit !(l < 6) }" && return; sleep 10; done; }
run() {  # run DB INPUT THREADS LOCK REP
  quiet
  local db=$1 in=$2 t=$3 lock=$4 rep=$5 d=$W/run
  rm -rf $d; mkdir -p $d
  local load=$(cut -d' ' -f1 /proc/loadavg)
  PROTAL_EXPERIMENT_BATCH=32 nice -n 5 /usr/bin/time -f "%e %U %S" -o $d.time ${BIN[$lock]} --db ${DB[$db]} \
      -1 ${R1[$in]} -2 ${R2[$in]} --prefix s -o $d -t $t --no_profile > $d.log 2>&1
  local a=$(grep 'Aligning reads took' $d.log | sed 's/.*took //' | seconds)
  local reader=$(awk -F'\t' '$1 == "Sequence reader" { print $4 }' $d/misc/s_runtime.tsv)
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $db $in $t $lock $rep $a $(awk -v a=$a 'BEGIN { printf "%.0f", 4995812 / a }') \
      $(tail -1 $d.time) $load $reader | tee -a $W/mutex.tsv
}
for rep in 1 2 3; do
  cat ${R1[@]} ${R2[@]} > /dev/null
  for t in 6 4; do for in in bgzf gzip plain; do for l in omp mutex; do run mini $in $t $l $rep; done; done; done
  for in in bgzf gzip plain; do for l in omp mutex; do run full $in 6 $l $rep; done; done
done
echo MUTEX DONE
