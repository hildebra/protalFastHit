#!/usr/bin/env bash
# Input-path audit: alignment of the 5M-pair sample (no profiling) by the current code with the batch-size
# experiment patch (protal-inexp; PROTAL_EXPERIMENT_BATCH, default 32), from BGZF (as simulated), standard gzip
# (pigz -6, one member: the zlib-ng path) and uncompressed FASTQ, against the mini database (3 species: almost no
# read aligns, so a pair costs k-mer extraction and lookups only and the input path is the limit) and the full
# database (765 species). Alternated, niced; rows into ~/mt-work/input/runs.tsv.
set -uo pipefail
SP=$(cd "$(dirname "$0")" && pwd)  # this folder
W=$HOME/mt-work/input; mkdir -p $W
bash $SP/build_wt.sh inexp $SP/batch_experiment.patch a212559 protal > $W/build.out 2>&1 || { cat $W/build.out; exit 1; }
cp $HOME/mt-work/inexp/src/build/protal $HOME/mt-work/bin/protal-inexp
B=$HOME/bench071
declare -A R1=( [bgzf]=$B/samples_deep/points/rl150_p5000000/sim/reads/rl150_p5000000_s_1_R1.fq.gz
                [gzip]=$HOME/mt-audit/plain/R1.std.gz [plain]=$HOME/mt-audit/plain/R1.fq )
declare -A R2=( [bgzf]=$B/samples_deep/points/rl150_p5000000/sim/reads/rl150_p5000000_s_1_R2.fq.gz
                [gzip]=$HOME/mt-audit/plain/R2.std.gz [plain]=$HOME/mt-audit/plain/R2.fq )
declare -A DB=( [mini]=$HOME/mt-work/rbfull/mini_db/protal_db [full]=$B/V071/protal_db )
cat ${R1[@]} ${R2[@]} > /dev/null
[ -f $W/runs.tsv ] || printf "db\tinput\tthreads\tbatch\trep\taligning_s\tpairs_per_s\twall_s\tuser_s\tsys_s\tload\n" > $W/runs.tsv
seconds() { awk '{ s = 0; for (i = 1; i <= NF; i++) { v = $i; if (v ~ /ms$/) { sub(/ms/, "", v); s += v / 1000 } else if (v ~ /m$/) { sub(/m/, "", v); s += v * 60 } else if (v ~ /s$/) { sub(/s/, "", v); s += v } } printf "%.3f", s }'; }
# Before each run, wait up to a minute for the load average to drop below 6 (other sessions share the machine).
quiet() { local i; for i in $(seq 6); do awk -v l=$(cut -d" " -f1 /proc/loadavg) "BEGIN { exit !(l < 6) }" && return; sleep 10; done; }
run() {  # run DB INPUT THREADS BATCH REP
  quiet
  local db=$1 in=$2 t=$3 batch=$4 rep=$5 d=$W/run
  rm -rf $d; mkdir -p $d
  local load=$(cut -d' ' -f1 /proc/loadavg)
  PROTAL_EXPERIMENT_BATCH=$batch nice -n 5 /usr/bin/time -f "%e %U %S" -o $d.time $HOME/mt-work/bin/protal-inexp --db ${DB[$db]} \
      -1 ${R1[$in]} -2 ${R2[$in]} --prefix s -o $d -t $t --no_profile > $d.log 2>&1
  local a=$(grep 'Aligning reads took' $d.log | sed 's/.*took //' | seconds)
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $db $in $t $batch $rep $a $(awk -v a=$a 'BEGIN { printf "%.0f", 4995812 / a }') \
      $(tail -1 $d.time) $load | tee -a $W/runs.tsv
}
for rep in 1 2; do
  for t in 1 2 4 6; do for in in bgzf gzip plain; do run mini $in $t 32 $rep; done; done
  for batch in 128 512; do for in in bgzf gzip plain; do run mini $in 6 $batch $rep; done; done
  for in in bgzf gzip plain; do run full $in 6 32 $rep; done
  for in in gzip plain; do run full $in 6 128 $rep; done
done
# The input path alone, on one reader thread (readbench of the audit): CPU seconds, which the load disturbs less
# than the wall time; compressed minus plain is the inflating threads' CPU.
[ -f $W/readbench.tsv ] || printf "input\trep\twall_s\tuser_s\tsys_s\tload\tline\n" > $W/readbench.tsv
for rep in 1 2 3; do
  for in in bgzf gzip plain; do
    quiet
    load=$(cut -d' ' -f1 /proc/loadavg)
    line=$(nice -n 5 /usr/bin/time -f "%e %U %S" -o $W/rb.time $HOME/mt-audit/readbench/readbench ${R1[$in]} ${R2[$in]} 1 2> /dev/null)
    printf "%s\t%s\t%s\t%s\t%s\n" $in $rep "$(tail -1 $W/rb.time | tr ' ' '\t')" $load "$line" | tee -a $W/readbench.tsv
  done
  for in in bgzf gzip plain; do  # the second file twice: the slower file's inflating cost alone
    quiet
    load=$(cut -d' ' -f1 /proc/loadavg)
    line=$(nice -n 5 /usr/bin/time -f "%e %U %S" -o $W/rb.time $HOME/mt-audit/readbench/readbench ${R2[$in]} ${R2[$in]} 1 2> /dev/null)
    printf "%s\t%s\t%s\t%s\t%s\n" ${in}_R2R2 $rep "$(tail -1 $W/rb.time | tr ' ' '\t')" $load "$line" | tee -a $W/readbench.tsv
  done
done
echo INPUT DONE
