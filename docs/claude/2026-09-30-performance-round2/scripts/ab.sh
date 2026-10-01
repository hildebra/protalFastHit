#!/bin/bash
# Alternating A/B of two binaries on one read set, alignment stage only:
#   [PIN=cpu] ab.sh BIN_A BIN_B DB READS_DIR THREADS REPS
# One line per run: name, wall, user, loop time, seeding time (seconds); then the minimum and median of each
# column. PIN pins every run to one vCPU (taskset), which removes part of this laptop's core-to-core noise.
source "$(dirname "$0")/env.sh"
A=$1; B=$2; db=$3; reads=$4; t=$5; reps=$6
r1=$(ls $reads/*_R1.fq* | head -1); r2=$(ls $reads/*_R2.fq* | head -1)
secs() { sed -E 's/.* took //; s/ms//; s/s / /' | awk '{ if (NF == 2) print $1 + $2 / 1000; else print $1 / 1000 }'; }
run() { # name bin
  local d=$PERF_DIR/runs2/ab_$1; rm -rf $d; mkdir -p $d
  /usr/bin/time -f "%e %U" -o $d/time.txt ${PIN:+taskset -c $PIN} $2 --db $db -1 $r1 -2 $r2 -o $d/out -t $t --no_qcmsa --no_profile --verbose > $d/stdout.log 2> $d/stderr.log
  local loop=$(grep "^OMP Loop handler took" $d/stdout.log | secs)
  local seed=$(grep -P "^\t\tSeeding took" $d/stdout.log | secs)
  echo "$1 $(cat $d/time.txt) $loop $seed"; rm -rf $d/out*
}
out=$(mktemp)
for i in $(seq 1 $reps); do echo "A_$i $(run a $A | cut -d' ' -f2-)"; echo "B_$i $(run b $B | cut -d' ' -f2-)"; done | tee $out
echo "column:        wall      user      loop      seeding"
for x in A B; do
  for stat in min median; do
    printf "%s %-8s" $x $stat
    for col in 2 3 4 5; do
      grep "^${x}_" $out | awk -v c=$col '{print $c}' | sort -n | awk -v s=$stat '{v[NR]=$1} END { printf " %9.3f", (s == "min" ? v[1] : v[int((NR+1)/2)]) }'
    done; echo
  done
done
rm -f $out
