#!/bin/bash
# As ab.sh, for any number of binaries, the order rotated every round (ab.sh always runs A first):
#   [PIN=cpu] abn.sh DB READS_DIR THREADS REPS NAME=BIN...
# One line per run: name, wall, user, loop time, seeding time (seconds); then the minimum and median of each column.
source "$(dirname "$0")/env.sh"
db=$1; reads=$2; t=$3; reps=$4; shift 4
names=(); bins=()
for nb in "$@"; do names+=("${nb%%=*}"); bins+=("${nb#*=}"); done
k=${#names[@]}
r1=$(ls $reads/*_R1.fq* | head -1); r2=$(ls $reads/*_R2.fq* | head -1)
secs() { sed -E 's/.* took //; s/ms//; s/s / /' | awk '{ if (NF == 2) print $1 + $2 / 1000; else print $1 / 1000 }'; }
run() { # name bin
  local d=$PERF_DIR/runs2/abn_$1; rm -rf $d; mkdir -p $d
  /usr/bin/time -f "%e %U" -o $d/time.txt ${PIN:+taskset -c $PIN} $2 --db $db -1 $r1 -2 $r2 -o $d/out -t $t --no_qcmsa --no_profile --verbose > $d/stdout.log 2> $d/stderr.log
  local loop=$(grep "^OMP Loop handler took" $d/stdout.log | secs)
  local seed=$(grep -P "^\t\tSeeding took" $d/stdout.log | secs)
  echo "$1 $(cat $d/time.txt) $loop $seed"; rm -rf $d/out*
}
out=$(mktemp)
for i in $(seq 1 $reps); do
  for j in $(seq 0 $((k - 1))); do
    x=$(( (i + j) % k ))
    run ${names[$x]} ${bins[$x]}
  done
done | tee $out
echo "column:        wall      user      loop      seeding"
for x in "${names[@]}"; do
  for stat in min median; do
    printf "%-6s %-8s" $x $stat
    for col in 2 3 4 5; do
      grep "^$x " $out | awk -v c=$col '{print $c}' | sort -n | awk -v s=$stat '{v[NR]=$1} END { printf " %9.3f", (s == "min" ? v[1] : v[int((NR+1)/2)]) }'
    done; echo
  done
done
rm -f $out
