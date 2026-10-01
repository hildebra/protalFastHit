#!/usr/bin/env bash
# Thread scaling of a 5M-pair sample, instrumented build (lock and hand-off times on stderr, MTAUDIT
# lines) and the plain 0.7.1 build, niced. Outputs in ~/mt-audit/runs/<name>/, logs next to them.
set -uo pipefail
M=${MT_AUDIT:-$HOME/mt-audit}
B=$HOME/bench071
DB=$B/V071/protal_db
R1=$B/samples_deep/points/rl150_p5000000/sim/reads/rl150_p5000000_s_1_R1.fq.gz
R2=${R1%_R1.fq.gz}_R2.fq.gz
mkdir -p $M/runs
run() {
  local name=$1 bin=$2 t=$3
  rm -rf $M/runs/$name; mkdir -p $M/runs/$name
  cat /proc/loadavg > $M/runs/$name.load
  nice -n 5 /usr/bin/time -v -o $M/runs/$name.time $M/bin/$bin --db $DB -1 $R1 -2 $R2 --prefix s1 --no_qcmsa \
      -o $M/runs/$name -t $t > $M/runs/$name.log 2> $M/runs/$name.err
  echo "$(date +%T) $name $(grep Elapsed $M/runs/$name.time | awk '{print $NF}') load $(cut -d' ' -f1-3 $M/runs/$name.load)"
}
for spec in ${SPECS:-"i1 protal-instr 1" "i2 protal-instr 2" "i3 protal-instr 3" "i4 protal-instr 4" "i6 protal-instr 6" "b6 protal-base 6" "i6b protal-instr 6" "b6b protal-base 6"}; do
  set -- $spec
  run $1 $2 $3
done
echo done
