#!/usr/bin/env bash
# Profile the scenarios with protal 0.6.0a and 0.7, timed (/usr/bin/time -v in runs/<run>.time).
#   v06.<sample>     0.6.0a, its own build of the database, the shipped model (random_forest.xml)
#   v07.<sample>     0.7, its build, the models trained on this world (build_gtdb_database.py, ~/tune/V2)
#   v07m06.<sample>  0.7 with 0.6's model: what the binary changes, apart from the model
#   v07.se|pb|ont.*  single-end, PacBio and Nanopore reads of the same communities (0.6 has none)
#   strainA.v06|v07  run A of the strain audit (42 samples, one strain per species, 1-50x), with qcmsa
# A run that completed is not run again (runs/<run>.done).
set -uo pipefail
B=${BENCH:-$HOME/bench07}
SRC6=${SRC6:-$HOME/protal-0.6.0a/src} SRC7=${SRC7:-$HOME/fix-build/src}
P6=${P6:-$SRC6/build/protal} P7=${P7:-$HOME/fix-build/bin/protal}
TEST=${TEST:-$HOME/tune/V2/test/points}
ACC=${ACC:-$HOME/audit5/accuracy}
T=${T:-6}
R=$B/runs
mkdir -p $R

run() {
  local name=$1; shift
  [ -f $R/$name.done ] && return
  rm -rf $R/$name
  mkdir -p $R/$name
  echo "$(date +%T) $name"
  if /usr/bin/time -v -o $R/$name.time "$@" -o $R/$name -t $T > $R/$name.log 2>&1; then
    touch $R/$name.done
  else
    echo "  $name failed ($?), see $R/$name.log"
  fi
}

for point in rl150_p500000 rl150_p10000 rl150_p500 rl100_p500000 rl250_p500000; do
  for s in 1 2; do
    id=${point}_s_$s
    r=$TEST/$point/sim/reads/$id
    pair=(-1 ${r}_R1.fq.gz -2 ${r}_R2.fq.gz --prefix $id --no_qcmsa)
    run v06.$id $P6 --db $B/tune/db06 "${pair[@]}"
    run v07.$id $P7 --db $B/tune/db07 "${pair[@]}"
    run v07m06.$id $P7 --db $B/tune/db07 --model $SRC7/scripts/random_forest.xml "${pair[@]}"
  done
done

for s in 1 2; do
  for point in rl150_p500000 rl150_p10000; do
    id=${point}_s_$s
    run v07.se.$id $P7 --db $B/tune/db07 -1 $TEST/$point/sim/reads/${id}_R1.fq.gz --prefix se_$id --no_qcmsa
  done
  for point in pb_b90000000 ont_b90000000 ont_b3000000; do
    id=${point}_s_$s
    type=${point%%_*}
    run v07.$type.$id $P7 --db $B/tune/db07 -1 $TEST/$point/sim/reads/$id.fq.gz --read_type $type --prefix $id --no_qcmsa
  done
done

run strainA.v06 $P6 --db $B/strainw/db06 --map $ACC/sim_A/protal.meta --qcmsa_script $SRC6/scripts/qcmsa.py
run strainA.v07 $P7 --db $B/strainw/db07 --map $ACC/sim_A/protal.meta --qcmsa_script $SRC7/scripts/qcmsa.py
echo "$(date +%T) done"
