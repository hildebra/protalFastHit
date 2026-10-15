#!/bin/bash
# The reduced r226 database's rl150 test points simulated again (the build removed their reads after profiling them)
# with the build's own simulate_metagenomes (5369c99) and the commands in their run_params.tsv, into
# OUT/<point>/sim; then their manifests and abundance matrices compared with the originals (the same seeds, the same
# genomes and read counts per sample).
# usage: resim.sh <out folder> <point> [...]     (WSL; one thread per point)
set -u
OUT=$1; shift
P=$HOME/c4test/db/work/test/points
mkdir -p $OUT
run_point() {
  local point=$1
  local cmd
  cmd=$(awk -F'\t' '$1=="command"{print $2}' $P/$point/sim/run_params.tsv)
  cmd=${cmd//"-o $P/$point/sim"/"-o $OUT/$point/sim"}
  cmd=${cmd//"--protal_metafile $P/$point/protal"/"--protal_metafile $OUT/$point/protal"}
  mkdir -p $OUT/$point
  echo "$cmd" > $OUT/$point/command.txt
  taskset -c 0-3 nice -n 5 $cmd > $OUT/$point/simulate.log 2>&1
  echo "$point: exit $?"
  for f in manifest.tsv abundance_matrix.tsv; do
    if cmp -s <(sed "s#$OUT/$point#X#g" $OUT/$point/sim/$f) <(sed "s#$P/$point#X#g" $P/$point/sim/$f); then
      echo "$point $f: same as the original"
    else
      echo "$point $f: DIFFERS"
    fi
  done
}
for point in "$@"; do run_point $point & done
wait
