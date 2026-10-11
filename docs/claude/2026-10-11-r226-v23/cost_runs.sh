#!/bin/bash
# What each feature group costs the profiler: the cost-probe build (~/b23p: 8d8c89b + cost_probe.pl) profiles existing
# SAMs of the small real world (~/realworld_sp_s1: its training database and samples, 170 species) under callgrind on
# one thread (core 3), with each group's evidence skipped in turn (PROTAL_SKIP; adjacency by --no_gene_neighbours).
# Each line of costs.tsv: sample, variant, instructions of the whole run and of ProfileSample (the profiling proper,
# without the start-up's model and table loading).
set -u
P=~/b23p/build/protal
DB=~/realworld_sp_s1/work/training_db
T=~/realworld_sp_s1/work/training/points
W=~/fcost/runs
mkdir -p $W
OUT=$W/costs.tsv
[ -f $OUT ] || echo -e "sample\tvariant\tskip\tother\ttotal_Ir\tprofile_Ir" > $OUT
declare -A SAMS=(
  [pe1M]=$T/rl150_p1000000/protal/alignments/rl150_p1000000_s_10.sam.zst
  [pb100M]=$T/pb_b100000000/protal/alignments/pb_b100000000_s_1.sam.gz
  [pe1Mx5]=$T/rl150_p1000000/protal/alignments/rl150_p1000000_s_1.sam.zst,$T/rl150_p1000000/protal/alignments/rl150_p1000000_s_2.sam.zst,$T/rl150_p1000000/protal/alignments/rl150_p1000000_s_3.sam.zst,$T/rl150_p1000000/protal/alignments/rl150_p1000000_s_4.sam.zst,$T/rl150_p1000000/protal/alignments/rl150_p1000000_s_5.sam.zst
)
run() {  # sample variant skip [extra options]
  local s=$1 v=$2 k=$3; shift 3
  grep -q "^$s	$v	" $OUT && return
  local d=$W/$s.$v; rm -rf $d; mkdir -p $d
  ( cd $d && PROTAL_SKIP=$k taskset -c 3 valgrind --tool=callgrind --callgrind-out-file=cg.out $P --profile_only ${SAMS[$s]} \
      --db $DB -o $d/out -t 1 --no_strains "$@" > run.log 2> vg.log )
  local tot=$(grep -oP "Collected : \K[0-9]+" $d/vg.log)
  local prof=$(callgrind_annotate --inclusive=yes $d/cg.out 2>/dev/null | grep -m1 "protal::ProfileSample(" | awk '{gsub(",", "", $1); print $1}')
  echo -e "$s\t$v\t$k\t$*\t$tot\t$prof" >> $OUT
  rm -f $d/cg.out
  echo "$s $v done"
}
for s in pe1M pb100M; do
  run $s all_on 0
  run $s no_ancestry 1
  run $s no_weights 2
  run $s no_alleles 4
  run $s no_gaps 8
  run $s no_distance 16
  run $s no_untried 32
  run $s no_consistency 64
  run $s no_adjacency 0 --no_gene_neighbours
  run $s no_anc_weights 3
  run $s none_of_these 127
done
run pe1Mx5 all_on 0
run pe1Mx5 no_ancestry 1
run pe1Mx5 none_of_these 127
echo COST RUNS DONE
