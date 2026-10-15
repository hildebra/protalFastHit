#!/bin/bash
# The sparse-first gate's loss of short and divergent reads (F1 audit, B2): reads from 40 species' own marker genes
# mutated at 0, 3, 5 and 10% (b2_reads.py), single-end 75/100/150 bp and paired-end 2x75/2x150, each aligned ungated
# (--sparse_first 0) and gated as GTDB builds gate pe and se (--sparse_first 4 --sparse_first_pick stride); per sample
# the fragments aligned and those whose primary record is on their own species (b2_count.py), into OUT/b2.tsv.
# usage: b2_gate_drop.sh <protal> <database.protal> <out folder> [divergences] [setups] [gates: 0 4, or 4m6 for
#        --sparse_first 4 --sparse_first_min_picks 6]     (WSL; 4 cores)
set -u
PROTAL=$1; DB=$2; OUT=$3
DIVERGENCES=${4:-0 0.03 0.05 0.10}; SETUPS=${5:-se75 se100 se150 pe75 pe150}; GATES=${6:-0 4}
HERE=$(cd "$(dirname "$0")" && pwd)
PY="python3 -I"
mkdir -p $OUT/reads $OUT/runs
if [ ! -f $OUT/ref/reference.fna ]; then
  $PROTAL --unpack_db --db $DB --unpack_dir $OUT/ref -t 4 > $OUT/unpack.log 2>&1 || { echo "unpack failed"; exit 1; }
fi
printf 'setup\tdivergence\tgate\tfragments\taligned\town_species\n' > $OUT/b2.tsv
for d in $DIVERGENCES; do
  for setup in $SETUPS; do
    len=${setup:2}; m=1; [ ${setup:0:2} = pe ] && m=2
    s=$OUT/reads/${setup}_d$d
    $PY $HERE/b2_reads.py $OUT/ref/reference.fna $s --divergence $d --length $len --mates $m 2> /dev/null
    frags=$(( $(wc -l < ${s}_1.fq) / 4 ))
    for gate in $GATES; do
      name=${setup}_d${d}_c$gate
      if [ $m = 2 ]; then reads="-1 ${s}_1.fq -2 ${s}_2.fq --read_type pe"; else reads="-1 ${s}_1.fq --read_type se"; fi
      # 4m6: C 4 with --sparse_first_min_picks 6 (a short read's stride halved until it has 6 picks).
      picks=""; c=$gate
      case $gate in *m*) c=${gate%m*}; picks="--sparse_first_min_picks ${gate#*m}" ;; esac
      taskset -c 0-3 nice -n 5 $PROTAL --db $DB $reads --prefix $name -o $OUT/runs/$name -t 4 --no_profile --no_strains \
          --force --sparse_first $c --sparse_first_pick stride $picks > $OUT/runs/$name.log 2>&1
      sam=$(ls $OUT/runs/$name/$name.sam* 2>/dev/null | head -1)
      if [ -z "$sam" ]; then echo "no SAM for $name (see $OUT/runs/$name.log)"; continue; fi
      printf '%s\t%s\tC%s\t%s\n' $setup $d $gate "$($PY $HERE/b2_count.py $sam $frags)" >> $OUT/b2.tsv
    done
  done
done
column -t $OUT/b2.tsv
