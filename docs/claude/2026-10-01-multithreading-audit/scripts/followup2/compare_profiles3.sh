#!/usr/bin/env bash
# Profiles existing SAMs (--profile_only, read only) with the 0.7.1 binary (base) and the parallel-profiling
# build (par), each sample alone and each database's samples together, at several thread counts, and compares
# every output file byte for byte (diff -r) and the logs without their timings. Niced.
set -uo pipefail
W=$HOME/mt-work/compare3
BASE=$HOME/mt-work/bin/protal-hb
PAR=$HOME/mt-work/bin/protal-hn
B=$HOME/bench071/runs
rm -rf $W; mkdir -p $W
declare -A FULL=( [pe5M]=$HOME/mt-audit/runs/b6/s1.sam.zst
                  [pe500k150]=$B/v071.full.pe.rl150_p500000_s_1/rl150_p500000_s_1.sam.zst
                  [pe500k100]=$B/v071.full.pe.rl100_p500000_s_1/rl100_p500000_s_1.sam.zst
                  [pe10k]=$B/v071.full.pe.rl150_p10000_s_1/rl150_p10000_s_1.sam.zst
                  [pe1k]=$B/v071.full.pe.rl150_p1000_s_1/rl150_p1000_s_1.sam.zst
                  [se500k]=$B/v071.full.se.rl150_p500000_s_1/rl150_p500000_s_1.sam.zst
                  [pb90M]=$B/v071.full.pb.pb_b90000000_s_1/pb_b90000000_s_1.sam.zst
                  [ont90M]=$B/v071.full.ont.ont_b90000000_s_1/ont_b90000000_s_1.sam.zst )
declare -A MISSING=( [pe500k150]=$B/v071.missing.pe.rl150_p500000_s_1/rl150_p500000_s_1.sam.zst
                     [se500k]=$B/v071.missing.se.rl150_p500000_s_1/rl150_p500000_s_1.sam.zst
                     [ont90M]=$B/v071.missing.ont.ont_b90000000_s_1/ont_b90000000_s_1.sam.zst )
DBFULL=$HOME/bench071/V071/protal_db
DBMISSING=$HOME/bench071/V071/training_db

run() {  # run NAME BIN DB THREADS SAMS PREFIXES
  local name=$1 bin=$2 db=$3 t=$4 sams=$5 prefixes=$6
  mkdir -p $W/$name
  nice -n 5 $bin --db $db --profile_only $sams --prefix $prefixes -o $W/$name/out -t $t --no_qcmsa > $W/$name/log 2>&1
  echo $? > $W/$name/exit
  grep -vE 'took|Took|[0-9]+(\.[0-9]+)?(ms|s)\b|Start parallel|threads' $W/$name/log > $W/$name/log.cmp
}
compare() {  # compare A B
  if diff -r -q $W/$1/out $W/$2/out > $W/diff_$1_$2 2>&1 && cmp -s $W/$1/exit $W/$2/exit; then
    echo "SAME $1 $2 ($(find $W/$1/out -type f | wc -l) files$(cmp -s $W/$1/log.cmp $W/$2/log.cmp && echo ', same log' || echo ', LOGS DIFFER'))"
  else
    echo "DIFFER $1 $2: $(head -3 $W/diff_$1_$2 | tr '\n' ' ')"
  fi
}
for set in FULL MISSING; do
  declare -n S=$set
  db=$([ $set = FULL ] && echo $DBFULL || echo $DBMISSING)
  all_sams=""; all_prefixes=""
  for sample in "${!S[@]}"; do
    sam=${S[$sample]}
    all_sams+="${all_sams:+,}$sam"; all_prefixes+="${all_prefixes:+,}$sample"
    run $set.$sample.base $BASE $db 6 $sam $sample
    for t in 1 6 8; do run $set.$sample.par$t $PAR $db $t $sam $sample; compare $set.$sample.base $set.$sample.par$t; done
  done
  run $set.all.base $BASE $db 6 $all_sams $all_prefixes
  for t in 6 16; do run $set.all.par$t $PAR $db $t $all_sams $all_prefixes; compare $set.all.base $set.all.par$t; done
done
echo done
