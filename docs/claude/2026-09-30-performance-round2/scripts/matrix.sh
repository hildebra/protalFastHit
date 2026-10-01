#!/bin/bash
# Alternating A/B wall-clock runs of the alignment stage (--no_profile --verbose):
#   matrix.sh REPS
# OLD_BIN/OLD_DB: the binary and database of the previous round (bin-old, db900);
# NEW_BIN/NEW_DB: the current build and a database built by it (db900n).
# Datasets: mix and w900, gzipped and plain FASTQ, 1 and $T threads. Results are appended to
# $PERF_DIR/runs2/matrix.tsv: label, wall, user, sys, maxrss, load, and the stage timers.
set -u
source "$(dirname "$0")/env.sh"
reps=${1:-3}; T=${T:-$(nproc)}
OLD_BIN=${OLD_BIN:-$PERF_DIR/bin-old/protal_avx2}; OLD_DB=${OLD_DB:-$PERF_DIR/db900}
NEW_BIN=${NEW_BIN:-$PERF_DIR/build-rel/protal_avx2}; NEW_DB=${NEW_DB:-$PERF_DIR/db900n}
out=$PERF_DIR/runs2; mkdir -p $out
tsv=$out/matrix.tsv
gz_dir() { if [ -d $PERF_DIR/reads/$1/reads ]; then echo $PERF_DIR/reads/$1/reads; else echo $PERF_DIR/reads/$1; fi; }
for ds in mix w900; do
  plain=$PERF_DIR/reads/${ds}_plain; mkdir -p $plain
  for m in R1 R2; do [ -s $plain/${ds}_$m.fq ] || zcat $(gz_dir $ds)/*_$m.fq.gz > $plain/${ds}_$m.fq; done
done
stage() { grep -E "^$1 took" $2 | head -1 | sed -E 's/.* took ([^ ]+( [^ ]+)?)( mean.*)?$/\1/' ; }
one() { # label bin db reads threads
  local label=$1 bin=$2 db=$3 reads=$4 t=$5
  local d=$out/$label; rm -rf $d; mkdir -p $d
  local r1=$(ls $reads/*_R1.fq* | head -1) r2=$(ls $reads/*_R2.fq* | head -1)
  local l0=$(cut -d' ' -f1-3 /proc/loadavg)
  /usr/bin/time -v -o $d/time.txt $bin --db $db -1 $r1 -2 $r2 -o $d/out -t $t --no_qcmsa --no_profile --verbose > $d/stdout.log 2> $d/stderr.log
  f() { grep "$1" $d/time.txt | awk '{print $NF}'; }
  printf "%s\t%s\t%s\t%s\t%s\tload=%s\tloop=%s\treader=%s\tseed=%s\talign=%s\tloadidx=%s\n" "$label" "$(f 'Elapsed (wall')" "$(f 'User time')" "$(f 'System time')" "$(f 'Maximum resident')" "$l0" \
    "$(stage 'OMP Loop handler' $d/stdout.log)" "$(stage 'Sequence reader' $d/stdout.log)" "$(stage 'Seed- and Anchor-finding' $d/stdout.log)" "$(stage 'Alignment handler' $d/stdout.log)" "$(stage 'Load Index' $d/stdout.log)" | tee -a $tsv
  rm -rf $d/out*
}
for rep in $(seq 1 $reps); do
  for ds in mix w900; do
    for cfg in "gz_t1:$(gz_dir $ds):1" "gz_t$T:$(gz_dir $ds):$T" "plain_t$T:$PERF_DIR/reads/${ds}_plain:$T"; do
      IFS=: read name reads t <<< "$cfg"
      [ "$t" = 1 ] && [ $rep -gt 2 ] && continue
      one old_${ds}_${name}_$rep $OLD_BIN $OLD_DB $reads $t
      one new_${ds}_${name}_$rep $NEW_BIN $NEW_DB $reads $t
    done
  done
done
