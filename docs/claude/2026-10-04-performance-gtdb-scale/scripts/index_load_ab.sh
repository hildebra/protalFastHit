#!/usr/bin/env bash
# The index load without the chunk value buffers: tests; e2e pe and pb against the previous build (fdad375, ~/tiebench's
# protal-after): outputs, "Load Index took" and the memory lines, alternated 3 times at 6 threads (and once at 1).
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/37e6327e-2743-4622-91c7-90dfd0315ee2/scratchpad
bash $S/build_work.sh "${TESTS:-*Seedmap*:*Index*:*Codec*:*Pack*}" 8 || exit 1
W=$HOME/perf-gtdb; DB=$W/e2e/db/database.protal; X=$W/idx; rm -rf $X; mkdir -p $X
cp $W/work/build/protal $X/protal.new; OLD=$HOME/tiebench/bin/protal-after
P=$HOME/bench071/samples/points; R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
PB=$HOME/bench071/samples_lr073/points/pb_b90000000/sim/reads/pb_b90000000_s_1.fq.gz
for t in ${THREADS:-6 6 6 1}; do for v in old new; do for s in pe pb; do
  bin=$X/protal.new; [ $v = old ] && bin=$OLD
  if [ $s = pe ]; then reads=(-1 $R1 -2 $R2); else reads=(-1 $PB); fi
  o=$X/$s.$v.t$t.$RANDOM; /usr/bin/time -f "%M" -o $o.rss $bin --db $DB "${reads[@]}" --read_type $s --prefix $s -o $o -t $t --no_qcmsa > $o.log 2>&1
  echo "$s $v t$t: $(grep -h '^Load Index took' $o.log); $(grep -h '^Memory after loading the index' $o.log | sed 's/Memory after loading the index: //'); max RSS $(awk '{printf "%.2f GB", $1/1048576}' $o.rss)"
  [ -e $X/$s.ref ] || ln -s $o $X/$s.ref
  diff -rq -x '*_runtime.tsv' -x '*.statistics.tsv' $X/$s.ref $o > /dev/null || echo "   OUTPUTS DIFFER from the first run: $(diff -rq -x '*_runtime.tsv' -x '*.statistics.tsv' $X/$s.ref $o | head -2)"
done; done; done
grep -h "Index in memory" $X/pe.ref.log
