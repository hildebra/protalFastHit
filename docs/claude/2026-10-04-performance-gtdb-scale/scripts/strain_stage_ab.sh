# The local cohort (pe and pb from their SAMs, qcMSA on) at 6 threads with 186c8ed and the current build, alternated twice.
W=$HOME/perf-gtdb; DB=$W/e2e/db/database.protal; X=$W/avx
for r in 1 2; do for v in before now; do
  bin=$W/work/build/protal; [ $v = before ] && { cp -u $HOME/tiebench/bin/protal-before $W/work/build/protal.before; bin=$W/work/build/protal.before; }
  o=$X/cohort.$v.$r; rm -rf $o; sed "s|^#OUTPUT_DIR\t.*|#OUTPUT_DIR\t$o|" $X/cohort.t6.map > $o.map
  $bin --db $DB --map $o.map -t 6 > $o.log 2>&1
  echo "$v $r: $(grep -hE '^Strain-level MSAs took|^Building the strain MSAs took|^qcMSA took' $o.log | tr '\n' ' ')"
done; done
diff -rq -x '*_runtime.tsv' -x '*.statistics.tsv' $X/cohort.before.1/strains $X/cohort.now.1/strains > /dev/null && echo "strain outputs the same" || echo "strain outputs differ"
