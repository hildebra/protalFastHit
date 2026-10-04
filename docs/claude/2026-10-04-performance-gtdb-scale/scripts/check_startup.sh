#!/usr/bin/env bash
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/37e6327e-2743-4622-91c7-90dfd0315ee2/scratchpad
bash $S/build_work.sh 'ProfileThreads.*:InputValidation.*:Parsing.*' 6 || exit 1
W=$HOME/perf-gtdb; N=$W/work/build/protal; DB=$HOME/bench071/V073/protal_db
P=$HOME/bench071/samples/points; R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
for v in default taxon_statistics; do
  o=$W/st.$v; rm -rf $o; extra=""; [ $v = taxon_statistics ] && extra="--taxon_statistics"
  $N --db $DB -1 $R1 -2 $R2 --read_type pe --prefix s -o $o -t 6 --no_qcmsa $extra > $o.log 2>&1
  echo "== $v: $(ls $o/misc/*.statistics.tsv 2>/dev/null | wc -l) statistics files; $(grep -cE '^Taxon statistics files took' $o.log) 'Taxon statistics' line"
  grep -E 'took' $o.log | grep -vE 'Profiling sample' | tr '\n' ';' | cut -c1-400; echo
done
diff <(grep -v 'took\|Freeing' $W/st.default.log) <(grep -v 'took\|Freeing' $W/st.taxon_statistics.log) > /dev/null && echo "logs otherwise identical" || echo "LOGS DIFFER"
diff -r -x '*_runtime.tsv' -x '*.statistics.tsv' $W/st.default $W/st.taxon_statistics > /dev/null && echo "outputs identical but for the statistics files" || echo "OUTPUTS DIFFER"
