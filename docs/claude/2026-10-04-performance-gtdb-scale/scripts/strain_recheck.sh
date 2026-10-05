#!/usr/bin/env bash
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/37e6327e-2743-4622-91c7-90dfd0315ee2/scratchpad
bash $S/build_work.sh 'StrainOutput*:Haplotypes*:FlexScan*' 6 || exit 1
W=$HOME/perf-gtdb; N=$W/work/build/protal; DB=$W/e2e/db/database.protal; X=$W/avx
for t in 6 1; do
  o=$X/cohort.t$t; rm -rf $o
  /usr/bin/time -f "%e s wall" $N --db $DB --map $X/cohort.t$t.map -t $t > $o.log 2> $o.time
  echo "threads $t: $(tail -1 $o.time)"
done
diff -rq -x '*_runtime.tsv' -x '*.statistics.tsv' $X/cohort.t6/strains $X/cohort.t1/strains > /dev/null && echo "strain outputs identical at 6 and 1 threads" || echo "STRAIN OUTPUTS DIFFER"
diff <(grep -E '^s__' $X/cohort.t6.log) <(grep -E '^s__' $X/cohort.t1.log) > /dev/null && echo "species lines in the same order ($(grep -cE '^s__' $X/cohort.t6.log) of $(grep -c . <(cut -f1 $X/cohort.t6/strains/*.tsv | grep '^s__' | sort -u)) species at line starts)" || { echo "ORDER DIFFERS"; diff <(grep -E '^s__' $X/cohort.t6.log) <(grep -E '^s__' $X/cohort.t1.log) | head -6; }
# the script's timer parser on both logs
took() {
    grep -m1 -oE "$1 took .*" "$2" | sed -e "s/^$1 took //" -e 's/ (.*//' -e 's/ mean over.*//' | awk '
        { s = 0; for (i = 1; i <= NF; i++) { v = $i; n = v + 0
              if (v ~ /^[0-9]+ms$/) s += n / 1000; else if (v ~ /^[0-9]+h$/) s += n * 3600
              else if (v ~ /^[0-9]+m$/) s += n * 60; else if (v ~ /^[0-9]+s$/) s += n }
          printf "%.3f", s; found = 1 }
        END { if (!found) printf "NA" }'
}
for t in 6 1; do echo "t$t: strains $(took 'Strain-level MSAs' $X/cohort.t$t.log) build $(took 'Building the strain MSAs' $X/cohort.t$t.log) qcmsa $(took 'qcMSA' $X/cohort.t$t.log) profiling $(took 'Profiling' $X/cohort.t$t.log)"; done
grep -c '^Strain-level MSAs took' $X/cohort.t1.log $X/cohort.t6.log
echo "== log sizes and qcmsa lines"; wc -c $X/cohort.t6.log $X/cohort.t1.log; grep -c '\[qcmsa\]' $X/cohort.t6.log; grep -E '^qcMSA:|\[qcmsa\] WARNING' $X/cohort.t6.log | head -3 | cut -c1-220
