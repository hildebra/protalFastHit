#!/usr/bin/env bash
# Parallel BGZF: full unit suite, then protal on the 500k-pair BGZF sample with 1 and 3 inflating threads per file
# (and the 5M-pair member-gzip pair's first 500k? no: the BGZF one only), outputs compared. 4 cores, niced.
set -uo pipefail
SC=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/cdfd1fc3-105f-41f2-9438-012dde2a0acc/scratchpad/gz
G=$HOME/gzpar; K=$G/work/build; T="taskset -c 0-3 nice -n 5"
DB=$HOME/bench071/V075/protal_db; P=$HOME/bench071/samples/points
bash $SC/build.sh work || exit 1
echo "== unit tests ($(date +%T))"
$T $K/tests/protal_tests --gtest_brief=1 > $G/unit_all.log 2>&1; echo "exit $?"
grep -E "^\[  (PASSED|FAILED|SKIPPED)|FAILED  \]" $G/unit_all.log | head
echo "== whole runs ($(date +%T))"
R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=${R1%_R1.fq.gz}_R2.fq.gz
head -c 18 $R1 | od -An -tx1
mkdir -p $G/cmp
for n in 1 3; do
  o=$G/cmp/pe.t$n; rm -rf $o
  PROTAL_INFLATE_THREADS=$n $T $K/protal --db $DB -1 $R1 -2 $R2 --read_type pe --prefix s -o $o -t 4 --no_qcmsa --verbose > $o.log 2>&1
  echo "  threads $n: exit $?; $(grep -E '^Input ' $o.log | tr '\n' ' ')"
done
a=$G/cmp/pe.t1; b=$G/cmp/pe.t3
sa=$(zstd -dc $a/s.sam.zst | grep -v '^@PG\|^@CO' | LC_ALL=C sort | md5sum | cut -c1-12)
sb=$(zstd -dc $b/s.sam.zst | grep -v '^@PG\|^@CO' | LC_ALL=C sort | md5sum | cut -c1-12)
diffs=$(cd $a && find . -type f ! -name '*.sam.zst' ! -name '*.err' ! -name '*_runtime.tsv' | while read -r f; do cmp -s "$f" "$b/$f" || echo "$f"; done | tr '\n' ' ')
echo "  SAM records $([ "$sa" = "$sb" ] && echo identical || echo "DIFFER $sa $sb") ($(zstd -dc $a/s.sam.zst | grep -vc '^@')); other files differing: ${diffs:-none}"
echo "== done ($(date +%T))"
