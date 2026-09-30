#!/bin/bash
# Anchored alignment against an older build, on one read set, 1 thread:
#   anchored_compare.sh OLD_BIN NEW_BIN DB READS_DIR
# 1. NEW --whole_read_alignment must reproduce every output of OLD (runtime.tsv aside).
# 2. NEW's SAM against OLD's, field by field; for records with another CIGAR or position, the WFA
#    penalty (mismatch 4, gap 6 + 2 per base) of each.
# 3. The profile: species and abundances.
# 4. The best and second-best pair bitscores behind MAPQ (--mapq_debug_output).
source "$(dirname "$0")/env.sh"
old=$1; new=$2; db=$3; src=$4
r1=$(ls $src/*_R1.fq* | head -1); r2=$(ls $src/*_R2.fq* | head -1)
out=$PERF_DIR/runs/anchored/$(basename $src); rm -rf $out; mkdir -p $out
run() { local bin=$1 dir=$2; shift 2; $bin --db $db -1 $r1 -2 $r2 -o $out/$dir -t 1 --no_qcmsa "$@" > $out/$dir.log 2> $out/$dir.err; }
run $old old --mapq_debug_output
run $new whole --whole_read_alignment --mapq_debug_output
run $new new --mapq_debug_output --verbose
grep -h 'Anchors aligned from' $out/new.log

d=0
for f in $(cd $out/old && find . -type f ! -name '*_runtime.tsv' | sort); do cmp -s $out/old/$f $out/whole/$f || { d=1; echo "  differs: $f"; }; done
echo "--whole_read_alignment vs OLD: $([ $d = 0 ] && echo 'every output identical' || echo 'OUTPUTS DIFFER')"

sam=$(cd $out/old && ls *.sam | head -1)
prof=${sam%.sam}.profile
# Records matched by (QNAME, mate bits, occurrence): key, FLAG, gene, position, MAPQ, CIGAR.
key() { grep -v '^@' $1 | awk -F'\t' -v OFS='\t' '{ m = and($2, 192); k = $1 "|" m; n[k]++; print k "|" n[k], $2, $3, $4, $5, $6 }' | sort -t$'\t' -k1,1; }
join -t$'\t' -a1 -a2 -e MISSING -o 0,1.2,1.3,1.4,1.5,1.6,2.2,2.3,2.4,2.5,2.6 <(key $out/old/$sam) <(key $out/new/$sam) | awk -F'\t' '
function penalty(c,   s, p, len, op) {
  s = c; p = 0
  while (match(s, /^[0-9]+[MXIDS=]/)) { len = substr(s, 1, RLENGTH - 1) + 0; op = substr(s, RLENGTH, 1)
    if (op == "X") p += 4 * len; else if (op == "I" || op == "D") p += 6 + 2 * len
    s = substr(s, RLENGTH + 1) }
  return p }
{ n++
  if ($2 == "MISSING") { onlynew++; next }
  if ($7 == "MISSING") { onlyold++; next }
  same = 1
  if ($3 != $8) { gene++; same = 0 }
  else {
    if ($4 != $9) { pos++; same = 0 }
    if ($6 != $11) { cig++; same = 0 }
    if ($4 != $9 || $6 != $11) { pa = penalty($6); pb = penalty($11)
      if (pa == pb) tie++; else if (pb < pa) better++; else worse++
      if (++shown <= 5) printf "  %s pos %s -> %s  %s (%d) -> %s (%d)\n", $1, $4, $9, $6, pa, $11, pb }
  }
  if ($5 != $10) { mapq++; same = 0; d = $10 - $5; if (d < 0) d = -d; if (d >= 10) mapq10++ }
  if ($2 != $7) { flag++; same = 0 }
  identical += same }
END { printf "SAM, %d records: identical %d (%.3f%%); other gene %d; other CIGAR %d, other position %d (penalty the same %d, lower %d, higher %d); other MAPQ %d (by >=10: %d), other flag %d; only in OLD %d, only in NEW %d\n",
      n, identical, 100 * identical / n, gene, cig, pos, tie, better, worse, mapq, mapq10, flag, onlyold, onlynew }'

cmp -s <(cut -f1,2 $out/old/$prof) <(cut -f1,2 $out/new/$prof) && echo "profile: same species" || echo "profile: SPECIES DIFFER"
paste <(cut -f3 $out/old/$prof) <(cut -f3 $out/new/$prof) | awk -F'\t' '$1 != $2 { n++; d = $1 - $2; if (d < 0) d = -d; if ($1 > 0 && d / $1 > m) m = d / $1 }
  END { printf "abundances: %d differ, largest by %.2e relative\n", n, m }'

# --mapq_debug_output: correct, MAPQ, best and second-best pair bitscore, read name.
scores() { grep -E '^[01]	' $1 | awk -F'\t' -v OFS='\t' '{ print $5, $2, $3, $4 }' | sort -k1,1; }
join -t$'\t' <(scores $out/old.err) <(scores $out/new.err) | awk -F'\t' '
  { n++; if ($2 != $5) mq++
    if ($3 != $6) { best++; if ($6 > $3) bu++; else bd++ }
    if ($4 != $7) { second++; if ($7 > $4) su++; else sd++ } }
  END { printf "pairs %d: MAPQ differs %d; best bitscore differs %d (higher %d, lower %d); second best %d (higher %d, lower %d)\n", n, mq, best, bu, bd, second, su, sd }'
