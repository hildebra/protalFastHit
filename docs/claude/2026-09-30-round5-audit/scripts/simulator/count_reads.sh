#!/bin/bash
# count_reads.sh SIMDIR: per (sample, genome) read pairs in reads/<sample>_R1.fq.gz vs manifest read_pairs.
# Genome = read name up to "_contig" (ART names reads <contig id>-<n>/1).
set -u
D=$1
for r1 in $D/reads/*_R1.fq.gz; do
  s=$(basename $r1 _R1.fq.gz)
  zcat $r1 | awk 'NR%4==1' | sed 's/^@//; s/_contig.*//; s/-[0-9]*\/1$//' | sort | uniq -c | awk -v s=$s '{print s"\t"$2"\t"$1}'
done > $D/written.tsv
awk -F'\t' 'FNR==NR{w[$1"|"$2]=$3; next}
  FNR==1{print "sample\tgenome\tgenome_len\trequested\twritten\tdiff\trel_diff\tvcov_manifest\tvcov_written"; next}
  {k=$1"|"$2; n=(k in w)?w[k]:0; tr+=$6; tw+=n; printf "%s\t%s\t%s\t%d\t%d\t%d\t%.4f\t%s\t%.4f\n", $1,$2,$5,$6,n,n-$6,(n-$6)/$6,$7,n*300/$5}
  END{printf "TOTAL requested %d written %d (%.4f)\n", tr, tw, (tw-tr)/tr}' $D/written.tsv $D/manifest.tsv
