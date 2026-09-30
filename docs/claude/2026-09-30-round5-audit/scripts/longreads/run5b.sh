#!/bin/bash
# FASTA input: which records differ from the FASTQ run's, and how.
S=$(dirname "$0")
W=~/audit6/longreads
G=~/strain-build/mini_db/gtdb_r226
D=$W/mini_db
OUT=$S/out_e5b.txt
cd $W/e1
{
python3 $S/eval_sam.py outfa/hifi.sam hifi.fa hifi.truth.pkl --gtdb $G --db $D --fasta_qual '?' --show 3 | grep -A3 -E "^(clip_sum|seq_mismatch)" | cut -c1-300
echo "== SEQ/QUAL of the same records, FASTQ vs FASTA"
paste <(grep -v '^@' out/hifi.sam | cut -f1-6,10 | sort) <(grep -v '^@' outfa/hifi.sam | cut -f1-6,10 | sort) | awk -F'\t' '$7 != $14 {print $1, $2, $3, $6, length($7), length($14)}' | cut -c1-200 | head
echo "== reads with problems: lengths in fq and fa"
for r in $(python3 $S/eval_sam.py outfa/hifi.sam hifi.fa hifi.truth.pkl --gtdb $G --db $D --fasta_qual '?' --show 100 | grep -A100 '^clip_sum' | grep "^    \[" | sed "s/^    \['\([^']*\)'.*/\1/" | sort -u | head -5); do
  echo "$r fq=$(awk -v n="@$r" 'NR%4==1{k=($1==n)} NR%4==2&&k{print length($0)}' hifi.fq) fa=$(awk -v n=">$r" '/^>/{k=($1==n); next} k{print length($0)}' hifi.fa)"
done
} > $OUT 2>&1
echo done
