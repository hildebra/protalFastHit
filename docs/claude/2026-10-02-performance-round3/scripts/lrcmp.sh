#!/usr/bin/env bash
# lrcmp.sh REF.sam.zst EXP.sam.zst: long-read SAM records of two builds compared by (read, gene, strand):
# records only in one, and for the common ones the same CIGAR, the same MAPQ, and the WFA penalty from the
# CIGAR (4 per X, 6 + 2 per base of each I or D run; S and H free) of the experiment against the reference.
pen() { zstdcat "$1" | grep -v '^@' | awk -F'\t' '{
  c = $6; p = 0; n = ""
  for (i = 1; i <= length(c); i++) { ch = substr(c, i, 1)
    if (ch ~ /[0-9]/) { n = n ch; continue }
    if (ch == "X") p += 4 * n; else if (ch == "I" || ch == "D") p += 6 + 2 * n
    n = "" }
  printf "%s|%s|%d\t%s\t%s\t%d\t%s\n", $1, $3, and($2, 16), c, $5, p, $4 }'; }
pen "$1" | sort -k1,1 > /tmp/lr_ref.tsv; pen "$2" | sort -k1,1 > /tmp/lr_exp.tsv
join -t$'\t' -j1 /tmp/lr_ref.tsv /tmp/lr_exp.tsv > /tmp/lr_both.tsv
echo "records: ref $(wc -l < /tmp/lr_ref.tsv), exp $(wc -l < /tmp/lr_exp.tsv), both $(wc -l < /tmp/lr_both.tsv), only ref $(join -t$'\t' -v1 /tmp/lr_ref.tsv /tmp/lr_exp.tsv | wc -l), only exp $(join -t$'\t' -v2 /tmp/lr_ref.tsv /tmp/lr_exp.tsv | wc -l)"
awk -F'\t' '{ n++; same += $2 == $6; mq += $3 == $7; pos += $5 == $9; d = $8 - $4; if (d < 0) better++; else if (d > 0) worse++; sum += d; if (d > 0) wsum += d; ref += $4 }
  END { printf "common: same CIGAR %.1f%%, same MAPQ %.1f%%, same POS %.1f%%; penalty exp-ref: lower %d, higher %d (mean +%.1f where higher), total %+d of %d (%+.2f%%)\n",
        100 * same / n, 100 * mq / n, 100 * pos / n, better, worse, worse ? wsum / worse : 0, sum, ref, 100 * sum / ref }' /tmp/lr_both.tsv
