#!/bin/bash
# CIGAR composition of a protal SAM (protal writes M for matches and X for mismatches):
#   sam_stats.sh SAM
grep -v '^@' $1 | awk -F'\t' '
{ n++; cig = $6
  if (cig ~ /[ID]/) indel++
  x = 0; s = cig
  while (match(s, /[0-9]+[=XIDMS]/)) { tok = substr(s, RSTART, RLENGTH); op = substr(tok, length(tok))
    if (op == "X") x += substr(tok, 1, length(tok)-1) + 0
    if (op == "S") soft++
    s = substr(s, RSTART + RLENGTH) }
  if (cig !~ /[ID]/) { if (x == 0) m0++; else if (x == 1) m1++; else if (x == 2) m2++; else if (x <= 5) m5++; else mbig++ }
}
END { printf "%d records; without indels: %d exact, %d with 1 mismatch, %d with 2, %d with 3-5, %d with >5; %d with I/D; %d soft clips\n", n, m0, m1, m2, m5, mbig, indel, soft }'
