#!/usr/bin/env bash
# A SAM of the r226 mix: the local pe SAM's records, and after each read's group N unmapped records (flag 4, ZF of 1-6 taxids).
set -euo pipefail
W=$HOME/perf-gtdb; X=$W/samread; mkdir -p $X
N=${1:-20}
zstd -dc $W/avx/pe.avx2.1/pe.sam.zst > $X/mapped.sam
awk -F'\t' -v N=$N 'BEGIN { srand(7); OFS = "\t" }
  /^@/ { print; next }
  { split($3, a, "_"); if (a[1] ~ /^[0-9]+$/ && !(a[1] in seen)) { seen[a[1]] = 1; tax[nt++] = a[1] } }
  $1 != last && NR > 1 && nt > 0 { for (i = 0; i < N; i++) { k = 1 + int(rand() * 6); zf = tax[int(rand() * nt)]
        for (j = 1; j < k; j++) zf = zf "," tax[int(rand() * nt)]
        printf "LH00123:45:22ABCDLT3:1:%d:%d:%d\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tZU:i:0\tZT:i:0\tZF:Z:%s\n", 1101 + int(u / 9000000), u % 30000, int(u / 30000) % 300, zf; u++ } }
  { last = $1; print }' $X/mapped.sam > $X/mixed.sam
ls -la $X/*.sam
