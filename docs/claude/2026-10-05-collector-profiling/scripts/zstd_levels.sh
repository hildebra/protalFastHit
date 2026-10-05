#!/usr/bin/env bash
# Size and single-thread speed of each read type's FASTQ as written (BGZF level 6 or gzip level 1) against gzip and
# zstd levels. Usage: bash zstd_levels.sh DIR (with pe150_R1.fq.gz, ultima.fq.gz, hifi.fq.gz, ont.fq.gz)
set -euo pipefail
cd "$1"
printf "%-10s %-9s %12s %7s %10s %10s\n" file codec bytes ratio comp_MB/s dec_MB/s
for f in pe150_R1 ultima hifi ont; do
  plain=/tmp/zl_$f.fq
  zcat $f.fq.gz > $plain
  n=$(stat -c %s $plain)
  printf "%-10s %-9s %12d %7.3f\n" $f written $(stat -c %s $f.fq.gz) $(echo "$(stat -c %s $f.fq.gz) / $n" | bc -l)
  for codec in gzip:1 gzip:6 zstd:1 zstd:3 zstd:6 zstd:9 zstd:12 zstd:19; do
    tool=${codec%:*}; level=${codec#*:}
    t0=$(date +%s.%N)
    if [ $tool = zstd ]; then zstd -$level -q -T1 -c $plain > /tmp/zl_out; else gzip -$level -c $plain > /tmp/zl_out; fi
    t1=$(date +%s.%N)
    $tool -d -c /tmp/zl_out > /dev/null
    t2=$(date +%s.%N)
    size=$(stat -c %s /tmp/zl_out)
    printf "%-10s %-9s %12d %7.3f %10.1f %10.1f\n" $f "$tool-$level" $size $(echo "$size / $n" | bc -l) \
      $(echo "$n / 1000000 / ($t1 - $t0)" | bc -l) $(echo "$n / 1000000 / ($t2 - $t1)" | bc -l)
  done
  rm -f $plain /tmp/zl_out
done
