#!/usr/bin/env bash
# Each full-database SAM of the 0.7.1 benchmark: its header compressed (zstd -3) against its read files' bytes.
P=$HOME/bench071/samples/points; D=$HOME/bench071/samples_deep/points
printf "sam\theader_zstd3\tread_bytes\tshare\n"
for f in $HOME/bench071/runs/v071.full.*/*.sam.zst; do
  run=$(basename $(dirname $f)); name=$(basename $f .sam.zst); type=$(echo $run | cut -d. -f3)
  point=${name%_s_*}
  base=$P/$point/sim/reads; [ -d $base ] || base=$D/$point/sim/reads
  case $type in
    pe) files="$base/${name}_R1.fq.gz $base/${name}_R2.fq.gz" ;;
    se) files="$base/${name}_R1.fq.gz" ;;
    *)  files="$base/${name}.fq.gz" ;;
  esac
  bytes=0; for x in $files; do [ -f $x ] && bytes=$((bytes + $(stat -c %s $x))); done
  h=$(zstdcat "$f" 2>/dev/null | awk '/^@/ { print; next } { exit }' | zstd -3 -c | wc -c)
  awk -v s="$run" -v h=$h -v b=$bytes 'BEGIN { printf "%s\t%d\t%d\t%.4f\n", s, h, b, b ? h / b : -1 }'
done | sort -t$'\t' -k4,4gr
