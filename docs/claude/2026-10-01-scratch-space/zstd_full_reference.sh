#!/bin/bash
# zstd settings for full_reference.fna: ratio, compression and decompression speed (WSL, 2026-10-01).
# Data: the synthetic GTDB-like tuning world's full reference (~/tune/V3/protal_db/full_reference.fna,
# species with simulated strains); real GTDB copies of a gene are more alike within species and further
# apart in the file, so the ratio is only a hint.
F=${1:-~/tune/V3/protal_db/full_reference.fna}
T=/tmp/zstd_full_ref_test.zst
ls -la "$F"
grep -c '>' "$F" | awk '{print $1, "records"}'
for opts in "-3" "-3 --long=27" "-6 --long=27" "-9 --long=27" "-12 --long=27" "-19 --long=27"; do
    start=$(date +%s.%N)
    nice -n 10 zstd -q -f -T4 $opts "$F" -o "$T"
    mid=$(date +%s.%N)
    zstd -q -dc "$T" > /dev/null
    end=$(date +%s.%N)
    awk -v o="$opts" -v a=$(stat -c %s "$F") -v b=$(stat -c %s "$T") -v c="$start" -v m="$mid" -v e="$end" \
        'BEGIN{printf "%-14s ratio %5.2f  compress %6.0f MB/s (4 threads)  decompress %6.0f MB/s (1 thread)\n", o, a/b, a/1e6/(m-c), a/1e6/(e-m)}'
done
rm -f "$T"
