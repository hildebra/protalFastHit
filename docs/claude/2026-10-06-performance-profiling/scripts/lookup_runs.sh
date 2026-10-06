#!/usr/bin/env bash
# The lookup bench in cache (a 4 MB arena) and in DRAM (a 6 GB arena), variants alternated, 3 rounds each.
B=$HOME/perf6/bench/lookup_bench
echo "== in cache, $(date +%T), load $(cut -d' ' -f1 /proc/loadavg)"
$B 0.004 200000 3
echo "== DRAM 6 GB, $(date +%T), load $(cut -d' ' -f1 /proc/loadavg)"
$B 6 200000 3
grep -i AnonHugePages /proc/meminfo
echo "LOOKUPDONE $(date +%T)"
