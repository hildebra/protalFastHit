#!/usr/bin/env bash
# The lookup bench in DRAM (6 GB arena), variants A, B, C alternated, 6 rounds; then in cache, 6 rounds.
B=$HOME/perf6/bench/lookup_bench
echo "== DRAM 6 GB, $(date +%T), load $(cut -d' ' -f1 /proc/loadavg)"
$B 6 150000 6 ABC
echo "== in cache, $(date +%T), load $(cut -d' ' -f1 /proc/loadavg)"
$B 0.004 150000 6 ABC
echo "ABCDONE $(date +%T)"
