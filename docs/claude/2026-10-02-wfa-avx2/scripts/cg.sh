#!/usr/bin/env bash
# Where WFA2 spends its instructions in wfabench (one round, the mix): callgrind, inclusive, WFA2's functions.
W=$HOME/mt-work/wfabench; cd $W
valgrind --tool=callgrind --callgrind-out-file=$W/cg.out $W/wfabench 1 > $W/cg.log 2>&1
callgrind_annotate --inclusive=yes $W/cg.out 2>/dev/null | grep -E "PROGRAM TOTALS|wavefront_|protal_wfa|WFA2Wrapper2::Alignment|alignEndsFree" | head -30
