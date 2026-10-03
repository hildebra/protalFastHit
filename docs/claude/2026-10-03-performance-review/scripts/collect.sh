#!/usr/bin/env bash
# results/ of the report: callgrind tables, stage timer lines of every run.
set -uo pipefail
W=$HOME/mt-work/perf4
R=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-03-performance-review/results
head -150 $W/cg.align.incl.txt | cut -c1-220 > $R/callgrind_align_pe100k_inclusive.txt
head -150 $W/cg.align.self.txt | cut -c1-220 > $R/callgrind_align_pe100k_self.txt
head -120 $W/cg.prof.incl.txt | cut -c1-220 > $R/callgrind_profile_pe500k_inclusive.txt
head -120 $W/cg.prof.self.txt | cut -c1-220 > $R/callgrind_profile_pe500k_self.txt
{ for n in pe500k_t1 pe500k_t6 pe5M_t6 ont90M_t6 pb90M_t6 pe500k_t1_np pe500k_t3_np pe500k_t6_np pe500k_t6_plain map4; do
    echo "=== $n"; grep -E "Elapsed|User time|Maximum resident" $W/o.$n.time | tr -s ' '; grep -E "took" $W/o.$n.log | sed 's/.*left    //' | grep -vE "^\s*$"; echo; done; } > $R/stage_times.txt
{ for n in pe500k_t1 pe500k_t6 pe5M_t6 ont90M_t6 pb90M_t6; do echo "=== $n"; cat $W/o.$n/misc/s_runtime.tsv; echo; done; } > $R/runtime_tsv.txt
cp $W/map4.tsv $R/map4.tsv
echo collected; wc -l $R/*.txt
