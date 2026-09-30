#!/bin/bash
# Peak and end RSS of the stages of a mem_trace.sh run: stage_peaks.sh RUN_DIR
# Stages are cut at protal's own timing lines ("... took", seconds since start in stdout.ts).
d=$1
t_align_end=$(awk -F'\t' '/Processing all samples took/{print $1}' $d/stdout.ts | head -1)
t_index_end=$(awk -F'\t' '/Load Index took/{print $1}' $d/stdout.ts | head -1)
t_prof_end=$(awk -F'\t' '/Profiling took/{print $1}' $d/stdout.ts | head -1)
awk -F'\t' -v a=$t_index_end -v b=$t_align_end -v c=$t_prof_end '
 { rss=$2/1024; anon=$3/1024; if ($1<a) {s="load"} else if ($1<b) {s="align"} else if ($1<c) {s="profile"} else {s="end"}
   if (rss>peak[s]) peak[s]=rss; last[s]=rss; if (!(s in seen)) {seen[s]=1; order[++k]=s} }
 END { for (i=1;i<=k;i++) printf "%-8s peak %7.0f MB  last %7.0f MB\n", order[i], peak[order[i]], last[order[i]] }' $d/mem.tsv
