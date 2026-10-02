#!/bin/bash
# The stages of a build_real_size.sh run: the build script's lines with their times, each collection's
# simulation and protal lines (the collector's log, and its simulation log when the build simulated in the
# background), the ART and pbsim3 CPU per collection (timed_tool.sh), the trainer's study times, and the
# whole run's CPU (/usr/bin/time).
# usage: build_summary.sh RUN_DIR [CONSOLE_LOG]
run=$1 db=$1/db console=${2:-$1/build.log}
here=$(cd "$(dirname "$0")" && pwd)
echo "== console (stages)"
grep -E "Converted|Found the genes|Wrote the files|Built|Collected|Trained|  (pe|se|pb|ont): |scores the|Adding|Ready|[0-9]/[0-9] |^\[[^]]*\]     " \
  "$console" | cut -c1-160
for c in training test; do
  echo "== $c collection"
  for log in "$db/${c}_data_simulation.log" "$db/${c}_data.log"; do
    [ -f "$log" ] && grep -E "simulating|simulated .* in|simulated \(|profiling|profiled" "$log" | grep -v "^simulating [a-z]*_b" | cut -c1-120
  done
  echo "-- ART and pbsim3 CPU seconds (all calls of the $c collection)"
  awk -F'\t' -v c="/$c/" 'index($5, c) { cpu[$1] += $3 + $4; wall[$1] += $2; n[$1]++
      if ($1 == "pbsim") { t = ($5 ~ /qshmm/) ? "pbsim ont" : "pbsim pb"; cpu[t] += $3 + $4; n[t]++ } }
    END { for (k in cpu) printf "  %-14s %6d calls, %8.0f CPU s, %8.0f wall s\n", k, n[k], cpu[k], wall[k] }' "$run/tool_times.tsv" | sort
  [ -f "$db/$c/profile_all/protal.log" ] && bash "$here/protal_log_summary.sh" "$db/$c/profile_all/protal.log" | grep -v "^  align"
done
echo "== trainer"
grep -h "^time:" "$db"/classifier_training*.log
echo "== whole run"
grep -hE "Elapsed|User time|System time|Maximum resident" "$run"/build*.time
