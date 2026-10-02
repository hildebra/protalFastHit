#!/bin/bash
# A stand-in for art_illumina or pbsim that runs the real one and appends its time to $TOOL_TIMES:
# tool, wall seconds, user seconds, system seconds, then the arguments (genome, coverage or depth).
#
# usage: a two-line script per tool that sets TOOL_NAME, TOOL_REAL (the real binary) and TOOL_TIMES and
#        execs this one (make_wrappers.sh); give it to simulate_metagenomes --art_path or the collector's --pbsim
name=${TOOL_NAME:-$(basename "$0")}
real=${TOOL_REAL:?TOOL_REAL: the real tool}
times=${TOOL_TIMES:?TOOL_TIMES: where to append the times}
start=$(date +%s.%N)
/usr/bin/time -f '%U\t%S' -o "$times.$$" "$real" "$@"
rc=$?
end=$(date +%s.%N)
printf '%s\t%s\t%s\t%s\n' "$name" "$(echo "$end - $start" | bc)" "$(tail -1 "$times.$$")" "$*" >> "$times"
rm -f "$times.$$"
exit $rc
