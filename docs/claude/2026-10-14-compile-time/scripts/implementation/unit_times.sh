#!/bin/bash
# Each translation unit compiled alone on one core (core 2), with -ftime-report: the split's units and Options.cpp
# (~/split-build) and the baseline's main.cpp (~/split-base), alternating, twice each; their own compile commands.
OUT=~/split-unit-times.tsv
cmd_of() {  # build dir, source path suffix -> its compile command (without -o)
    local c
    c=$(grep -o "\"command\": \"[^\"]*/$2\"" $1/compile_commands.json | head -1 | sed 's/^"command": "//; s/"$//')
    echo "${c%% -o *}"
}
dir_of() { grep -B1 "\"command\": \"[^\"]*/$2\"" $1/compile_commands.json | grep -o '"directory": "[^"]*"' | head -1 | sed 's/"directory": "//; s/"$//'; }
one() {  # label, build dir, source suffix, source path
    local label=$1 b=$2 suffix=$3 src=$4 cmd dir
    cmd=$(cmd_of $b $suffix); dir=$(dir_of $b $suffix)
    (cd ${dir:-$b} && /usr/bin/time -f "%e\t%U" -o /tmp/ut_time.txt taskset -c 2 nice -n 5 $cmd -ftime-report -c $src -o /tmp/ut.o 2> /tmp/ut_report.txt)
    local fe opt ggc text
    fe=$(awk -F: '/phase parsing|phase lang. deferred/ {split($2, a, " "); s += a[1]} END {print s}' /tmp/ut_report.txt)
    opt=$(awk -F: '/phase opt and generate/ {split($2, a, " "); print a[1]}' /tmp/ut_report.txt)
    ggc=$(grep -E '^ TOTAL' /tmp/ut_report.txt | awk '{print $NF}')
    text=$(size /tmp/ut.o | awk 'NR==2 {print $1}')
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$label" "$(cut -f1 /tmp/ut_time.txt)" "$(cut -f2 /tmp/ut_time.txt)" "$fe" "$opt" "$ggc" "$text" | tee -a $OUT
}
printf 'unit\twall_s\tuser_s\tfront_end_user_s\topt_user_s\tggc\ttext_bytes\n' > $OUT
uptime
for round in 1 2; do
    one "before: main.cpp" ~/split-base/build src/main.cpp ~/split-base/src/src/main.cpp
    for u in Run/Run.cpp Run/Load.cpp Run/Align.cpp Run/BuildStage.cpp Run/Profile.cpp Run/Strains.cpp Options.cpp main.cpp; do
        one "after: $u" ~/split-build/build src/$u ~/split-build/src/src/$u
    done
done
uptime
rm -f /tmp/ut.o /tmp/ut_time.txt /tmp/ut_report.txt
