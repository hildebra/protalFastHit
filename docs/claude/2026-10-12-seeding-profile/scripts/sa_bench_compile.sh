#!/usr/bin/env bash
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
taskset -c 0-3 nice -n 5 bash "$A/seedbench_build.sh" "$A/seedbench.cpp" $HOME/seed/head/src $HOME/seed/head/build 2 > $HOME/seed/bench_build.log 2>&1
echo "exit $?"; grep -E "error|warning: unused" -A3 $HOME/seed/bench_build.log | head -40; tail -2 $HOME/seed/bench_build.log
