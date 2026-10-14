#!/bin/bash
# Time main.cpp.o: preprocess, front end only, -O0, -O3 (the real one) with -ftime-report.
W=~/ctime
B=$W/build
cd $B
CMD=$(grep -o '"command": "[^"]*protal.dir/src/main.cpp[^"]*"' compile_commands.json | head -1 | sed 's/^"command": "//; s/"$//')
BASE=${CMD% -o *}
SRC=$W/src/src/main.cpp
R=$W/res
mkdir -p $R
T="taskset -c 0-3 nice -n 5"
uptime
# Preprocessed size
$T $BASE -E $SRC -o $R/main.ii
echo "preprocessed lines: $(wc -l < $R/main.ii), bytes: $(wc -c < $R/main.ii)"
# Lines per origin file (by # line markers), top 40
awk '/^# [0-9]+ "/ {f=$3; next} {n[f]++} END {for (k in n) print n[k], k}' $R/main.ii | sort -rn > $R/lines_by_file.txt
head -40 $R/lines_by_file.txt
# Front end only
/usr/bin/time -f "syntax-only: %e s, %M KB" $T $BASE -fsyntax-only $SRC
# -O0
/usr/bin/time -f "O0: %e s, %M KB" $T ${BASE/-O3/-O0} -c $SRC -o $R/main_O0.o
# The real one, with the time report
/usr/bin/time -f "O3: %e s, %M KB" $T $BASE -ftime-report -c $SRC -o $R/main_O3.o 2> $R/ftime_O3.txt
tail -1 $R/ftime_O3.txt
ls -la $R/main_O0.o $R/main_O3.o
uptime
