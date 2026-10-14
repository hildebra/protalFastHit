#!/bin/bash
# Rebuild times of protal after an edit, ab6bf25 (one main.cpp) against the split (src/Run/ units): the file is
# touched in the build's source copy and `ninja -j4 protal` timed (wall, and CPU of all its processes), the two trees
# alternating per scenario so that both see the same load. Cores 0-3.
BASE=~/split-base
SPLIT=~/split-build
OUT=~/split-times.tsv
T="taskset -c 0-3 nice -n 5"
build() {  # tree, label, file
    local tree=$1 label=$2 file=$3
    touch $tree/src/$file
    local log=$tree/rebuild_$(echo $label | tr ' /' '__').log
    /usr/bin/time -f "%e\t%U\t%S" -o $tree/time.txt $T ninja -C $tree/build -j4 protal > $log 2>&1
    local steps
    steps=$(grep -c '^\[' $log)
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$label" "$(basename $tree)" "$file" "$(cat $tree/time.txt)" "$steps" "$(uptime | sed 's/.*average: //')" | tee -a $OUT
}
# Both up to date first.
$T ninja -C $BASE/build -j4 protal > /dev/null 2>&1
$T ninja -C $SPLIT/build -j4 protal > /dev/null 2>&1
printf 'scenario\ttree\tfile\twall_s\tuser_s\tsys_s\tsteps\tload\n' > $OUT
scenario() {  # label, file in base, file in split
    build $BASE "$1" "$2"
    build $SPLIT "$1" "$3"
}
clean_build() {  # tree, label: everything of protal compiled again (ninja -t clean), the libraries too
    local tree=$1 label=$2
    ninja -C $tree/build -t clean protal > /dev/null 2>&1
    /usr/bin/time -f "%e\t%U\t%S" -o $tree/time.txt $T ninja -C $tree/build -j4 protal > $tree/clean.log 2>&1
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$label" "$(basename $tree)" "-" "$(cat $tree/time.txt)" "$(grep -c '^\[' $tree/clean.log)" "$(uptime | sed 's/.*average: //')" | tee -a $OUT
}
for round in 1 2; do
    clean_build $BASE "from scratch (libraries too)"
    clean_build $SPLIT "from scratch (libraries too)"
    scenario "profiler header" src/Profiling/Profiler.h src/Profiling/Profiler.h
    scenario "alignment header" src/Core/ChainAnchorFinder.h src/Core/ChainAnchorFinder.h
    scenario "build header" src/Build.h src/Build.h
    scenario "strain stage code" src/RunProtal.h src/Run/Strains.cpp
    scenario "run start-up code" src/RunProtal.h src/Run/Run.cpp
    scenario "version history" src/main.cpp src/main.cpp
    scenario "options header (all)" src/Options.h src/Options.h
    scenario "genome loader (all)" src/SequenceUtils/GenomeLoader.h src/SequenceUtils/GenomeLoader.h
done
