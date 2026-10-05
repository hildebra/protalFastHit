#!/usr/bin/env bash
# Builds the 900-species world with the baseline and the change, traced (trace_build.sh), and compares every
# file they leave in the database folder byte for byte, and their logs without timings, memory lines and paths
# as sorted lines (the change moves phases). compare_builds.sh THREADS [protal args]
here=$(cd "$(dirname "$0")" && pwd)
t=$1; shift
BIN=$HOME/buildmem/base/build/protal bash $here/trace_build.sh base_t$t $t "$@" > /dev/null
BIN=$HOME/buildmem/new/build/protal bash $here/trace_build.sh new_t$t $t "$@" > /dev/null
a=$HOME/buildmem/runs/base_t$t; b=$HOME/buildmem/runs/new_t$t
for run in $a $b; do echo "== $(basename $run)"; cat $run/phases.txt; done
echo "== files"
(cd $a/db && find . -type f | sort) > /tmp/files_a; (cd $b/db && find . -type f | sort) > /tmp/files_b
diff /tmp/files_a /tmp/files_b && echo "same $(wc -l < /tmp/files_a) files"
while read -r f; do cmp -s "$a/db/$f" "$b/db/$f" && echo "identical  $f" || echo "DIFFERENT  $f"; done < /tmp/files_a
echo "== logs (sorted, without timings, memory lines and paths)"
clean() { cut -f2 $1/stdout.ts | grep -v -E " took |^Memory after |^Total available memory|^Gene tables: " | sed -e "s#$1/db#DB#g" -e "s#buildmem/[a-z]*/build/protal#BIN#g" | sort; }
diff <(clean $a) <(clean $b) && echo "same lines"
