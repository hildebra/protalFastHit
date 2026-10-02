#!/usr/bin/env bash
# The stage timers with the final patch (lookups of the next read timed as seeding) against HEAD's build:
# 500k pairs at -t 1, alternated, three times each; the seeding rows of s_runtime.tsv and the alignment time.
W=$HOME/mt-work/pfrun; DB=$HOME/bench071/V071/protal_db
R5=$HOME/bench071/samples/points/rl150_p500000/sim/reads/rl150_p500000_s_1
declare -A BIN=( [ref]=$HOME/mt-work/pfref/src/build/protal [pf2]=$HOME/mt-work/pfcheck2/src/build/protal )
for rep in 1 2 3; do
  order="ref pf2"; [ $((rep % 2)) = 0 ] && order="pf2 ref"
  for b in $order; do
    d=$W/t.$b; rm -rf $d
    ${BIN[$b]} --db $DB -1 ${R5}_R1.fq.gz -2 ${R5}_R2.fq.gz --no_profile --prefix s -o $d -t 1 --no_qcmsa > $W/tlog 2>&1 || echo "FAIL $b"
    al=$(grep 'Aligning reads took' $W/tlog | sed 's/.*took //')
    rows=$(awk -F'\t' '$1 == "Seeding" || $1 == "Seed-finding operator" || $1 == "Seed- and Anchor-finding" { printf "%s=%s  ", $1, $2 }' $d/misc/s_runtime.tsv)
    echo -e "$rep\t$b\talign $al\t$rows"
  done
done
