#!/bin/bash
set -u
R=~/strain-build/mini_db/protal_db/full_reference.fna
grep '>' $R | head -5
grep '>' $R | cut -d_ -f1 | sort | uniq -c
grep '>' $R | sort | uniq -d | head -3
grep '>' $R | awk -F_ '{print $1}' | sort | uniq -c
cat ~/strain-build/mini_db/protal_db/genome2tiid.tsv | wc -l
head -3 ~/strain-build/mini_db/gtdb_r226/simulation/divergence.tsv
wc -l ~/strain-build/mini_db/gtdb_r226/simulation/divergence.tsv
ls ~/strain-build/src/scripts ~/strain-build/src/scripts/mini_db
ls ~/strain-build/src/examples/mini_db
cat ~/strain-build/src/examples/mini_db/community.tsv
~/strain-build/bin/protal --help 2>&1 | head -150
