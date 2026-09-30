#!/bin/bash
set -u
ls ~/strain-build/mini_db/gtdb_r226/ ~/strain-build/mini_db/gtdb_r226/simulation | head -30
wc -l ~/strain-build/mini_db/gtdb_r226/simulation/genomes.tsv
cut -f2 ~/strain-build/mini_db/gtdb_r226/simulation/genomes.tsv | sort | uniq -c
head -3 ~/strain-build/mini_db/gtdb_r226/simulation/genomes.tsv
cat ~/strain-build/mini_db/protal_db/genome2tiid.tsv | head -40
cmp ~/strain-build/mini_db/gtdb_r226/bac120_taxonomy_r226.tsv ~/audit5/world/gtdb_r226/bac120_taxonomy_r226.tsv && echo SAMEWORLD
tail -30 ~/strain-build/mini_db/build.log
which Rscript python3 samtools
ls ~/strain-build/src/scripts | head -40
