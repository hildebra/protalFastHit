#!/bin/bash
set -u
mkdir -p ~/audit6/simulator
ls ~/strain-build/bin/
ls ~/audit5/world/gtdb_r226/simulation/ | head
wc -l ~/audit5/world/gtdb_r226/simulation/genomes.tsv
head -3 ~/audit5/world/gtdb_r226/simulation/genomes.tsv
ls ~/strain-build/mini_db/protal_db | head -30
ls ~/audit5/accuracy/prot_A_P/profiles/ | head
ls ~/genes_study/cong/out/profiles/ | head
which art_illumina; art_illumina 2>&1 | head -5
nproc; uptime
