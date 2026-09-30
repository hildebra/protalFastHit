#!/bin/bash
# Make a synthetic GTDB-like world, a fake GTDB mirror of it, and the download folder (fake NCBI datasets).
exec > ~/audit6/gtdb/world.out 2>&1
set -eux
PY=~/protal-train/bin/python
S=~/audit6/gtdb/src/scripts
W=~/audit6/gtdb/w
rm -rf $W; mkdir -p $W
cd $W
$PY $S/mini_db/gtdb_like_lineages.py --species 120 --archaea 0.1 --seed 1 > lineages.txt
wc -l lineages.txt
cut -d';' -f1 lineages.txt | sort | uniq -c
time $PY $S/mini_db/simulate_gtdb_release.py --outdir $W/gtdb --lineages lineages.txt --genomes_per_species 3 \
   --genome_length 60000 --strain_divergence 0.002-0.015 --species_divergence 0.015-0.04 --seed 3
du -sh $W/gtdb
# fake mirror, with the helpers of the repository's own test
cd $S/mini_db
$PY - <<EOF
import sys, os
sys.path.insert(0, "$S/mini_db")
import test_mini_db as t
t.gtdb_mirror("$W/gtdb", "$W/mirror")
open("$W/datasets", "w").write(t.FAKE_DATASETS.replace("#!/usr/bin/env python3", "#!$PY"))
os.chmod("$W/datasets", 0o755)
EOF
find $W/mirror | head -30
cat $W/mirror/release226/226.0/MD5SUM.txt
