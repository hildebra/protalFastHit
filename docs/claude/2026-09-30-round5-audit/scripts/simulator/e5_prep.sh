#!/bin/bash
# E5 prep: unpacked copy of the mini DB; a 4000-pair read subset of E1 s_1 with its truth file.
set -u
W=~/audit6/simulator
mkdir -p $W/e5; cd $W/e5
rm -rf dbsrc; mkdir dbsrc
cp ~/strain-build/mini_db/protal_db/database.protal dbsrc/
~/strain-build/bin/protal --unpack_db --db dbsrc/database.protal > unpack.log 2>&1; echo "unpack exit $?"
rm dbsrc/database.protal
ls -la dbsrc
cat dbsrc/internal_taxonomy.dmp
zcat $W/e1/sim/reads/s_1_R1.fq.gz | head -16000 > r1.fq
zcat $W/e1/sim/reads/s_1_R2.fq.gz | head -16000 > r2.fq
cp $W/e1/sim/protal_goldstd/s_1.profile_truth truth.txt
