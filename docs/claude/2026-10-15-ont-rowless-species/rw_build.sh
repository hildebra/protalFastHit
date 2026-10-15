#!/bin/bash
# Small ONT-only build of the real small world (~/realworld, 6 GTDB families) with protal at audit-fixes 55f850d,
# for the rowless-ONT-species study. Runs on cores 0-3, niced.
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
df -h ~ | tail -1
cd ~/ontu
rm -rf rw rw_scratch
mkdir -p rw_scratch
taskset -c 0-3 nice -n 5 python3 src/scripts/build_gtdb_database.py \
  --inputs ~/realworld --release 226 --outdir ~/ontu/rw \
  --protal ~/ontu/build/protal --simulator ~/ontu/build/simulate_metagenomes \
  -t 4 --scratch ~/ontu/rw_scratch --keep-free 2 \
  --read-types ont --model gbm --model-gates none --scenarios none \
  --long-read-bases 300000,1500000,6000000,30000000,150000000 --long-read-samples 12 \
  --test-long-read-bases 150000,1000000,5000000,25000000,250000000 --test-long-read-samples 4 \
  --share-logs
echo "exit $?"
