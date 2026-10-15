#!/bin/bash
# Small ONT-only builds of the real small world that differ only in protal's identity floor: 0.85 (the ONT default; no
# -a), 0.80, 0.75. protal at 55f850d does not record -a in the database (run settings in build_settings.tsv came with
# 0.8.5), so a wrapper gives -a to every protal call of the build: the samples are all ONT, so it sets the ONT floor of
# the collections, the parity check and the reports alike. Same seed, same samples. Cores 0-3, niced, one at a time.
# Usage: exp_builds.sh FLOOR... (default: 0.80 0.75; the baseline 0.85 ran with an earlier version whose wrapper gave -a
# to --build only, which is the ONT default anyway)
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
cd ~/ontu
for A in ${@:-0.80 0.75}; do
  W=~/ontu/protal_a$A
  printf '#!/bin/bash\nexec /home/falk/ontu/build/protal "$@" -a %s\n' "$A" > $W
  chmod +x $W
  O=~/ontu/rwA$A
  rm -rf $O ${O}_scratch && mkdir -p ${O}_scratch
  echo "=== floor $A $(date +%T)"
  taskset -c 0-3 nice -n 5 python3 src/scripts/build_gtdb_database.py \
    --inputs ~/realworld --release 226 --outdir $O \
    --protal $W --simulator ~/ontu/build/simulate_metagenomes \
    -t 4 --scratch ${O}_scratch --keep-free 2 \
    --read-types ont --model gbm --model-gates none --scenarios none \
    --long-read-bases 300000,1500000,6000000,30000000,150000000 --long-read-samples 24 \
    --test-long-read-bases 150000,1000000,5000000,25000000,250000000 --test-long-read-samples 8 \
    > $O.log 2>&1
  echo "exit $? $(date +%T)"
  grep -E "F1 with|marker genes got reads ->" $O.log | cut -c1-250
  grep -m1 "max score ani" ${O}_scratch/training/profile_all/run1/protal.log
  md5sum ${O}_scratch/training/training_data_ont.tsv | cut -c1-12
done
echo done
