#!/bin/bash
# Alignment stage (--no_profile) at 1 and N threads, with gzip and with plain FASTQ input:
#   scaling.sh DB READS_DIR LABEL [THREADS]
# The gzip/plain difference is the cost of inflating the input inside the reader's critical section.
source "$(dirname "$0")/env.sh"
db=$1; src=$2; ds=$3; t=${4:-$(nproc)}
plain=$PERF_DIR/reads/${ds}_plain; mkdir -p $plain
for m in R1 R2; do [ -s $plain/${ds}_$m.fq ] || zcat $src/*_$m.fq.gz > $plain/${ds}_$m.fq; done
for m in R1 R2; do /usr/bin/time -f "zcat $m alone: %e s wall" zcat $src/*_$m.fq.gz > /dev/null; done
bash $here/run_one.sh align_${ds}_t$t $db $src $t --no_profile --verbose
bash $here/run_one.sh align_${ds}_plain_t$t $db $plain $t --no_profile --verbose
bash $here/run_one.sh align_${ds}_t1 $db $src 1 --no_profile --verbose
