#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH -o /hpc-home/hildebra/dev/protalFastHit/runs/ancestry_truth_test.sh.otxt
#SBATCH -e /hpc-home/hildebra/dev/protalFastHit/runs/ancestry_truth_test.sh.etxt
#SBATCH --mem=64000
#SBATCH --export=ALL
#SBATCH -p "qib-medium,ei-medium"
#SBATCH --time=8:00:00
#SBATCH -J protal_atp
#SBATCH --localscratch=ssd:100
#
# The true-positive test of the ancestry sites, strain alleles and polymorphic sites at full size: two worlds of 36
# test genera, their databases, reads, protal runs and models, then the report
# (scripts/ancestry_truth_test.py; docs/claude/2026-10-09-ancestry-true-positive-test/README.md). It works on the node's
# local SSD and copies the report, the models' outputs, the tables and the logs back to $KEEP as one tarball
# (tens of MB); the rest (worlds, reads, SAMs: a few GB) stays on the scratch.
#
#   sbatch scripts/mini_db/run_ancestry_truth_test.sh
#
# Settings by environment: PDIR (the protal checkout with build/), KEEP (where the tarball goes), THREADS, WORLD_ARGS
# (more simulate_ancestry_world.py options, e.g. "--replicates 2").

echo "SLURM job ID: $SLURM_JOB_ID"
echo $HOSTNAME;
set -eo pipefail

eval "$(micromamba shell hook -s bash)";if [[ $CONDA_DEFAULT_ENV != protal-db-build ]]; then micromamba activate protal-db-build; fi

PDIR=${PDIR:-/hpc-home/hildebra/dev/protalFastHit/}
KEEP=${KEEP:-/hpc-home/hildebra/dev/protalFastHit/runs/ancestry_truth_test_$SLURM_JOB_ID}
THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-32}}
OUT=${SLURM_LOCAL_SCRATCH:-/tmp}/atp

echo "checkout $PDIR at $(git -C $PDIR rev-parse --short HEAD), $THREADS threads, work in $OUT, results to $KEEP"
$PDIR/build/protal --version
mkdir -p $KEEP

python3 $PDIR/scripts/ancestry_truth_test.py run --outdir $OUT -t $THREADS \
    --protal $PDIR/build/protal --simulate $PDIR/build/simulate_metagenomes --train-python python3 \
    --world-args "${WORLD_ARGS:-}" 2>&1 | tee $KEEP/console.log

tar -C $OUT -czf $KEEP/ancestry_truth_report.tar.gz report models test/tables train/tables test/logs train/logs \
    test/world/simulation/genera.tsv test/world/simulation/roles.tsv test/runs/default/protal.log
echo "report: $KEEP/ancestry_truth_report.tar.gz (report/summary.md inside)"
