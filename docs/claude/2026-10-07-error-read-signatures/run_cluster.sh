#!/bin/bash
# The paired-end FP/FN error reads of the r226 v15 build, with their alignment features (error_read_features.py),
# from the samples the build left on q512n10's local SSD; packed into one archive to copy back.
#
#   sbatch docs/claude/2026-10-07-error-read-signatures/run_cluster.sh           # from the protal checkout
#   RT=se sbatch ...                                                             # another read type
#
# Settings by environment: BUILD (the build's OUTDIR), SCRATCH (its --scratch), RT (read type), THREADS, MAXFRAG.
#SBATCH -J v15_error_reads
#SBATCH --nodelist=q512n10
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 8:00:00
#SBATCH -o v15_error_reads_%j.log
set -euo pipefail

BUILD=${BUILD:-/hpc-home/hildebra/DB/protal/protal0.7.8_r226_v15}
SCRATCH=${SCRATCH:-/nbi/local/ssd/24027266/protalDBbuild_1}
RT=${RT:-pe}
THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-16}}
MAXFRAG=${MAXFRAG:-20}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
[ -f "$HERE/error_read_features.py" ] || HERE=${SLURM_SUBMIT_DIR:-$PWD}/docs/claude/2026-10-07-error-read-signatures
OUT=$BUILD/error_read_features/$RT

echo "node $(hostname), build $BUILD, scratch $SCRATCH, read type $RT, $THREADS threads"
for need in "$SCRATCH/training/points" "$SCRATCH/training_db/genome2tiid.tsv"; do
    if [ ! -e "$need" ]; then
        echo "missing: $need -- the build's samples are not on this node (any more)." >&2
        echo "Run on the node of the build (q512n10 for job 24027266), or rerun the build with the same --scratch." >&2
        exit 1
    fi
done
[ -d "$SCRATCH/test/points" ] || echo "no $SCRATCH/test/points: training samples only"
du -sh "$SCRATCH/training/points" "$SCRATCH/test/points" 2>/dev/null || true

python3 "$HERE/error_read_features.py" --build "$BUILD" --scratch "$SCRATCH" --read-type "$RT" \
    --max-fragments "$MAXFRAG" --threads "$THREADS" --out "$OUT"

tar -C "$BUILD/error_read_features" -czf "$BUILD/error_read_features_${RT}.tar.gz" "$RT"
ls -l "$BUILD/error_read_features_${RT}.tar.gz"
echo "copy $BUILD/error_read_features_${RT}.tar.gz to local/v15/ and unpack it there"
