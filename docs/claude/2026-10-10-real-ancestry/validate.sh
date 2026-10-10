#!/bin/bash
# The column weights' mappings of the <set> databases (build_db.sh) against GTDB's own alignments of the
# representatives' proteins; usage: validate.sh <set> [alignment ...] (default: aa nt). Outputs <set>_<align>_validation.*.
set -eu
SET=$1; shift; ALIGNS=${*:-aa nt}
HERE=$(cd "$(dirname "$0")" && pwd)
PY=~/micromamba/envs/protal-db-build/bin/python
for ALIGN in $ALIGNS; do
  DB=~/real_ancestry/${SET}_db_$ALIGN
  [ -s $DB/unpacked/column_weights.tsv ] || ~/bidx/build/protal --db $DB/database.protal --unpack_db --unpack_dir $DB/unpacked > /dev/null
  cp $DB/genome2tiid.tsv $DB/gene2geneid.tsv $DB/unpacked/
  taskset -c 0-3 nice -n 5 $PY $HERE/validate_alignments.py $DB/unpacked \
      ~/GTDB/r226/genomic_files_reps/${SET}_msa_marker_genes_reps_r226.tar.gz $HERE/${SET}_${ALIGN}_validation
done
