#!/bin/bash
# snapshot of a build output folder: file names, sizes, mtimes, md5 of the small ones
OUT=~/audit6/gtdb/${1:-b1}
TAG=${2:-snap}
{
for d in $OUT $OUT/protal_db $OUT/training_db; do
  echo "== $d"; ls -la --time-style=+%T $d
done
echo "== points"; ls $OUT/training/points 2>/dev/null
for p in $OUT/training/points/*; do echo "$p: $(find $p -name '*.truth_annotated' | wc -l) dumps, $(find $p -name '*.sam.zst' | wc -l) sams"; done
md5sum $OUT/heldout_species.txt $OUT/training_db/database.protal $OUT/protal_db/database.protal 2>/dev/null
} > ~/audit6/gtdb/$TAG.txt 2>&1
cat ~/audit6/gtdb/$TAG.txt
