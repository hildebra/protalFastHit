#!/bin/bash
# The three e1 ONT reads in the x_drop runs of e8 (-t 2, whole sample), and the x-drop setting logged.
S=$(dirname "$0")
cd ~/audit6/longreads/e8
{
for d in out_xd1000 out_xd0; do
  echo $d
  grep -v '^@' $d/ont.sam | awk '($1=="ont_124"||$1=="ont_152"||$1=="ont_300") && $3 ~ /_(38|103|82)$/ {print "   ", $1,$2,$3,$5,length($10)}'
done
grep -E 'x-drop' ont.xd0.log ont.xd1000.log
python3 $S/eval_sam.py out_xd0/ont.sam ../e1/ont.fq ../e1/ont.truth.pkl --gtdb ~/strain-build/mini_db/gtdb_r226 --db ../mini_db --exact 0.95 --show 5 | grep -A5 "^wrong_taxon	" | cut -c1-120
grep -B2 'Invalid after' ont.xd0.err | grep Record
} > $S/out_e8b.txt 2>&1
