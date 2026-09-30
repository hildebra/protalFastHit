#!/bin/bash
# Experiment 5: FASTA input (no qualities) of the e1 HiFi and ONT reads: constant QUAL, records as
# for FASTQ; profiles. Also e1's profiles against the community truth.
S=$(dirname "$0")
W=~/audit6/longreads
P=${P:-~/strain-build/bin/protal}
G=~/strain-build/mini_db/gtdb_r226
D=$W/mini_db
OUT=$S/out_e5.txt
cd $W/e1
{
$P --db $D -1 hifi.fa --read_type pb --model_pb $D/model_pe.xml -o outfa --prefix hifi -t 2 --sam_format sam --no_qcmsa --no_strains > hifi_fa.log 2>&1; echo "hifi fasta rc=$?"
$P --db $D -1 ont.fa --read_type ont --model_ont $D/model_pe.xml -o outfa --prefix ont -t 2 --sam_format sam --no_qcmsa --no_strains > ont_fa.log 2>&1; echo "ont fasta rc=$?"
python3 $S/eval_sam.py outfa/hifi.sam hifi.fa hifi.truth.pkl --gtdb $G --db $D --fasta_qual '?' --show 2 | grep -vE "^    \["
python3 $S/eval_sam.py outfa/ont.sam ont.fa ont.truth.pkl --gtdb $G --db $D --fasta_qual '3' --exact 0.95 --show 2 | grep -vE "^    \["
echo "== same records FASTQ vs FASTA (fields 1-6)?"
for s in hifi ont; do
  diff <(grep -v '^@' out/$s.sam | cut -f1-6 | sort) <(grep -v '^@' outfa/$s.sam | cut -f1-6 | sort) > /dev/null && echo "$s: identical" || echo "$s: differ ($(diff <(grep -v '^@' out/$s.sam | cut -f1-6 | sort) <(grep -v '^@' outfa/$s.sam | cut -f1-6 | sort) | grep -c '^[<>]') lines)"
done
echo "== profiles (FASTQ)"
for s in hifi ont; do echo "-- $s"; cat out/$s.profile; done
echo "== profiles (FASTA)"
for s in hifi ont; do echo "-- $s"; cat outfa/$s.profile; done
} > $OUT 2>&1
echo done
