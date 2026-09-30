#!/bin/bash
# Experiment 11: a mini database of close relatives (species 0.5% from their genus, ~1% apart;
# strains 0.2%), HiFi and ONT reads of the community: read-level consensus (ZR:i:1) and MAPQ.
S=$(dirname "$0")
W=~/audit6/longreads
P=~/strain-build/bin/protal
OUT=$S/out_e11.txt
cd ~/strain-build/src
{
if [ ! -s $W/mini_db_close/protal_db/database.protal ]; then
  PROTAL=$P bash scripts/mini_db/build_mini_db.sh $W/mini_db_close --species_divergence 0.005 --strain_divergence 0.002 > $W/minidb_close.log 2>&1 || { tail -20 $W/minidb_close.log; exit 1; }
  $P --unpack_db --db $W/mini_db_close/protal_db/database.protal -t 2 > $W/minidb_close_unpack.log 2>&1
fi
tail -2 $W/minidb_close.log
G=$W/mini_db_close/gtdb_r226
D=$W/mini_db_close/protal_db
ls $D
mkdir -p $W/e11 && cd $W/e11
python3 $S/gen_reads.py community hifi --gtdb $G --db $D --reads 400 --platform hifi --seed 31 --name "m64001_c/{}/ccs"
python3 $S/gen_reads.py community ont --gtdb $G --db $D --reads 300 --platform ont --seed 32 --median 20000 --sigma 0.6 --name "ontc_{}"
for s in hifi:pb ont:ont; do
  n=${s%%:*}; t=${s#*:}
  for xd in 1000 0; do
    rm -rf out_xd$xd/$n.*
    $P --db $D -1 $n.fq --read_type $t --model_$t $D/model_pe.xml -o out_xd$xd --prefix $n -t 2 --sam_format sam --no_qcmsa --no_strains --x_drop $xd > $n.xd$xd.log 2> $n.xd$xd.err
    echo "######## $n --x_drop $xd rc=$? $(grep settled $n.xd$xd.log) invalid=$(grep -c 'Invalid after alignment' $n.xd$xd.err)"
    python3 $S/eval_sam.py out_xd$xd/$n.sam $n.fq $n.truth.pkl --gtdb $G --db $D --exact 0.95 --show 3 | cut -c1-220
    cut -f1,3 out_xd$xd/$n.profile
  done
done
python3 - <<'EOF'
import pickle, collections
for f in ("hifi", "ont"):
    T = pickle.load(open(f + ".truth.pkl", "rb"))
    b = collections.Counter()
    for t in T["reads"]:
        b[t["acc"]] += t["length"]
    print(f, "bases per genome", dict(b))
EOF
} > $OUT 2>&1
echo done
