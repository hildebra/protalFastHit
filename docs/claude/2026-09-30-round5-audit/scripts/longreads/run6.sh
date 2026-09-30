#!/bin/bash
# Experiment 6: one map with pe, pb and ont samples of two Mockella alpha strains, profiled with
# the pe model standing in for all (--model_pb/--model_ont); strain MSA rows against the strains'
# true genes. Then the same with a placeholder pb model (scores 0): the pb samples' taxa must not
# pass in the statistics or enter the MSA.
S=$(dirname "$0")
W=~/audit6/longreads
P=${P:-~/strain-build/bin/protal}
G=~/strain-build/mini_db/gtdb_r226
D=$W/mini_db
OUT=$S/out_e6.txt
mkdir -p $W/e6 && cd $W/e6
{
A=GCA_999001002.1; B=GCA_999001003.1
for x in a:$A b:$B; do
  s=${x%%:*}; acc=${x#*:}
  [ -f pe_${s}_R1.fq ] || { printf "accession\trelative_abundance\n$acc\t1\n" > comm_$s.tsv; python3 $W/src/scripts/mini_db/simulate_reads.py --genomes $G/simulation/genomes.tsv --community comm_$s.tsv --out_prefix pe_$s --pairs 20000 --seed 3 > /dev/null; }
  [ -f pb_$s.fq ] || python3 $S/gen_reads.py community pb_$s --gtdb $G --db $D --community $acc:1 --reads 380 --platform hifi --median 15000 --sigma 0.2 --seed 11 --name "pb_${s}_{}"
  [ -f ont_$s.fq ] || python3 $S/gen_reads.py community ont_$s --gtdb $G --db $D --community $acc:1 --reads 280 --platform ont --median 20000 --sigma 0.5 --seed 12 --name "ont_${s}_{}"
done
for run in pe_model placeholder; do
  mkdir -p $run
  printf "#OUTPUT_DIR\t$W/e6/$run\n#SAMPLEID\tPREFIX\tFIRST\tSECOND\tREAD_TYPE\n" > $run/samples.map
  for s in a b; do
    printf "pe_$s\tpe_$s\t$W/e6/pe_${s}_R1.fq\t$W/e6/pe_${s}_R2.fq\tpe\n" >> $run/samples.map
    printf "pb_$s\tpb_$s\t$W/e6/pb_$s.fq\t-\tpb\n" >> $run/samples.map
    printf "ont_$s\tont_$s\t$W/e6/ont_$s.fq\t-\tont\n" >> $run/samples.map
  done
done
python3 $W/src/scripts/placeholder_models.py -o models --read_types pb > /dev/null
$P --db $D --map pe_model/samples.map --model_pb $D/model_pe.xml --model_ont $D/model_pe.xml -t 2 --sam_format sam --no_qcmsa > pe_model.log 2>&1; echo "pe_model rc=$?"
$P --db $D --map placeholder/samples.map --model_pb models/model_PB.xml --model_ont $D/model_pe.xml -t 2 --sam_format sam --no_qcmsa > placeholder.log 2>&1; echo "placeholder rc=$?"
for run in pe_model placeholder; do
  echo "=========== $run"
  grep -E "read types|Model of|Align the|settled|minimum allele|Warning|placeholder" $run.log | cut -c1-200
  for f in $run/profiles/*.profile $run/*.profile; do [ -f $f ] && { echo "-- $f"; cat $f; }; done 2>/dev/null
  ls $run $run/misc $run/strains 2>/dev/null | head -40
  for f in $run/misc/*alpha*.statistics.tsv; do echo "-- $f"; cat $f; done
  msa=$(ls $run/strains/*alpha*.msa.fna 2>/dev/null | grep -v '\.msa\.fna\.' | head -1)
  part=$(ls $run/strains/*alpha*partition* 2>/dev/null | head -1)
  echo "msa=$msa part=$part"
  [ -n "$msa" ] && grep '>' $msa
  [ -n "$msa" ] && python3 $S/eval_msa.py $msa $part --gtdb $G --db $D --taxid 2 --sample pe_a=$A --sample pb_a=$A --sample ont_a=$A --sample pe_b=$B --sample pb_b=$B --sample ont_b=$B
done
} > $OUT 2>&1
echo done
