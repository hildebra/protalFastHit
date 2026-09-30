#!/bin/bash
# Experiment 8: the default --x_drop 1000 against --x_drop 0 on the e1 community reads (HiFi, ONT)
# and on ONT reads at 5% errors: wrong-species records whose own species has the gene, genes found,
# failed alignments (debug build), and time.
S=$(dirname "$0")
W=~/audit6/longreads
P=~/strain-build/bin/protal
PA=$W/build-asan/protal
G=~/strain-build/mini_db/gtdb_r226
D=$W/mini_db
OUT=$S/out_e8.txt
mkdir -p $W/e8 && cd $W/e8
{
[ -f ont5.fq ] || python3 $S/gen_reads.py community ont5 --gtdb $G --db $D --reads 300 --platform ont5 --seed 21 --median 20000 --sigma 0.6 --name "ont5_{}"
for s in hifi:pb ont:ont ont5:ont; do
  n=${s%%:*}; t=${s#*:}
  fq=$W/e1/$n.fq; tp=$W/e1/$n.truth.pkl
  [ -f $fq ] || { fq=$W/e8/$n.fq; tp=$W/e8/$n.truth.pkl; }
  for xd in 1000 0; do
    rm -rf out_xd$xd/$n.*
    /usr/bin/time -f "%e s" $P --db $D -1 $fq --read_type $t --model_$t $D/model_pe.xml -o out_xd$xd --prefix $n -t 2 --sam_format sam --no_qcmsa --no_strains --x_drop $xd > $n.xd$xd.log 2> $n.xd$xd.err
    echo "######## $n --x_drop $xd rc=$? time $(tail -1 $n.xd$xd.err) invalid=$(grep -c 'Invalid after alignment' $n.xd$xd.err) logged $(grep 'x-drop' $n.xd$xd.log)"
    python3 $S/eval_sam.py out_xd$xd/$n.sam $fq $tp --gtdb $G --db $D --exact 0.95 --show 0 | grep -E "^(full_|partial|records|wrong|primary|supplementary|ZR|mapq)" | tr '\n' ' ' | sed 's/mapq (correct/\n    mapq (correct/'
    echo
    grep -A0 -E "abundance|^GCF" out_xd$xd/$n.profile | awk -F'\t' '{printf "    %s %s\n", $1, $3}'
  done
done
echo "== AlignAnchor failures in the debug build (ONT e1 reads)"
for xd in 1000 0; do
  PROTAL_LR_DEBUG=1 ASAN_OPTIONS=detect_leaks=0 $PA --db $D -1 $W/e1/ont.fq --read_type ont --model_ont $D/model_pe.xml -o dbg_xd$xd --prefix ont -t 2 --sam_format sam --no_qcmsa --no_strains --no_profile --x_drop $xd > dbg$xd.log 2>&1
  echo "--x_drop $xd rc=$? AlignAnchor failures: $(grep -c 'align-fail AlignAnchor' dbg$xd.log), proxy-ani failures: $(grep -c 'align-fail proxy-ani' dbg$xd.log), link-outside: $(grep -c 'align-fail link-outside' dbg$xd.log), invalid: $(grep -c '^Invalid alignment' dbg$xd.log), candidates: $(grep -c '^LRDBG cand' dbg$xd.log)"
  grep -E "runtime error|AddressSanitizer" dbg$xd.log | sort | uniq -c | head -5
done
} > $OUT 2>&1
echo done
