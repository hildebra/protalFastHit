set -e
cd ~/protal-mem
r=~/protal-mem/reads/w900s
m=~/protal-mem/cohort_small.map
{ printf '#OUTPUT_DIR\t%s/runs/massif2/out\n#INPUT_DIR\t%s\n#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tPROFILE\n' $HOME/protal-mem $r
  for i in 1 2 3 4; do printf 's%s\tw900s_R1.fq.gz\tw900s_R2.fq.gz\ts%s.sam.zst\ts%s\ts%s.profile\n' $i $i $i $i; done; } > $m
rm -rf runs/massif2; mkdir -p runs/massif2
valgrind --tool=massif --massif-out-file=runs/massif2/massif.out --threshold=0.5 --depth=16 --detailed-freq=1 --max-snapshots=200 --time-unit=B \
  ./build/protal --db ~/protal-perf/db900 --map $m -t 1 --no_qcmsa > runs/massif2/stdout.log 2> runs/massif2/stderr.log
ms_print runs/massif2/massif.out > runs/massif2/ms_print.txt
echo done
