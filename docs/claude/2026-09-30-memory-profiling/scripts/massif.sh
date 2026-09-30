set -e
cd ~/protal-mem
r=~/protal-perf/reads/w900/reads
mkdir -p reads/w900s
for m in 1 2; do zcat $r/w900_1_R$m.fq.gz | head -n 400000 | gzip -1 > reads/w900s/w900s_R$m.fq.gz; done
rm -rf runs/massif; mkdir -p runs/massif
valgrind --tool=massif --massif-out-file=runs/massif/massif.out --threshold=0.3 --depth=14 --detailed-freq=1000000 --time-unit=B \
  ./build/protal --db ~/protal-perf/db900 -1 reads/w900s/w900s_R1.fq.gz -2 reads/w900s/w900s_R2.fq.gz -o runs/massif/out -t 1 --no_qcmsa > runs/massif/stdout.log 2> runs/massif/stderr.log
ms_print runs/massif/massif.out > runs/massif/ms_print.txt
echo done
