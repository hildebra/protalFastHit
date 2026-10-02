#!/usr/bin/env bash
# Results of the performance round 3 into its report folder.
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad
R=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-02-performance-round3/results
W=$HOME/mt-work/perf3
# instructions: totals and the alignment parts, per build and workload
{ for n in pe100k.lr4 pe100k.id ont_b3000000.ref ont_b3000000.lr ont_b3000000.reseed ont_b3000000.noreseed2 ont_b3000000.reseed2 ont_b3000000.noreseed3 ont_b3000000.reseed3 ont_b3000000.lr4 \
            pb_b3000000.ref pb_b3000000.lr pb_b3000000.reseed pb_b3000000.noreseed2 pb_b3000000.reseed2 pb_b3000000.noreseed3 pb_b3000000.reseed3 pb_b3000000.lr4 hifi3M.ref hifi3M.lr4; do
    f=$W/cg.$n.incl.txt; [ -f $f ] || f=$W/cg.$n.reseed.incl.txt; [ -f $f ] || { callgrind_annotate --inclusive=yes --threshold=99.5 $W/cg.$n 2>/dev/null | c++filt > $W/cg.$n.incl.txt; f=$W/cg.$n.incl.txt; }
    echo "== $n: $(grep 'PROGRAM TOTALS' $f | awk '{print $1}')"
    grep -E 'RunLongReads|RunPairedEnd|AnchoredAligner::(Align|Flank|Gap|Reseed|Piece)|Between|alignEndsFree|alignEnd2End|slab_reap_repurpose|DecodeChunk|Seedmap::Load ' $f | sed -E 's/\[\/home[^]]*\]//; s/\(([^()]|\([^()]*\))*\)//g' | awk '{ $1 = $1; print }' | cut -c1-120 | awk '!seen[$2]++' | head -12
  done; } > $R/callgrind_parts.txt
cp $HOME/mt-work/onebcg/pe100k.cl2.txt $R/callgrind_pe100k_head_self.txt 2>/dev/null; head -80 $R/callgrind_pe100k_head_self.txt > $R/tmp && mv $R/tmp $R/callgrind_pe100k_head_self.txt
bash $S/perf3/top.sh ont3M cl2 40 > $R/callgrind_ont3M_head.txt 2>&1
bash $S/perf3/calls_all.sh > $R/wfa2_calls_and_costs_ont3M.txt 2>&1
bash $S/perf3/stages.sh > $R/stage_times_1thread.txt 2>&1
cat $S/perf3/lr_run.out $S/perf3/lr_run2.out $S/perf3/lr_run3.out $S/perf3/lr_run4.out 2>/dev/null | grep -E "^==|^records|^common|same SAM|SAM text|sorted SAM|outputs" > $R/long_read_records.txt
grep -h -E "links per|longer flank|anchors:|read\+gene" $W/../perf3/*.out 2>/dev/null > /dev/null
cp $W/lrbench/times.tsv $R/long_read_benchmark_times.tsv
perl $S/perf3/lr_score.pl > $R/long_read_benchmark.txt; perl $S/perf3/f1_v072.pl >> $R/long_read_benchmark.txt
cp $S/perf3/id_check.out $R/identical_changes_checks.txt
ls -la $R
