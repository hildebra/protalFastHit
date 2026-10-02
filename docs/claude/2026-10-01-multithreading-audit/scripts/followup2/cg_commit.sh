#!/usr/bin/env bash
# Callgrind, one thread, --profile_only of the 500k-pair SAM: f359d49 (fbase2) and f359d49 + 85e53e6 + a212559
# (final2, the committed src). Summaries into ~/mt-work/cgc/summary_<name>.txt.
set -uo pipefail
SP=$(cd "$(dirname "$0")" && pwd)  # this folder
W=$HOME/mt-work/cgc; mkdir -p $W
SAM=$HOME/bench071/runs/v071.full.pe.rl150_p500000_s_1/rl150_p500000_s_1.sam.zst
for name in fbase2 final2; do
  cp $HOME/mt-work/$name/src/build/protal $W/protal-$name
  rm -rf $W/out_$name
  nice valgrind --tool=callgrind --callgrind-out-file=$W/cg_$name.out $W/protal-$name --db $HOME/bench071/V071/protal_db \
    --profile_only $SAM --prefix p -o $W/out_$name -t 1 --no_qcmsa > $W/log_$name 2> $W/err_$name
  callgrind_annotate --inclusive=yes $W/cg_$name.out 2>/dev/null | grep -E "PROGRAM TOTALS|SamReader::Advance|PrepareMAPQ\(|CollectGroups|ZSTD_decompressStream|Taxon::AddSam|WriteSparseProfile\(|CoveredPortion|NormalizeCigar|SamFromTokens|LineSplitter::Split" \
    | sed -E 's/\(protal::profiler::MicrobialProfile const&.*//; s/<protal::profiler::Profiler<protal::ScoreAlignments>::ProfileSam.*//; s/\[\/home.*//' | head -16 > $W/summary_$name.txt
  echo "== $name"; cat $W/summary_$name.txt
done
diff -r -q $W/out_fbase2 $W/out_final2 > /dev/null && echo "CG SAME outputs" || echo "CG DIFFER outputs"
echo CG DONE
