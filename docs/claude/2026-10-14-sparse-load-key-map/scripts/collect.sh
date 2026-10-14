#!/bin/bash
# Copy the small result files of ~/sparseprof/runs into the report's results folder.
R=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-14-sparse-load-key-map/results
mkdir -p $R
cd ~/sparseprof/runs
cp load_base/loads.tsv $R/load_base.tsv
cp load_ab/loads.tsv $R/load_ab1.tsv
cp load_ab2/loads.tsv $R/load_ab2.tsv
grep -v "^Load\|^Read\|^Pre\|^ " load_ab2/phases.txt > $R/load_ab2_phases.txt
for d in keys_v1 keys2 keys3; do grep -E "^==|key maps|bucket lines|huge pages|same slots|reads,|^(D|D2|C|Cp|Cf|Ca|Cl|M|Ma|B)_|synthetic" $d/log.txt > $R/$d.txt; done
cd cg2
for f in sparse2c_old1 sparse2c_old0 sparse2d_old1 sparse2d_old0 sparse2d_fill wholec_old1 wholec_old0 wholed_old0; do
  echo "#### $f"; callgrind_annotate $f.out 2>/dev/null | sed -n '/^Ir/,$p' | grep -v "^-\|^$\|file:function\|^Ir" | head -14 | cut -c1-200
done > $R/callgrind_self.txt
bash /mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/cf462e00-4c32-4801-b980-f96668bbdc3a/scratchpad/summ.sh base > $R/load_base_medians.tsv
bash /mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/cf462e00-4c32-4801-b980-f96668bbdc3a/scratchpad/summ_ab.sh > $R/load_ab1_medians.tsv
bash /mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/cf462e00-4c32-4801-b980-f96668bbdc3a/scratchpad/summ_ab2.sh > $R/load_ab2_medians.tsv
bash /mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/cf462e00-4c32-4801-b980-f96668bbdc3a/scratchpad/popcnt_callers.sh > $R/popcountdi2_callers.txt
ls -la $R; du -sh $R
