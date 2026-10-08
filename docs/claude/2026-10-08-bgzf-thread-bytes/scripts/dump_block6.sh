#!/bin/bash
# Dumps block 6 of both s2 files (deflate payload: after the 18-byte header, before the 8-byte footer) and diffs them.
set -u
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/99a2c9fe-d722-41e1-bc32-3fc6baabdeb2/scratchpad
cd ~/lrdet/fail39
perl $SP/deflate_dump.pl one_s2.fq.gz $((141815 + 18)) $((1905 - 26)) > one_b6.txt
perl $SP/deflate_dump.pl four_s2.fq.gz $((141815 + 18)) $((1901 - 26)) > four_b6.txt
wc -l one_b6.txt four_b6.txt
tail -1 one_b6.txt; tail -1 four_b6.txt
grep -c " M " one_b6.txt four_b6.txt
diff <(head -6 one_b6.txt) <(head -6 four_b6.txt) | cut -c1-400
echo "--- token diff"
diff one_b6.txt four_b6.txt | grep -v -E "^[<>] (LL|DL|CL|HLIT)" | head -60
