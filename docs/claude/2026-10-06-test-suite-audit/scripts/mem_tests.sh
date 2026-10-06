#!/usr/bin/env bash
# Peak RSS and wall time of the slowest unit tests (one process each).
cd ~/testaudit/src
OUT=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/57f6bd73-f81f-460d-bd11-1b9ef875fd6f/scratchpad/results/mem_tests.tsv
printf 'test\twall_s\tmax_rss_mb\n' > "$OUT"
for t in PackedIndex.ABuildIntoThePackedLayoutWritesTheSameIndex \
         PackedIndex.TaxidAndGeneAsOneNumberComeBackForEveryGeneCount \
         PackedIndex.TheColumnFormatPacksChunkByChunkToTheSameValues \
         PackedIndex.AValueOutsideTheLayoutStopsTheColumnLoad \
         PackedIndex.HoldsEveryFlexCellAndEntryOfTheStoredLayout \
         ProfileSam.AZstdSamOfManyFramesIsDecompressedOnThreadsOfItsOwn \
         IndexCodec.ColumnsCompressBetterThanRawCells \
         ProfileSam.OnSeveralThreadsTheProfileIsTheSame; do
    /usr/bin/time -f "%e\t%M" -o /tmp/ta_time.txt nice -n 10 build/tests/protal_tests --gtest_filter="$t" > /dev/null 2>&1
    printf '%s\t%s\n' "$t" "$(awk -F'\t' '{printf "%s\t%.0f", $1, $2/1024}' /tmp/ta_time.txt)" >> "$OUT"
done
free -g | head -2 >> "$OUT"
cat "$OUT"
