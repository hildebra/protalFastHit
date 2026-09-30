#!/bin/bash
# callgrind of --profile_only on one SAM (1 thread): parsing against the rest of profiling.
#   cg_profile.sh DB SAM
source "$(dirname "$0")/env.sh"
db=$1; sam=$2
rm -rf $OUT/po; mkdir -p $OUT/po
valgrind --tool=callgrind --callgrind-out-file=$OUT/cg_prof.out $BIN --db $db --profile_only $sam --prefix p -o $OUT/po -t 1 --no_qcmsa > /dev/null 2>&1
callgrind_annotate --inclusive=yes --threshold=99.5 $OUT/cg_prof.out 2>/dev/null > $OUT/cg_prof.txt
grep "PROGRAM TOTALS" $OUT/cg_prof.txt
echo "records: $(grep -vc '^@' $sam)"
grep -E "ProfileWrapper|ProfileSam|SamReader::Next|LineSplitter::Split|strtol|strtoul|SamEntry::SamEntry|MicrobialProfile::AddSam|Taxon::AddSam|WriteSparseProfile|CheckReference|LoadModel|LoadAllGenomes|ProtalDB::ProtalDB" \
  $OUT/cg_prof.txt | cg_lines | awk '!seen[$2]++'
