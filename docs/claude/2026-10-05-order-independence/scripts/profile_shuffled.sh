#!/usr/bin/env bash
# The profiling stage alone on a SAM and on the same SAM with its reads in another order (each read's records kept
# together and in their order, the reads shuffled with a fixed seed): build <name>'s --profile_only on both, at 6 threads.
#   profile_shuffled.sh <name> <run dir with s.sam.zst> <pe|pb> <out dir>
set -uo pipefail
NAME=$1; RUN=$2; TYPE=$3; O=$4
W=$HOME/det-order; BIN=$W/$NAME/tree/build/protal; DB=$HOME/bench071/V073/protal_db/database.protal
mkdir -p $O
zstd -dcq $RUN/s.sam.zst > $O/in.sam
grep '^@' $O/in.sam > $O/shuffled.sam
grep -v '^@' $O/in.sam | awk -F'\t' '$1 != last { if (NR > 1) printf "\n"; last = $1; printf "%s", $0; next } { printf "\001%s", $0 } END { printf "\n" }' \
  | shuf --random-source=<(yes 2026) | tr '\001' '\n' >> $O/shuffled.sam
echo "records: $(grep -vc '^@' $O/in.sam) in, $(grep -vc '^@' $O/shuffled.sam) shuffled; same set: $(cmp -s <(grep -v '^@' $O/in.sam | LC_ALL=C sort) <(grep -v '^@' $O/shuffled.sam | LC_ALL=C sort) && echo yes || echo NO)"
for v in in shuffled; do
  rm -rf $O/p.$v
  $BIN --db $DB --profile_only $O/$v.sam --read_type $TYPE --prefix s -o $O/p.$v -t 6 --no_qcmsa > $O/p.$v.log 2>&1 || echo "profile $v failed"
done
d=$(diff -r -x '*_runtime.tsv' -x '*.statistics.tsv' $O/p.in $O/p.shuffled)
[ -z "$d" ] && echo "$NAME $TYPE: profile of the shuffled SAM same" || { echo "$NAME $TYPE: profile of the shuffled SAM DIFFERS"; echo "$d" | head -${LINES_SHOWN:-20}; }
rm -f $O/in.sam $O/shuffled.sam
