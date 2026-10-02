#!/usr/bin/env bash
# The SAM headers' sizes: @SQ lines, plain bytes and zstd -3 bytes, for the full database's SAMs of the 0.7.1
# benchmark and the 5M-pair run; and the full database's gene count.
for f in $HOME/mt-audit/runs/b6/s1.sam.zst $HOME/bench071/runs/v071.full.*/*.sam.zst; do
  h=$(zstdcat "$f" 2>/dev/null | awk '/^@/ { print; next } { exit }')
  printf "%s\tSQ\t%s\tplain\t%s\tzstd3\t%s\tsam\t%s\n" "$(basename $(dirname $f))/$(basename $f)" "$(printf '%s\n' "$h" | grep -c '^@SQ')" \
    "$(printf '%s\n' "$h" | wc -c)" "$(printf '%s\n' "$h" | zstd -3 -c | wc -c)" "$(stat -c %s "$f")"
done
echo "db genes: $(wc -l < $HOME/bench071/V071/protal_db/gene2geneid.tsv 2>/dev/null || echo '?')"
ls $HOME/bench071/V071/protal_db | head
