#!/usr/bin/env bash
# unquoted_glob.sh - what protal does when the shell expands an unquoted --profile_only pattern into several arguments
# (docs/claude/2026-10-07-sam-combine). Needs run.sh's work folder.
set -u
PROTAL=${PROTAL:-$HOME/samcombine/bin/protal}
DB=${DB:-$HOME/samcombine/db/database.protal}
W=${W:-$HOME/samcombine/work}
cd "$W" || exit 1
rm -rf unquoted
taskset -c 0-3 nice -n 5 "$PROTAL" --db "$DB" -t 4 --no_qcmsa -o unquoted --profile_only run1/alignments/*.sam.zst run2/alignments/*.sam.zst \
  > unquoted.log 2>&1
echo "exit $?; profiles: $(ls unquoted/*.profile 2>/dev/null | xargs -n1 basename | tr '\n' ' ')"
grep -i "warning\|error\|unmatched" unquoted.log | head -5
