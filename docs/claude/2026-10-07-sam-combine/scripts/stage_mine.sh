#!/usr/bin/env bash
# stage_mine.sh STAGE FILE[:DROP,DROP...]... - copies this session's version of each FILE into STAGE (same relative path):
# the checkout's file, or for FILE:DROPS the committed file (HEAD) with the checkout's changes except the hunks starting
# at the old lines DROPS (another session's uncommitted hunks in a file this session changes too), applied with context
# (docs/claude/2026-10-07-sam-combine). Run from the checkout, in Git Bash.
set -eu
STAGE=$1; shift
HERE=$(cd "$(dirname "$0")" && pwd)
for spec in "$@"; do
  file=${spec%%:*}
  mkdir -p "$STAGE/$(dirname "$file")"
  if [ "$spec" = "$file" ]; then
    tr -d '\r' < "$file" > "$STAGE/$file"
    continue
  fi
  drops=${spec#*:}
  index=$(mktemp)
  git diff -U3 HEAD -- "$file" | awk -v drop="${drops//,/ }" -f "$HERE/drop_hunks.awk" > "$index.patch"
  GIT_INDEX_FILE=$index git read-tree HEAD
  GIT_INDEX_FILE=$index git apply --cached "$index.patch"
  git cat-file -p "$(GIT_INDEX_FILE=$index git ls-files -s "$file" | awk '{print $2}')" > "$STAGE/$file"
  rm -f "$index" "$index.patch"
  echo "$file: HEAD + $(grep -c '^@@' <(git diff -U3 HEAD -- "$file") ) hunks less ${drops}"
done
