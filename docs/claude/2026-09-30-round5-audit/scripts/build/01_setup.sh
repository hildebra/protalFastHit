#!/bin/bash
# Copy the source (commit 39a8585) into ~/audit6/build/src and check it matches the worktree.
set -u
A=~/audit6/build
mkdir -p $A
rsync -a --delete --exclude /build ~/strain-build/src/ $A/src/
WT=/mnt/c/Users/hildebra/Documents/locDev/protal/.claude/worktrees/strain-fixes
echo "== diff vs worktree (excluding .git .idea .vscode .claude build dirs)"
diff -rq --exclude=.git --exclude=.idea --exclude=.vscode --exclude=.claude --exclude=build --exclude='build-*' --exclude=__pycache__ --exclude=strain_test_out --exclude=work $WT $A/src 2>&1 | head -30
echo "diff rc=${PIPESTATUS[0]}"
ls $A/src
ls -la $A/src/src/protal_config.h 2>&1
