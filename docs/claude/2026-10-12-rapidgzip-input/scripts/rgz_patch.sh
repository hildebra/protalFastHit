#!/bin/bash
# Try the fork's commits on upstream ISA-L 2.32.1: which apply, which conflict.
set -uo pipefail
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad
W=~/rgz
cd $W
rm -rf split isal
mkdir split
cd split
git mailsplit -o. $SP/rapidgzip-dl/patches/fork.patch >/dev/null
ls
cd $W
mkdir isal
tar -xzf ~/isal/isa-l-2.32.1.tar.gz -C isal --strip-components=1
cd isal
git init -q
git add -A
git -c user.email=x@y -c user.name=x commit -qm upstream-2.32.1
for n in 0001 0002 0004 0005 0006 0007; do
  f=$W/split/$n
  subj=$(grep -m1 '^Subject:' $f)
  echo "== $n $subj"
  git diff --stat=200 HEAD >/dev/null
  if git apply --check $f 2>/tmp/err_$n; then
    git apply $f && git -c user.email=x@y -c user.name=x commit -qam "$n" && echo "applied cleanly"
  else
    echo "conflicts:"; head -20 /tmp/err_$n
    # files it touches
    grep '^diff --git' $f | sed 's/^diff --git a\///; s/ b\/.*//'
  fi
done
git log --oneline | head
