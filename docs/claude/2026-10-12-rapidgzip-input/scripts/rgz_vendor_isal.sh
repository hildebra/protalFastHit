#!/bin/bash
# lib/isa-l: ISA-L 2.32.1 (the patched tree in ~/rgz/isal) without doc/, examples/, programs/, tests/, tools/.
set -euo pipefail
P=/mnt/c/Users/hildebra/Documents/locDev/protal
cd ~/rgz/isal
rm -f include/igzip_lib.h.rej
test ! -e $P/lib/isa-l
mkdir $P/lib/isa-l
tar -cf - --exclude=./.git --exclude=./doc --exclude=./examples --exclude=./programs --exclude=./tests --exclude=./tools . | tar -xf - -C $P/lib/isa-l
ls $P/lib/isa-l
du -sh $P/lib/isa-l
find $P/lib/isa-l -type f | wc -l
# the upstream tree, for diffs later
rm -rf ~/rgz/isal-upstream
mkdir ~/rgz/isal-upstream
tar -xzf ~/isal/isa-l-2.32.1.tar.gz -C ~/rgz/isal-upstream --strip-components=1
