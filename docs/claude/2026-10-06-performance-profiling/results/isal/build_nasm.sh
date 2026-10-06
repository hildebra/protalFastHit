#!/bin/bash
# Fresh nasm 2.16.03 build from the release tarball, timed.
set -e
cd ~/isal
rm -rf nasm-build && mkdir nasm-build && tar xf nasm-2.16.03.tar.xz -C nasm-build
cd nasm-build/nasm-2.16.03
t0=$(date +%s.%N)
./configure --prefix=$HOME/isal/prefix > ../configure.log 2>&1
t1=$(date +%s.%N)
make -j6 > ../make.log 2>&1
t2=$(date +%s.%N)
make install > ../install.log 2>&1
t3=$(date +%s.%N)
awk -v a=$t0 -v b=$t1 -v c=$t2 -v d=$t3 'BEGIN{printf "configure %.1f s, make -j6 %.1f s, install %.1f s, total %.1f s\n", b-a, c-b, d-c, d-a}'
grep -iE 'perl|asciidoc|xmlto|autoconf|automake' ../configure.log || true
"$HOME/isal/prefix/bin/nasm" -v
