#!/bin/bash
# Configure-level experiments on a scratch copy of the source (never the repository):
# stale-header deletion, ignored -DCMAKE_CXX_FLAGS_RELEASE, $GCC-derived ar, OpenMP off,
# include order into cPMML and a simulated case-insensitive options.h/Options.h clash,
# the -L directories on the link line.
set -u
A=~/audit6/build
E=$A/src_exp
rm -rf $E; rsync -a --exclude '/build*' $A/src/ $E/
cd $E
echo "== (a) a stale src/protal_config.h is deleted from the source tree at configure time"
echo '// stale' > src/protal_config.h
cmake -S . -B b1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS_RELEASE="-O1 -g" > $A/exp_b1.log 2>&1; echo "rc=$?"
ls -la src/protal_config.h 2>&1; grep -i 'stale' $A/exp_b1.log
echo "== (b) -DCMAKE_CXX_FLAGS_RELEASE='-O1 -g' on the command line vs what targets get"
grep -E '^CXX_FLAGS' b1/CMakeFiles/protal.dir/flags.make; grep CMAKE_CXX_FLAGS_RELEASE b1/CMakeCache.txt
echo "== (f) include order for cPMML (protal's src/ first?)"
grep -E '^CXX_INCLUDES' b1/src/cPMML/CMakeFiles/cPMML.dir/flags.make | tr ' ' '\n' | grep -n '^-I' | head -8
grep -rn '#include "options.h"' lib/cPMML/src lib/cPMML/include | head -3
echo "== (g) -L directories on protal's link line, and whether they exist"
for d in $(tr ' ' '\n' < b1/CMakeFiles/protal.dir/link.txt | grep '^-L' | sed 's/^-L//' | sort -u); do [ -d "$d" ] && echo "exists  $d" || echo "MISSING $d"; done
echo "== (e) simulated case-insensitive file system: src/options.h resolving to src/Options.h"
ln -s Options.h src/options.h
cmake --build b1 --target cPMML -j 2 > $A/exp_clash.log 2>&1; echo "cPMML build rc=$?"; grep -m3 -E 'error' $A/exp_clash.log
rm src/options.h
echo "== (c) GCC set in the environment to a compiler name whose '-ar' does not exist"
ls /usr/bin/ | grep -E 'gcc.*-ar|-gcc-ar' | head
GCC=gcc-13 cmake -S . -B b2 -DCMAKE_BUILD_TYPE=Release > $A/exp_b2.log 2>&1; echo "configure rc=$?"; grep -E 'GCC-prefixed|Ranlib' $A/exp_b2.log
cmake --build b2 --target gzstream_lib -j 2 > $A/exp_b2_build.log 2>&1; echo "gzstream_lib build rc=$?"; grep -m3 -iE 'error|not found|No such' $A/exp_b2_build.log
echo "== (d) OpenMP disabled: which directories still get -fopenmp"
cmake -S . -B b3 -DCMAKE_BUILD_TYPE=Release -DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=ON > $A/exp_b3.log 2>&1; echo "configure rc=$?"
for f in b3/CMakeFiles/protal.dir/flags.make b3/src/gzstream/CMakeFiles/gzstream_lib.dir/flags.make b3/src/cPMML/CMakeFiles/cPMML.dir/flags.make; do echo "-- $f"; grep -E '^CXX_FLAGS' $f; done
echo "== default configure output noise (message() lines without STATUS)"
grep -vE '^-- ' $A/rel_configure.log | head -20
echo DONE
