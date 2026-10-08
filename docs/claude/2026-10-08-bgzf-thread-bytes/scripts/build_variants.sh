#!/bin/bash
# In ~/bgzfdet (401c4f5): bin/tests_newtest_oldcode = the new test file with the old Bgzf.h; bin/tests_fixed = both
# changed files. Built on cores 0-3.
set -u
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=~/bgzfdet
cd $W
cp $REPO/tests/test_RunStatusAndBgzf.cpp src/tests/test_RunStatusAndBgzf.cpp
taskset -c 0-3 nice -n 5 ninja -C build -j4 protal_tests > ninja_newtest_oldcode.log 2>&1
rc=$?
echo "newtest_oldcode ninja rc $rc $(date +%T)"
[ $rc -eq 0 ] && cp build/tests/protal_tests bin/tests_newtest_oldcode || tail -30 ninja_newtest_oldcode.log
cp $REPO/src/IO/Bgzf.h src/src/IO/Bgzf.h
taskset -c 0-3 nice -n 5 ninja -C build -j4 protal_tests > ninja_fixed.log 2>&1
rc=$?
echo "fixed ninja rc $rc $(date +%T)"
[ $rc -eq 0 ] && cp build/tests/protal_tests bin/tests_fixed || tail -30 ninja_fixed.log
grep -c "Building CXX" ninja_fixed.log
grep -i -E "warning" ninja_fixed.log | grep -v lto | head
cmp -s $REPO/src/IO/Bgzf.h src/src/IO/Bgzf.h && cmp -s $REPO/tests/test_RunStatusAndBgzf.cpp src/tests/test_RunStatusAndBgzf.cpp && echo "tree = working files"
ls -la bin
