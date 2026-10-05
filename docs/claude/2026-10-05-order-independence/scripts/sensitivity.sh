#!/usr/bin/env bash
# Does the new unit test notice order dependence? The change's tree with ExactSum summing in a plain double (the sums as
# before the change), tests built and run: ProfileSam.TheProfileDoesNotDependOnTheOrderOfTheReads should fail.
set -uo pipefail
W=$HOME/det-order; E=$W/sensitivity; mkdir -p $E
rsync -a --delete $W/work/tree/ $E/tree/ --exclude '/build/'
perl -0pi -e 's/__int128 m_units = 0;/double m_units = 0;/; s/m_units \+= static_cast<__int128>\(value \* kUnit\);/m_units += value * kUnit;/' $E/tree/src/Utilities/ExactSum.h
grep -n 'double m_units' $E/tree/src/Utilities/ExactSum.h || exit 1
cd $E/tree && nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > $W/sensitivity.configure.log 2>&1 && \
  nice cmake --build build --target protal_tests -j 6 > $W/sensitivity.build.log 2>&1 || { echo FAIL build; exit 1; }
$E/tree/build/tests/protal_tests --gtest_filter='ProfileSam.*:ExactSum.*' 2>&1 | grep -E '^\[ +(OK|FAILED) +\]|PASSED|FAILED'
