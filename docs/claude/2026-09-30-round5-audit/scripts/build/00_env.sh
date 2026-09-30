#!/bin/bash
# Probe the WSL toolchain and the prebuilt state for the build/packaging audit.
set -u
echo "== tools"; for t in cmake ninja g++ gcc just python3 zstd rsync nproc shellcheck clang++ ld.gold ld.lld dash; do printf '%s: ' $t; command -v $t || echo none; done
cmake --version | head -1; g++ --version | head -1; ninja --version 2>/dev/null
echo "== static libs"; ls -la /usr/lib/x86_64-linux-gnu/{libzstd.a,libdeflate.a,libc.a,libgtest.a} 2>&1; ls /usr/lib/gcc/x86_64-linux-gnu/*/libstdc++.a /usr/lib/gcc/x86_64-linux-gnu/*/libgomp.a 2>&1
echo "== gtest"; ls /usr/include/gtest 2>&1 | head -3; ls /usr/lib/x86_64-linux-gnu/cmake/GTest 2>&1
echo "== cpu"; grep -m1 '^flags' /proc/cpuinfo | tr ' ' '\n' | grep -E '^(avx|avx2|bmi1|bmi2|fma|f16c|movbe|abm|lzcnt|avx512f)$' | tr '\n' ' '; echo; grep -m1 'model name' /proc/cpuinfo
echo "== strain-build"; ls -la ~/strain-build ~/strain-build/bin ~/strain-build/mini_db ~/strain-build/mini_db/protal_db 2>&1 | head -60
echo "== src head"; (cd ~/strain-build/src && git log --oneline -1 2>&1; git status --short 2>&1 | head)
echo "== load"; uptime; free -g; nproc
echo "== audit6"; ls -la ~/audit6 ~/audit6/build 2>&1
