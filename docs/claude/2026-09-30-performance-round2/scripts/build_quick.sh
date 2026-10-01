#!/bin/bash
# Experiment build of the cheap fixes the profile points at, applied to a copy of the tree:
#   1. KmerUtils::ReverseComplement takes a const reference and fills a pre-sized string through a 256-entry
#      table (was: by value, appended one character at a time through a switch);
#   2. SimpleAlignmentHandler copies neither the read (fwd) nor, per anchor, the reverse complement (rev);
#   3. the per-read stage timers are compiled out (Benchmark::Start/Stop empty).
# The repository is not changed.
set -e
source "$(dirname "$0")/env.sh"
rm -rf $PERF_DIR/src-q; mkdir -p $PERF_DIR/src-q
rsync -a --exclude '/data' --exclude '/build*/' --exclude '/.git' $PERF_DIR/src/ $PERF_DIR/src-q/
cd $PERF_DIR/src-q/src
perl -0pi -e 's/    inline std::string ReverseComplement\(std::string forward\) \{\n        std::string reverse = "";\n        const char \* seq = forward.c_str\(\);\n        for \(int i = forward.length\(\)-1; i >= 0; i--\) \{\n            reverse \+= Complement\(forward\[i\]\);\n        \}\n        return reverse;\n    \}/    inline std::string ReverseComplement(std::string const& forward) {\n        static const std::array<char, 256> table = [] { std::array<char, 256> t{}; for (int i = 0; i < 256; i++) t[i] = Complement(static_cast<char>(i)); return t; }();\n        size_t const n = forward.size();\n        std::string reverse(n, \x27N\x27);\n        for (size_t i = 0; i < n; i++) reverse[i] = table[static_cast<unsigned char>(forward[n - 1 - i])];\n        return reverse;\n    }/' SequenceUtils/KmerUtils.h
grep -n "ReverseComplement(std::string" SequenceUtils/KmerUtils.h
perl -0pi -e 's/std::string& fwd, std::string rev, bool allow_heuristic_alignment/std::string const& fwd, std::string const& rev, bool allow_heuristic_alignment/; s/            auto fwd = sequence;\n            auto rev = KmerUtils::ReverseComplement\(fwd\);/            std::string const& fwd = sequence;\n            auto const rev = KmerUtils::ReverseComplement(fwd);/' Core/AlignmentStrategy.h
grep -n "std::string const& fwd" Core/AlignmentStrategy.h
perl -0pi -e 's/void Start\(bool new_sample=true\) \{\s*samplings \+= new_sample;\s*start_time = std::chrono::steady_clock::now\(\);\s*\}/void Start(bool new_sample=true) { samplings += new_sample; }/; s/void Stop\(\) \{\s*auto const stop_time = std::chrono::steady_clock::now\(\);\s*time_sum \+= duration_cast<std::chrono::nanoseconds>\(stop_time - start_time\)\.count\(\);\s*\}/void Stop() {}/' Utilities/Benchmark.h
grep -n "#include <array>" SequenceUtils/KmerUtils.h || sed -i '0,/#include/s//#include <array>\n#include/' SequenceUtils/KmerUtils.h
cd $PERF_DIR
cmake -S src-q -B build-q -G Ninja -DCMAKE_BUILD_TYPE=Release > build-q.cmake.log 2>&1
cmake --build build-q --target protal_avx2 -j "$(nproc)" > build-q.log 2>&1
ls -la build-q/protal_avx2
