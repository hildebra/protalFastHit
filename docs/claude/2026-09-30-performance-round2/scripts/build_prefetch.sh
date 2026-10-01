#!/bin/bash
# Experiment build: the current source plus software prefetching of the key map (all of a read's
# control blocks first, then the value blocks of the smallest lookups), to see what the seeding
# stage's memory latency costs. Applied to a copy of the tree; the repository is not changed.
set -e
source "$(dirname "$0")/env.sh"
rm -rf $PERF_DIR/src-pf; mkdir -p $PERF_DIR/src-pf
rsync -a --exclude '/data' --exclude '/build*/' --exclude '/.git' $PERF_DIR/src/ $PERF_DIR/src-pf/
cd $PERF_DIR/src-pf/src
perl -0pi -e 's/(        uint64_t FlexKey\(uint64_t key\) \{)/        void PrefetchKey(uint64_t key) {\n            uint64_t const b = ControlBlockIndex(MainKey(key));\n            __builtin_prefetch(m_keymap + b);\n            __builtin_prefetch(m_keymap + b + 12);\n        }\n\n$1/' Hash/Seedmap.h
perl -0pi -e 's/(        inline void Get\(std::vector<LookupPointer>& result, size_t &kmer, uint32_t readpos\) \{)/        inline void Prefetch(size_t kmer) { m_sm.PrefetchKey(kmer); }\n\n$1/' Hash/KmerLookup.h
perl -0pi -e 's/(            m_lookups.clear\(\);\n\n            \/\/ Get ranges in values \(no seeds yet\)\n)/            m_lookups.clear();\n            for (auto [mmer, pos] : kmer_list) m_kmer_lookup.Prefetch(mmer);\n\n            \/\/ Get ranges in values (no seeds yet)\n/; s/(            uint32_t previous_size = 0;\n            m_successful_lookups = 0;)/            for (size_t i = 0, n = std::min<size_t>(m_lookups.size(), 12); i < n; i++) {\n                __builtin_prefetch(m_lookups[i].values_begin);\n                if (m_lookups[i].flex_begin) __builtin_prefetch(m_lookups[i].flex_begin);\n            }\n\n$1/' Core/ChainAnchorFinder.h
grep -n "Prefetch\|prefetch" Hash/Seedmap.h Hash/KmerLookup.h Core/ChainAnchorFinder.h
cd $PERF_DIR
cmake -S src-pf -B build-pf -G Ninja -DCMAKE_BUILD_TYPE=Release > build-pf.cmake.log 2>&1
cmake --build build-pf --target protal_avx2 -j "$(nproc)" > build-pf.log 2>&1
ls -la build-pf/protal_avx2
