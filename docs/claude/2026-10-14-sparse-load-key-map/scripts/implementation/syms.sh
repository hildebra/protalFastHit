#!/bin/bash
# The symbols of ChainAnchorFinder's lookup functions and operator() in a protal binary, demangled.
nm -C ${1:-~/keysfix/new/build/protal} | grep -E "ChainAnchorFinder<protal::KmerLookupSM>::(LookupKmers|SparseHits|SparsePrepare|operator\(\))|KmerLookupSM::ScoreAtLevel|CompactKeys::CountRanks" | cut -c1-40,80-260 | sed 's/std::vector<std::pair<unsigned long, unsigned long>, std::allocator<std::pair<unsigned long, unsigned long> > >/KmerList/g' | cut -c1-200
