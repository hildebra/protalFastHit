#pragma once

// A SamSink for the unit tests of the alignment output: the records into a stream as they come, the genes they name into
// a set.

#include <algorithm>
#include <cstdint>
#include <mutex>
#include <ostream>
#include <unordered_set>
#include <vector>
#include "IO/SamFile.h"

namespace protal::test {

    class SamStreamSink : public SamSink {
    public:
        explicit SamStreamSink(std::ostream& os) : m_os(os) {}

        void Write(char const* data, size_t size, std::vector<uint64_t>& genes) override {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_os.write(data, static_cast<std::streamsize>(size));
            m_genes.insert(genes.begin(), genes.end());
            genes.clear();
        }

        // The genes named so far, sorted.
        std::vector<uint64_t> Genes() const {
            std::lock_guard<std::mutex> lock(m_mutex);
            std::vector<uint64_t> genes(m_genes.begin(), m_genes.end());
            std::sort(genes.begin(), genes.end());
            return genes;
        }

    private:
        std::ostream& m_os;
        mutable std::mutex m_mutex;
        std::unordered_set<uint64_t> m_genes;
    };
}
