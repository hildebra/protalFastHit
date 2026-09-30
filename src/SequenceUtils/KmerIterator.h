//
// Created by fritsche on 14/08/22.
//

#pragma once

#include <algorithm>
#include <array>
#include <cstring>
#include <cassert>
#include <type_traits>
#include <vector>
#include "FastxReader.h"
#include "KmerUtils.h"
#include "Minimizer.h"
#include "Constants.h"

namespace protal {



    // make templated
//    template<typename ContextFreeMinimizer=Syncmer15>
    class CanonicalKmerIterator {
    private:
//        FastxRecord m_record;
        uint32_t m_key_bytes;
        uint32_t m_last_pos_shift;
        uint32_t m_k;
        size_t m_mask;
        uint32_t m_pos = 0;
//        std::string m_sequence;
        char* m_seq = nullptr;
        size_t m_seq_length;
        bool m_first = true;

        size_t m_canonical_kmer_fwd = 0;
        size_t m_canonical_kmer_rev = 0;

        void GetKey(size_t& key) {
            for (int i = 0; i < m_k; i++) {
                key |= KmerUtils::BaseToInt(m_seq[i], 0) << (2 * (m_k - i - 1));
            }
        }

        void GetKeyC(size_t& key) {
            for (int i = 0; i < m_k; i++) {
                key |= KmerUtils::BaseToIntC(m_seq[i], 0) << (2llu * i);
            }
        }

        void RollKey(size_t &key) {
            key <<= 2;
            key &= m_mask;
            key |= KmerUtils::BaseToInt(m_seq[m_pos + m_k - 1], 0);
        }

        void RollKeyC(size_t &key) {
            key >>= 2;
            key |= KmerUtils::BaseToIntC(m_seq[m_pos + m_k - 1], 0) << 2llu * (m_k - 1);
        }

    public:
        int64_t GetPos() const {
            return m_pos - 1;
        }

        void Reset() {
            m_pos = 0;
            m_first = true;
            m_canonical_kmer_fwd = 0;
            m_canonical_kmer_rev = 0;
        }
//
//        void SetRecord(FastxRecord& record) {
//            m_record = record;
//            m_sequence = record.sequence;
//            Reset();
//        }

        // TODO this is not ideal -  should not have to convert const char* to char*
        void SetSequence(char* seq, size_t len) {
            m_seq = seq;
            m_seq_length = len;
            Reset();
        }
//
//        void SetSequence(std::string const& sequence) {
//            m_sequence = sequence;
//            Reset();
//        }

        bool HasNext() {
            return (m_pos < m_seq_length - m_k + 1);
        }


        void Init(size_t k) {
            m_k = k;
            m_key_bytes = (2 * m_k + 8 - 1) / 8;
            m_last_pos_shift = (8 - ((2 * m_k) % 8)) % 8;

            m_mask = (1llu << (m_k*2)) -1;
        }

        void Init(CanonicalKmerIterator const& other) {
            m_k = other.m_k;
            Init(m_k);
        }

        std::string ToString() const {
            std::string str = "";
            str += "k: " + std::to_string(m_k);
            return str;
        }

        CanonicalKmerIterator() {};

        CanonicalKmerIterator(size_t k) : m_k(k) {
            Init(m_k);
        }


        ~CanonicalKmerIterator() {}

        inline void KmerList(KmerList& list) {
            size_t kmer;
            list.clear();
            while ((*this)(kmer)) {
                list.emplace_back( KmerElement(kmer, GetPos()) );
            }
        }

        inline bool operator () (size_t& key) {
            if (!HasNext()) return false;
            if (m_first) {
                GetKey(m_canonical_kmer_fwd);
                GetKeyC(m_canonical_kmer_rev);
                m_first = false;
            } else {
                RollKey(m_canonical_kmer_fwd);
                RollKeyC(m_canonical_kmer_rev);
            }
            bool take_normal = m_canonical_kmer_fwd < m_canonical_kmer_rev;
            key = (take_normal * m_canonical_kmer_fwd) + (!take_normal * m_canonical_kmer_rev);

            m_pos++;
            return true;
        }
    };

    template<typename CFMinimizer>
//    requires ContextFreeMinimizer<CFMinimizer>
    class SimpleKmerHandler {
    private:
        CFMinimizer m_minimizer{};

        const uint32_t m_k{};
        const uint32_t m_m{};

        size_t m_kmer{};
        size_t m_kmer_fwd{};
        size_t m_kmer_rev{};

        size_t m_mmer{};
        size_t m_mmer_fwd{};
        size_t m_mmer_rev{};

        size_t m_mask{};
        size_t m_mmask{};
        size_t m_mshift{};

        uint32_t m_pos = 0;

//        char* m_seq = nullptr;
        std::string_view m_seq{};
        bool m_first = true;

        size_t m_total_kmers = 0;
        size_t m_total_minimizers = 0;

        size_t m_canonical_kmer_fwd = 0;
        size_t m_canonical_kmer_rev = 0;

        void GetKey(size_t& key) {
            for (int i = 0; i < m_k; i++) {
                key |= KmerUtils::BaseToInt(m_seq[i], 0) << (2 * (m_k - i - 1));
            }
        }

        void GetKeyC(size_t& key) {
            for (int i = 0; i < m_k; i++) {
                key |= KmerUtils::BaseToIntC(m_seq[i], 0) << (2llu * i);
            }
        }

        void RollKey(size_t &key) {
            key <<= 2;
            key &= m_mask;
            key |= KmerUtils::BaseToInt(m_seq[m_pos + m_k - 1], 0);
        }

        void RollKeyC(size_t &key) {
            key >>= 2;
            key |= KmerUtils::BaseToIntC(m_seq[m_pos + m_k - 1], 0) << 2llu * (m_k - 1);
        }

        // Base codes of the forward and the reverse (complemented) k-mer, as BaseToInt(c, 0) and
        // BaseToIntC(c, 0): anything but ACGT (e.g. N) is 0 on both strands.
        static constexpr std::array<uint8_t, 256> kForwardCode = [] {
            std::array<uint8_t, 256> code{};
            code['C'] = code['c'] = 1; code['G'] = code['g'] = 2; code['T'] = code['t'] = 3;
            return code;
        }();
        static constexpr std::array<uint8_t, 256> kReverseCode = [] {
            std::array<uint8_t, 256> code{};
            code['A'] = code['a'] = 3; code['C'] = code['c'] = 2; code['G'] = code['g'] = 1;
            return code;
        }();

        // Closed syncmers of a whole sequence (ScanClosedSyncmers), for the layout protal uses: the
        // core starts at a whole base of the k-mer, is the syncmer's k, and has at most 16 s-mers of
        // at most 14 bases (an s-mer and its index fit 32 bits).
        bool m_scan = false;
        std::vector<uint32_t> m_smers_fwd;  // (s-mer at i << 4), forward strand
        std::vector<uint32_t> m_smers_rev;  // (reverse complement of that s-mer << 4), coded as the reverse k-mer

        bool ScanApplies() const {
            if constexpr (std::is_same_v<CFMinimizer, ClosedSyncmer>) {
                return m_mshift % 2 == 0 && m_m <= m_k && m_minimizer.CoreLength() == m_m &&
                       m_minimizer.SmerCount() >= 1 && m_minimizer.SmerCount() <= 16 && m_minimizer.SmerLength() <= 14;
            }
            return false;
        }

        // The same k-mers, in the same order, as the window-by-window loop in operator(): for each
        // k-mer window, the canonical core (of the forward or the reverse k-mer, the reverse on a
        // tie) is a closed syncmer if the first minimum of its s-mers is at t or at count-1-t. The
        // loop recomputed every s-mer of every core and found the minimum with data-dependent
        // branches; here each s-mer is computed once per strand, and each window's first minimum
        // is the smallest (s-mer << 4 | index), both strands evaluated without branching.
        void ScanClosedSyncmers(KmerList& list) {
            size_t const n = m_seq.length();
            uint32_t const s = m_minimizer.SmerLength(), count = m_minimizer.SmerCount(), t = m_minimizer.T();
            uint32_t const smask = m_minimizer.Mask(), sfull = (1u << (2 * s)) - 1;
            size_t const core = m_k - m_m - m_mshift / 2;  // the core's first base in the window

            m_smers_fwd.resize(n - s + 1);
            m_smers_rev.resize(n - s + 1);
            uint32_t f = 0, r = 0;
            for (size_t i = 0; i < n; i++) {
                auto const c = static_cast<unsigned char>(m_seq[i]);
                f = ((f << 2) | kForwardCode[c]) & sfull;
                r = (r >> 2) | (uint32_t{kReverseCode[c]} << (2 * (s - 1)));
                if (i + 1 >= s) {
                    m_smers_fwd[i + 1 - s] = (f & smask) << 4;
                    m_smers_rev[i + 1 - s] = (r & smask) << 4;
                }
            }

            uint64_t kf = 0, kr = 0;
            for (size_t i = 0; i < n; i++) {
                auto const c = static_cast<unsigned char>(m_seq[i]);
                kf = ((kf << 2) | kForwardCode[c]) & m_mask;
                kr = (kr >> 2) | (uint64_t{kReverseCode[c]} << (2 * (m_k - 1)));
                if (i + 1 < m_k) continue;
                size_t const p = i + 1 - m_k;
                m_total_kmers++;
                uint64_t const core_f = (kf >> m_mshift) & m_mmask, core_r = (kr >> m_mshift) & m_mmask;
                // The core's s-mers in its own order: forward at p+core+j; the reverse core's j-th is
                // the reverse complement of the forward one at p+core+count-1-j.
                uint32_t const* fwd = &m_smers_fwd[p + core];
                uint32_t const* rev = &m_smers_rev[p + core];
                uint32_t fmin = fwd[0], rmin = rev[count - 1];
                for (uint32_t j = 1; j < count; j++) {
                    fmin = std::min(fmin, fwd[j] | j);
                    rmin = std::min(rmin, rev[count - 1 - j] | j);
                }
                bool const forward = core_f < core_r;
                uint32_t const first = (forward ? fmin : rmin) & 0xF;
                if (first == t || first == count - 1 - t) {
                    m_total_minimizers++;
                    list.emplace_back(KmerElement(forward ? kf : kr, p));
                }
            }
            m_pos = static_cast<uint32_t>(n - m_k + 1);
        }

    public:
        /**
         * Implementing KmerStatisticsConcept
         * @returns Total Kmers processed
         */
        size_t TotalKmers() {
            return m_total_kmers;
        }

        /**
         * Implementing MinimizerStatisticsConcept
         * @returns Total Kmers processed
         */
        size_t TotalMinimizers() {
            return m_total_minimizers;
        }

        SimpleKmerHandler(size_t k, size_t m, CFMinimizer& minimizer) :
                m_k(k), m_m(m), m_mshift((k-m)), m_mmask((1llu << (m*2)) -1), m_mask((1llu << (k*2)) -1), m_minimizer(minimizer) {
            m_scan = ScanApplies();
        }
        SimpleKmerHandler(SimpleKmerHandler const& other) :
                m_k(other.m_k), m_m(other.m_m), m_mshift((other.m_k-other.m_m)), m_mmask((1llu << (other.m_m*2)) -1), m_mask((1llu << (other.m_k*2)) -1), m_minimizer(other.m_minimizer) {
            m_scan = ScanApplies();
        }

        // Whether operator() scans whole sequences (ScanClosedSyncmers) rather than window by window.
        bool Scans() const {
            return m_scan;
        }


        int64_t GetPos() const {
            return m_pos - 1;
        }

        void Reset() {
            m_pos = 0;
            m_first = true;
            m_canonical_kmer_fwd = 0;
            m_canonical_kmer_rev = 0;

            m_total_kmers = 0;
            m_total_minimizers = 0;
        }

        void SetSequence(std::string_view const& sequence) {
            m_seq = sequence;
            Reset();
        }
        void SetSequence(std::string_view const&& sequence) {
            m_seq = sequence;
            Reset();
        }

        bool HasNext() {
            return (m_pos < m_seq.length() - m_k + 1);
        }


        [[nodiscard]] std::string ToString() const {
            std::string str = "";
            str += "k: " + std::to_string(m_k);
            return str;
        }

        // The k-mers of sequence whose canonical core passes the minimizer, with their window start:
        // for ClosedSyncmer by scanning the whole sequence (ScanClosedSyncmers), else window by window.
        inline void operator () (std::string_view const& sequence, KmerList& list) {
            if (m_scan) {
                list.clear();
                SetSequence(sequence);
                if (m_seq.length() >= m_k) ScanClosedSyncmers(list);
                return;
            }
            WindowByWindow(sequence, list);
        }

        inline void operator () (std::string_view const&& sequence, KmerList& list) {
            (*this)(sequence, list);
        }

        // The definition: each k-mer window in turn, its canonical core tested by the minimizer.
        inline void WindowByWindow(std::string_view const& sequence, KmerList& list) {
            list.clear();
            SetSequence(sequence);

            if (m_seq.length() < m_k) {
                return;
            }

            while (NextWrapper(m_kmer_fwd, m_kmer_rev)) {
                m_mmer_fwd = (m_kmer_fwd >> m_mshift) & m_mmask;
                m_mmer_rev = (m_kmer_rev >> m_mshift) & m_mmask;

                m_kmer = m_mmer_fwd < m_mmer_rev ? m_kmer_fwd : m_kmer_rev;
                m_mmer = m_mmer_fwd < m_mmer_rev ? m_mmer_fwd : m_mmer_rev;

                if (m_minimizer(m_mmer)) {
                    m_total_minimizers++;
                    list.emplace_back(KmerElement(m_kmer, GetPos()));
                }
            }
        }

        inline bool Next(size_t& key) {
            exit(127);
            while (NextWrapper(m_kmer)) {
                if (m_minimizer(m_kmer)) {
                    key = m_kmer;
                    return true;
                }
            }
            return false;
        }

        inline bool NextWrapper(size_t& key) {
            exit(127);
            if (!HasNext()) return false;
            m_total_kmers++;
            if (m_first) {
                GetKey(m_canonical_kmer_fwd);
                GetKeyC(m_canonical_kmer_rev);
                m_first = false;
            } else {
                RollKey(m_canonical_kmer_fwd);
                RollKeyC(m_canonical_kmer_rev);
            }
            bool take_normal = m_canonical_kmer_fwd < m_canonical_kmer_rev;
            key = (take_normal * m_canonical_kmer_fwd) + (!take_normal * m_canonical_kmer_rev);

            m_pos++;
            return true;
        }

        inline bool NextWrapper(size_t& key_fwd, size_t& key_rev) {
            if (!HasNext()) return false;
            m_total_kmers++;
            if (m_first) {
                GetKey(m_canonical_kmer_fwd);
                GetKeyC(m_canonical_kmer_rev);
                m_first = false;
            } else {
                RollKey(m_canonical_kmer_fwd);
                RollKeyC(m_canonical_kmer_rev);
            }
            key_fwd = m_canonical_kmer_fwd;
            key_rev = m_canonical_kmer_rev;

            m_pos++;
            return true;
        }
    };
}