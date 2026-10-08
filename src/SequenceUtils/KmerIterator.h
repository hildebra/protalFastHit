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
// x86 with GCC or Clang: the closed-syncmer scan can evaluate windows with AVX2 or AVX-512, compiled
// for them with function attributes and chosen at run time (SimpleKmerHandler::Level, Utilities/SimdLevel.h).
#if (defined(__x86_64__) || defined(__i386__)) && (defined(__GNUC__) || defined(__clang__))
#define PROTAL_SYNCMER_AVX2 1
#include <immintrin.h>
#endif
#include "SimdLevel.h"
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
        // core starts at a whole base of the k-mer, is the syncmer's k, is at most 15 bases (30 bits),
        // and has at most 16 s-mers of at most 14 bases (an s-mer and its index fit 32 bits).
        bool m_scan = false;
        simd::Level m_level = simd::Level::scalar;  // windows one by one, 8 at a time (AVX2) or 16 (AVX-512)
        std::vector<uint32_t> m_smers_fwd;  // (s-mer at i << 4), forward strand
        std::vector<uint32_t> m_smers_rev;  // (reverse complement of that s-mer << 4), coded as the reverse k-mer

        bool ScanApplies() const {
            if constexpr (std::is_same_v<CFMinimizer, ClosedSyncmer>) {
                return m_mshift % 2 == 0 && m_m <= m_k && m_k <= 31 && m_m <= 15 && m_minimizer.CoreLength() == m_m &&
                       m_minimizer.SmerCount() >= 1 && m_minimizer.SmerCount() <= 16 && m_minimizer.SmerLength() <= 14;
            }
            return false;
        }

        // The same k-mers, in the same order, as WindowByWindow: for each k-mer window, the
        // canonical core (of the forward or the reverse k-mer, the reverse on a tie) is a closed
        // syncmer if the first minimum of its s-mers is at t or at count-1-t. The window loop
        // recomputed every s-mer of every core and found the minimum with data-dependent branches;
        // here each s-mer is computed once per strand, and a window's first minimum is the smallest
        // (s-mer << 4 | index), both strands evaluated and the canonical one selected, without
        // branches: 16 windows at a time with AVX-512 or 8 with AVX2 where the CPU has it, else one by one.
        void ScanClosedSyncmers(KmerList& list) {
            size_t const windows = m_seq.length() - m_k + 1;
#ifdef PROTAL_SYNCMER_AVX2
            if (m_level == simd::Level::avx512) ScanWindowsAvx512(list);
            else if (m_level == simd::Level::avx2) ScanWindowsAvx2(list);
            else
#endif
            ScanWindows(list);
            m_total_kmers = windows;
            m_total_minimizers = list.size();
            m_pos = static_cast<uint32_t>(windows);
        }

        // The index (low 4 bits) of the first minimum of the canonical core's s-mers. In the core's
        // own order the forward s-mers are fwd[0..count); the reverse core's j-th is the reverse
        // complement of the forward one at count-1-j, rev[count-1-j].
        static uint32_t FirstMinimum(uint32_t const* fwd, uint32_t const* rev, uint32_t count, bool forward) {
            uint32_t fmin = fwd[0], rmin = rev[count - 1];
            for (uint32_t j = 1; j < count; j++) {
                fmin = std::min(fmin, fwd[j] | j);
                rmin = std::min(rmin, rev[count - 1 - j] | j);
            }
            return (forward ? fmin : rmin) & 0xF;
        }

        // The s-mers of the sequence on both strands (m_smers_fwd, m_smers_rev).
        void FillSmers() {
            size_t const n = m_seq.length();
            uint32_t const s = m_minimizer.SmerLength();
            uint32_t const smask = m_minimizer.Mask(), sfull = (1u << (2 * s)) - 1;
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
        }

        // One window at a time, the k-mers rolled along.
        void ScanWindows(KmerList& list) {
            FillSmers();
            size_t const n = m_seq.length();
            uint32_t const count = m_minimizer.SmerCount(), t = m_minimizer.T();
            size_t const core = m_k - m_m - m_mshift / 2;  // the core's first base in the window
            uint64_t kf = 0, kr = 0;
            for (size_t i = 0; i < n; i++) {
                auto const c = static_cast<unsigned char>(m_seq[i]);
                kf = ((kf << 2) | kForwardCode[c]) & m_mask;
                kr = (kr >> 2) | (uint64_t{kReverseCode[c]} << (2 * (m_k - 1)));
                if (i + 1 < m_k) continue;
                size_t const p = i + 1 - m_k;
                bool const forward = ((kf >> m_mshift) & m_mmask) < ((kr >> m_mshift) & m_mmask);
                uint32_t const first = FirstMinimum(&m_smers_fwd[p + core], &m_smers_rev[p + core], count, forward);
                if (first == t || first == count - 1 - t) list.emplace_back(KmerElement(forward ? kf : kr, p));
            }
        }

#ifdef PROTAL_SYNCMER_AVX2
        // Per base, the forward and the reverse code (kForwardCode, kReverseCode) and 64 zeros after the sequence; the
        // codes of 2, 4 and 8 consecutive bases (Bases16's steps), and of 16 on each strand (from base i: the forward
        // ones big-endian, the reverse ones little-endian; n + 16 of each, zeros past the sequence).
        std::vector<uint8_t> m_code_fwd, m_code_rev, m_bases2, m_bases4;
        std::vector<uint16_t> m_bases8;
        std::vector<uint32_t> m_b16_fwd, m_b16_rev;
        // The s-mer arrays hold this many zeros past their last s-mer, for the vector loops' last windows.
        static constexpr size_t kSmerPadding = 16;

        // The codes of 16 consecutive bases from the per-base codes `code` (n + 32 of them at least, zeros after the
        // sequence), by doubling: 2 bases from 1, 4 from 2, 8 from 4, 16 from 8, each step a loop the compiler
        // vectorises for the instruction set of the function this is inlined into (always inlined: it has no target
        // of its own). Big-endian (the first base highest, as the forward k-mer is rolled) or little-endian (the
        // first base lowest, as the reverse one is). The 16 at i are valid for i < n + 16.
        __attribute__((always_inline)) inline void Bases16(uint8_t const* code, size_t n, bool big_endian, std::vector<uint32_t>& out) {
            m_bases2.resize(n + 30);
            m_bases4.resize(n + 28);
            m_bases8.resize(n + 24);
            out.resize(n + 16);
            uint8_t* b2 = m_bases2.data();
            uint8_t* b4 = m_bases4.data();
            uint16_t* b8 = m_bases8.data();
            uint32_t* b16 = out.data();
            if (big_endian) {
                for (size_t i = 0; i < n + 30; i++) b2[i] = static_cast<uint8_t>((code[i] << 2) | code[i + 1]);
                for (size_t i = 0; i < n + 28; i++) b4[i] = static_cast<uint8_t>((b2[i] << 4) | b2[i + 2]);
                for (size_t i = 0; i < n + 24; i++) b8[i] = static_cast<uint16_t>((uint32_t{b4[i]} << 8) | b4[i + 4]);
                for (size_t i = 0; i < n + 16; i++) b16[i] = (uint32_t{b8[i]} << 16) | b8[i + 8];
            } else {
                for (size_t i = 0; i < n + 30; i++) b2[i] = static_cast<uint8_t>(code[i] | (code[i + 1] << 2));
                for (size_t i = 0; i < n + 28; i++) b4[i] = static_cast<uint8_t>(b2[i] | (b2[i + 2] << 4));
                for (size_t i = 0; i < n + 24; i++) b8[i] = static_cast<uint16_t>(b4[i] | (uint32_t{b4[i + 4]} << 8));
                for (size_t i = 0; i < n + 16; i++) b16[i] = b8[i] | (uint32_t{b8[i + 8]} << 16);
            }
        }

        // From the per-base codes (m_code_fwd, m_code_rev): the 16-base codes of both strands, and the s-mers of both
        // strands from them (s <= 14: the forward s-mer at q is the top 2s bits of the 16 bases at q, the reverse
        // one their low 2s bits), as FillSmers gives them, with kSmerPadding zeros after. Always inlined, as Bases16.
        __attribute__((always_inline)) inline void Bases16AndSmers() {
            size_t const n = m_seq.length();
            uint32_t const s = m_minimizer.SmerLength();
            uint32_t const smask = m_minimizer.Mask(), sfull = (1u << (2 * s)) - 1;
            unsigned const sshift = 32 - 2 * s;
            Bases16(m_code_fwd.data(), n, true, m_b16_fwd);
            Bases16(m_code_rev.data(), n, false, m_b16_rev);
            size_t const smers = n - s + 1;
            m_smers_fwd.resize(smers + kSmerPadding);
            m_smers_rev.resize(smers + kSmerPadding);
            std::fill(m_smers_fwd.begin() + static_cast<std::ptrdiff_t>(smers), m_smers_fwd.end(), 0u);
            std::fill(m_smers_rev.begin() + static_cast<std::ptrdiff_t>(smers), m_smers_rev.end(), 0u);
            uint32_t const* bf = m_b16_fwd.data();
            uint32_t const* br = m_b16_rev.data();
            uint32_t* sf = m_smers_fwd.data();
            uint32_t* sr = m_smers_rev.data();
            for (size_t q = 0; q < smers; q++) sf[q] = ((bf[q] >> sshift) & smask) << 4;
            for (size_t q = 0; q < smers; q++) sr[q] = (br[q] & sfull & smask) << 4;
        }

        // A window's k-mers and cores from the 16-base codes (k <= 31, m <= 15). Forward: the 32 bases from p
        // big-endian, cut to the top 2k bits; its core (bits mshift to mshift + 2m) starts at base p + k - m - mshift/2,
        // the top 2m bits of the 16 there. Reverse: the 32 bases from p little-endian, cut to the low 2k bits; its core
        // starts at base p + mshift/2, the low 2m bits of the 16 there.
        uint64_t KmerFwd(size_t p) const {
            return ((uint64_t{m_b16_fwd[p]} << 32) | m_b16_fwd[p + 16]) >> (64 - 2 * m_k);
        }
        uint64_t KmerRev(size_t p) const {
            return (uint64_t{m_b16_rev[p]} | (uint64_t{m_b16_rev[p + 16]} << 32)) & m_mask;
        }
        uint32_t CoreFwd(size_t p) const {
            return m_b16_fwd[p + m_k - m_m - m_mshift / 2] >> (32 - 2 * m_m);
        }
        uint32_t CoreRev(size_t p) const {
            return m_b16_rev[p + m_mshift / 2] & static_cast<uint32_t>(m_mmask);
        }

        // Window p one by one (the vector loops' last windows), from the arrays the vector fill made.
        void ScanWindowFromArrays(size_t p, KmerList& list) const {
            uint32_t const count = m_minimizer.SmerCount(), t = m_minimizer.T();
            size_t const core = m_k - m_m - m_mshift / 2;
            bool const forward = CoreFwd(p) < CoreRev(p);
            uint32_t const first = FirstMinimum(&m_smers_fwd[p + core], &m_smers_rev[p + core], count, forward);
            if (first == t || first == count - 1 - t) list.emplace_back(KmerElement(forward ? KmerFwd(p) : KmerRev(p), p));
        }

        // The per-base codes, 32 at a time, the rest one by one: the low nibble of the upper-case letter tells A (1),
        // C (3), G (7) and T (4) apart, the table holds their forward codes; the reverse one is 3 minus that;
        // anything else is 0 on both strands. Then the 16-base codes and the s-mers (Bases16AndSmers, for AVX2).
        __attribute__((target("avx2"))) void FillAvx2() {
            size_t const n = m_seq.length();
            m_code_fwd.resize(n + 64);
            m_code_rev.resize(n + 64);
            std::fill(m_code_fwd.begin() + static_cast<std::ptrdiff_t>(n), m_code_fwd.end(), uint8_t{ 0 });
            std::fill(m_code_rev.begin() + static_cast<std::ptrdiff_t>(n), m_code_rev.end(), uint8_t{ 0 });
            uint8_t* cf = m_code_fwd.data();
            uint8_t* cr = m_code_rev.data();
            __m256i const table = _mm256_setr_epi8(0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
                                                   0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0);
            __m256i const upper = _mm256_set1_epi8(static_cast<char>(0xDF)), low_nibble = _mm256_set1_epi8(0x0F);
            __m256i const a = _mm256_set1_epi8('A'), c = _mm256_set1_epi8('C'), g = _mm256_set1_epi8('G'), t = _mm256_set1_epi8('T');
            __m256i const three = _mm256_set1_epi8(3);
            size_t i = 0;
            for (; i + 32 <= n; i += 32) {
                __m256i const u = _mm256_and_si256(_mm256_loadu_si256(reinterpret_cast<__m256i const*>(m_seq.data() + i)), upper);
                __m256i const base = _mm256_or_si256(_mm256_or_si256(_mm256_cmpeq_epi8(u, a), _mm256_cmpeq_epi8(u, c)),
                                                     _mm256_or_si256(_mm256_cmpeq_epi8(u, g), _mm256_cmpeq_epi8(u, t)));
                __m256i const code = _mm256_shuffle_epi8(table, _mm256_and_si256(u, low_nibble));
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(cf + i), _mm256_and_si256(code, base));
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(cr + i), _mm256_and_si256(_mm256_sub_epi8(three, code), base));
            }
            for (; i < n; i++) {
                auto const ch = static_cast<unsigned char>(m_seq[i]);
                cf[i] = kForwardCode[ch];
                cr[i] = kReverseCode[ch];
            }
            Bases16AndSmers();
        }

        // As ScanWindows, 8 windows at a time (FillAvx2's arrays; a window's cores from the 16-base codes, its k-mer
        // only if it is a syncmer); the last windows (fewer than 8) one by one. Compiled for AVX2 whatever the build
        // targets, and only called where the CPU has it.
        __attribute__((target("avx2"))) void ScanWindowsAvx2(KmerList& list) {
            FillAvx2();
            size_t const windows = m_seq.length() - m_k + 1;
            uint32_t const count = m_minimizer.SmerCount(), t = m_minimizer.T();
            size_t const core = m_k - m_m - m_mshift / 2, reverse_core = m_mshift / 2;
            uint32_t const* sf = m_smers_fwd.data();
            uint32_t const* sr = m_smers_rev.data();
            uint32_t const* bf = m_b16_fwd.data();
            uint32_t const* br = m_b16_rev.data();
            __m128i const core_shift = _mm_cvtsi32_si128(static_cast<int>(32 - 2 * m_m));
            __m256i const core_mask = _mm256_set1_epi32(static_cast<int>(m_mmask));
            __m256i const t_first = _mm256_set1_epi32(static_cast<int>(t));
            __m256i const t_last = _mm256_set1_epi32(static_cast<int>(count - 1 - t));
            __m256i const index_bits = _mm256_set1_epi32(0xF);
            size_t p = 0;
            for (; p + 8 <= windows; p += 8) {
                // Lane l is window p+l: its j-th forward s-mer is at sf[p + core + j + l].
                __m256i fmin = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(sf + p + core));
                __m256i rmin = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(sr + p + core + count - 1));
                for (uint32_t j = 1; j < count; j++) {
                    __m256i const index = _mm256_set1_epi32(static_cast<int>(j));
                    __m256i const f = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(sf + p + core + j));
                    __m256i const r = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(sr + p + core + count - 1 - j));
                    fmin = _mm256_min_epu32(fmin, _mm256_or_si256(f, index));
                    rmin = _mm256_min_epu32(rmin, _mm256_or_si256(r, index));
                }
                __m256i const core_f = _mm256_srl_epi32(_mm256_loadu_si256(reinterpret_cast<__m256i const*>(bf + p + core)), core_shift);
                __m256i const core_r = _mm256_and_si256(_mm256_loadu_si256(reinterpret_cast<__m256i const*>(br + p + reverse_core)), core_mask);
                // Cores are below 2^30, so the signed comparison orders them.
                __m256i const forward = _mm256_cmpgt_epi32(core_r, core_f);
                __m256i const first = _mm256_and_si256(_mm256_blendv_epi8(rmin, fmin, forward), index_bits);
                __m256i const hit = _mm256_or_si256(_mm256_cmpeq_epi32(first, t_first), _mm256_cmpeq_epi32(first, t_last));
                auto hits = static_cast<unsigned>(_mm256_movemask_ps(_mm256_castsi256_ps(hit)));
                auto const forwards = static_cast<unsigned>(_mm256_movemask_ps(_mm256_castsi256_ps(forward)));
                while (hits) {
                    auto const lane = static_cast<unsigned>(__builtin_ctz(hits));
                    list.emplace_back(KmerElement((forwards >> lane) & 1 ? KmerFwd(p + lane) : KmerRev(p + lane), p + lane));
                    hits &= hits - 1;
                }
            }
            for (; p < windows; p++) ScanWindowFromArrays(p, list);
        }

        // As FillAvx2 with AVX-512 (simd::Level::avx512): the per-base codes 64 at a time, the last ones by a masked
        // load (nothing read past the sequence), and the compiler's loops of Bases16AndSmers vectorised for AVX-512.
#if defined(__GNUC__) && !defined(__clang__)
#define PROTAL_SYNCMER_AVX512 \
    __attribute__((target("avx2,bmi,bmi2,lzcnt,popcnt,fma,avx512f,avx512bw,avx512vl,avx512dq,avx512vbmi,avx512vbmi2,avx512vpopcntdq,prefer-vector-width=512")))
#else
#define PROTAL_SYNCMER_AVX512 PROTAL_TARGET_AVX512
#endif
        PROTAL_SYNCMER_AVX512 void FillAvx512() {
            size_t const n = m_seq.length();
            m_code_fwd.resize(n + 64);
            m_code_rev.resize(n + 64);
            uint8_t* cf = m_code_fwd.data();
            uint8_t* cr = m_code_rev.data();
            __m512i const table = _mm512_broadcast_i32x4(_mm_setr_epi8(0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0));
            __m512i const upper = _mm512_set1_epi8(static_cast<char>(0xDF)), low_nibble = _mm512_set1_epi8(0x0F);
            __m512i const a = _mm512_set1_epi8('A'), c = _mm512_set1_epi8('C'), g = _mm512_set1_epi8('G'), t = _mm512_set1_epi8('T');
            __m512i const three = _mm512_set1_epi8(3);
            // Every 64 bytes, the last step's lanes past the sequence loaded as 0 (not a base: code 0); the zeros up to
            // n + 64 come from that step and the one after it.
            for (size_t i = 0; i < n + 64; i += 64) {
                size_t const left = i < n ? n - i : 0;
                __mmask64 const in = left >= 64 ? ~__mmask64{ 0 } : static_cast<__mmask64>(_bzhi_u64(~uint64_t{ 0 }, static_cast<unsigned>(left)));
                __m512i const u = _mm512_and_si512(_mm512_maskz_loadu_epi8(in, m_seq.data() + i), upper);
                __mmask64 const base = _mm512_cmpeq_epi8_mask(u, a) | _mm512_cmpeq_epi8_mask(u, c) |
                                       _mm512_cmpeq_epi8_mask(u, g) | _mm512_cmpeq_epi8_mask(u, t);
                __m512i const code = _mm512_shuffle_epi8(table, _mm512_and_si512(u, low_nibble));
                // Stores of whole 64 bytes up to the array's n + 64 (masked where the array ends before).
                size_t const room = n + 64 - i;
                __mmask64 const store = room >= 64 ? ~__mmask64{ 0 } : static_cast<__mmask64>(_bzhi_u64(~uint64_t{ 0 }, static_cast<unsigned>(room)));
                _mm512_mask_storeu_epi8(cf + i, store, _mm512_maskz_mov_epi8(base, code));
                _mm512_mask_storeu_epi8(cr + i, store, _mm512_maskz_sub_epi8(base, three, code));
            }
            Bases16AndSmers();
        }

        // As ScanWindowsAvx2, 16 windows at a time with AVX-512, the last step masked to the windows left (the arrays
        // are padded for its loads). Only called where the CPU has the AVX-512 level.
        PROTAL_SYNCMER_AVX512 void ScanWindowsAvx512(KmerList& list) {
            FillAvx512();
            size_t const windows = m_seq.length() - m_k + 1;
            uint32_t const count = m_minimizer.SmerCount(), t = m_minimizer.T();
            size_t const core = m_k - m_m - m_mshift / 2, reverse_core = m_mshift / 2;
            uint32_t const* sf = m_smers_fwd.data();
            uint32_t const* sr = m_smers_rev.data();
            uint32_t const* bf = m_b16_fwd.data();
            uint32_t const* br = m_b16_rev.data();
            __m128i const core_shift = _mm_cvtsi32_si128(static_cast<int>(32 - 2 * m_m));
            __m512i const core_mask = _mm512_set1_epi32(static_cast<int>(m_mmask));
            __m512i const t_first = _mm512_set1_epi32(static_cast<int>(t));
            __m512i const t_last = _mm512_set1_epi32(static_cast<int>(count - 1 - t));
            __m512i const index_bits = _mm512_set1_epi32(0xF);
            for (size_t p = 0; p < windows; p += 16) {
                size_t const left = windows - p;
                __mmask16 const valid = left >= 16 ? __mmask16{ 0xFFFF } : static_cast<__mmask16>((1u << left) - 1);
                // Lane l is window p+l: its j-th forward s-mer is at sf[p + core + j + l]. Lanes past the windows read
                // the arrays' padding (kSmerPadding; the 16-base codes reach n + 16) and are masked off below.
                __m512i fmin = _mm512_loadu_si512(sf + p + core);
                __m512i rmin = _mm512_loadu_si512(sr + p + core + count - 1);
                for (uint32_t j = 1; j < count; j++) {
                    __m512i const index = _mm512_set1_epi32(static_cast<int>(j));
                    fmin = _mm512_min_epu32(fmin, _mm512_or_si512(_mm512_loadu_si512(sf + p + core + j), index));
                    rmin = _mm512_min_epu32(rmin, _mm512_or_si512(_mm512_loadu_si512(sr + p + core + count - 1 - j), index));
                }
                __m512i const core_f = _mm512_srl_epi32(_mm512_loadu_si512(bf + p + core), core_shift);
                __m512i const core_r = _mm512_and_si512(_mm512_loadu_si512(br + p + reverse_core), core_mask);
                __mmask16 const forward = _mm512_cmplt_epu32_mask(core_f, core_r);
                __m512i const first = _mm512_and_si512(_mm512_mask_blend_epi32(forward, rmin, fmin), index_bits);
                unsigned hits = static_cast<unsigned>(valid & (_mm512_cmpeq_epi32_mask(first, t_first) | _mm512_cmpeq_epi32_mask(first, t_last)));
                auto const forwards = static_cast<unsigned>(forward);
                while (hits) {
                    auto const lane = static_cast<unsigned>(__builtin_ctz(hits));
                    list.emplace_back(KmerElement((forwards >> lane) & 1 ? KmerFwd(p + lane) : KmerRev(p + lane), p + lane));
                    hits &= hits - 1;
                }
            }
        }
#endif

        // The highest level of the scan this CPU supports (AVX2 or AVX-512 on x86 with GCC or Clang).
        static simd::Level CpuLevel() {
#ifdef PROTAL_SYNCMER_AVX2
            return simd::CpuLevel();
#else
            return simd::Level::scalar;
#endif
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
            m_level = m_scan ? simd::Min(simd::Default(), CpuLevel()) : simd::Level::scalar;
        }
        SimpleKmerHandler(SimpleKmerHandler const& other) :
                m_k(other.m_k), m_m(other.m_m), m_mshift((other.m_k-other.m_m)), m_mmask((1llu << (other.m_m*2)) -1), m_mask((1llu << (other.m_k*2)) -1), m_minimizer(other.m_minimizer) {
            m_scan = ScanApplies();
            m_level = other.m_level;
        }

        // Whether operator() scans whole sequences (ScanClosedSyncmers) rather than window by window,
        // and how the scan evaluates windows (Level: the CPU's highest, capped by PROTAL_SIMD; simd::Default).
        bool Scans() const {
            return m_scan;
        }

        simd::Level Level() const {
            return m_level;
        }

        // Evaluates windows at `wanted`, or the highest level below it the CPU supports (scalar if the handler does
        // not scan), e.g. to test each; returns the level used.
        simd::Level Use(simd::Level wanted) {
            m_level = m_scan ? simd::Min(wanted, CpuLevel()) : simd::Level::scalar;
            return m_level;
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