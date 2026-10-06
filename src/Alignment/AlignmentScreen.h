// AlignmentScreen.h - a k-mer screen before WFA2. A candidate alignment of a read into its gene
// window is refused without aligning when the read and the window share too few k-mers for any
// alignment within the score budget to exist: the q-gram lemma. The candidate fails exactly as
// WFA2 would have failed it (no alignment under the budget exists, so none can be found), and its
// taxon stays one of the read's failed candidates.
//
// The bound: an alignment of score at most B touches (changes a base of, or inserts into) at most
// B / c of the read's k-mers, where c is the smallest score an operation pays per k-mer it touches:
// a mismatch (score `mismatch`) touches at most k k-mers; an insertion run of g read bases
// (`gap_opening` + g * `gap_extension`) at most k - 1 + g; a deletion run at most k - 1. Every
// other k-mer of the read's aligned part is intact: k read bases aligned to k equal, consecutive
// window bases, so it occurs in the window. With free ends, the alignment may leave up to
// `begin_free` read bases before it and `end_free` after it unaligned; the read bases between
// are aligned by every alignment, and only their k-mers count.
//
// On GTDB-scale databases most candidates of real reads fail (docs/claude/2026-10-04-performance-gtdb-scale):
// a failing WFA2 alignment explores its whole budget, where the screen costs a pass over the read
// and the window. The window's k-mers come from the gene's packed bytes, and the read's from its
// strands packed once for all its candidates (ReadKmers; docs/claude/2026-10-06-performance-profiling).
#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "SequenceUtils/KmerUtils.h"
#include "SequenceUtils/PackedSequence.h"

namespace protal {
    class AlignmentScreen {
    public:
        static constexpr size_t kMaxK = 10;    // the stamps are indexed by 2-bit codes of up to 20 bits
        static constexpr size_t kDefaultK = 6;  // for windows of up to kShortWindow bases; longer windows take kLongK
        static constexpr size_t kShortWindow = 512;
        static constexpr size_t kLongK = 7;

        // The k for a window: measured on random reads of 5-25% divergence (test AlignmentScreen.RefusesOnlyWhatWFA2Fails),
        // k = 6 refused most of the failures of 150 bp reads (47%; k = 8: 19%) and k = 7 most of those of 500-1500
        // base windows (39%): a smaller k leaves more k-mers intact per edit but matches more by chance in a long window.
        static size_t KFor(size_t window_length) { return window_length <= kShortWindow ? kDefaultK : kLongK; }

        // The most k-mers an alignment of score at most `max_score - 1` (an alignment whose score reaches max_score
        // is abandoned, WFA2Wrapper2::Alignment) can touch, with the gap-affine penalties given.
        static int64_t MaxTouched(int max_score, int mismatch, int gap_opening, int gap_extension, size_t k_mer = kDefaultK) {
            int64_t const budget = static_cast<int64_t>(max_score) - 1;
            if (budget < 0) return -1;  // nothing aligns
            int64_t const k = static_cast<int64_t>(k_mer);
            int64_t touched = 0;
            if (mismatch > 0) touched = std::max(touched, budget * k / mismatch);
            int64_t const open = static_cast<int64_t>(gap_opening) + static_cast<int64_t>(gap_extension);
            if (open > 0) touched = std::max(touched, budget * k / open);
            if (gap_extension > 0) touched = std::max(touched, budget / gap_extension);
            return touched;
        }

        // The fewest k-mers of `aligned` read bases that occur in the window if an alignment within the budget
        // exists; 0 or less: the screen cannot refuse.
        static int64_t Required(size_t aligned, int max_score, int mismatch, int gap_opening, int gap_extension, size_t k = kDefaultK) {
            if (aligned < k) return 0;
            int64_t const touched = MaxTouched(max_score, mismatch, gap_opening, gap_extension, k);
            if (touched < 0) return static_cast<int64_t>(aligned);  // more than there are: refused
            return static_cast<int64_t>(aligned - k + 1) - touched;
        }

        // A read's k-mers for the screen, made once for all its candidates: a mate's (2.2 at r226) read one of its two
        // strands, and a long read's candidate windows (19.5 per HiFi read) are stretches of the read or of its reverse
        // complement, mostly over the same genes. Each strand is packed 2 bits a base, the first base lowest (as the
        // genes are, packed::Pack, AVX2), when a candidate first reads it, with a bitmap of its bases other than A, C, G
        // or T (either case, as KmerUtils::BaseToInt takes them): a k-mer with one of those counts as shared. The k-mer
        // at base p of a strand is then the 2k bits at bit 2 (p % 4) of the 8 bytes from byte p / 4, for every k, four
        // k-mers to a load, as the window's are (StampPacked).
        class ReadKmers {
        public:
            static uint32_t Any(size_t k) { return uint32_t{ 1 } << (2 * k); }

            // A strand, packed: `packed` with 8 bytes to spare past its bases, `other` with a bit for each base other
            // than A, C, G or T (and a word to spare), `any_other` whether it has one.
            struct Strand {
                std::vector<uint8_t> packed;
                std::vector<uint64_t> other;
                bool any_other = false;

                // Whether one of the k bases from base p is not A, C, G or T.
                bool OtherIn(size_t p, size_t k) const {
                    size_t const word = p >> 6, bit = p & 63;
                    uint64_t bits = other[word] >> bit;
                    if (bit + k > 64) bits |= other[word + 1] << (64 - bit);
                    return (bits & ((uint64_t{ 1 } << k) - 1)) != 0;
                }
            };

            // The read whose k-mers come next, and its reverse complement as KmerUtils::ReverseComplementInto writes it
            // (none: made from the read if a candidate needs it); the last read's are forgotten. Both must outlive
            // their use.
            void Set(std::string_view read, std::string_view reverse = {}) {
                m_read = read;
                m_reverse = reverse;
                m_made[0] = m_made[1] = false;
            }

            size_t Length() const { return m_read.size(); }

            // A strand of the read, packed when first asked for.
            Strand const& Get(bool forward) {
                size_t const s = forward ? 0 : 1;
                if (!m_made[s]) {
                    if (!forward && m_reverse.size() != m_read.size()) {
                        KmerUtils::ReverseComplementInto(m_read, m_reverse_buffer);
                        m_reverse = m_reverse_buffer;
                    }
                    Make(m_strands[s], forward ? m_read : m_reverse);
                    m_made[s] = true;
                }
                return m_strands[s];
            }

        private:
            static void Make(Strand& strand, std::string_view bases) {
                size_t const n = bases.size();
                strand.packed.assign(packed::Bytes(n) + 8, 0);
                packed::Pack(bases.data(), n, strand.packed.data());
                strand.other.assign(n / 64 + 2, 0);
                strand.any_other = MarkOther(bases, strand.other.data());
            }

            // Sets a bit in `other` for each base other than A, C, G or T (either case); whether there was one.
            static bool MarkOther(std::string_view bases, uint64_t* other) {
                size_t i = 0;
#ifdef PROTAL_PACKED_AVX2
                if (packed::Avx2Enabled().load(std::memory_order_relaxed)) i = MarkOtherAvx2(bases, other);
#endif
                uint64_t any = 0;
                for (size_t w = 0; w < i / 64; w++) any |= other[w];
                for (; i < bases.size(); i++) {
                    if (KmerUtils::BaseToInt(bases[i]) > 3) {
                        other[i >> 6] |= uint64_t{ 1 } << (i & 63);
                        any = 1;
                    }
                }
                return any != 0;
            }

#ifdef PROTAL_PACKED_AVX2
            // A bit for each of the 32 characters at s that is A, C, G or T (either case). (A function of its own: a lambda
            // in an AVX2 function is not compiled for AVX2.)
            __attribute__((target("avx2"))) static uint32_t AcgtMaskAvx2(char const* s) {
                __m256i const u = _mm256_and_si256(_mm256_loadu_si256(reinterpret_cast<__m256i const*>(s)),
                                                   _mm256_set1_epi8(static_cast<char>(0xDF)));
                __m256i const is_base = _mm256_or_si256(_mm256_or_si256(_mm256_cmpeq_epi8(u, _mm256_set1_epi8('A')),
                                                                        _mm256_cmpeq_epi8(u, _mm256_set1_epi8('C'))),
                                                        _mm256_or_si256(_mm256_cmpeq_epi8(u, _mm256_set1_epi8('G')),
                                                                        _mm256_cmpeq_epi8(u, _mm256_set1_epi8('T'))));
                return static_cast<uint32_t>(_mm256_movemask_epi8(is_base));
            }

            // MarkOther's bits for the first bases, 64 at a time; how many bases it marked.
            __attribute__((target("avx2"))) static size_t MarkOtherAvx2(std::string_view bases, uint64_t* other) {
                size_t i = 0;
                for (; i + 64 <= bases.size(); i += 64) {
                    uint64_t const is_base = uint64_t{ AcgtMaskAvx2(bases.data() + i) } | (uint64_t{ AcgtMaskAvx2(bases.data() + i + 32) } << 32);
                    other[i >> 6] = ~is_base;
                }
                return i;
            }
#endif

            std::string_view m_read, m_reverse;
            std::string m_reverse_buffer;
            std::array<Strand, 2> m_strands;
            std::array<bool, 2> m_made{};
        };

        // Where a candidate's read lies in a read whose packed strands its candidates share (ReadKmers): the
        // candidate's read from base fwd_offset of the forward strand, or for a reverse candidate its reverse
        // complement from base rev_offset of the reverse strand (a short read: both 0; a long read's window [s, e) of a
        // read of length n: s and n - e).
        struct ReadStretch {
            ReadKmers* kmers = nullptr;  // none: the candidate's read is packed for it alone
            size_t fwd_offset = 0, rev_offset = 0;
        };

        // Whether an alignment of `read` into `window` of score below max_score may exist, with up to begin_free
        // read bases before and end_free after it left unaligned: true unless the read's aligned part shares
        // fewer k-mers with the window than any such alignment leaves intact. A read k-mer with a base other
        // than ACGT counts as shared (WFA2 matches equal characters, whatever they are); window k-mers with such
        // a base match nothing.
        bool MayAlign(std::string_view read, size_t begin_free, size_t end_free, std::string_view window,
                      int max_score, int mismatch, int gap_opening, int gap_extension, size_t k = 0) {
            if (k == 0) k = KFor(window.size());
            k = std::min(k, kMaxK);
            int64_t required = 0;
            if (auto const answer = Bound(read.size(), begin_free, end_free, max_score, mismatch, gap_opening, gap_extension, k, required)) {
                return *answer;
            }
            Stamp(window, k);
            m_own.Set(read);
            return SharesEnough(m_own, true, 0, read.size(), begin_free, end_free, k, required);
        }

        // As MayAlign, with the window as the bases [begin, end) of a gene's 2-bit packed bytes (PackedSequence.h: base i
        // at bits 2 (i % 4) of byte i / 4; `bytes` bytes in all, as packed::Bytes gives them for the gene's length). Every
        // packed base is A, C, G or T, as every base of the decoded window is, so the window has the same k-mers and the
        // answer is MayAlign's on the decoded window; but a refused candidate (most of them at GTDB scale) is never
        // decoded, and the k-mers come four to an 8-byte load (docs/claude/2026-10-06-performance-profiling).
        bool MayAlignPacked(std::string_view read, size_t begin_free, size_t end_free, uint8_t const* packed, size_t bytes,
                            size_t begin, size_t end, int max_score, int mismatch, int gap_opening, int gap_extension, size_t k = 0) {
            m_own.Set(read);
            return MayAlignPacked(m_own, true, 0, read.size(), begin_free, end_free, packed, bytes, begin, end, max_score,
                                  mismatch, gap_opening, gap_extension, k);
        }

        // As MayAlignPacked, for a candidate's read of `length` bases from base `offset` of a strand of `kmers`' read
        // (ReadStretch): the strand is packed once for all the read's candidates, and the exits tested every kExitEvery
        // k-mers (the shared k-mers only grow and shared plus remaining only shrinks, so the answer is the same).
        bool MayAlignPacked(ReadKmers& kmers, bool forward, size_t offset, size_t length, size_t begin_free, size_t end_free,
                            uint8_t const* packed, size_t bytes, size_t begin, size_t end, int max_score, int mismatch,
                            int gap_opening, int gap_extension, size_t k = 0) {
            if (k == 0) k = KFor(end > begin ? end - begin : 0);
            k = std::min(k, kMaxK);
            int64_t required = 0;
            if (auto const answer = Bound(length, begin_free, end_free, max_score, mismatch, gap_opening, gap_extension, k, required)) {
                return *answer;
            }
            StampPacked(packed, bytes, begin, end, k);
            return SharesEnough(kmers, forward, offset, length, begin_free, end_free, k, required);
        }

        // Candidates refused so far.
        size_t Refused() const { return m_refused; }
        void CountRefused() { m_refused++; }

    private:
        // The answer when the bound alone gives it (nothing to align, the screen cannot refuse, or more shared k-mers
        // required than the read has); else nullopt, with the k-mers required.
        static std::optional<bool> Bound(size_t read_length, size_t begin_free, size_t end_free, int max_score, int mismatch,
                                         int gap_opening, int gap_extension, size_t k, int64_t& required) {
            if (begin_free + end_free >= read_length) return true;
            size_t const aligned = read_length - begin_free - end_free;
            required = Required(aligned, max_score, mismatch, gap_opening, gap_extension, k);
            if (required <= 0) return true;
            if (required > static_cast<int64_t>(aligned - k + 1)) return false;
            return std::nullopt;
        }

        // A fresh generation stamp for a window's k-mers, the stamps sized for k and ReadKmers::Any(k), which is always
        // stamped: a read k-mer with a base other than ACGT counts as shared.
        void NextGeneration(size_t k) {
            size_t const any = ReadKmers::Any(k);
            if (m_stamps.size() < any + 1) { m_stamps.assign(any + 1, 0); m_generation = 0; }
            if (++m_generation == 0) {  // wrapped: every old stamp reads as 0, so start over at 1
                std::fill(m_stamps.begin(), m_stamps.end(), 0);
                m_generation = 1;
            }
            m_stamps[any] = m_generation;
        }

        // K-mers are coded with their first base in the lowest bits (the packed bytes' order), the window's and the read's
        // alike: the next base goes in at the top.
        static uint32_t Roll(uint32_t code, uint32_t base, size_t k) {
            return (code >> 2) | (base << (2 * (k - 1)));
        }

        // Marks the window's k-mers with a fresh generation stamp.
        void Stamp(std::string_view window, size_t k) {
            NextGeneration(k);
            uint32_t code = 0;
            size_t valid = 0;
            for (char const c : window) {
                uint32_t const base = static_cast<uint32_t>(KmerUtils::BaseToInt(c));
                if (base > 3) { valid = 0; continue; }
                valid++;
                code = Roll(code, base, k);
                if (valid >= k) m_stamps[code] = m_generation;
            }
        }

        // As Stamp, for the bases [begin, end) of packed bytes: the k-mer at base p is the 2k bits at bit 2 (p % 4) of the
        // 8 bytes from byte p / 4, so an aligned load gives four k-mers (2 (3 + k) <= 26 bits); the k-mers whose load
        // would pass the gene's last byte are put together base by base.
        void StampPacked(uint8_t const* packed, size_t bytes, size_t begin, size_t end, size_t k) {
            NextGeneration(k);
            if (end < begin + k) return;
            uint32_t const mask = (1u << (2 * k)) - 1;
            uint16_t const g = m_generation;
            size_t const last = end - k;  // the first base of the window's last k-mer
            auto code_at = [&](size_t p) {
                uint32_t code = 0;
                for (size_t j = 0; j < k; j++) code |= static_cast<uint32_t>((packed[(p + j) >> 2] >> (2 * ((p + j) & 3))) & 3u) << (2 * j);
                return code;
            };
            size_t p = begin;
            for (; p <= last && (p & 3) != 0; p++) m_stamps[code_at(p)] = g;
            for (; p + 3 <= last && (p >> 2) + 8 <= bytes; p += 4) {
                uint64_t w;
                std::memcpy(&w, packed + (p >> 2), 8);
                m_stamps[static_cast<uint32_t>(w) & mask] = g;
                m_stamps[static_cast<uint32_t>(w >> 2) & mask] = g;
                m_stamps[static_cast<uint32_t>(w >> 4) & mask] = g;
                m_stamps[static_cast<uint32_t>(w >> 6) & mask] = g;
            }
            for (; p <= last; p++) m_stamps[code_at(p)] = g;
        }

        static constexpr size_t kExitEvery = 16;  // a multiple of 4

        // Whether the candidate's aligned part, the bases [offset + begin_free, offset + length - end_free) of a strand of
        // `kmers`' read, has `required` k-mers stamped in the window (one with a base other than A, C, G or T counts as
        // shared: ReadKmers::Any, always stamped). Its k-mers start at the bases [first, end) (Bound left at least one),
        // four to a load from base 4i on; the exits are tested every kExitEvery k-mers, which stops at the same answer.
        bool SharesEnough(ReadKmers& kmers, bool forward, size_t offset, size_t length, size_t begin_free, size_t end_free,
                          size_t k, int64_t required) const {
            assert(offset + length <= kmers.Length());
            auto const& strand = kmers.Get(forward);
            uint8_t const* const packed = strand.packed.data();
            uint64_t const mask = (uint64_t{ 1 } << (2 * k)) - 1;
            uint32_t const any = ReadKmers::Any(k);
            uint16_t const g = m_generation;
            uint16_t const* const stamps = m_stamps.data();
            auto code_at = [&](size_t p) -> uint32_t {
                uint64_t w;
                std::memcpy(&w, packed + (p >> 2), 8);
                return static_cast<uint32_t>((w >> (2 * (p & 3))) & mask);
            };
            size_t p = offset + begin_free;
            size_t const end = offset + length - end_free - k + 1;
            int64_t shared = 0;
            int64_t remaining = static_cast<int64_t>(end - p);  // k-mers not yet counted
            while (p < end) {
                // The first group ends at a multiple of 4, so that the others start at one.
                size_t const group_end = std::min(end, (p + kExitEvery) & ~size_t{ 3 });
                remaining -= static_cast<int64_t>(group_end - p);
                if (strand.any_other) {
                    for (; p < group_end; p++) shared += stamps[strand.OtherIn(p, k) ? any : code_at(p)] == g;
                } else {
                    for (; p < group_end && (p & 3) != 0; p++) shared += stamps[code_at(p)] == g;
                    for (; p + 4 <= group_end; p += 4) {
                        uint64_t w;
                        std::memcpy(&w, packed + (p >> 2), 8);
                        shared += (stamps[w & mask] == g) + (stamps[(w >> 2) & mask] == g) + (stamps[(w >> 4) & mask] == g) +
                                  (stamps[(w >> 6) & mask] == g);
                    }
                    for (; p < group_end; p++) shared += stamps[code_at(p)] == g;
                }
                if (shared >= required) return true;
                if (shared + remaining < required) return false;
            }
            return shared >= required;
        }

        std::vector<uint16_t> m_stamps;  // per k-mer code, the generation that last saw it in a window; sized on first use
        uint16_t m_generation = 0;
        size_t m_refused = 0;
        ReadKmers m_own;  // the k-mers of a candidate's read screened without a ReadStretch
    };
}
