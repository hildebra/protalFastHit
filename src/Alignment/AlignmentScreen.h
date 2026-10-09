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
//
// IndelBound, at the end: an exact bound for ONT reads in place of the k-mer screen, from the indel distance, which
// refuses the failing candidates of ONT reads that the k-mer screen cannot (docs/claude/2026-10-09-ont-wfa2-skipping.md).
#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "SequenceUtils/KmerUtils.h"
#include "SequenceUtils/PackedSequence.h"
#include "TargetClones.h"

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#include <immintrin.h>  // _addcarry_u64 (IndelBound)
#endif

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

    // An exact bound before WFA2 for ONT reads, in place of the k-mer screen: the indel distance
    // (docs/claude/2026-10-09-ont-wfa2-skipping.md). At ONT's identity floor of 0.85 an alignment within the budget may
    // touch every k-mer of its window, so the k-mer screen refuses next to nothing (none of 145,184 candidates on the
    // bench sample, 1,191 of 17.9M on r226's host sample), while a failing candidate (a gene window against unrelated read
    // sequence: 99.8% of them on host reads at r226) is still far from aligning. PacBio reads keep the k-mer screen: the
    // candidates it passes mostly align, and the bound on them cost 20% more instructions than it saved.
    //
    // The bound: a mismatch is a deletion and an insertion, and a gap of g bases costs gap_opening + g gap_extension, so
    // an alignment's score is at least q / 2 per base of the indel distance of what it aligns (the fewest insertions and
    // deletions that turn one into the other), q = min(mismatch, 2 gap_extension). With free ends it aligns the read
    // bases [a, n - b), a <= begin_free and b <= end_free, to the window bases [c, m - d), c <= ref_begin_free and
    // d <= ref_end_free: at least n_req = n - begin_free - end_free and m_req = m - ref_begin_free - ref_end_free bases.
    // Their indel distance is their lengths less twice their longest common subsequence (LCS), which is at most the whole
    // read's and window's. An alignment that succeeds scores at most max_score - 1 (WFA2Wrapper2::Alignment), so
    //     q (n_req + m_req - 2 LCS(read, window)) > 2 (max_score - 1)
    // means that none exists, and the candidate fails as WFA2 would have failed it. The anchored aligner's path through
    // the links is one such alignment (its pieces share the budget, and its free bases are the window's), so it fails
    // too.
    //
    // The LCS is computed bit-parallel over the window's bases, one read base at a time (Allison and Dix 1986; Crochemore
    // et al. 2001), from the gene's packed bytes, so a refused window is never decoded. A read base other than A, C, G or
    // T matches every window base, as does a window base other than those: the LCS can only grow, and the bound stays
    // exact. It costs up to m n / 64 word steps: ~345k instructions per candidate on the bench ONT sample, where a refused
    // candidate saved ~6.4M of WFA2 (docs/claude/2026-10-09-ont-wfa2-skipping.md, section 7), so it pays from ~5% refused.
    // WFA2's cost grows with the budget only (wf-adaptive keeps its band narrow), so windows of more than kMaxCells read x
    // window bases are left to WFA2.
    class IndelBound {
    public:
        static constexpr uint64_t kMaxCells = uint64_t{ 1 } << 27;

        // The largest LCS of a candidate with `required` read and window bases (n_req + m_req) that the bound refuses;
        // -1: it refuses none (or the penalties give no bound).
        static int64_t MaxRefusedLcs(size_t required, int max_score, int mismatch, int gap_opening, int gap_extension) {
            int64_t const q = std::min<int64_t>(mismatch, 2 * static_cast<int64_t>(gap_extension));
            if (q <= 0 || gap_opening < 0) return -1;
            // q (required - 2 L) >= 2 max_score - 1
            int64_t const numerator = q * static_cast<int64_t>(required) - 2 * static_cast<int64_t>(max_score) + 1;
            return numerator < 0 ? -1 : numerator / (2 * q);
        }

        // Whether an alignment of `read` into the bases [begin, end) of a gene's 2-bit packed bytes (PackedSequence.h) of
        // score below max_score may exist, with up to begin_free and end_free read bases and ref_begin_free and
        // ref_end_free window bases free at the ends: false only if the bound above refuses it.
        bool MayAlignPacked(std::string_view read, size_t begin_free, size_t end_free, uint8_t const* packed, size_t begin, size_t end,
                            size_t ref_begin_free, size_t ref_end_free, int max_score, int mismatch, int gap_opening, int gap_extension) {
            size_t const m = end > begin ? end - begin : 0;
            int64_t max_lcs = 0;
            if (auto const answer = Bound(read.size(), begin_free, end_free, m, ref_begin_free, ref_end_free, max_score, mismatch,
                                          gap_opening, gap_extension, max_lcs)) {
                return *answer;
            }
            MasksPacked(packed, begin, m);
            return LcsAbove(read, m, max_lcs);
        }

        // As MayAlignPacked, with the window as text.
        bool MayAlign(std::string_view read, size_t begin_free, size_t end_free, std::string_view window, size_t ref_begin_free,
                      size_t ref_end_free, int max_score, int mismatch, int gap_opening, int gap_extension) {
            int64_t max_lcs = 0;
            if (auto const answer = Bound(read.size(), begin_free, end_free, window.size(), ref_begin_free, ref_end_free, max_score,
                                          mismatch, gap_opening, gap_extension, max_lcs)) {
                return *answer;
            }
            Masks(window);
            return LcsAbove(read, window.size(), max_lcs);
        }

        // The LCS of `read` and `window` as the bound takes it (bases other than A, C, G or T match every base), without
        // its exits: for tests.
        size_t Lcs(std::string_view read, std::string_view window) {
            Masks(window);
            LcsAbove(read, window.size(), std::numeric_limits<int64_t>::max(), false);
            return static_cast<size_t>(m_lcs);
        }

    private:
        static constexpr size_t kExitEvery = 64;  // read bases between the exits' tests

        // The answer when the lengths alone give it (the bound refuses nothing, the window is too large to be worth it,
        // or not even an LCS of every base of the shorter one would be enough); else nullopt, with the largest LCS
        // refused.
        static std::optional<bool> Bound(size_t n, size_t begin_free, size_t end_free, size_t m, size_t ref_begin_free,
                                         size_t ref_end_free, int max_score, int mismatch, int gap_opening, int gap_extension,
                                         int64_t& max_lcs) {
            size_t const n_req = begin_free + end_free < n ? n - begin_free - end_free : 0;
            size_t const m_req = ref_begin_free + ref_end_free < m ? m - ref_begin_free - ref_end_free : 0;
            max_lcs = MaxRefusedLcs(n_req + m_req, max_score, mismatch, gap_opening, gap_extension);
            if (max_lcs < 0) return true;
            if (static_cast<int64_t>(std::min(n, m)) <= max_lcs) return false;
            if (static_cast<uint64_t>(n) * m > kMaxCells) return true;
            return std::nullopt;
        }

        // The match masks of a window of m bases: m_masks[c] has bit i for each window base i that read base code c
        // (KmerUtils::BaseToInt: A, C, G, T 0-3, anything else 4) matches; code 4 matches every base.
        void Clear(size_t m) {
            size_t const words = (m + 63) / 64;
            for (size_t c = 0; c < 4; c++) m_masks[c].assign(words, 0);
            m_masks[4].assign(words, ~uint64_t{ 0 });
            if (m % 64 != 0) m_masks[4][words - 1] = (uint64_t{ 1 } << (m % 64)) - 1;
        }

        void MasksPacked(uint8_t const* packed, size_t begin, size_t m) {
            Clear(m);
            for (size_t i = 0; i < m; i++) {
                size_t const p = begin + i;
                m_masks[(packed[p >> 2] >> (2 * (p & 3))) & 3u][i >> 6] |= uint64_t{ 1 } << (i & 63);
            }
        }

        void Masks(std::string_view window) {
            Clear(window.size());
            for (size_t i = 0; i < window.size(); i++) {
                uint64_t const code = KmerUtils::BaseToInt(window[i]);
                uint64_t const bit = uint64_t{ 1 } << (i & 63);
                if (code <= 3) {
                    m_masks[code][i >> 6] |= bit;
                } else {  // matches every read base
                    for (size_t c = 0; c < 4; c++) m_masks[c][i >> 6] |= bit;
                }
            }
        }

        // The zeros of V's lowest p bits: the LCS of the read bases so far and the window's first p bases (each bit of V is
        // a column of the dynamic programming table, a zero where its value grows).
        int64_t ZerosBelow(size_t p) const {
            size_t const full = p / 64;
            int64_t ones = 0;
            for (size_t w = 0; w < full; w++) ones += __builtin_popcountll(m_v[w]);
            if (p % 64 != 0) ones += __builtin_popcountll(m_v[full] & ((uint64_t{ 1 } << (p % 64)) - 1));
            return static_cast<int64_t>(p) - ones;
        }

        // Whether the LCS of `read` and the window of m bases in m_masks is above max_lcs. Per read base with match mask
        // M: V = (V + (V & M)) | (V & ~M), from V all ones; the LCS is the number of zeros. Every kExitEvery read bases, two
        // exits, which give the final answer: the LCS so far only grows; and a common subsequence of the whole read is one
        // of the read so far with the window's first p bases followed by one of the r read bases left with the rest, so the
        // LCS is at most max over p of ZerosBelow(p) + min(r, m - p), which is ZerosBelow(m - r) + r (r < m; each bit adds
        // at most one zero): an unrelated window's matches so far are spread over all of it, which this counts.
        // The step's addition carries from word to word, which bounds its speed: on x86-64 one add-with-carry per word
        // (the carry stays in the flag), elsewhere the carry from two comparisons. Cloned for x86-64-v3 (TargetClones.h),
        // where x & ~match is one instruction.
        PROTAL_CLONE_V3 bool LcsAbove(std::string_view read, size_t m, int64_t max_lcs, bool exits = true) {
            size_t const words = (m + 63) / 64;
            m_v.assign(words, ~uint64_t{ 0 });
            uint64_t* const v = m_v.data();
            size_t const n = read.size();
            m_lcs = 0;
            size_t i = 0;
            while (i < n) {
                size_t const block_end = std::min(n, i + kExitEvery);
                for (; i < block_end; i++) {
                    uint64_t const* const mask = m_masks[KmerUtils::BaseToInt(read[i])].data();
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
                    unsigned char carry = 0;
                    for (size_t w = 0; w < words; w++) {
                        uint64_t const x = v[w], match = mask[w];
                        unsigned long long sum;
                        carry = _addcarry_u64(carry, x, x & match, &sum);
                        v[w] = static_cast<uint64_t>(sum) | (x & ~match);
                    }
#else
                    uint64_t carry = 0;
                    for (size_t w = 0; w < words; w++) {
                        uint64_t const x = v[w], match = mask[w];
                        uint64_t const sum = x + (x & match);
                        uint64_t const with_carry = sum + carry;
                        carry = static_cast<uint64_t>(sum < x) | static_cast<uint64_t>(with_carry < sum);
                        v[w] = with_carry | (x & ~match);
                    }
#endif
                }
                m_lcs = ZerosBelow(m);
                if (!exits) continue;
                if (m_lcs > max_lcs) return true;
                size_t const rest = n - i;
                int64_t const most = rest >= m ? static_cast<int64_t>(m) : ZerosBelow(m - rest) + static_cast<int64_t>(rest);
                if (most <= max_lcs) return false;
            }
            return m_lcs > max_lcs;
        }

        std::array<std::vector<uint64_t>, 5> m_masks;
        std::vector<uint64_t> m_v;
        int64_t m_lcs = 0;
    };
}
