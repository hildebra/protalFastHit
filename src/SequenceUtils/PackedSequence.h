//
// Gene sequences at two bits per base.
//
// The reference genes are the second largest part of protal's memory (~17 GB of the 59 of the full
// GTDB r226 database at one byte per base). They are held packed, four bases to a byte, and a
// Gene::Sequence() call decodes the gene into a GeneSequence (a buffer on the stack, for genes up to
// 4096 bases).
//
// Coding: A 0, C 1, G 2, T 3; base i of a gene is in byte i/4 at bits 2*(i%4), every gene starts on
// a byte. Anything that is not A, C, G or T has no code of its own and is stored as a base:
// N as A, and an IUPAC ambiguity code as the first base it can stand for in the order A C G T
// (R, W, M, D, H, V as A; Y, S, B as C; K as G). Lowercase letters are read as uppercase; any other
// character is stored as A. So decoding a gene gives back its sequence up to these codes.
//
// Packing and unpacking use AVX2 where the CPU has it (chosen at run time, like the syncmer scan) and
// a table otherwise; both give the same bytes.
//
#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <ostream>
#include <string>
#include <string_view>

// x86 with GCC or Clang: packing and unpacking can use AVX2, compiled for it by function attribute.
#if (defined(__x86_64__) || defined(__i386__)) && (defined(__GNUC__) || defined(__clang__))
#define PROTAL_PACKED_AVX2 1
#include <immintrin.h>
#endif

namespace protal::packed {

    // Bytes that n bases take.
    constexpr size_t Bytes(size_t bases) { return (bases + 3) / 4; }

    // The 2-bit code of each input character (see the file comment).
    inline constexpr std::array<uint8_t, 256> kCode = [] {
        std::array<uint8_t, 256> code{};  // everything else is A
        auto set = [&](char upper, uint8_t value) {
            code[static_cast<uint8_t>(upper)] = value;
            code[static_cast<uint8_t>(upper | 0x20)] = value;
        };
        set('C', 1); set('G', 2); set('T', 3);
        set('Y', 1); set('S', 1); set('B', 1);  // C or T; C or G; C, G or T
        set('K', 2);                            // G or T
        return code;
    }();

    inline constexpr char kBases[4] = { 'A', 'C', 'G', 'T' };

    // The four characters of a packed byte, as a little-endian word.
    inline constexpr std::array<uint32_t, 256> kWord = [] {
        std::array<uint32_t, 256> word{};
        for (uint32_t b = 0; b < 256; b++) {
            for (uint32_t j = 0; j < 4; j++) word[b] |= static_cast<uint32_t>(kBases[(b >> (2 * j)) & 3]) << (8 * j);
        }
        return word;
    }();

    inline uint8_t PackByte(const char* s, size_t bases) {
        uint8_t b = 0;
        for (size_t j = 0; j < bases; j++) b |= static_cast<uint8_t>(kCode[static_cast<uint8_t>(s[j])] << (2 * j));
        return b;
    }

    // Packs n bases into Bytes(n) bytes at d (the last byte's unused bits are 0).
    inline void PackScalar(const char* s, size_t n, uint8_t* d) {
        size_t i = 0, o = 0;
        for (; i + 4 <= n; i += 4, o++) d[o] = PackByte(s + i, 4);
        if (i < n) d[o] = PackByte(s + i, n - i);
    }

    // Unpacks n bases (from Bytes(n) bytes) into n characters at d.
    inline void UnpackScalar(const uint8_t* s, size_t n, char* d) {
        size_t i = 0, o = 0;
        for (; i + 4 <= n; i += 4, o++) std::memcpy(d + i, &kWord[s[o]], 4);
        if (i < n) std::memcpy(d + i, &kWord[s[o]], n - i);
    }

#ifdef PROTAL_PACKED_AVX2
    // 32 characters to 8 bytes. A block of only A, C, G, T (either case) is coded in vector
    // registers; a block with anything else goes through the table.
    __attribute__((target("avx2"))) inline void PackAvx2(const char* s, size_t n, uint8_t* d) {
        // The low nibble tells the four letters apart: A 0x41 -> 1, C 0x43 -> 3, G 0x47 -> 7, T 0x54 -> 4;
        // the table holds their codes at those positions.
        const __m256i table = _mm256_setr_epi8(0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 1, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0);
        const __m256i upper = _mm256_set1_epi8(static_cast<char>(0xDF));
        const __m256i low_nibble = _mm256_set1_epi8(0x0F);
        const __m256i a = _mm256_set1_epi8('A'), c = _mm256_set1_epi8('C'), g = _mm256_set1_epi8('G'), t = _mm256_set1_epi8('T');
        const __m256i pairs = _mm256_set1_epi16(0x0401);      // byte 0 * 1 + byte 1 * 4
        const __m256i quads = _mm256_set1_epi32(0x00100001);  // word 0 * 1 + word 1 * 16
        const __m256i first_bytes = _mm256_setr_epi8(0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                                     0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
        size_t i = 0, o = 0;
        for (; i + 32 <= n; i += 32, o += 8) {
            __m256i v = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(s + i));
            __m256i u = _mm256_and_si256(v, upper);
            __m256i is_base = _mm256_or_si256(_mm256_or_si256(_mm256_cmpeq_epi8(u, a), _mm256_cmpeq_epi8(u, c)),
                                              _mm256_or_si256(_mm256_cmpeq_epi8(u, g), _mm256_cmpeq_epi8(u, t)));
            if (_mm256_movemask_epi8(is_base) != -1) {
                PackScalar(s + i, 32, d + o);
                continue;
            }
            __m256i code = _mm256_shuffle_epi8(table, _mm256_and_si256(u, low_nibble));
            __m256i p1 = _mm256_maddubs_epi16(code, pairs);   // c0 + 4 c1, per pair
            __m256i p2 = _mm256_madd_epi16(p1, quads);        // c0 + 4 c1 + 16 c2 + 64 c3, per 4 bases
            __m256i bytes = _mm256_shuffle_epi8(p2, first_bytes);
            __m128i lo = _mm256_castsi256_si128(bytes), hi = _mm256_extracti128_si256(bytes, 1);
            _mm_storel_epi64(reinterpret_cast<__m128i*>(d + o), _mm_unpacklo_epi32(lo, hi));
        }
        if (i < n) PackScalar(s + i, n - i, d + o);
    }

    // 32 bytes to 128 characters.
    __attribute__((target("avx2"))) inline void Unpack128Avx2(const uint8_t* s, char* d) {
        const __m256i letters = _mm256_setr_epi8('A', 'C', 'G', 'T', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 'A', 'C', 'G', 'T', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        const __m256i three = _mm256_set1_epi8(3);
        __m256i x = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(s));
        // The bases at 2-bit position 0, 1, 2, 3 of every byte, as letters.
        __m256i c0 = _mm256_shuffle_epi8(letters, _mm256_and_si256(x, three));
        __m256i c1 = _mm256_shuffle_epi8(letters, _mm256_and_si256(_mm256_srli_epi16(x, 2), three));
        __m256i c2 = _mm256_shuffle_epi8(letters, _mm256_and_si256(_mm256_srli_epi16(x, 4), three));
        __m256i c3 = _mm256_shuffle_epi8(letters, _mm256_and_si256(_mm256_srli_epi16(x, 6), three));
        // Interleave them: byte b gives c0[b] c1[b] c2[b] c3[b]. Within a 128-bit lane, r0..r3 hold the
        // bases of its bytes 0-3, 4-7, 8-11 and 12-15.
        __m256i u01l = _mm256_unpacklo_epi8(c0, c1), u01h = _mm256_unpackhi_epi8(c0, c1);
        __m256i u23l = _mm256_unpacklo_epi8(c2, c3), u23h = _mm256_unpackhi_epi8(c2, c3);
        __m256i r0 = _mm256_unpacklo_epi16(u01l, u23l), r1 = _mm256_unpackhi_epi16(u01l, u23l);
        __m256i r2 = _mm256_unpacklo_epi16(u01h, u23h), r3 = _mm256_unpackhi_epi16(u01h, u23h);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(d), _mm256_permute2x128_si256(r0, r1, 0x20));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(d + 32), _mm256_permute2x128_si256(r2, r3, 0x20));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(d + 64), _mm256_permute2x128_si256(r0, r1, 0x31));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(d + 96), _mm256_permute2x128_si256(r2, r3, 0x31));
    }

    __attribute__((target("avx2"))) inline void UnpackAvx2(const uint8_t* s, size_t n, char* d) {
        size_t i = 0;
        for (; i + 128 <= n; i += 128) Unpack128Avx2(s + i / 4, d + i);
        if (i < n) {
            // Fewer than 128 bases left: through a padded copy, so that no byte past the gene is read
            // or written.
            alignas(32) uint8_t in[32] = {};
            alignas(32) char out[128];
            std::memcpy(in, s + i / 4, Bytes(n - i));
            Unpack128Avx2(in, out);
            std::memcpy(d + i, out, n - i);
        }
    }
#endif

    inline bool CpuHasAvx2() {
#ifdef PROTAL_PACKED_AVX2
        return __builtin_cpu_supports("avx2");
#else
        return false;
#endif
    }

    // Whether packing and unpacking use AVX2. On where the CPU has it; tests switch it to compare both.
    inline std::atomic<bool>& Avx2Enabled() {
        static std::atomic<bool> enabled{ CpuHasAvx2() };
        return enabled;
    }

    inline void UseAvx2(bool use) { Avx2Enabled().store(use && CpuHasAvx2(), std::memory_order_relaxed); }

    inline void Pack(const char* s, size_t n, uint8_t* d) {
#ifdef PROTAL_PACKED_AVX2
        if (Avx2Enabled().load(std::memory_order_relaxed)) return PackAvx2(s, n, d);
#endif
        PackScalar(s, n, d);
    }

    inline void Unpack(const uint8_t* s, size_t n, char* d) {
#ifdef PROTAL_PACKED_AVX2
        if (Avx2Enabled().load(std::memory_order_relaxed)) return UnpackAvx2(s, n, d);
#endif
        UnpackScalar(s, n, d);
    }

    // Packs `n` bases, which are bases first .. first + n - 1 of the gene whose packed bytes start at
    // `gene`, into it. The bytes must be zero before (a fresh calloc): the byte at either end of the
    // range can be shared with the neighbouring range of another thread and is set with an atomic or;
    // the bytes inside are written plainly. So threads can fill disjoint ranges of one gene.
    inline void PackInto(uint8_t* gene, size_t first, const char* s, size_t n) {
        if (n == 0) return;
        uint8_t* p = gene + first / 4;
        if (size_t const head = first % 4) {
            size_t const take = std::min(n, 4 - head);
            std::atomic_ref<uint8_t>(*p).fetch_or(static_cast<uint8_t>(PackByte(s, take) << (2 * head)), std::memory_order_relaxed);
            s += take; n -= take; p++;
            if (n == 0) return;
        }
        size_t const whole = n / 4;
        Pack(s, whole * 4, p);
        if (size_t const tail = n % 4) {
            std::atomic_ref<uint8_t>(p[whole]).fetch_or(PackByte(s + whole * 4, tail), std::memory_order_relaxed);
        }
    }
}

namespace protal {
    // A decoded gene: its bases as characters in a buffer of its own (on the stack for genes up to
    // kInline bases), or a view of a sequence that lives elsewhere (a test's string). It lives as long
    // as the object, which is neither copied nor moved: keep it in a variable while its view is used
    // (`auto const reference = gene.Sequence();`), a view taken from a temporary dangles.
    class GeneSequence {
    public:
        static constexpr size_t kInline = 4096;

        GeneSequence(const uint8_t* packed, size_t bases) {
            char* out = m_buffer;
            if (bases > kInline) {
                m_heap = std::make_unique_for_overwrite<char[]>(bases);
                out = m_heap.get();
            }
            if (bases > 0) packed::Unpack(packed, bases, out);
            m_view = std::string_view(out, bases);
        }

        explicit GeneSequence(std::string_view external) : m_view(external) {}

        GeneSequence(GeneSequence const&) = delete;
        GeneSequence& operator=(GeneSequence const&) = delete;

        operator std::string_view() const { return m_view; }
        std::string_view View() const { return m_view; }
        const char* data() const { return m_view.data(); }
        size_t size() const { return m_view.size(); }
        size_t length() const { return m_view.size(); }
        bool empty() const { return m_view.empty(); }
        char operator[](size_t i) const { return m_view[i]; }
        auto begin() const { return m_view.begin(); }
        auto end() const { return m_view.end(); }
        std::string_view substr(size_t pos, size_t n = std::string_view::npos) const { return m_view.substr(pos, n); }

        friend bool operator==(GeneSequence const& a, std::string_view b) { return a.m_view == b; }
        friend bool operator==(std::string_view a, GeneSequence const& b) { return a == b.m_view; }
        friend std::ostream& operator<<(std::ostream& os, GeneSequence const& s) { return os << s.m_view; }

    private:
        std::string_view m_view;
        std::unique_ptr<char[]> m_heap;
        char m_buffer[kInline];  // not initialised
    };
}
