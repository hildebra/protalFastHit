// Syncmer extraction: protal's SimpleKmerHandler<ClosedSyncmer> against a version that computes each
// read's 7-mers once per strand. Checks the k-mer lists are identical read by read and times both.
//   bench_syncmer READS.fq [READS2.fq ...]
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>
#include "SequenceUtils/KmerIterator.h"

using namespace protal;

// Same syncmers as SimpleKmerHandler{31, 15, ClosedSyncmer{15, 7, 2, full}}: for each 31-mer window,
// the canonical 15-mer core (the middle 15 bases of the 31-mer or of its reverse complement,
// whichever is smaller; the reverse on a tie) is a closed syncmer if the first minimum of its nine
// 7-mers is the 3rd or the 7th.
class FastSyncmers {
    static constexpr int k = 31, s = 7;
    static constexpr uint64_t kmask = (uint64_t{1} << (2 * k)) - 1, cmask = (uint64_t{1} << 30) - 1;
    std::array<uint8_t, 256> m_fwd{}, m_rev{};
    std::vector<uint16_t> m_f7, m_r7;
    uint32_t m_smask;

public:
    explicit FastSyncmers(bool full_smer_mask) : m_smask(ClosedSyncmer::SmerMask(s, full_smer_mask)) {
        for (int c = 0; c < 256; c++) {
            m_fwd[c] = static_cast<uint8_t>(KmerUtils::BaseToInt(static_cast<char>(c), 0));
            m_rev[c] = static_cast<uint8_t>(KmerUtils::BaseToIntC(static_cast<char>(c), 0));
        }
    }

    void operator()(std::string_view seq, KmerList& list) {
        list.clear();
        size_t const n = seq.size();
        if (n < static_cast<size_t>(k)) return;
        // f7[j]: the 7-mer at j, first base in the high bits; r7[j]: its reverse complement, encoded
        // as the reverse k-mer is (so that N reads as 0 on both strands, as in SimpleKmerHandler).
        m_f7.resize(n - s + 1);
        m_r7.resize(n - s + 1);
        uint32_t f = 0, r = 0;
        for (size_t i = 0; i < n; i++) {
            auto const c = static_cast<unsigned char>(seq[i]);
            f = ((f << 2) | m_fwd[c]) & 0x3FFF;
            r = (r >> 2) | (uint32_t{m_rev[c]} << 12);
            if (i + 1 >= static_cast<size_t>(s)) {
                m_f7[i + 1 - s] = static_cast<uint16_t>(f & m_smask);
                m_r7[i + 1 - s] = static_cast<uint16_t>(r & m_smask);
            }
        }
        uint64_t kf = 0, kr = 0;
        for (size_t i = 0; i < n; i++) {
            auto const c = static_cast<unsigned char>(seq[i]);
            kf = ((kf << 2) | m_fwd[c]) & kmask;
            kr = (kr >> 2) | (uint64_t{m_rev[c]} << (2 * (k - 1)));
            if (i + 1 < static_cast<size_t>(k)) continue;
            size_t const p = i + 1 - k;  // window start
            uint64_t const core_f = (kf >> 16) & cmask, core_r = (kr >> 16) & cmask;
            bool const forward = core_f < core_r;
            // The core's 7-mers in its own order: forward f7[p+8+j]; reverse r7[p+16-j].
            uint16_t const* v = forward ? &m_f7[p + 8] : &m_r7[p + 8];
            int best = 0;
            uint16_t min = forward ? v[0] : v[8];
            for (int j = 1; j < 9; j++) {
                uint16_t const x = forward ? v[j] : v[8 - j];
                if (x < min) { min = x; best = j; }
            }
            if (best == 2 || best == 6) list.emplace_back(forward ? kf : kr, p);
        }
    }
};

// As FastSyncmers, without data-dependent branches: each 7-mer is packed with its index in the core
// (value << 4 | index), so the smallest packed value is the first minimum; both strands are
// evaluated and the canonical one selected.
class BranchFreeSyncmers {
    static constexpr int k = 31, s = 7;
    static constexpr uint64_t kmask = (uint64_t{1} << (2 * k)) - 1, cmask = (uint64_t{1} << 30) - 1;
    std::array<uint8_t, 256> m_fwd{}, m_rev{};
    std::vector<uint32_t> m_f7, m_r7;
    uint32_t m_smask;

public:
    explicit BranchFreeSyncmers(bool full_smer_mask) : m_smask(ClosedSyncmer::SmerMask(s, full_smer_mask)) {
        for (int c = 0; c < 256; c++) {
            m_fwd[c] = static_cast<uint8_t>(KmerUtils::BaseToInt(static_cast<char>(c), 0));
            m_rev[c] = static_cast<uint8_t>(KmerUtils::BaseToIntC(static_cast<char>(c), 0));
        }
    }

    void operator()(std::string_view seq, KmerList& list) {
        list.clear();
        size_t const n = seq.size();
        if (n < static_cast<size_t>(k)) return;
        m_f7.resize(n - s + 1);
        m_r7.resize(n - s + 1);
        uint32_t f = 0, r = 0;
        for (size_t i = 0; i < n; i++) {
            auto const c = static_cast<unsigned char>(seq[i]);
            f = ((f << 2) | m_fwd[c]) & 0x3FFF;
            r = (r >> 2) | (uint32_t{m_rev[c]} << 12);
            if (i + 1 >= static_cast<size_t>(s)) {
                m_f7[i + 1 - s] = (f & m_smask) << 4;
                m_r7[i + 1 - s] = (r & m_smask) << 4;
            }
        }
        uint64_t kf = 0, kr = 0;
        for (size_t i = 0; i < n; i++) {
            auto const c = static_cast<unsigned char>(seq[i]);
            kf = ((kf << 2) | m_fwd[c]) & kmask;
            kr = (kr >> 2) | (uint64_t{m_rev[c]} << (2 * (k - 1)));
            if (i + 1 < static_cast<size_t>(k)) continue;
            size_t const p = i + 1 - k;
            uint64_t const core_f = (kf >> 16) & cmask, core_r = (kr >> 16) & cmask;
            uint32_t const* fv = &m_f7[p + 8];
            uint32_t const* rv = &m_r7[p + 8];
            uint32_t fmin = fv[0], rmin = rv[8];
            for (uint32_t j = 1; j < 9; j++) {
                fmin = std::min(fmin, fv[j] | j);
                rmin = std::min(rmin, rv[8 - j] | j);
            }
            uint32_t const best = (core_f < core_r ? fmin : rmin) & 0xF;
            if (best == 2 || best == 6) list.emplace_back(core_f < core_r ? kf : kr, p);
        }
    }
};

#ifdef __AVX2__
#include <immintrin.h>
// As BranchFreeSyncmers, 8 windows at a time: a scalar pass stores the s-mers and each window's
// cores and k-mers, then AVX2 takes the minima of 8 windows per step.
class Avx2Syncmers {
    static constexpr int k = 31, s = 7, count = 9, t = 2, core = 8;
    static constexpr uint64_t kmask = (uint64_t{1} << (2 * k)) - 1, cmask = (uint64_t{1} << 30) - 1;
    std::array<uint8_t, 256> m_fwd{}, m_rev{};
    std::vector<uint32_t> m_f7, m_r7, m_cf, m_cr;
    std::vector<uint64_t> m_kf, m_kr;
    uint32_t m_smask;

public:
    explicit Avx2Syncmers(bool full_smer_mask) : m_smask(ClosedSyncmer::SmerMask(s, full_smer_mask)) {
        for (int c = 0; c < 256; c++) {
            m_fwd[c] = static_cast<uint8_t>(KmerUtils::BaseToInt(static_cast<char>(c), 0));
            m_rev[c] = static_cast<uint8_t>(KmerUtils::BaseToIntC(static_cast<char>(c), 0));
        }
    }

    void operator()(std::string_view seq, KmerList& list) {
        list.clear();
        size_t const n = seq.size();
        if (n < static_cast<size_t>(k)) return;
        size_t const windows = n - k + 1;
        m_f7.resize(n - s + 1); m_r7.resize(n - s + 1);
        m_cf.resize(windows); m_cr.resize(windows); m_kf.resize(windows); m_kr.resize(windows);
        uint32_t f = 0, r = 0;
        uint64_t kf = 0, kr = 0;
        for (size_t i = 0; i < n; i++) {
            auto const c = static_cast<unsigned char>(seq[i]);
            f = ((f << 2) | m_fwd[c]) & 0x3FFF;
            r = (r >> 2) | (uint32_t{m_rev[c]} << 12);
            kf = ((kf << 2) | m_fwd[c]) & kmask;
            kr = (kr >> 2) | (uint64_t{m_rev[c]} << (2 * (k - 1)));
            if (i + 1 >= static_cast<size_t>(s)) {
                m_f7[i + 1 - s] = (f & m_smask) << 4;
                m_r7[i + 1 - s] = (r & m_smask) << 4;
            }
            if (i + 1 >= static_cast<size_t>(k)) {
                size_t const p = i + 1 - k;
                m_kf[p] = kf; m_kr[p] = kr;
                m_cf[p] = static_cast<uint32_t>((kf >> 16) & cmask);
                m_cr[p] = static_cast<uint32_t>((kr >> 16) & cmask);
            }
        }
        __m256i const t1 = _mm256_set1_epi32(t), t2 = _mm256_set1_epi32(count - 1 - t), low = _mm256_set1_epi32(0xF);
        size_t p = 0;
        for (; p + 8 <= windows; p += 8) {
            __m256i fmin = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(&m_f7[p + core]));
            __m256i rmin = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(&m_r7[p + core + count - 1]));
            for (int j = 1; j < count; j++) {
                __m256i const idx = _mm256_set1_epi32(j);
                fmin = _mm256_min_epu32(fmin, _mm256_or_si256(_mm256_loadu_si256(reinterpret_cast<__m256i const*>(&m_f7[p + core + j])), idx));
                rmin = _mm256_min_epu32(rmin, _mm256_or_si256(_mm256_loadu_si256(reinterpret_cast<__m256i const*>(&m_r7[p + core + count - 1 - j])), idx));
            }
            __m256i const cf = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(&m_cf[p]));
            __m256i const cr = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(&m_cr[p]));
            __m256i const forward = _mm256_cmpgt_epi32(cr, cf);           // core_f < core_r
            __m256i const first = _mm256_and_si256(_mm256_blendv_epi8(rmin, fmin, forward), low);
            __m256i const hit = _mm256_or_si256(_mm256_cmpeq_epi32(first, t1), _mm256_cmpeq_epi32(first, t2));
            unsigned bits = static_cast<unsigned>(_mm256_movemask_ps(_mm256_castsi256_ps(hit)));
            unsigned const fbits = static_cast<unsigned>(_mm256_movemask_ps(_mm256_castsi256_ps(forward)));
            while (bits) {
                unsigned const l = static_cast<unsigned>(__builtin_ctz(bits));
                list.emplace_back((fbits >> l) & 1 ? m_kf[p + l] : m_kr[p + l], p + l);
                bits &= bits - 1;
            }
        }
        for (; p < windows; p++) {  // the last < 8 windows
            uint32_t fmin = m_f7[p + core], rmin = m_r7[p + core + count - 1];
            for (uint32_t j = 1; j < count; j++) {
                fmin = std::min(fmin, m_f7[p + core + j] | j);
                rmin = std::min(rmin, m_r7[p + core + count - 1 - j] | j);
            }
            bool const fw = m_cf[p] < m_cr[p];
            uint32_t const first = (fw ? fmin : rmin) & 0xF;
            if (first == t || first == count - 1 - t) list.emplace_back(fw ? m_kf[p] : m_kr[p], p);
        }
    }
};
#endif

int main(int argc, char** argv) {
    std::vector<std::string> reads;
    for (int a = 1; a < argc; a++) {
        std::ifstream in(argv[a]);
        std::string line;
        for (size_t l = 0; std::getline(in, line); l++) if (l % 4 == 1) reads.push_back(line);
    }
    std::printf("%zu reads\n", reads.size());

    ClosedSyncmer minimizer{15, 7, 2, true};
    SimpleKmerHandler<ClosedSyncmer> current{31, 15, minimizer};
    FastSyncmers fast{true};
    BranchFreeSyncmers branchfree{true};
    KmerList a, b, c;

    // "window": the window-by-window definition; "handler": SimpleKmerHandler as built (the scan
    // where the source has it); "avx2": 8 windows at a time.
    auto window = [&current](std::string_view seq, KmerList& list) { current.WindowByWindow(seq, list); };
#ifdef __AVX2__
    Avx2Syncmers avx2{true};
    KmerList d;
#endif

    // Equality, read by read (and the reads reversed-complemented, and with Ns).
    size_t mismatched = 0, mismatched_bf = 0, mismatched_avx2 = 0, total = 0;
    for (auto const& read : reads) {
        for (int variant = 0; variant < 3; variant++) {
            std::string seq = read;
            if (variant == 1) seq = KmerUtils::ReverseComplement(seq);
            if (variant == 2) for (size_t i = 0; i < seq.size(); i += 37) seq[i] = 'N';
            window(std::string_view(seq), a);
            current(std::string_view(seq), b);
            branchfree(std::string_view(seq), c);
            mismatched += a != b;
            mismatched_bf += a != c;
#ifdef __AVX2__
            avx2(std::string_view(seq), d);
            mismatched_avx2 += a != d;
#endif
            total += a.size();
        }
    }
    std::printf("syncmers %zu, reads with lists other than the definition's: handler %zu, branch-free %zu, avx2 %zu\n",
                total, mismatched, mismatched_bf, mismatched_avx2);

    auto time = [&](auto& handler, KmerList& list) {
        size_t sum = 0;
        auto t0 = std::chrono::steady_clock::now();
        for (auto const& read : reads) { handler(std::string_view(read), list); sum += list.size(); }
        return std::pair{ std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(), sum };
    };
    for (int rep = 0; rep < 3; rep++) {
        auto [win, s1] = time(window, a);
        auto [cur, s2] = time(current, b);
        auto [bf, s3] = time(branchfree, c);
        double av = 0; size_t s4 = 0;
#ifdef __AVX2__
        std::tie(av, s4) = time(avx2, d);
#endif
        std::printf("ns/read: window %.0f, handler %.0f (%.1fx), branch-free %.0f (%.1fx), avx2 %.0f (%.1fx)  [%zu %zu %zu %zu]\n",
                    1e9 * win / reads.size(), 1e9 * cur / reads.size(), win / cur, 1e9 * bf / reads.size(), win / bf,
                    1e9 * av / reads.size(), av > 0 ? win / av : 0.0, s1, s2, s3, s4);
    }
#ifdef HAVE_USE_AVX2
    // The handler's own evaluations: one by one, and 8 at a time with AVX2 (runtime dispatch).
    SimpleKmerHandler<ClosedSyncmer> scalar{31, 15, minimizer}, vector{31, 15, minimizer};
    scalar.UseAvx2(false);
    std::printf("handler uses AVX2 by default: %d\n", vector.UsesAvx2());
    size_t diff = 0;
    for (auto const& read : reads) { scalar(std::string_view(read), a); vector(std::string_view(read), b); window(std::string_view(read), c); diff += (a != c) + (b != c); }
    std::printf("handler lists other than the definition's: %zu\n", diff);
    for (int rep = 0; rep < 3; rep++) {
        auto [win, s1] = time(window, a);
        auto [sc, s2] = time(scalar, b);
        auto [ve, s3] = time(vector, c);
        std::printf("ns/read: window %.0f, handler one by one %.0f (%.1fx), handler AVX2 %.0f (%.1fx)  [%zu %zu %zu]\n",
                    1e9 * win / reads.size(), 1e9 * sc / reads.size(), win / sc, 1e9 * ve / reads.size(), win / ve, s1, s2, s3);
    }
#endif
}
