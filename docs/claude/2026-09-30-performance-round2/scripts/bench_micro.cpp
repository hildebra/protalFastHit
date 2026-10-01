// Micro-benchmarks for the SIMD and other candidates of the round-2 report. Each section times
// protal's current code (copied or included) against a candidate on read-sized inputs and checks
// that both give the same result.
//   g++ -O3 -march=x86-64-v3 -I<src>/SequenceUtils -I<src> bench_micro.cpp -o bench_micro && ./bench_micro
#include <immintrin.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

using Clock = std::chrono::steady_clock;

template <class F>
static double NsPerCall(F&& f, size_t calls) {
    auto t0 = Clock::now();
    for (size_t i = 0; i < calls; i++) f(i);
    auto t1 = Clock::now();
    return std::chrono::duration<double, std::nano>(t1 - t0).count() / static_cast<double>(calls);
}

static volatile uint64_t g_sink;

// ---------------------------------------------------------------- 1. clock reads (Benchmark::Start/Stop)
static void ClockCost() {
    double ns = NsPerCall([](size_t) { g_sink += Clock::now().time_since_epoch().count(); }, 20'000'000);
    std::printf("steady_clock::now(): %.1f ns per call\n", ns);
}

// ---------------------------------------------------------------- 2. reverse complement
static char Complement(char c) {  // as KmerUtils::Complement: A<->T, C<->G, anything else stays
    switch (c) {
        case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A';
        case 'a': return 't'; case 'c': return 'g'; case 'g': return 'c'; case 't': return 'a';
        default: return c;
    }
}
// As KmerUtils::ReverseComplement: by value, appended one character at a time.
static std::string RcCurrent(std::string forward) {
    std::string reverse = "";
    for (int i = static_cast<int>(forward.length()) - 1; i >= 0; i--) reverse += Complement(forward[i]);
    return reverse;
}
static void RcLut(char const* s, size_t n, char* d) {
    static const auto table = [] {
        struct T { char c[256]; } t;
        for (int i = 0; i < 256; i++) t.c[i] = Complement(static_cast<char>(i));
        return t;
    }();
    for (size_t i = 0; i < n; i++) d[i] = table.c[static_cast<unsigned char>(s[n - 1 - i])];
}
// 32 bases at a time: complement ACGT (upper case) with a byte shuffle keyed on bits 1..3 of the
// letter, reverse the bytes of the vector. Other characters (N, lower case, IUPAC) fall to the scalar path.
__attribute__((target("avx2"))) static void RcAvx2(char const* s, size_t n, char* d) {
    // A=0x41 C=0x43 G=0x47 T=0x54 ; index = (c & 0x0F) maps A->1 C->3 G->7 T->4 (all distinct)
    alignas(16) static const char lut16[16] = {0, 'T', 0, 'G', 'A', 0, 0, 'C', 0, 0, 0, 0, 0, 0, 0, 0};
    __m256i const lut = _mm256_broadcastsi128_si256(_mm_load_si128(reinterpret_cast<__m128i const*>(lut16)));
    __m256i const rev = _mm256_setr_epi8(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
                                         15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    __m256i const low = _mm256_set1_epi8(0x0F);
    size_t i = 0;
    for (; i + 32 <= n; i += 32) {
        __m256i v = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(s + n - 32 - i));
        __m256i c = _mm256_shuffle_epi8(lut, _mm256_and_si256(v, low));
        // ACGT map to themselves+complement; anything else would come out 0 or wrong: patch below
        __m256i ok = _mm256_or_si256(_mm256_or_si256(_mm256_cmpeq_epi8(v, _mm256_set1_epi8('A')), _mm256_cmpeq_epi8(v, _mm256_set1_epi8('C'))),
                                     _mm256_or_si256(_mm256_cmpeq_epi8(v, _mm256_set1_epi8('G')), _mm256_cmpeq_epi8(v, _mm256_set1_epi8('T'))));
        __m256i r = _mm256_blendv_epi8(v, c, ok);  // non-ACGT stay as they are (N -> N)
        r = _mm256_shuffle_epi8(r, rev);
        r = _mm256_permute2x128_si256(r, r, 0x01);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(d + i), r);
    }
    for (; i < n; i++) d[i] = Complement(s[n - 1 - i]);
}
static void ReverseComplementBench() {
    std::mt19937 rng(1);
    std::vector<std::string> reads(1024);
    for (auto& r : reads) {
        r.resize(150);
        for (auto& c : r) c = "ACGT"[rng() & 3];
        if (rng() % 8 == 0) r[rng() % 150] = 'N';
    }
    std::string out(150, ' '), out2(150, ' ');
    for (auto const& r : reads) {
        RcLut(r.data(), r.size(), out.data());
        RcAvx2(r.data(), r.size(), out2.data());
        if (out != RcCurrent(r) || out2 != out) { std::puts("reverse complement MISMATCH"); return; }
    }
    size_t const calls = 3'000'000;
    std::printf("reverse complement of 150 bp: current %.1f ns, table %.1f ns, AVX2 %.1f ns\n",
                NsPerCall([&](size_t i) { g_sink += RcCurrent(reads[i & 1023])[0]; }, calls),
                NsPerCall([&](size_t i) { RcLut(reads[i & 1023].data(), 150, out.data()); g_sink += out[0]; }, calls),
                NsPerCall([&](size_t i) { RcAvx2(reads[i & 1023].data(), 150, out.data()); g_sink += out[0]; }, calls));
}

// ---------------------------------------------------------------- 3. flex-block scan (KmerLookupSM::GetFromLookup)
static uint32_t Similarity(uint32_t a, uint32_t b) {
    return __builtin_popcount(((~(a ^ b) >> 1) & ~(a ^ b)) & 0x55555555u);
}
struct FlexResult { int max; int max_count; std::vector<uint32_t> hits; };
// As the loop in KmerLookupSM: similarity per cell into a vector, the maximum and its count, then the
// indices at the maximum.
static void FlexCurrent(uint32_t const* cells, size_t n, uint32_t key, std::vector<uint16_t>& flex_vector, FlexResult& r) {
    flex_vector.clear();
    int max = 0, max_count = 0;
    for (size_t i = 0; i < n; i++) {
        auto sim = Similarity(cells[i], key);
        flex_vector.emplace_back(sim);
        if (sim > max) { max = sim; max_count = 0; }
        max_count += (sim == static_cast<uint32_t>(max));
    }
    r.max = max; r.max_count = max_count; r.hits.clear();
    for (size_t i = 0; i < flex_vector.size(); i++)
        if (flex_vector[i] == max) r.hits.push_back(static_cast<uint32_t>(i));
}
// Byte-wise popcount of the equal-base mask, 8 cells at a time, then a second pass for the hits.
__attribute__((target("avx2"))) static inline __m256i SimilarityAvx2(__m256i cells, __m256i key) {
    __m256i x = _mm256_xor_si256(cells, key);
    __m256i eq = _mm256_andnot_si256(_mm256_or_si256(x, _mm256_srli_epi32(x, 1)), _mm256_set1_epi32(0x55555555));
    // popcount of eq (bits only at even positions): pairs of bits -> 0/1 each; sum with shuffles
    __m256i lut = _mm256_setr_epi8(0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4);
    __m256i low = _mm256_set1_epi8(0x0F);
    __m256i cnt = _mm256_add_epi8(_mm256_shuffle_epi8(lut, _mm256_and_si256(eq, low)),
                                  _mm256_shuffle_epi8(lut, _mm256_and_si256(_mm256_srli_epi16(eq, 4), low)));
    return _mm256_madd_epi16(_mm256_maddubs_epi16(cnt, _mm256_set1_epi8(1)), _mm256_set1_epi16(1));
}
__attribute__((target("avx2"))) static void FlexAvx2(uint32_t const* cells, size_t n, uint32_t key, FlexResult& r) {
    __m256i const k = _mm256_set1_epi32(static_cast<int>(key));
    __m256i vmax = _mm256_setzero_si256();
    size_t i = 0;
    int smax = 0;
    for (; i + 8 <= n; i += 8)
        vmax = _mm256_max_epi32(vmax, SimilarityAvx2(_mm256_loadu_si256(reinterpret_cast<__m256i const*>(cells + i)), k));
    alignas(32) int lanes[8];
    _mm256_store_si256(reinterpret_cast<__m256i*>(lanes), vmax);
    for (int l = 0; l < 8; l++) smax = std::max(smax, lanes[l]);
    for (size_t j = i; j < n; j++) smax = std::max<int>(smax, static_cast<int>(Similarity(cells[j], key)));
    r.max = smax; r.max_count = 0; r.hits.clear();
    __m256i const vm = _mm256_set1_epi32(smax);
    for (i = 0; i + 8 <= n; i += 8) {
        unsigned m = static_cast<unsigned>(_mm256_movemask_ps(_mm256_castsi256_ps(
            _mm256_cmpeq_epi32(SimilarityAvx2(_mm256_loadu_si256(reinterpret_cast<__m256i const*>(cells + i)), k), vm))));
        while (m) { r.hits.push_back(static_cast<uint32_t>(i + __builtin_ctz(m))); m &= m - 1; }
    }
    for (size_t j = i; j < n; j++) if (static_cast<int>(Similarity(cells[j], key)) == smax) r.hits.push_back(static_cast<uint32_t>(j));
    r.max_count = static_cast<int>(r.hits.size());
}
static void FlexBench() {
    std::mt19937 rng(2);
    for (size_t n : {3, 8, 16, 64, 256, 1024}) {
        std::vector<uint32_t> cells(n + 64);
        for (auto& c : cells) c = static_cast<uint32_t>(rng());
        std::vector<uint16_t> fv;
        FlexResult a, b;
        uint32_t key = cells[n / 2] ^ 0x00040010u;  // close to one cell
        FlexCurrent(cells.data(), n, key, fv, a);
        FlexAvx2(cells.data(), n, key, b);
        if (a.max != b.max || a.max_count != b.max_count || a.hits != b.hits) { std::printf("flex MISMATCH n=%zu\n", n); continue; }
        size_t calls = 20'000'000 / n + 100000;
        double cur = NsPerCall([&](size_t i) { FlexCurrent(cells.data(), n, key + static_cast<uint32_t>(i & 7), fv, a); g_sink += a.max; }, calls);
        double avx = NsPerCall([&](size_t i) { FlexAvx2(cells.data(), n, key + static_cast<uint32_t>(i & 7), b); g_sink += b.max; }, calls);
        std::printf("flex block of %4zu cells: current %8.1f ns (%.2f/cell), AVX2 %8.1f ns (%.2f/cell)\n", n, cur, cur / n, avx, avx / n);
    }
}

// ---------------------------------------------------------------- 4. exact-match run length (seed extension, link verification)
static size_t MatchRunScalar(char const* a, char const* b, size_t n) {
    size_t i = 0;
    while (i < n && a[i] == b[i]) i++;
    return i;
}
__attribute__((target("avx2"))) static size_t MatchRunAvx2(char const* a, char const* b, size_t n) {
    size_t i = 0;
    for (; i + 32 <= n; i += 32) {
        unsigned m = static_cast<unsigned>(_mm256_movemask_epi8(_mm256_cmpeq_epi8(
            _mm256_loadu_si256(reinterpret_cast<__m256i const*>(a + i)), _mm256_loadu_si256(reinterpret_cast<__m256i const*>(b + i)))));
        if (m != 0xFFFFFFFFu) return i + static_cast<size_t>(__builtin_ctz(~m));
    }
    while (i < n && a[i] == b[i]) i++;
    return i;
}
static void MatchRunBench() {
    std::mt19937 rng(3);
    std::string a(300, 'A'), b;
    for (auto& c : a) c = "ACGT"[rng() & 3];
    for (size_t run : {5, 30, 100, 150}) {
        b = a;
        if (run < 150) b[run] = b[run] == 'A' ? 'C' : 'A';
        if (MatchRunScalar(a.data(), b.data(), 150) != MatchRunAvx2(a.data(), b.data(), 150)) { std::puts("match run MISMATCH"); continue; }
        size_t calls = 20'000'000;
        double s = NsPerCall([&](size_t) { g_sink += MatchRunScalar(a.data(), b.data(), 150); }, calls);
        double v = NsPerCall([&](size_t) { g_sink += MatchRunAvx2(a.data(), b.data(), 150); }, calls);
        std::printf("match run %3zu of 150: scalar %.1f ns, AVX2 %.1f ns\n", run, s, v);
    }
}

int main() {
    ClockCost();
    ReverseComplementBench();
    FlexBench();
    MatchRunBench();
    return 0;
}
