// Feasibility prototype: parallel inflate of ONE gzip member with ISA-L alone, by decoding pieces that
// start at deflate block boundaries with "labelled" fake 32 KiB windows and resolving the labels later.
//
// Modes:
//   isalpar bounds <in.gz> <out.tsv>                    sequential pass recording every deflate block start
//   isalpar seq    <in.gz> [chunkKiB]                    sequential ISA-L baselines
//   isalpar cost   <in.gz> <bounds.tsv> <pieceMiB> [fq] [stats.tsv]   per-step CPU costs (1 thread) + markers
//   isalpar par    <in.gz> <bounds.tsv> <pieceMiB> <threads> <sched> <variant> [gap] [fq]
//        sched: spec (all threads speculate, main thread resolves in order)
//               lead (main thread decodes in order with the true window, threads-1 speculate ahead)
//        variant: 2 (two passes, injective pair labels) or 3 (three passes, labels from the brief)
#include <isa-l/igzip_lib.h>
#include <isa-l/crc.h>

#include <immintrin.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

static constexpr uint32_t WIN = 32768;

static double now()
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

static double cpuSeconds()
{
    rusage ru{};
    getrusage(RUSAGE_SELF, &ru);
    return ru.ru_utime.tv_sec + ru.ru_utime.tv_usec * 1e-6 + ru.ru_stime.tv_sec + ru.ru_stime.tv_usec * 1e-6;
}

static long maxRssMiB()
{
    rusage ru{};
    getrusage(RUSAGE_SELF, &ru);
    return ru.ru_maxrss / 1024;
}

struct Mapped
{
    const uint8_t* p = nullptr;
    size_t n = 0;
};

static Mapped mapFile(const char* path, bool populate)
{
    int fd = open(path, O_RDONLY);
    if (fd < 0) { perror(path); exit(1); }
    struct stat sb{};
    fstat(fd, &sb);
    void* m = mmap(nullptr, sb.st_size, PROT_READ, MAP_PRIVATE | (populate ? MAP_POPULATE : 0), fd, 0);
    if (m == MAP_FAILED) { perror("mmap"); exit(1); }
    close(fd);
    return { static_cast<const uint8_t*>(m), static_cast<size_t>(sb.st_size) };
}

struct Buf
{
    uint8_t* p = nullptr;
    size_t size = 0;
    size_t cap = 0;
    void reserve(size_t n)
    {
        if (n <= cap) return;
        p = static_cast<uint8_t*>(realloc(p, n));
        if (!p) { fprintf(stderr, "out of memory\n"); exit(1); }
        cap = n;
    }
    ~Buf() { free(p); }
};

// ---------------------------------------------------------------------------------------------
// Labels.
// Variant 3 (the brief): L1(p) = 0x80 | (p & 0x7f), L2(p) = 0x80 | ((p >> 7) & 0x7f), L3(p) = p >> 14.
// Variant 2 (injective pair): A(p) = p & 0xff, B(p) = (A(p) + 1 + (p >> 8)) & 0xff, so A != B always and
//   p = A | (((B - A - 1) & 0xff) << 8).
static uint8_t L1[WIN], L2[WIN], L3[WIN], LA[WIN], LB[WIN];

static void initLabels()
{
    for (uint32_t p = 0; p < WIN; ++p) {
        L1[p] = uint8_t(0x80 | (p & 0x7f));
        L2[p] = uint8_t(0x80 | ((p >> 7) & 0x7f));
        L3[p] = uint8_t(p >> 14);
        LA[p] = uint8_t(p & 0xff);
        LB[p] = uint8_t((LA[p] + 1 + (p >> 8)) & 0xff);
    }
}

// ---------------------------------------------------------------------------------------------
static inline uint64_t bitPos(const inflate_state& st, const uint8_t* base)
{
    return uint64_t(st.next_in - base) * 8 - uint64_t(st.read_in_length);
}

static std::runtime_error err(const std::string& what) { return std::runtime_error(what); }

// Decodes raw deflate from block start startBit up to block start endBit (UINT64_MAX: through the final
// block) into out (out.size = decoded bytes). dict: 32 KiB window or nullptr.
static void inflateRange(const Mapped& in, uint64_t startBit, uint64_t endBit, const uint8_t* dict, Buf& out,
                         inflate_state& st, size_t sizeHint)
{
    isal_inflate_init(&st);
    st.crc_flag = ISAL_DEFLATE;
    st.points_to_stop_at = ISAL_STOPPING_POINT_NONE;
    if (dict && isal_inflate_set_dict(&st, dict, WIN) != COMP_OK) throw err("isal_inflate_set_dict failed");

    size_t byte = startBit >> 3;
    const unsigned sh = startBit & 7;
    if (sh) {
        st.read_in = uint64_t(in.p[byte] >> sh);
        st.read_in_length = int32_t(8 - sh);
        ++byte;
    }
    st.next_in = const_cast<uint8_t*>(in.p + byte);
    const bool toEnd = endBit == UINT64_MAX;
    const size_t endByte = toEnd ? in.n : size_t(endBit >> 3);
    constexpr size_t MARGIN = 64;
    const size_t p1End = toEnd ? in.n : (endByte > byte + MARGIN ? endByte - MARGIN : byte);
    st.avail_in = uint32_t(p1End - byte);

    out.size = 0;
    out.reserve(sizeHint);
    auto prep = [&] {
        if (out.cap - out.size < (size_t(1) << 20)) out.reserve(std::max(out.cap * 2, out.size + (size_t(4) << 20)));
        st.next_out = out.p + out.size;
        st.avail_out = uint32_t(std::min<size_t>(out.cap - out.size, size_t(1) << 30));
    };
    auto call = [&](bool& progressed) {
        prep();
        const auto oldIn = st.avail_in;
        const auto oldOut = st.next_out;
        const auto oldRil = st.read_in_length;
        const int ret = isal_inflate(&st);
        out.size = size_t(st.next_out - out.p);
        if (ret < 0) throw err("isal_inflate error " + std::to_string(ret) + " at bit " + std::to_string(startBit));
        progressed = oldIn != st.avail_in || oldOut != st.next_out || oldRil != st.read_in_length;
    };

    // Phase 1: no stopping points; the input ends MARGIN bytes before the piece end.
    bool progressed = false;
    if (toEnd) {
        for (;;) {
            call(progressed);
            if (st.block_state == ISAL_BLOCK_FINISH) return;
            if (!progressed) throw err("no progress before the end of the stream");
        }
    }
    while (st.avail_in > 0) {
        call(progressed);
        if (st.block_state == ISAL_BLOCK_FINISH) throw err("stream ended before the piece end");
        if (!progressed) break;
    }
    // Phase 2: stop at every block end until the piece end.
    const size_t at = size_t(st.next_in - in.p);
    st.avail_in = uint32_t(std::min(in.n, endByte + 4096) - at);
    st.points_to_stop_at = ISAL_STOPPING_POINT_END_OF_BLOCK;
    for (;;) {
        call(progressed);
        if (st.stopped_at == ISAL_STOPPING_POINT_END_OF_BLOCK) {
            const uint64_t pos = bitPos(st, in.p);
            if (pos == endBit) {
                if (st.tmp_out_valid != st.tmp_out_processed) throw err("output left in tmp_out_buffer");
                return;
            }
            if (pos > endBit) throw err("passed the piece end: " + std::to_string(pos) + " > " + std::to_string(endBit));
            continue;
        }
        if (st.block_state == ISAL_BLOCK_FINISH || st.block_state == ISAL_BLOCK_INPUT_DONE)
            throw err("stream ended before the piece end (phase 2)");
        if (!progressed) throw err("no progress in phase 2");
    }
}

// ---------------------------------------------------------------------------------------------
// Merges: markers are packed as (output offset << 16) | window position.
static void merge3(const uint8_t* o1, const uint8_t* o2, const uint8_t* o3, size_t n, std::vector<uint64_t>& mk)
{
    mk.clear();
    size_t i = 0;
    for (; i + 32 <= n; i += 32) {
        const __m256i a = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(o1 + i));
        const __m256i c = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(o3 + i));
        const uint32_t eq = uint32_t(_mm256_movemask_epi8(_mm256_cmpeq_epi8(a, c)));
        if (eq == 0xFFFFFFFFu) continue;
        uint32_t ne = ~eq;
        while (ne) {
            const size_t x = i + size_t(__builtin_ctz(ne));
            ne &= ne - 1;
            const uint32_t p = (o1[x] & 0x7fu) | ((o2[x] & 0x7fu) << 7) | (uint32_t(o3[x]) << 14);
            mk.push_back((uint64_t(x) << 16) | p);
        }
    }
    for (; i < n; ++i) {
        if (o1[i] != o3[i]) {
            const uint32_t p = (o1[i] & 0x7fu) | ((o2[i] & 0x7fu) << 7) | (uint32_t(o3[i]) << 14);
            mk.push_back((uint64_t(i) << 16) | p);
        }
    }
}

static void merge2(const uint8_t* oa, const uint8_t* ob, size_t n, std::vector<uint64_t>& mk)
{
    mk.clear();
    size_t i = 0;
    for (; i + 32 <= n; i += 32) {
        const __m256i a = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(oa + i));
        const __m256i b = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(ob + i));
        const uint32_t eq = uint32_t(_mm256_movemask_epi8(_mm256_cmpeq_epi8(a, b)));
        if (eq == 0xFFFFFFFFu) continue;
        uint32_t ne = ~eq;
        while (ne) {
            const size_t x = i + size_t(__builtin_ctz(ne));
            ne &= ne - 1;
            const uint32_t p = oa[x] | (((uint32_t(ob[x]) - oa[x] - 1) & 0xffu) << 8);
            mk.push_back((uint64_t(x) << 16) | p);
        }
    }
    for (; i < n; ++i) {
        if (oa[i] != ob[i]) {
            const uint32_t p = oa[i] | (((uint32_t(ob[i]) - oa[i] - 1) & 0xffu) << 8);
            mk.push_back((uint64_t(i) << 16) | p);
        }
    }
}

static void resolve(uint8_t* out, const std::vector<uint64_t>& mk, const uint8_t* window)
{
    for (const uint64_t m : mk) out[m >> 16] = window[m & 0xffff];
}

// Branchless variant: for every 16-byte half-block with a marker, compute 16-bit indices into
// T = window ++ identity(256) with AVX2 and look all 16 bytes up.
static size_t mergeResolveTable(int variant, uint8_t* a, const uint8_t* b, const uint8_t* c, size_t n, const uint8_t* window)
{
    alignas(64) uint8_t T[WIN + 256];
    memcpy(T, window, WIN);
    for (uint32_t v = 0; v < 256; ++v) T[WIN + v] = uint8_t(v);
    const uint8_t* cmp = variant == 3 ? c : b;
    const __m256i ones16 = _mm256_set1_epi16(1), ff16 = _mm256_set1_epi16(0xff), x7f = _mm256_set1_epi16(0x7f);
    const __m256i win16 = _mm256_set1_epi16(int16_t(WIN)), mask15 = _mm256_set1_epi16(0x7fff);
    alignas(32) uint16_t idx[16];
    size_t i = 0, markers = 0;
    for (; i + 16 <= n; i += 16) {
        const __m128i a8 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(a + i));
        const __m128i c8 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(cmp + i));
        const uint32_t eq = uint32_t(_mm_movemask_epi8(_mm_cmpeq_epi8(a8, c8)));
        if (eq == 0xFFFFu) continue;
        markers += size_t(__builtin_popcount(~eq & 0xFFFFu));
        const __m256i A = _mm256_cvtepu8_epi16(a8);
        const __m256i B = _mm256_cvtepu8_epi16(_mm_loadu_si128(reinterpret_cast<const __m128i*>(b + i)));
        __m256i p;
        if (variant == 3) {
            const __m256i C = _mm256_cvtepu8_epi16(c8);
            p = _mm256_or_si256(_mm256_or_si256(_mm256_and_si256(A, x7f), _mm256_slli_epi16(_mm256_and_si256(B, x7f), 7)),
                                _mm256_slli_epi16(C, 14));
            p = _mm256_and_si256(p, mask15);
        } else {
            const __m256i hi = _mm256_and_si256(_mm256_sub_epi16(_mm256_sub_epi16(B, A), ones16), ff16);
            p = _mm256_and_si256(_mm256_or_si256(A, _mm256_slli_epi16(hi, 8)), mask15);
        }
        const __m256i same = _mm256_cvtepi8_epi16(_mm_cmpeq_epi8(a8, c8));  // 0xffff where not a marker
        const __m256i ix = _mm256_blendv_epi8(p, _mm256_add_epi16(A, win16), same);
        _mm256_store_si256(reinterpret_cast<__m256i*>(idx), ix);
        for (int j = 0; j < 16; ++j) a[i + j] = T[idx[j]];
    }
    for (; i < n; ++i) {
        if (a[i] != cmp[i]) {
            ++markers;
            uint32_t p;
            if (variant == 3) p = (a[i] & 0x7fu) | ((b[i] & 0x7fu) << 7) | (uint32_t(c[i]) << 14);
            else p = a[i] | (((uint32_t(b[i]) - a[i] - 1) & 0xffu) << 8);
            a[i] = window[p & (WIN - 1)];
        }
    }
    return markers;
}

static int g_mrImpl = 1;  // 0: bit loop, 1: branchless table (env MR_IMPL)

static size_t mergeResolveBits(int variant, uint8_t* a, const uint8_t* b, const uint8_t* c, size_t n, const uint8_t* window);

// Fused merge + resolve without a marker list: a holds pass-1 (variant 3) or pass-A (variant 2) bytes and
// is resolved in place with the true window; b (and c) are the other passes. Returns the marker count.
static size_t mergeResolve(int variant, uint8_t* a, const uint8_t* b, const uint8_t* c, size_t n, const uint8_t* window)
{
    return g_mrImpl == 1 ? mergeResolveTable(variant, a, b, c, n, window) : mergeResolveBits(variant, a, b, c, n, window);
}

static size_t mergeResolveBits(int variant, uint8_t* a, const uint8_t* b, const uint8_t* c, size_t n, const uint8_t* window)
{
    const uint8_t* cmp = variant == 3 ? c : b;
    size_t i = 0, markers = 0;
    auto one = [&](size_t x) {
        uint32_t p;
        if (variant == 3) p = (a[x] & 0x7fu) | ((b[x] & 0x7fu) << 7) | (uint32_t(c[x]) << 14);
        else p = a[x] | (((uint32_t(b[x]) - a[x] - 1) & 0xffu) << 8);
        a[x] = window[p & (WIN - 1)];
    };
    for (; i + 32 <= n; i += 32) {
        const __m256i va = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(a + i));
        const __m256i vc = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(cmp + i));
        const uint32_t eq = uint32_t(_mm256_movemask_epi8(_mm256_cmpeq_epi8(va, vc)));
        if (eq == 0xFFFFFFFFu) continue;
        uint32_t ne = ~eq;
        markers += size_t(__builtin_popcount(ne));
        while (ne) {
            one(i + size_t(__builtin_ctz(ne)));
            ne &= ne - 1;
        }
    }
    for (; i < n; ++i) {
        if (a[i] != cmp[i]) {
            ++markers;
            one(i);
        }
    }
    return markers;
}

// Keeps the last 32 KiB of all output so far.
static void pushWindow(uint8_t* win, const uint8_t* data, size_t n)
{
    if (n >= WIN) {
        memcpy(win, data + n - WIN, WIN);
    } else {
        memmove(win, win + n, WIN - n);
        memcpy(win + WIN - n, data, n);
    }
}

// ---------------------------------------------------------------------------------------------
struct Bounds
{
    uint64_t deflateStartBit = 0, deflateEndBit = 0;
    uint64_t totalOut = 0;
    uint32_t crc = 0, isize = 0;
    std::vector<uint64_t> bit;  // block starts (the first is deflateStartBit)
    std::vector<uint64_t> out;  // output offset at each block start
};

static size_t gzipHeaderBytes(const Mapped& in)
{
    auto st = std::make_unique<inflate_state>();
    isal_inflate_init(st.get());
    st->next_in = const_cast<uint8_t*>(in.p);
    st->avail_in = uint32_t(std::min<size_t>(in.n, 1 << 20));
    isal_gzip_header hdr{};
    isal_gzip_header_init(&hdr);
    if (isal_read_gzip_header(st.get(), &hdr) != ISAL_DECOMP_OK) throw err("bad gzip header");
    if (st->read_in_length != 0) throw err("unexpected bits after the gzip header");
    return size_t(st->next_in - in.p);
}

static int modeBounds(const char* gz, const char* outTsv)
{
    const Mapped in = mapFile(gz, true);
    const size_t hdr = gzipHeaderBytes(in);
    auto stp = std::make_unique<inflate_state>();
    inflate_state& st = *stp;

    std::vector<uint64_t> bits, outs;
    Buf buf;
    buf.reserve(size_t(8) << 20);
    const double t0 = now(), c0 = cpuSeconds();
    isal_inflate_init(&st);
    st.crc_flag = ISAL_DEFLATE;
    st.points_to_stop_at = ISAL_STOPPING_POINT_END_OF_BLOCK;
    st.next_in = const_cast<uint8_t*>(in.p + hdr);
    st.avail_in = uint32_t(in.n - hdr);
    uint64_t total = 0;
    uint32_t crc = 0;
    bits.push_back(uint64_t(hdr) * 8);
    outs.push_back(0);
    uint64_t endBit = 0;
    for (;;) {
        if (buf.cap - buf.size < (size_t(1) << 20)) buf.size = 0;
        st.next_out = buf.p + buf.size;
        st.avail_out = uint32_t(buf.cap - buf.size);
        uint8_t* before = st.next_out;
        const int ret = isal_inflate(&st);
        if (ret < 0) throw err("isal_inflate error in bounds pass");
        const size_t got = size_t(st.next_out - before);
        crc = crc32_gzip_refl(crc, before, got);
        total += got;
        buf.size += got;
        if (st.block_state == ISAL_BLOCK_FINISH) {
            endBit = bitPos(st, in.p);
            break;
        }
        if (st.stopped_at == ISAL_STOPPING_POINT_END_OF_BLOCK) {
            const uint64_t pos = bitPos(st, in.p);
            if (st.bfinal && (st.block_state == ISAL_BLOCK_INPUT_DONE)) continue;
            if (pos != bits.back()) {
                bits.push_back(pos);
                outs.push_back(total);
            }
        }
    }
    const double t1 = now(), c1 = cpuSeconds();
    // Drop a recorded "start" that is the end of the final block.
    while (!bits.empty() && bits.back() >= endBit) { bits.pop_back(); outs.pop_back(); }
    uint32_t fcrc = 0, fsize = 0;
    memcpy(&fcrc, in.p + in.n - 8, 4);
    memcpy(&fsize, in.p + in.n - 4, 4);

    FILE* f = fopen(outTsv, "w");
    fprintf(f, "#deflate_start_bit\t%llu\tdeflate_end_bit\t%llu\ttotal_out\t%llu\tcrc\t%u\tisize\t%u\n",
            (unsigned long long)bits.front(), (unsigned long long)endBit, (unsigned long long)total, fcrc, fsize);
    for (size_t i = 0; i < bits.size(); ++i)
        fprintf(f, "%llu\t%llu\n", (unsigned long long)bits[i], (unsigned long long)outs[i]);
    fclose(f);
    printf("bounds\tblocks\t%zu\ttotal_out\t%llu\tcrc_ok\t%d\tisize_ok\t%d\twall_s\t%.3f\tcpu_s\t%.3f\tGBps\t%.3f\t"
           "mean_block_out_KiB\t%.1f\tmean_block_in_KiB\t%.1f\n",
           bits.size(), (unsigned long long)total, int(crc == fcrc), int(uint32_t(total) == fsize), t1 - t0, c1 - c0,
           total / 1e9 / (t1 - t0), total / 1024.0 / bits.size(), (endBit - bits.front()) / 8.0 / 1024.0 / bits.size());
    return 0;
}

static Bounds readBounds(const char* path)
{
    Bounds b;
    FILE* f = fopen(path, "r");
    if (!f) { perror(path); exit(1); }
    unsigned long long a, e, t;
    unsigned c, s;
    if (fscanf(f, "#deflate_start_bit\t%llu\tdeflate_end_bit\t%llu\ttotal_out\t%llu\tcrc\t%u\tisize\t%u\n", &a, &e, &t, &c, &s) != 5)
        throw err("bad bounds header");
    b.deflateStartBit = a;
    b.deflateEndBit = e;
    b.totalOut = t;
    b.crc = c;
    b.isize = s;
    unsigned long long x, y;
    while (fscanf(f, "%llu\t%llu\n", &x, &y) == 2) {
        b.bit.push_back(x);
        b.out.push_back(y);
    }
    fclose(f);
    return b;
}

struct PieceDef
{
    uint64_t startBit, endBit, outOff;
};

static std::vector<PieceDef> choosePieces(const Bounds& b, double pieceMiB)
{
    const uint64_t step = uint64_t(pieceMiB * 1048576.0) * 8;
    std::vector<PieceDef> v;
    uint64_t target = b.bit[0] + step;
    v.push_back({ b.bit[0], 0, 0 });
    for (size_t i = 1; i < b.bit.size(); ++i) {
        if (b.bit[i] >= target) {
            v.back().endBit = b.bit[i];
            v.push_back({ b.bit[i], 0, b.out[i] });
            target = b.bit[i] + step;
        }
    }
    v.back().endBit = UINT64_MAX;
    return v;
}

static size_t hintFor(const PieceDef& p, const Mapped& in)
{
    const uint64_t endBit = p.endBit == UINT64_MAX ? uint64_t(in.n) * 8 : p.endBit;
    return size_t((endBit - p.startBit) / 8) * 5 + (size_t(1) << 20);
}

// ---------------------------------------------------------------------------------------------
static int modeSeq(const char* gz, size_t chunkKiB)
{
    const double tm0 = now();
    const Mapped in = mapFile(gz, true);
    const double tm1 = now();
    auto stp = std::make_unique<inflate_state>();
    inflate_state& st = *stp;
    Buf buf;
    buf.reserve(size_t(1) << 20);
    for (int mode = 0; mode < 2; ++mode) {
        // mode 0: gzip with CRC check (isal_inflate crc_flag ISAL_GZIP); mode 1: raw deflate, no CRC
        const double t0 = now(), c0 = cpuSeconds();
        isal_inflate_init(&st);
        size_t pos = 0;
        if (mode == 0) {
            st.crc_flag = ISAL_GZIP;
        } else {
            st.crc_flag = ISAL_DEFLATE;
            pos = gzipHeaderBytes(in);
        }
        st.next_in = const_cast<uint8_t*>(in.p + pos);
        st.avail_in = 0;
        uint64_t total = 0;
        for (;;) {
            if (st.avail_in == 0) {
                const size_t at = size_t(st.next_in - in.p);
                const size_t chunk = chunkKiB ? chunkKiB * 1024 : in.n;
                st.avail_in = uint32_t(std::min(chunk, in.n - at));
            }
            st.next_out = buf.p;
            st.avail_out = uint32_t(buf.cap);
            const int ret = isal_inflate(&st);
            if (ret < 0) throw err("isal_inflate error (seq) " + std::to_string(ret));
            total += size_t(st.next_out - buf.p);
            if (st.block_state == ISAL_BLOCK_FINISH) break;
        }
        const double t1 = now(), c1 = cpuSeconds();
        printf("seq\t%s\tchunk_KiB\t%zu\ttotal_out\t%llu\twall_s\t%.3f\tcpu_s\t%.3f\tGBps\t%.3f\tcpu_s_per_GB\t%.3f\t"
               "mmap_s\t%.3f\tmaxrss_MiB\t%ld\n",
               mode == 0 ? "gzip_crc" : "raw_nocrc", chunkKiB, (unsigned long long)total, t1 - t0, c1 - c0,
               total / 1e9 / (t1 - t0), (c1 - c0) / (total / 1e9), tm1 - tm0, maxRssMiB());
    }
    return 0;
}

// ---------------------------------------------------------------------------------------------
struct MarkerStats
{
    size_t n = 0, markers = 0, runs = 0, farP = 0, ambig12 = 0, lastMarker = 0, inLast32K = 0;
    size_t firstMiB = 0, lastMiB = 0;
};

static MarkerStats markerStats(const std::vector<uint64_t>& mk, size_t n)
{
    MarkerStats s;
    s.n = n;
    s.markers = mk.size();
    uint64_t prevX = UINT64_MAX, prevP = UINT64_MAX;
    for (const uint64_t m : mk) {
        const uint64_t x = m >> 16, p = m & 0xffff;
        if (!(x == prevX + 1 && p == prevP + 1)) ++s.runs;
        prevX = x;
        prevP = p;
        if (p >= 16384) ++s.farP;
        if ((p & 0x7f) == ((p >> 7) & 0x7f)) ++s.ambig12;
        if (x + WIN >= n) ++s.inLast32K;
        if (x < (1u << 20)) ++s.firstMiB;
        if (x + (1u << 20) >= n) ++s.lastMiB;
        s.lastMarker = x;
    }
    return s;
}

static int modeCost(const char* gz, const char* boundsPath, double pieceMiB, const char* fq, const char* statsPath)
{
    initLabels();
    const Mapped in = mapFile(gz, true);
    const Bounds b = readBounds(boundsPath);
    const auto pieces = choosePieces(b, pieceMiB);
    Mapped ref;
    if (fq && *fq && strcmp(fq, "-") != 0) ref = mapFile(fq, false);
    auto stp = std::make_unique<inflate_state>();
    inflate_state& st = *stp;
    Buf T, B1, B2, B3;
    std::vector<uint64_t> mk3, mk2;
    uint8_t win[WIN];
    memset(win, 0, WIN);
    double tF2 = 0, tF3 = 0, tTrue = 0, tP1 = 0, tP2 = 0, tP3 = 0, tM3 = 0, tR3 = 0, tA = 0, tB = 0, tM2 = 0, tR2 = 0, tFirst = 0;
    uint64_t outSpec = 0, outAll = 0, totalMarkers = 0, totalRuns = 0, totalFar = 0, totalAmbig = 0;
    size_t piecesDying = 0, mismatches = 0, refMismatches = 0;
    uint32_t crc = 0;
    int lineMod = 0;
    uint64_t typeBytes[4] = {}, typeMarkers[4] = {};
    FILE* sf = statsPath && *statsPath ? fopen(statsPath, "w") : nullptr;
    if (sf) fprintf(sf, "piece\tout_bytes\tmarkers\tmarker_share\truns\tp_ge_16384\tambig_L1eqL2\tlast_marker_frac\tin_last_32K\tfirst_MiB\tlast_MiB\n");
    for (size_t k = 0; k < pieces.size(); ++k) {
        const auto& pc = pieces[k];
        const size_t hint = hintFor(pc, in);
        double t = now();
        inflateRange(in, pc.startBit, pc.endBit, k ? win : nullptr, T, st, hint);
        const double dTrue = now() - t;
        if (k == 0) {
            tFirst += dTrue;
        } else {
            tTrue += dTrue;
            outSpec += T.size;
            // Variant 3
            t = now(); inflateRange(in, pc.startBit, pc.endBit, L1, B1, st, hint); tP1 += now() - t;
            t = now(); inflateRange(in, pc.startBit, pc.endBit, L2, B2, st, hint); tP2 += now() - t;
            t = now(); inflateRange(in, pc.startBit, pc.endBit, L3, B3, st, hint); tP3 += now() - t;
            if (B1.size != T.size || B2.size != T.size || B3.size != T.size) throw err("pass sizes differ");
            t = now(); merge3(B1.p, B2.p, B3.p, T.size, mk3); tM3 += now() - t;
            t = now(); resolve(B1.p, mk3, win); tR3 += now() - t;
            if (memcmp(B1.p, T.p, T.size) != 0) ++mismatches;
            // Variant 2
            t = now(); inflateRange(in, pc.startBit, pc.endBit, LA, B2, st, hint); tA += now() - t;
            t = now(); inflateRange(in, pc.startBit, pc.endBit, LB, B3, st, hint); tB += now() - t;
            if (B2.size != T.size || B3.size != T.size) throw err("pass sizes differ (v2)");
            t = now(); merge2(B2.p, B3.p, T.size, mk2); tM2 += now() - t;
            if (mk2 != mk3) ++mismatches;
            t = now(); resolve(B2.p, mk2, win); tR2 += now() - t;
            if (memcmp(B2.p, T.p, T.size) != 0) ++mismatches;
            // Variant 2 again, fused merge + resolve (no marker list).
            inflateRange(in, pc.startBit, pc.endBit, LA, B2, st, hint);
            t = now(); const size_t fm = mergeResolve(2, B2.p, B3.p, nullptr, T.size, win); tF2 += now() - t;
            if (fm != mk3.size() || memcmp(B2.p, T.p, T.size) != 0) ++mismatches;
            // Variant 3 fused: B1 (pass 1) is resolved already, so decode pass 1 again.
            inflateRange(in, pc.startBit, pc.endBit, L1, B1, st, hint);
            inflateRange(in, pc.startBit, pc.endBit, L2, B2, st, hint);
            inflateRange(in, pc.startBit, pc.endBit, L3, B3, st, hint);
            t = now(); const size_t f3 = mergeResolve(3, B1.p, B2.p, B3.p, T.size, win); tF3 += now() - t;
            if (f3 != mk3.size() || memcmp(B1.p, T.p, T.size) != 0) ++mismatches;
            const MarkerStats s = markerStats(mk3, T.size);
            totalMarkers += s.markers;
            totalRuns += s.runs;
            totalFar += s.farP;
            totalAmbig += s.ambig12;
            if (s.markers == 0 || s.lastMarker + WIN < s.n) ++piecesDying;
            if (sf)
                fprintf(sf, "%zu\t%zu\t%zu\t%.6f\t%zu\t%zu\t%zu\t%.4f\t%zu\t%zu\t%zu\n", k, s.n, s.markers,
                        double(s.markers) / s.n, s.runs, s.farP, s.ambig12, s.n ? double(s.lastMarker) / s.n : 0.0,
                        s.inLast32K, s.firstMiB, s.lastMiB);
        }
        // Markers by FASTQ line type (0 header, 1 sequence, 2 '+', 3 quality), newline counted with its line.
        {
            size_t pos = 0, mi = 0;
            int lt = lineMod;
            while (pos < T.size) {
                const void* nl = memchr(T.p + pos, '\n', T.size - pos);
                const size_t end = nl ? size_t(static_cast<const uint8_t*>(nl) - T.p) + 1 : T.size;
                if (k > 0) {
                    typeBytes[lt] += end - pos;
                    while (mi < mk3.size() && (mk3[mi] >> 16) < end) {
                        ++typeMarkers[lt];
                        ++mi;
                    }
                }
                pos = end;
                if (nl) lt = (lt + 1) & 3;
            }
            lineMod = lt;
        }
        if (pc.outOff != outAll) throw err("output offset mismatch");
        if (ref.p && memcmp(ref.p + outAll, T.p, T.size) != 0) ++refMismatches;
        crc = crc32_gzip_refl(crc, T.p, T.size);
        outAll += T.size;
        pushWindow(win, T.p, T.size);
    }
    if (sf) fclose(sf);
    const double g = outSpec / 1e9;
    printf("cost\tpiece_MiB\t%.2f\tpieces\t%zu\tout_total\t%llu\tcrc_ok\t%d\tsize_ok\t%d\tref_checked\t%d\tref_mismatch_pieces\t%zu\t"
           "spec_mismatch_pieces\t%zu\n",
           pieceMiB, pieces.size(), (unsigned long long)outAll, int(crc == b.crc), int(outAll == b.totalOut), int(ref.p != nullptr),
           refMismatches, mismatches);
    printf("cost_s_per_GB\ttrue_window_1pass\t%.4f\tv3_pass1\t%.4f\tv3_pass2\t%.4f\tv3_pass3\t%.4f\tv3_merge\t%.4f\tv3_resolve\t%.4f\t"
           "v2_passA\t%.4f\tv2_passB\t%.4f\tv2_merge\t%.4f\tv2_resolve\t%.4f\tv2_fused_merge_resolve\t%.4f\t"
           "v3_fused_merge_resolve\t%.4f\tfirst_piece_s\t%.4f\n",
           tTrue / g, tP1 / g, tP2 / g, tP3 / g, tM3 / g, tR3 / g, tA / g, tB / g, tM2 / g, tR2 / g, tF2 / g, tF3 / g, tFirst);
    printf("markers\ttotal\t%llu\tshare_of_spec_output\t%.6f\truns\t%llu\tmean_run_len\t%.2f\tp_ge_16384\t%llu\tshare_p_ge_16384\t%.4f\t"
           "ambig_L1eqL2\t%llu\tpieces_where_markers_die\t%zu\tof\t%zu\n",
           (unsigned long long)totalMarkers, double(totalMarkers) / outSpec, (unsigned long long)totalRuns,
           totalRuns ? double(totalMarkers) / totalRuns : 0.0, (unsigned long long)totalFar,
           totalMarkers ? double(totalFar) / totalMarkers : 0.0, (unsigned long long)totalAmbig, piecesDying, pieces.size() - 1);
    static const char* names[4] = { "header", "sequence", "plus", "quality" };
    for (int i = 0; i < 4; ++i)
        printf("markers_by_line\t%s\tbytes_share\t%.4f\tmarkers_share_of_all_markers\t%.4f\tmarker_share_within_line_type\t%.4f\n",
               names[i], double(typeBytes[i]) / outSpec, totalMarkers ? double(typeMarkers[i]) / totalMarkers : 0.0,
               typeBytes[i] ? double(typeMarkers[i]) / typeBytes[i] : 0.0);
    return 0;
}

// ---------------------------------------------------------------------------------------------
struct Piece
{
    PieceDef def{};
    // 0 free, 1 claimed (decoding), 2 decoded, 3 chained (window after it known, resolve queued), 4 resolved
    int state = 0;  // guarded by the mutex
    bool spec = false;
    Buf* a = nullptr;  // output (pass 1 / pass A; resolved in place)
    Buf* b = nullptr;  // pass 2 / pass B
    Buf* c = nullptr;  // pass 3 (variant 3)
    std::unique_ptr<uint8_t[]> win;  // true window before the piece (spec pieces only)
    size_t markers = 0;
};

struct Pool
{
    std::mutex m;
    std::vector<Buf*> free;
    Buf* get()
    {
        std::lock_guard<std::mutex> g(m);
        if (free.empty()) return new Buf();
        Buf* b = free.back();
        free.pop_back();
        return b;
    }
    void put(Buf*& b)
    {
        if (!b) return;
        std::lock_guard<std::mutex> g(m);
        free.push_back(b);
        b = nullptr;
    }
    ~Pool() { for (Buf* b : free) delete b; }
};

static void decodeSpec(const Mapped& in, Piece& pc, int variant, inflate_state& st, Pool& pool)
{
    const size_t hint = hintFor(pc.def, in);
    pc.a = pool.get();
    pc.b = pool.get();
    if (variant == 3) {
        pc.c = pool.get();
        inflateRange(in, pc.def.startBit, pc.def.endBit, L1, *pc.a, st, hint);
        inflateRange(in, pc.def.startBit, pc.def.endBit, L2, *pc.b, st, hint);
        inflateRange(in, pc.def.startBit, pc.def.endBit, L3, *pc.c, st, hint);
    } else {
        inflateRange(in, pc.def.startBit, pc.def.endBit, LA, *pc.a, st, hint);
        inflateRange(in, pc.def.startBit, pc.def.endBit, LB, *pc.b, st, hint);
    }
    if (pc.b->size != pc.a->size || (pc.c && pc.c->size != pc.a->size)) throw err("pass sizes differ");
    pc.spec = true;
}

static int modePar(const char* gz, const char* boundsPath, double pieceMiB, int threads, const std::string& sched,
                   int variant, int gap, const char* fq)
{
    initLabels();
    const double tm0 = now();
    const Mapped in = mapFile(gz, true);
    const double tm1 = now();
    const Bounds bnd = readBounds(boundsPath);
    const auto defs = choosePieces(bnd, pieceMiB);
    const size_t n = defs.size();
    Mapped ref;
    if (fq && *fq && strcmp(fq, "-") != 0) ref = mapFile(fq, false);
    std::unique_ptr<Piece[]> pcs(new Piece[n]);
    for (size_t i = 0; i < n; ++i) pcs[i].def = defs[i];

    const bool lead = sched == "lead";
    const int nWorkers = lead ? threads - 1 : threads;
    if (gap <= 0) gap = 2 * threads + 4;
    const size_t maxAhead = lead ? size_t(gap + 2 * threads) : size_t(2 * threads + 2);

    Pool pool;
    std::mutex mtx;
    std::condition_variable cv;
    // All guarded by mtx:
    size_t chainPos = 0;    // next piece whose window is known (in lead mode: the lead's position)
    size_t consumePos = 0;  // next piece to consume
    size_t specNext = 0;
    std::vector<size_t> resolveQ;
    size_t resolveHead = 0;
    bool done = false, failed = false;
    std::string failure;
    double workerResolveS = 0, workerDecodeS = 0;

    auto runResolve = [&](size_t k) {
        Piece& pc = pcs[k];
        pc.markers = mergeResolve(variant, pc.a->p, pc.b->p, pc.c ? pc.c->p : nullptr, pc.a->size, pc.win.get());
        pool.put(pc.b);
        pool.put(pc.c);
        pc.win.reset();
    };

    const double t0 = now(), c0 = cpuSeconds();
    std::vector<std::thread> ths;
    for (int t = 0; t < nWorkers; ++t) {
        ths.emplace_back([&] {
            auto st = std::make_unique<inflate_state>();
            double myResolve = 0, myDecode = 0;
            try {
                for (;;) {
                    size_t k = 0;
                    bool isResolve = false;
                    {
                        std::unique_lock<std::mutex> g(mtx);
                        for (;;) {
                            if (done || failed) {
                                workerResolveS += myResolve;
                                workerDecodeS += myDecode;
                                return;
                            }
                            if (resolveHead < resolveQ.size()) {
                                k = resolveQ[resolveHead++];
                                isResolve = true;
                                break;
                            }
                            const size_t cand = lead ? std::max(specNext, chainPos + size_t(gap)) : specNext;
                            if (cand < n && cand < consumePos + maxAhead) {
                                specNext = cand + 1;
                                if (pcs[cand].state == 0) {
                                    pcs[cand].state = 1;
                                    k = cand;
                                    break;
                                }
                                continue;
                            }
                            cv.wait(g);
                        }
                    }
                    Piece& pc = pcs[k];
                    const double ts = now();
                    if (isResolve) {
                        runResolve(k);
                        myResolve += now() - ts;
                        std::lock_guard<std::mutex> g(mtx);
                        pc.state = 4;
                    } else {
                        if (k == 0) {
                            pc.a = pool.get();
                            inflateRange(in, pc.def.startBit, pc.def.endBit, nullptr, *pc.a, *st, hintFor(pc.def, in));
                            pc.spec = false;
                        } else {
                            decodeSpec(in, pc, variant, *st, pool);
                        }
                        myDecode += now() - ts;
                        std::lock_guard<std::mutex> g(mtx);
                        pc.state = 2;
                    }
                    cv.notify_all();
                }
            } catch (const std::exception& e) {
                std::lock_guard<std::mutex> g(mtx);
                if (failure.empty()) failure = e.what();
                failed = true;
                cv.notify_all();
            }
        });
    }

    // Main thread: chain (true windows in order), consumption (CRC, reference check), and in "lead" mode the
    // in-order decoding of pieces nobody speculated on with the true window, plus resolve tasks while idle.
    auto st = std::make_unique<inflate_state>();
    uint8_t cur[WIN];
    memset(cur, 0, WIN);
    std::vector<uint8_t> tailA(WIN);
    uint32_t crc = 0;
    uint64_t outAll = 0, markers = 0;
    size_t nSpec = 0, nLead = 0, refMismatch = 0, offMismatch = 0;
    double waitS = 0, chainS = 0, leadDecodeS = 0, leadResolveS = 0, crcS = 0;
    try {
        for (;;) {
            bool progress = false;
            size_t cp, cons;
            int chainState = -1;
            {
                std::lock_guard<std::mutex> g(mtx);
                if (failed) throw err(failure);
                cp = chainPos;
                cons = consumePos;
                if (cp < n) {
                    if (lead && pcs[cp].state == 0) {
                        pcs[cp].state = 1;
                        chainState = 100;  // claimed by the lead
                    } else {
                        chainState = pcs[cp].state;
                    }
                }
            }
            if (cons == n) break;
            if (cp < n && chainState == 100) {
                const double t = now();
                Piece& pc = pcs[cp];
                pc.a = pool.get();
                inflateRange(in, pc.def.startBit, pc.def.endBit, cp ? cur : nullptr, *pc.a, *st, hintFor(pc.def, in));
                pc.spec = false;
                ++nLead;
                leadDecodeS += now() - t;
                chainState = 2;
            }
            if (cp < n && chainState == 2) {
                const double t = now();
                Piece& pc = pcs[cp];
                bool queue = false;
                if (pc.spec) {
                    pc.win.reset(new uint8_t[WIN]);
                    memcpy(pc.win.get(), cur, WIN);
                    const size_t sz = pc.a->size, nt = std::min<size_t>(sz, WIN);
                    memcpy(tailA.data(), pc.a->p + sz - nt, nt);
                    mergeResolve(variant, tailA.data(), pc.b->p + sz - nt, pc.c ? pc.c->p + sz - nt : nullptr, nt, pc.win.get());
                    pushWindow(cur, tailA.data(), nt);
                    queue = true;
                    ++nSpec;
                } else {
                    pushWindow(cur, pc.a->p, pc.a->size);
                }
                chainS += now() - t;
                {
                    std::lock_guard<std::mutex> g(mtx);
                    pc.state = queue ? 3 : 4;
                    if (queue) resolveQ.push_back(cp);
                    chainPos = cp + 1;
                }
                cv.notify_all();
                progress = true;
            }
            for (;;) {
                {
                    std::lock_guard<std::mutex> g(mtx);
                    cons = consumePos;
                    if (cons >= chainPos || pcs[cons].state != 4) break;
                }
                Piece& pc = pcs[cons];
                const double tc = now();
                crc = crc32_gzip_refl(crc, pc.a->p, pc.a->size);
                crcS += now() - tc;
                markers += pc.markers;
                if (pc.def.outOff != outAll) ++offMismatch;
                if (ref.p && (outAll + pc.a->size > ref.n || memcmp(ref.p + outAll, pc.a->p, pc.a->size) != 0)) ++refMismatch;
                outAll += pc.a->size;
                pool.put(pc.a);
                {
                    std::lock_guard<std::mutex> g(mtx);
                    consumePos = cons + 1;
                }
                cv.notify_all();
                progress = true;
            }
            if (progress) continue;
            // Nothing to do in order: in lead mode (or without workers) take a resolve task, else wait.
            size_t task = SIZE_MAX;
            {
                std::unique_lock<std::mutex> g(mtx);
                if ((lead || nWorkers == 0) && resolveHead < resolveQ.size()) {
                    task = resolveQ[resolveHead++];
                } else {
                    const double t = now();
                    cv.wait(g, [&] {
                        if (failed) return true;
                        if (chainPos < n && (pcs[chainPos].state == 2 || (lead && pcs[chainPos].state == 0))) return true;
                        if (consumePos < chainPos && pcs[consumePos].state == 4) return true;
                        if ((lead || nWorkers == 0) && resolveHead < resolveQ.size()) return true;
                        return false;
                    });
                    waitS += now() - t;
                }
            }
            if (task != SIZE_MAX) {
                const double t = now();
                runResolve(task);
                leadResolveS += now() - t;
                {
                    std::lock_guard<std::mutex> g(mtx);
                    pcs[task].state = 4;
                }
                cv.notify_all();
            }
        }
    } catch (const std::exception& e) {
        {
            std::lock_guard<std::mutex> g(mtx);
            failed = true;
            if (failure.empty()) failure = e.what();
        }
        cv.notify_all();
        for (auto& th : ths) th.join();
        fprintf(stderr, "FAILED: %s\n", failure.c_str());
        return 1;
    }
    {
        std::lock_guard<std::mutex> g(mtx);
        done = true;
    }
    cv.notify_all();
    for (auto& th : ths) th.join();
    const double t1 = now(), c1 = cpuSeconds();
    const bool ok = crc == bnd.crc && outAll == bnd.totalOut && refMismatch == 0 && offMismatch == 0;
    printf("par\tsched\t%s\tvariant\t%d\tpiece_MiB\t%.2f\tthreads\t%d\tgap\t%d\tpieces\t%zu\tspec_pieces\t%zu\tlead_pieces\t%zu\t"
           "wall_s\t%.3f\tGBps\t%.3f\tcpu_s\t%.3f\tcpu_s_per_GB\t%.3f\tmaxrss_MiB\t%ld\tmmap_s\t%.3f\t"
           "main_wait_s\t%.3f\tchain_s\t%.3f\tlead_decode_s\t%.3f\tlead_resolve_s\t%.3f\tcrc_s\t%.3f\t"
           "worker_decode_s\t%.3f\tworker_resolve_s\t%.3f\tmarkers\t%llu\t"
           "crc_ok\t%d\tsize_ok\t%d\tref_checked\t%d\tref_mismatch\t%zu\tOK\t%d\n",
           sched.c_str(), variant, pieceMiB, threads, gap, n, nSpec, nLead, t1 - t0, outAll / 1e9 / (t1 - t0), c1 - c0,
           (c1 - c0) / (outAll / 1e9), maxRssMiB(), tm1 - tm0, waitS, chainS, leadDecodeS, leadResolveS, crcS, workerDecodeS,
           workerResolveS, (unsigned long long)markers, int(crc == bnd.crc), int(outAll == bnd.totalOut), int(ref.p != nullptr),
           refMismatch, int(ok));
    return ok ? 0 : 2;
}

int main(int argc, char** argv)
{
    if (argc < 3) {
        fprintf(stderr, "usage: see the header of isalpar.cpp\n");
        return 1;
    }
    const std::string mode = argv[1];
    if (const char* e = getenv("MR_IMPL")) g_mrImpl = atoi(e);
    try {
        if (mode == "bounds" && argc >= 4) return modeBounds(argv[2], argv[3]);
        if (mode == "seq") return modeSeq(argv[2], argc >= 4 ? size_t(atol(argv[3])) : 0);
        if (mode == "cost" && argc >= 5)
            return modeCost(argv[2], argv[3], atof(argv[4]), argc >= 6 ? argv[5] : nullptr, argc >= 7 ? argv[6] : nullptr);
        if (mode == "par" && argc >= 8)
            return modePar(argv[2], argv[3], atof(argv[4]), atoi(argv[5]), argv[6], atoi(argv[7]),
                           argc >= 9 ? atoi(argv[8]) : 0, argc >= 10 ? argv[9] : nullptr);
    } catch (const std::exception& e) {
        fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    fprintf(stderr, "bad arguments\n");
    return 1;
}
