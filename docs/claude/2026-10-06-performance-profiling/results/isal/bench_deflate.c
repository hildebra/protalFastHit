// BGZF-style block compression benchmark: the input text is cut into 65280-byte
// blocks (the BGZF payload size) and each block is compressed on its own into a
// raw deflate stream with an output limit of 65536 - 18 - 8 = 65510 bytes.
//   isal1/2/3 : isal_deflate_stateless, levels 1..3, gzip_flag IGZIP_DEFLATE
//   ld1/ld6   : libdeflate_deflate_compress, levels 1 and 6
// Every compressed block is round-trip checked once with libdeflate.
// Extra: per-block raw inflate of the ld6 blocks (isal_inflate_stateless vs
// libdeflate_deflate_decompress) and CRC32 throughput (ISA-L, libdeflate, zlib-ng).
//
// usage: bench_deflate FILE [limit_MB] [reps]
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <isa-l/igzip_lib.h>
#include <isa-l/crc.h>
#include <libdeflate.h>
#include <zlib-ng.h>

#define BLK 65280
#define OUT_LIMIT (65536 - 18 - 8)

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}
static double cpu(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static int cmpd(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : x > y;
}

enum { M_ISAL1, M_ISAL2, M_ISAL3, M_LD1, M_LD6, NM };
static const char *mname[NM] = { "isal_l1", "isal_l2", "isal_l3", "libdeflate_l1", "libdeflate_l6" };

static uint8_t *lvlbuf[4];
static uint32_t lvlsize[4] = { 0, ISAL_DEF_LVL1_DEFAULT, ISAL_DEF_LVL2_DEFAULT, ISAL_DEF_LVL3_DEFAULT };
static struct libdeflate_compressor *ldc[13];

// compress one block; returns compressed size (0 on failure)
static size_t comp(int m, const uint8_t *in, size_t len, uint8_t *out)
{
    if (m <= M_ISAL3) {
        int level = m - M_ISAL1 + 1;
        struct isal_zstream s;
        isal_deflate_stateless_init(&s);
        s.level = level;
        s.level_buf = lvlbuf[level];
        s.level_buf_size = lvlsize[level];
        s.gzip_flag = IGZIP_DEFLATE;
        s.flush = NO_FLUSH;
        s.end_of_stream = 1;
        s.next_in = (uint8_t *)in;
        s.avail_in = (uint32_t)len;
        s.next_out = out;
        s.avail_out = OUT_LIMIT;
        int ret = isal_deflate_stateless(&s);
        if (ret != COMP_OK) { fprintf(stderr, "isal_deflate_stateless ret %d\n", ret); exit(2); }
        return s.total_out;
    } else {
        int level = m == M_LD1 ? 1 : 6;
        size_t r = libdeflate_deflate_compress(ldc[level], in, len, out, OUT_LIMIT);
        if (r == 0) { fprintf(stderr, "libdeflate overflow\n"); exit(2); }
        return r;
    }
}

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: %s FILE [limit_MB] [reps]\n", argv[0]); return 1; }
    size_t limit = argc > 2 ? (size_t)atol(argv[2]) * 1000000 : 0;
    int reps = argc > 3 ? atoi(argv[3]) : 3;

    int fd = open(argv[1], O_RDONLY);
    if (fd < 0) { perror(argv[1]); return 1; }
    struct stat st;
    fstat(fd, &st);
    size_t n = st.st_size;
    if (limit && limit < n) n = limit;
    uint8_t *text = malloc(n);
    size_t got = 0;
    while (got < n) {
        ssize_t r = read(fd, text + got, n - got);
        if (r <= 0) { perror("read"); return 1; }
        got += r;
    }
    close(fd);
    size_t nblk = (n + BLK - 1) / BLK;

    for (int l = 1; l <= 3; l++) lvlbuf[l] = malloc(lvlsize[l]);
    ldc[1] = libdeflate_alloc_compressor(1);
    ldc[6] = libdeflate_alloc_compressor(6);
    uint8_t *cbuf = malloc(nblk * (size_t)OUT_LIMIT);   // all compressed blocks of one method
    size_t *csz = malloc(nblk * sizeof(size_t));
    uint8_t *dbuf = malloc(BLK + 64);
    printf("input %zu bytes, %zu blocks of %d; level_buf sizes L1 %u L2 %u L3 %u\n", n, nblk, BLK,
           lvlsize[1], lvlsize[2], lvlsize[3]);

    double t[NM][64], tc[NM][64];
    size_t total[NM];
    uint32_t chk[NM];
    struct libdeflate_decompressor *ldd = libdeflate_alloc_decompressor();

    for (int r = 0; r < reps && r < 64; r++) {
        for (int m = 0; m < NM; m++) {
            size_t tot = 0;
            double c0 = cpu(), w0 = now();
            for (size_t b = 0; b < nblk; b++) {
                size_t len = (b + 1) * BLK <= n ? BLK : n - b * BLK;
                csz[b] = comp(m, text + b * BLK, len, cbuf + b * (size_t)OUT_LIMIT);
                tot += csz[b];
            }
            double w1 = now(), c1 = cpu();
            t[m][r] = w1 - w0;
            tc[m][r] = c1 - c0;
            // checksum over the compressed bytes: identical output across runs?
            uint32_t c = 0;
            for (size_t b = 0; b < nblk; b++) c = zng_crc32(c, cbuf + b * (size_t)OUT_LIMIT, (uint32_t)csz[b]);
            if (r == 0) {
                total[m] = tot;
                chk[m] = c;
                // round trip with libdeflate
                for (size_t b = 0; b < nblk; b++) {
                    size_t len = (b + 1) * BLK <= n ? BLK : n - b * BLK, act = 0;
                    if (libdeflate_deflate_decompress(ldd, cbuf + b * (size_t)OUT_LIMIT, csz[b], dbuf, BLK, &act) != LIBDEFLATE_SUCCESS ||
                        act != len || memcmp(dbuf, text + b * BLK, len) != 0) {
                        fprintf(stderr, "%s: round trip failed at block %zu\n", mname[m], b);
                        return 2;
                    }
                }
            } else if (tot != total[m] || c != chk[m]) {
                fprintf(stderr, "%s: output differs between repetitions\n", mname[m]);
                return 2;
            }
            printf("rep %d %-14s %.3f s wall %.3f s cpu\n", r, mname[m], t[m][r], tc[m][r]);
            fflush(stdout);
        }
    }

    printf("\n%-14s %14s %8s %12s %12s %12s %12s  %s\n", "method", "compressed", "ratio", "MB/s(min t)", "MB/s(med t)",
           "cpuMB/s(min)", "cpuMB/s(med)", "crc32(compressed)");
    for (int m = 0; m < NM; m++) {
        double a[64], b[64];
        memcpy(a, t[m], reps * sizeof(double));
        memcpy(b, tc[m], reps * sizeof(double));
        qsort(a, reps, sizeof(double), cmpd);
        qsort(b, reps, sizeof(double), cmpd);
        printf("%-14s %14zu %8.4f %12.1f %12.1f %12.1f %12.1f  %08x\n", mname[m], total[m], (double)total[m] / n,
               n / 1e6 / a[0], n / 1e6 / a[reps / 2], n / 1e6 / b[0], n / 1e6 / b[reps / 2], chk[m]);
    }

    // ---- Extra: per-block raw inflate of libdeflate level-6 blocks ----
    size_t tot = 0;
    for (size_t b = 0; b < nblk; b++) {
        size_t len = (b + 1) * BLK <= n ? BLK : n - b * BLK;
        csz[b] = comp(M_LD6, text + b * BLK, len, cbuf + b * (size_t)OUT_LIMIT);
        tot += csz[b];
    }
    struct inflate_state *is = malloc(sizeof(*is));
    isal_inflate_init(is);
    double ti[64], tl[64], tci[64], tcl[64], tcz[64];
    for (int r = 0; r < reps && r < 64; r++) {
        double w0 = now();
        for (size_t b = 0; b < nblk; b++) {
            size_t len = (b + 1) * BLK <= n ? BLK : n - b * BLK;
            is->next_in = cbuf + b * (size_t)OUT_LIMIT;
            is->avail_in = (uint32_t)csz[b];
            is->next_out = dbuf;
            is->avail_out = BLK;
            is->crc_flag = ISAL_DEFLATE;
            int ret = isal_inflate_stateless(is);
            if (ret != ISAL_DECOMP_OK || is->total_out != len) { fprintf(stderr, "isal_inflate_stateless %d\n", ret); return 2; }
        }
        double w1 = now();
        for (size_t b = 0; b < nblk; b++) {
            size_t len = (b + 1) * BLK <= n ? BLK : n - b * BLK, act;
            if (libdeflate_deflate_decompress(ldd, cbuf + b * (size_t)OUT_LIMIT, csz[b], dbuf, BLK, &act) != LIBDEFLATE_SUCCESS || act != len) {
                fprintf(stderr, "libdeflate decompress failed\n");
                return 2;
            }
        }
        double w2 = now();
        uint32_t c1 = crc32_gzip_refl(0, text, n);
        double w3 = now();
        uint32_t c2 = libdeflate_crc32(0, text, n);
        double w4 = now();
        uint32_t c3 = 0;
        for (size_t off = 0; off < n; off += 1u << 30) c3 = zng_crc32(c3, text + off, (uint32_t)((n - off) < (1u << 30) ? n - off : (1u << 30)));
        double w5 = now();
        if (c1 != c2 || c1 != c3) { fprintf(stderr, "crc mismatch %08x %08x %08x\n", c1, c2, c3); return 2; }
        ti[r] = w1 - w0; tl[r] = w2 - w1; tci[r] = w3 - w2; tcl[r] = w4 - w3; tcz[r] = w5 - w4;
    }
    qsort(ti, reps, sizeof(double), cmpd);
    qsort(tl, reps, sizeof(double), cmpd);
    qsort(tci, reps, sizeof(double), cmpd);
    qsort(tcl, reps, sizeof(double), cmpd);
    qsort(tcz, reps, sizeof(double), cmpd);
    printf("\nper-block raw inflate of the libdeflate_l6 blocks (%zu compressed bytes), MB/s of output, min-time / median-time:\n", tot);
    printf("  isal_inflate_stateless        %8.1f %8.1f\n", n / 1e6 / ti[0], n / 1e6 / ti[reps / 2]);
    printf("  libdeflate_deflate_decompress %8.1f %8.1f\n", n / 1e6 / tl[0], n / 1e6 / tl[reps / 2]);
    printf("CRC32 over the whole text, MB/s min-time / median-time:\n");
    printf("  crc32_gzip_refl               %8.1f %8.1f\n", n / 1e6 / tci[0], n / 1e6 / tci[reps / 2]);
    printf("  libdeflate_crc32              %8.1f %8.1f\n", n / 1e6 / tcl[0], n / 1e6 / tcl[reps / 2]);
    printf("  zng_crc32                     %8.1f %8.1f\n", n / 1e6 / tcz[0], n / 1e6 / tcz[reps / 2]);
    return 0;
}
