// Streaming gzip decompression benchmark: ISA-L isal_inflate vs zlib-ng zng_inflate.
// The compressed file is loaded into memory once; each run feeds it in 128 KiB
// reads (memcpy into an input buffer, like read() from the page cache) and
// decompresses into a 1 MiB output buffer that is then discarded.
// Multi-member files (concatenated gzip, BGZF) are handled by resetting the
// decoder after each member and continuing in the same input buffer.
//
// usage: bench_inflate isal|zng FILE [reps] [verify]
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
#include <zlib-ng.h>

#define IN_CHUNK (128 * 1024)
#define OUT_SIZE (1024 * 1024)

static uint8_t *load(const char *path, size_t *n)
{
    int fd = open(path, O_RDONLY);
    if (fd < 0) { perror(path); exit(1); }
    struct stat st;
    fstat(fd, &st);
    uint8_t *b = malloc(st.st_size ? st.st_size : 1);
    size_t got = 0;
    while (got < (size_t)st.st_size) {
        ssize_t r = read(fd, b + got, st.st_size - got);
        if (r <= 0) { perror("read"); exit(1); }
        got += r;
    }
    close(fd);
    *n = got;
    return b;
}

static double now(clockid_t c)
{
    struct timespec ts;
    clock_gettime(c, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

struct result { uint64_t out; uint32_t crc; int members; };

// ISA-L: isal_inflate with crc_flag ISAL_GZIP (header parsed, CRC32 + ISIZE checked)
static struct result run_isal(const uint8_t *file, size_t n, uint8_t *in, uint8_t *out, int verify)
{
    struct inflate_state *st = malloc(sizeof(*st));
    struct result r = {0, 0, 0};
    size_t pos = 0;
    isal_inflate_init(st);
    st->crc_flag = ISAL_GZIP;
    st->avail_in = 0;
    for (;;) {
        if (st->avail_in == 0) {
            if (pos >= n) break;
            size_t c = n - pos < IN_CHUNK ? n - pos : IN_CHUNK;
            memcpy(in, file + pos, c);
            pos += c;
            st->next_in = in;
            st->avail_in = (uint32_t)c;
        }
        st->next_out = out;
        st->avail_out = OUT_SIZE;
        int ret = isal_inflate(st);
        if (ret != ISAL_DECOMP_OK) {
            fprintf(stderr, "isal_inflate error %d at input offset %zu\n", ret, pos - st->avail_in);
            exit(2);
        }
        size_t produced = st->next_out - out;
        r.out += produced;
        if (verify) r.crc = zng_crc32(r.crc, out, (uint32_t)produced);
        if (st->block_state == ISAL_BLOCK_FINISH) {
            r.members++;
            if (st->avail_in == 0 && pos >= n) break;
            // Next member starts at next_in (possibly in the next read)
            isal_inflate_reset(st);
        }
    }
    if (st->block_state != ISAL_BLOCK_FINISH) { fprintf(stderr, "isal: truncated input\n"); exit(2); }
    free(st);
    return r;
}

// zlib-ng: zng_inflate with windowBits 15+16 (gzip only)
static struct result run_zng(const uint8_t *file, size_t n, uint8_t *in, uint8_t *out, int verify)
{
    zng_stream zs;
    memset(&zs, 0, sizeof(zs));
    struct result r = {0, 0, 0};
    size_t pos = 0;
    int ret = 0, done = 0;
    if (zng_inflateInit2(&zs, 15 + 16) != Z_OK) { fprintf(stderr, "zng init failed\n"); exit(2); }
    zs.avail_in = 0;
    for (;;) {
        if (zs.avail_in == 0) {
            if (pos >= n) break;
            size_t c = n - pos < IN_CHUNK ? n - pos : IN_CHUNK;
            memcpy(in, file + pos, c);
            pos += c;
            zs.next_in = in;
            zs.avail_in = (uint32_t)c;
        }
        zs.next_out = out;
        zs.avail_out = OUT_SIZE;
        ret = zng_inflate(&zs, Z_NO_FLUSH);
        if (ret != Z_OK && ret != Z_STREAM_END && ret != Z_BUF_ERROR) {
            fprintf(stderr, "zng_inflate error %d\n", ret);
            exit(2);
        }
        size_t produced = zs.next_out - out;
        r.out += produced;
        if (verify) r.crc = zng_crc32(r.crc, out, (uint32_t)produced);
        done = 0;
        if (ret == Z_STREAM_END) {
            r.members++;
            done = 1;
            if (zs.avail_in == 0 && pos >= n) break;
            zng_inflateReset(&zs);
        }
    }
    if (!done) { fprintf(stderr, "zng: truncated input\n"); exit(2); }
    zng_inflateEnd(&zs);
    return r;
}

int main(int argc, char **argv)
{
    if (argc < 3) { fprintf(stderr, "usage: %s isal|zng FILE [reps] [verify]\n", argv[0]); return 1; }
    int use_isal = strcmp(argv[1], "isal") == 0;
    int reps = argc > 3 ? atoi(argv[3]) : 3;
    int verify = argc > 4 ? atoi(argv[4]) : 0;
    size_t n;
    uint8_t *file = load(argv[2], &n);
    uint8_t *in = malloc(IN_CHUNK), *out = malloc(OUT_SIZE);
    for (int i = 0; i < reps; i++) {
        double w0 = now(CLOCK_MONOTONIC), c0 = now(CLOCK_PROCESS_CPUTIME_ID);
        struct result r = use_isal ? run_isal(file, n, in, out, verify) : run_zng(file, n, in, out, verify);
        double w1 = now(CLOCK_MONOTONIC), c1 = now(CLOCK_PROCESS_CPUTIME_ID);
        printf("%s\trep=%d\tin=%zu\tout=%llu\tmembers=%d\twall_s=%.4f\tcpu_s=%.4f\twall_MBps=%.1f\tcpu_MBps=%.1f%s",
               argv[1], i, n, (unsigned long long)r.out, r.members, w1 - w0, c1 - c0,
               r.out / 1e6 / (w1 - w0), r.out / 1e6 / (c1 - c0), verify ? "" : "\n");
        if (verify) printf("\tcrc32=%08x\n", r.crc);
        fflush(stdout);
    }
    return 0;
}
