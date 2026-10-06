// Behaviour checks for the ISA-L igzip API (inflate wrappers, multi-member
// streams, return codes, stateless calls, deflate_stateless limits, CRC32).
// usage: test_semantics TEXTFILE
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <isa-l/igzip_lib.h>
#include <isa-l/crc.h>
#include <libdeflate.h>
#include <zlib-ng.h>

extern uint32_t crc32_gzip_refl_base(uint32_t, uint8_t *, uint64_t);
extern uint32_t crc32_gzip_refl_by8(uint32_t, const uint8_t *, uint64_t);
extern uint32_t crc32_gzip_refl_avx2(uint32_t, const uint8_t *, uint64_t);

static int fails = 0;
#define CHECK(cond, ...)                                                   \
    do {                                                                   \
        if (!(cond)) { printf("  FAIL: "); printf(__VA_ARGS__); printf("\n"); fails++; } \
    } while (0)

static uint8_t *text;
static size_t textlen;

// ---------- gzip member construction ----------
enum { HDR_PLAIN, HDR_NAME, HDR_BGZF, HDR_ALL };
static const char *hdrname[] = { "no flags", "FNAME", "FEXTRA (BGZF)", "FEXTRA+FNAME+FCOMMENT+FHCRC" };

static size_t make_member(uint8_t *out, int hdr, const uint8_t *data, size_t len)
{
    size_t p = 0;
    uint8_t flg = 0;
    if (hdr == HDR_NAME) flg = 0x08;
    if (hdr == HDR_BGZF) flg = 0x04;
    if (hdr == HDR_ALL) flg = 0x04 | 0x08 | 0x10 | 0x02;
    uint8_t base[10] = { 0x1f, 0x8b, 8, flg, 0x11, 0x22, 0x33, 0x44, 0, 3 };
    memcpy(out, base, 10);
    p = 10;
    size_t bsize_pos = 0;
    if (flg & 0x04) {
        if (hdr == HDR_BGZF) {
            uint8_t x[8] = { 6, 0, 'B', 'C', 2, 0, 0, 0 };
            memcpy(out + p, x, 8);
            bsize_pos = p + 6;
            p += 8;
        } else {
            uint8_t x[2 + 23] = { 23, 0 };
            for (int i = 0; i < 23; i++) x[2 + i] = (uint8_t)(i * 7 + 1);
            memcpy(out + p, x, sizeof(x));
            p += sizeof(x);
        }
    }
    if (flg & 0x08) { const char *nm = "some_reads_file_name.fq"; memcpy(out + p, nm, strlen(nm) + 1); p += strlen(nm) + 1; }
    if (flg & 0x10) { const char *cm = "a gzip comment that is longer than a few bytes"; memcpy(out + p, cm, strlen(cm) + 1); p += strlen(cm) + 1; }
    if (flg & 0x02) {
        uint32_t c = zng_crc32(0, out, (uint32_t)p);
        out[p++] = c & 0xff;
        out[p++] = (c >> 8) & 0xff;
    }
    struct libdeflate_compressor *c = libdeflate_alloc_compressor(6);
    size_t cl = libdeflate_deflate_compress(c, data, len, out + p, 1 << 20);
    libdeflate_free_compressor(c);
    p += cl;
    uint32_t crc = zng_crc32(0, data, (uint32_t)len);
    for (int i = 0; i < 4; i++) out[p++] = (crc >> (8 * i)) & 0xff;
    for (int i = 0; i < 4; i++) out[p++] = ((uint32_t)len >> (8 * i)) & 0xff;
    if (bsize_pos) { out[bsize_pos] = (p - 1) & 0xff; out[bsize_pos + 1] = ((p - 1) >> 8) & 0xff; }
    return p;
}

// Stream-decode buf with isal_inflate (crc_flag ISAL_GZIP), feeding `chunk` bytes per
// call and offering `outchunk` bytes of output per call (at a deliberately odd address).
// Returns 0 on success (output == expect), else the failing isal ret code or 99.
static int stream_decode(const uint8_t *buf, size_t n, size_t chunk, size_t outchunk,
                         const uint8_t *expect, size_t explen, int *members_out)
{
    struct inflate_state *st = malloc(sizeof(*st));
    uint8_t *outbuf = malloc(explen + outchunk + 64);
    uint8_t *inbuf = malloc(chunk + 16);
    size_t pos = 0, op = 0;
    int members = 0, ret = 0;
    isal_inflate_init(st);
    st->crc_flag = ISAL_GZIP;
    st->avail_in = 0;
    long guard = 0;
    for (;;) {
        if (++guard > 400000000L) { ret = 98; break; }
        if (st->avail_in == 0) {
            if (pos >= n) break;
            size_t c = n - pos < chunk ? n - pos : chunk;
            memcpy(inbuf + 1, buf + pos, c);   // odd input address too
            pos += c;
            st->next_in = inbuf + 1;
            st->avail_in = (uint32_t)c;
        }
        uint8_t *o = outbuf + 3 + op;          // unaligned output pointer
        size_t room = explen + outchunk - op;
        st->next_out = o;
        st->avail_out = (uint32_t)(room < outchunk ? room : outchunk);
        ret = isal_inflate(st);
        if (ret != ISAL_DECOMP_OK) break;
        op += st->next_out - o;
        if (st->block_state == ISAL_BLOCK_FINISH) {
            members++;
            if (st->avail_in == 0 && pos >= n) break;
            isal_inflate_reset(st);
        }
    }
    if (ret == 0 && st->block_state != ISAL_BLOCK_FINISH) ret = 97;   // truncated
    if (ret == 0 && (op != explen || memcmp(outbuf + 3, expect, explen) != 0)) ret = 99;
    *members_out = members;
    free(st);
    free(outbuf);
    free(inbuf);
    return ret;
}

// Same as stream_decode, but the caller parses every member header with
// isal_read_gzip_header into a persistent struct (crc_flag ISAL_GZIP_NO_HDR_VER).
static int stream_decode_hdr(const uint8_t *buf, size_t n, size_t chunk, const uint8_t *expect, size_t explen, int *members_out)
{
    struct inflate_state *st = malloc(sizeof(*st));
    struct isal_gzip_header gh;
    uint8_t *outbuf = malloc(explen + 64);
    uint8_t *inbuf = malloc(chunk + 16);
    size_t pos = 0, op = 0;
    int members = 0, ret = 0, in_header = 1;
    isal_inflate_init(st);
    st->crc_flag = ISAL_GZIP_NO_HDR_VER;
    isal_gzip_header_init(&gh);
    st->avail_in = 0;
    for (;;) {
        if (st->avail_in == 0) {
            if (pos >= n) break;
            size_t c = n - pos < chunk ? n - pos : chunk;
            memcpy(inbuf, buf + pos, c);
            pos += c;
            st->next_in = inbuf;
            st->avail_in = (uint32_t)c;
        }
        if (in_header) {
            ret = isal_read_gzip_header(st, &gh);
            if (ret == ISAL_END_INPUT) { ret = 0; continue; }   // header continues in the next read
            if (ret != 0) break;
            in_header = 0;
        }
        st->next_out = outbuf + op;
        st->avail_out = (uint32_t)(explen + 64 - op);
        ret = isal_inflate(st);
        if (ret != ISAL_DECOMP_OK) break;
        op = st->next_out - outbuf;
        if (st->block_state == ISAL_BLOCK_FINISH) {
            members++;
            if (st->avail_in == 0 && pos >= n) break;
            isal_inflate_reset(st);
            isal_gzip_header_init(&gh);
            in_header = 1;
        }
    }
    if (ret == 0 && st->block_state != ISAL_BLOCK_FINISH) ret = 97;
    if (ret == 0 && (op != explen || memcmp(outbuf, expect, explen) != 0)) ret = 99;
    *members_out = members;
    free(st); free(outbuf); free(inbuf);
    return ret;
}

static void test_streaming(void)
{
    {
        printf("== 3b: 5-member stream, headers parsed by the caller (isal_read_gzip_header + ISAL_GZIP_NO_HDR_VER)\n");
        uint8_t *mm = malloc(8 << 20), *expect = malloc(8 << 20);
        size_t n = 0, e = 0, off = 0;
        size_t sl[5] = { 70000, 1234, 65280, 300000, 17 };
        int hd[5] = { HDR_NAME, HDR_PLAIN, HDR_BGZF, HDR_ALL, HDR_NAME };
        for (int i = 0; i < 5; i++) {
            n += make_member(mm + n, hd[i], text + off, sl[i]);
            memcpy(expect + e, text + off, sl[i]);
            e += sl[i];
            off += sl[i];
        }
        size_t chunks[] = { 1, 2, 3, 5, 7, 11, 13, 64, 1000, 4096, 131072, 1 << 30 };
        int nbad = 0;
        for (size_t ci = 0; ci < sizeof(chunks) / sizeof(*chunks); ci++) {
            int members;
            int r = stream_decode_hdr(mm, n, chunks[ci], expect, e, &members);
            if (r != 0 || members != 5) { printf("  chunk %zu: ret %d members %d\n", chunks[ci], r, members); nbad++; }
        }
        printf("  %s\n", nbad ? "FAILURES above" : "ok for input chunks 1,2,3,5,7,11,13,64,1000,4096,128Ki,all");
        free(mm); free(expect);
    }
    printf("== 3a/3b: isal_inflate, crc_flag ISAL_GZIP, multi-member, byte-wise feeding\n");
    uint8_t *mm = malloc(8 << 20), *expect = malloc(8 << 20);
    size_t n = 0, e = 0;
    size_t sl[5] = { 70000, 1234, 65280, 300000, 17 };
    int hd[5] = { HDR_NAME, HDR_PLAIN, HDR_BGZF, HDR_ALL, HDR_NAME };
    size_t off = 0;
    for (int i = 0; i < 5; i++) {
        n += make_member(mm + n, hd[i], text + off, sl[i]);
        memcpy(expect + e, text + off, sl[i]);
        e += sl[i];
        off += sl[i];
    }
    size_t chunks[] = { 1, 2, 3, 5, 7, 11, 13, 64, 1000, 4096, 131072, 1 << 30 };
    size_t outs[] = { 1, 7, 1 << 20 };
    for (size_t ci = 0; ci < sizeof(chunks) / sizeof(*chunks); ci++)
        for (size_t oi = 0; oi < 3; oi++) {
            if (chunks[ci] < 64 && outs[oi] < 8) continue;   // too slow, covered by other combinations
            int members;
            int r = stream_decode(mm, n, chunks[ci], outs[oi], expect, e, &members);
            CHECK(r == 0 && members == 5, "5-member stream, in-chunk %zu out-chunk %zu: ret %d members %d", chunks[ci], outs[oi], r, members);
        }
    // Each header type alone, fed 1..16 bytes at a time
    for (int h = 0; h < 4; h++) {
        size_t m = make_member(mm, h, text, 5000);
        int bad = 0, firstbad = 0, badret = 0;
        for (size_t c = 1; c <= 16; c++) {
            int members;
            int r = stream_decode(mm, m, c, 1 << 20, text, 5000, &members);
            if (r != 0) { if (!bad) { firstbad = (int)c; badret = r; } bad++; }
        }
        printf("  single member, header %-28s: %s", hdrname[h], bad ? "" : "ok for input chunks of 1..16 bytes\n");
        if (bad) printf("FAILS for %d of 16 chunk sizes (first: chunk %d, ret %d)\n", bad, firstbad, badret);
    }
    // Same, but header split exactly at each offset of the HDR_ALL header (2 reads: k bytes, rest)
    {
        size_t m = make_member(mm, HDR_ALL, text, 5000);
        int hdrlen = 10 + 25 + 24 + 47 + 2;
        printf("  HDR_ALL member split into two reads at offset k (header is %d bytes):", hdrlen);
        int nbad = 0;
        for (int k = 1; k <= hdrlen + 2; k++) {
            struct inflate_state *st = malloc(sizeof(*st));
            uint8_t *ob = malloc(6000);
            isal_inflate_init(st);
            st->crc_flag = ISAL_GZIP;
            st->next_in = mm; st->avail_in = k; st->next_out = ob; st->avail_out = 6000;
            int r1 = isal_inflate(st);
            st->next_in = mm + k; st->avail_in = (uint32_t)(m - k);
            int r2 = r1 ? r1 : isal_inflate(st);
            int ok = r2 == 0 && st->block_state == ISAL_BLOCK_FINISH && st->total_out == 5000 && !memcmp(ob, text, 5000);
            if (!ok) { printf(" k=%d(ret %d)", k, r2); nbad++; }
            free(st); free(ob);
        }
        printf(nbad ? "  <- failing split points\n" : " all ok\n");
    }
    // Header supplied by the caller with isal_read_gzip_header, then crc_flag ISAL_GZIP_NO_HDR_VER
    {
        size_t m = make_member(mm, HDR_ALL, text, 5000);
        int nbad = 0;
        for (size_t c = 1; c <= 16; c++) {
            struct inflate_state *st = malloc(sizeof(*st));
            struct isal_gzip_header gh;
            uint8_t *ob = malloc(6000);
            isal_gzip_header_init(&gh);
            isal_inflate_init(st);
            st->crc_flag = ISAL_GZIP_NO_HDR_VER;
            size_t pos = 0;
            int r = ISAL_END_INPUT;
            st->avail_in = 0;
            while (r == ISAL_END_INPUT && pos < m) {
                size_t k = m - pos < c ? m - pos : c;
                st->next_in = mm + pos; st->avail_in = (uint32_t)k; pos += k;
                r = isal_read_gzip_header(st, &gh);
            }
            st->next_out = ob; st->avail_out = 6000;
            while (r == 0 && st->block_state != ISAL_BLOCK_FINISH) {
                if (st->avail_in == 0) {
                    if (pos >= m) break;
                    size_t k = m - pos < c ? m - pos : c;
                    st->next_in = mm + pos; st->avail_in = (uint32_t)k; pos += k;
                }
                r = isal_inflate(st);
            }
            if (!(r == 0 && st->block_state == ISAL_BLOCK_FINISH && st->total_out == 5000)) nbad++;
            free(st); free(ob);
        }
        printf("  HDR_ALL via persistent isal_read_gzip_header + ISAL_GZIP_NO_HDR_VER, chunks 1..16: %s\n", nbad ? "FAILS" : "ok");
    }

    printf("== 3a/3c: return codes\n");
    {
        size_t m = make_member(mm, HDR_NAME, text, 100000);
        struct inflate_state *st = malloc(sizeof(*st));
        uint8_t *ob = malloc(200000);
        int r;
        // need more input
        isal_inflate_init(st); st->crc_flag = ISAL_GZIP;
        st->next_in = mm; st->avail_in = (uint32_t)(m / 2); st->next_out = ob; st->avail_out = 200000;
        r = isal_inflate(st);
        printf("  half the input: ret %d, avail_in %u, block_state %d (FINISH=%d)\n", r, st->avail_in, st->block_state, ISAL_BLOCK_FINISH);
        // output full
        isal_inflate_init(st); st->crc_flag = ISAL_GZIP;
        st->next_in = mm; st->avail_in = (uint32_t)m; st->next_out = ob; st->avail_out = 1000;
        r = isal_inflate(st);
        printf("  1000-byte output buffer: ret %d, avail_out %u, avail_in %u, block_state %d\n", r, st->avail_out, st->avail_in, st->block_state);
        // missing trailer bytes
        isal_inflate_init(st); st->crc_flag = ISAL_GZIP;
        st->next_in = mm; st->avail_in = (uint32_t)(m - 3); st->next_out = ob; st->avail_out = 200000;
        r = isal_inflate(st);
        printf("  input ends 3 bytes before the end of the trailer: ret %d, block_state %d (CHECKSUM_CHECK=%d), total_out %u\n", r, st->block_state, ISAL_CHECKSUM_CHECK, st->total_out);
        st->next_in = mm + m - 3; st->avail_in = 3;
        r = isal_inflate(st);
        printf("    ... then the last 3 bytes: ret %d, block_state %d\n", r, st->block_state);
        struct { const char *what; size_t at; uint8_t x; } bad[] = {
            { "CRC32 byte flipped", m - 8, 0x01 }, { "ISIZE byte flipped", m - 2, 0x01 },
            { "ID2 wrong (0x8c)", 1, 0x07 }, { "CM = 7", 2, 0x0f }, { "deflate data byte flipped", 400, 0x5a } };
        for (size_t i = 0; i < sizeof(bad) / sizeof(*bad); i++) {
            mm[bad[i].at] ^= bad[i].x;
            isal_inflate_init(st); st->crc_flag = ISAL_GZIP;
            st->next_in = mm; st->avail_in = (uint32_t)m; st->next_out = ob; st->avail_out = 200000;
            r = isal_inflate(st);
            printf("  %-26s: isal_inflate ret %d\n", bad[i].what, r);
            mm[bad[i].at] ^= bad[i].x;
        }
        size_t m2 = make_member(mm, HDR_ALL, text, 1000);
        mm[m2 - 1000 / 2] ^= 0;   // no-op, keep layout
        // header CRC16 wrong: flip a comment byte
        mm[10 + 25 + 24 + 3] ^= 0x20;
        isal_inflate_init(st); st->crc_flag = ISAL_GZIP;
        st->next_in = mm; st->avail_in = (uint32_t)m2; st->next_out = ob; st->avail_out = 200000;
        r = isal_inflate(st);
        printf("  %-26s: isal_inflate ret %d\n", "FHCRC mismatch", r);
        // trailing zeros after a complete member: reset + continue
        m = make_member(mm, HDR_NAME, text, 1000);
        memset(mm + m, 0, 16);
        isal_inflate_init(st); st->crc_flag = ISAL_GZIP;
        st->next_in = mm; st->avail_in = (uint32_t)m + 16; st->next_out = ob; st->avail_out = 200000;
        r = isal_inflate(st);
        printf("  member + 16 zero bytes: ret %d, block_state %d, avail_in %u (expect 16: next_in is just past the trailer)\n", r, st->block_state, st->avail_in);
        isal_inflate_reset(st);
        r = isal_inflate(st);
        printf("    ... isal_inflate_reset and continue on the zeros: ret %d\n", r);
        free(st); free(ob);
    }
    free(mm); free(expect);
}

static void test_stateless(void)
{
    printf("== 3e: isal_inflate_stateless\n");
    uint8_t *mm = malloc(1 << 20), *ob = malloc(1 << 20);
    struct inflate_state *st = malloc(sizeof(*st));
    size_t len = 65280;
    struct libdeflate_compressor *c = libdeflate_alloc_compressor(6);
    size_t cl = libdeflate_deflate_compress(c, text, len, mm, 1 << 20);
    uint32_t crc = zng_crc32(0, text, (uint32_t)len);
    int r;
    const char *flagname[] = { "ISAL_DEFLATE", "ISAL_GZIP", "ISAL_GZIP_NO_HDR", "", "", "", "ISAL_GZIP_NO_HDR_VER" };
    int flags[] = { ISAL_DEFLATE, ISAL_GZIP_NO_HDR, ISAL_GZIP_NO_HDR_VER };
    for (int i = 0; i < 3; i++) {
        isal_inflate_init(st);
        st->crc_flag = flags[i];
        st->next_in = mm; st->avail_in = (uint32_t)cl; st->next_out = ob; st->avail_out = 1 << 20;
        r = isal_inflate_stateless(st);
        printf("  raw deflate only, %-20s: ret %d, total_out %u, crc %08x (crc32 %08x), avail_in left %u\n", flagname[flags[i]], r, st->total_out, st->crc, crc, st->avail_in);
    }
    // with trailer appended
    for (int k = 0; k < 4; k++) mm[cl + k] = (crc >> (8 * k)) & 0xff;
    for (int k = 0; k < 4; k++) mm[cl + 4 + k] = ((uint32_t)len >> (8 * k)) & 0xff;
    for (int i = 1; i < 3; i++) {
        isal_inflate_init(st);
        st->crc_flag = flags[i];
        st->next_in = mm; st->avail_in = (uint32_t)cl + 8 + 5; st->next_out = ob; st->avail_out = 1 << 20;
        r = isal_inflate_stateless(st);
        printf("  deflate + 8-byte trailer + 5 extra bytes, %-20s: ret %d, crc %08x, avail_in left %u\n", flagname[flags[i]], r, st->crc, st->avail_in);
    }
    mm[cl] ^= 1;
    isal_inflate_init(st);
    st->crc_flag = ISAL_GZIP_NO_HDR_VER;
    st->next_in = mm; st->avail_in = (uint32_t)cl + 8; st->next_out = ob; st->avail_out = 1 << 20;
    r = isal_inflate_stateless(st);
    printf("  wrong CRC in trailer, ISAL_GZIP_NO_HDR_VER: ret %d\n", r);
    mm[cl] ^= 1;
    isal_inflate_init(st);
    st->crc_flag = ISAL_DEFLATE;
    st->next_in = mm; st->avail_in = (uint32_t)cl; st->next_out = ob; st->avail_out = (uint32_t)len - 1;
    r = isal_inflate_stateless(st);
    printf("  avail_out one byte short: ret %d\n", r);
    isal_inflate_init(st);
    st->crc_flag = ISAL_DEFLATE;
    st->next_in = mm; st->avail_in = (uint32_t)cl / 2; st->next_out = ob; st->avail_out = 1 << 20;
    r = isal_inflate_stateless(st);
    printf("  half the deflate input: ret %d\n", r);
    // whole gzip member (BGZF header) with ISAL_GZIP
    size_t m = make_member(mm, HDR_BGZF, text, len);
    isal_inflate_init(st);
    st->crc_flag = ISAL_GZIP;
    st->next_in = mm; st->avail_in = (uint32_t)m + 28; st->next_out = ob; st->avail_out = 1 << 20;
    r = isal_inflate_stateless(st);
    printf("  BGZF member + 28 bytes, ISAL_GZIP: ret %d, total_out %u, crc %08x, avail_in left %u\n", r, st->total_out, st->crc, st->avail_in);
    // stateless call without isal_inflate_init (garbage state)
    memset(st, 0xa5, sizeof(*st));
    st->crc_flag = ISAL_GZIP;
    st->next_in = mm; st->avail_in = (uint32_t)m; st->next_out = ob; st->avail_out = 1 << 20;
    r = isal_inflate_stateless(st);
    printf("  same, state memset to 0xa5 instead of isal_inflate_init: ret %d, total_out %u\n", r, st->total_out);
    libdeflate_free_compressor(c);
    free(mm); free(ob); free(st);
}

static size_t isal_comp(int level, uint8_t *lb, uint32_t lbsize, const uint8_t *in, size_t len, uint8_t *out, uint32_t avail_out, int *ret, int gzip_flag)
{
    struct isal_zstream s;
    isal_deflate_stateless_init(&s);
    s.level = level;
    s.level_buf = lb;
    s.level_buf_size = lbsize;
    s.gzip_flag = gzip_flag;
    s.flush = NO_FLUSH;
    s.end_of_stream = 1;
    s.next_in = (uint8_t *)in;
    s.avail_in = (uint32_t)len;
    s.next_out = out;
    s.avail_out = avail_out;
    *ret = isal_deflate_stateless(&s);
    return s.total_out;
}

static void test_deflate(void)
{
    printf("== 3f: isal_deflate_stateless\n");
    uint32_t lbs[4] = { 0, ISAL_DEF_LVL1_DEFAULT, ISAL_DEF_LVL2_DEFAULT, ISAL_DEF_LVL3_DEFAULT };
    printf("  ISAL_DEF_LVLx_MIN: %u %u %u, _DEFAULT (=_LARGE): %u %u %u, _EXTRA_LARGE: %u %u %u\n",
           ISAL_DEF_LVL1_MIN, ISAL_DEF_LVL2_MIN, ISAL_DEF_LVL3_MIN, lbs[1], lbs[2], lbs[3],
           ISAL_DEF_LVL1_EXTRA_LARGE, ISAL_DEF_LVL2_EXTRA_LARGE, ISAL_DEF_LVL3_EXTRA_LARGE);
    uint8_t *rnd = malloc(65280), *out = malloc(1 << 20), *out2 = malloc(1 << 20), *dec = malloc(1 << 20);
    uint64_t x = 88172645463325252ULL;
    for (int i = 0; i < 65280; i++) { x ^= x << 13; x ^= x >> 7; x ^= x << 17; rnd[i] = (uint8_t)x; }
    struct libdeflate_decompressor *d = libdeflate_alloc_decompressor();
    for (int level = 0; level <= 3; level++) {
        uint8_t *lb = level ? malloc(lbs[level]) : NULL;
        int ret;
        size_t z = isal_comp(level, lb, lbs[level], rnd, 65280, out, 65510, &ret, IGZIP_DEFLATE);
        size_t act = 0;
        int dr = libdeflate_deflate_decompress(d, out, z, dec, 65280, &act);
        printf("  L%d random 65280 bytes, avail_out 65510: ret %d, size %zu, first byte 0x%02x (BTYPE %d), decodes %s\n", level, ret, z, out[0], (out[0] >> 1) & 3,
               dr == 0 && act == 65280 && !memcmp(dec, rnd, 65280) ? "ok" : "BAD");
        z = isal_comp(level, lb, lbs[level], rnd, 65280, out, 65285, &ret, IGZIP_DEFLATE);
        printf("     avail_out 65285 (= stored size): ret %d size %zu; ", ret, z);
        z = isal_comp(level, lb, lbs[level], rnd, 65280, out, 65284, &ret, IGZIP_DEFLATE);
        printf("avail_out 65284: ret %d; ", ret);
        z = isal_comp(level, lb, lbs[level], text, 65280, out, 1000, &ret, IGZIP_DEFLATE);
        printf("FASTQ block into 1000 bytes: ret %d\n", ret);
        if (level) {
            // determinism vs. level_buf history/garbage and vs. total_in
            size_t z1 = isal_comp(level, lb, lbs[level], text + 5 * 65280, 65280, out, 65510, &ret, IGZIP_DEFLATE);
            for (int b = 0; b < 20; b++) isal_comp(level, lb, lbs[level], text + b * 65280 + 17, 65280, out2, 65510, &ret, IGZIP_DEFLATE);
            size_t z2 = isal_comp(level, lb, lbs[level], text + 5 * 65280, 65280, out2, 65510, &ret, IGZIP_DEFLATE);
            int same_dirty = z1 == z2 && !memcmp(out, out2, z1);
            memset(lb, 0x5c, lbs[level]);
            z2 = isal_comp(level, lb, lbs[level], text + 5 * 65280, 65280, out2, 65510, &ret, IGZIP_DEFLATE);
            int same_garbage = z1 == z2 && !memcmp(out, out2, z1);
            uint8_t *lb2 = malloc(ISAL_DEF_LVL3_EXTRA_LARGE);
            uint32_t sz2 = level == 1 ? ISAL_DEF_LVL1_EXTRA_LARGE : level == 2 ? ISAL_DEF_LVL2_EXTRA_LARGE : ISAL_DEF_LVL3_EXTRA_LARGE;
            z2 = isal_comp(level, lb2, sz2, text + 5 * 65280, 65280, out2, 65510, &ret, IGZIP_DEFLATE);
            int same_bigger = z1 == z2 && !memcmp(out, out2, z1);
            uint32_t sz3 = level == 1 ? ISAL_DEF_LVL1_SMALL : level == 2 ? ISAL_DEF_LVL2_SMALL : ISAL_DEF_LVL3_SMALL;
            size_t z3 = isal_comp(level, lb2, sz3, text + 5 * 65280, 65280, out2, 65510, &ret, IGZIP_DEFLATE);
            int same_small = z1 == z3 && !memcmp(out, out2, z1);
            // stream reused without isal_deflate_stateless_init: total_in left at 98765
            struct isal_zstream s;
            isal_deflate_stateless_init(&s);
            s.level = level; s.level_buf = lb; s.level_buf_size = lbs[level]; s.gzip_flag = IGZIP_DEFLATE;
            s.flush = NO_FLUSH; s.end_of_stream = 1; s.total_in = 98765;
            s.next_in = text + 5 * 65280; s.avail_in = 65280; s.next_out = out2; s.avail_out = 65510;
            isal_deflate_stateless(&s);
            int same_totin = s.total_out == z1 && !memcmp(out, out2, z1);
            printf("     same bytes for the same block: after other blocks in the same level_buf %s, level_buf full of garbage %s,\n"
                   "       EXTRA_LARGE level_buf %s (size %zu), SMALL level_buf %s (size %zu vs %zu), total_in=98765 on entry %s (size %zu)\n",
                   same_dirty ? "yes" : "NO", same_garbage ? "yes" : "NO", same_bigger ? "yes" : "NO", z2,
                   same_small ? "yes" : "NO", z3, z1, same_totin ? "yes" : "NO", (size_t)s.total_out);
            free(lb2);
        }
        if (level == 1) {
            size_t z1 = isal_comp(1, lb, lbs[1], text, 65280, out, 65510, &ret, IGZIP_DEFLATE);
            size_t z0 = isal_comp(1, NULL, 0, text, 65280, out2, 65510, &ret, IGZIP_DEFLATE);
            printf("     L1 with level_buf NULL: ret %d size %zu (with DEFAULT level_buf: %zu)\n", ret, z0, z1);
        }
        free(lb);
    }
    // gzip_flag IGZIP_GZIP and IGZIP_GZIP_NO_HDR
    {
        uint8_t *lb = malloc(lbs[1]);
        int ret;
        size_t z = isal_comp(1, lb, lbs[1], text, 65280, out, 70000, &ret, IGZIP_GZIP);
        printf("  IGZIP_GZIP: ret %d, size %zu, header bytes:", ret, z);
        for (int i = 0; i < 10; i++) printf(" %02x", out[i]);
        size_t act = 0;
        int dr = libdeflate_gzip_decompress(d, out, z, dec, 65280, &act);
        printf(", libdeflate_gzip_decompress %d\n", dr);
        free(lb);
    }
    libdeflate_free_decompressor(d);
    free(rnd); free(out); free(out2); free(dec);
}

static void test_crc(void)
{
    printf("== 3g: crc32_gzip_refl vs zng_crc32 / libdeflate_crc32\n");
    uint8_t *b = malloc(1 << 20);
    uint64_t x = 1;
    for (int i = 0; i < (1 << 20); i++) { x = x * 6364136223846793005ULL + 1442695040888963407ULL; b[i] = (uint8_t)(x >> 56); }
    int bad = 0, n = 0;
    for (int len = 0; len < 3000; len += (len < 300 ? 1 : 97))
        for (int off = 0; off < 16; off += 5)
            for (int k = 0; k < 3; k++) {
                uint32_t init = k == 0 ? 0 : k == 1 ? 0xffffffffu : 0x12345678u;
                uint32_t a = crc32_gzip_refl(init, b + off, len);
                uint32_t z = zng_crc32(init, b + off, len);
                uint32_t l = libdeflate_crc32(init, b + off, len);
                uint32_t a0 = crc32_gzip_refl_base(init, b + off, len);
                uint32_t a1 = crc32_gzip_refl_by8(init, b + off, len);
                uint32_t a2 = crc32_gzip_refl_avx2(init, b + off, len);
                n++;
                if (a != z || a != l || a != a0 || a != a1 || a != a2) bad++;
            }
    uint32_t whole = crc32_gzip_refl(0, b, 1 << 20);
    uint32_t chained = crc32_gzip_refl(crc32_gzip_refl(crc32_gzip_refl(0, b, 12345), b + 12345, 1), b + 12346, (1 << 20) - 12346);
    printf("  %d (init, buffer, length) cases: %d mismatches (dispatched, _base, _by8, _avx2 vs zng_crc32 and libdeflate_crc32)\n", n, bad);
    printf("  chaining: whole %08x, three chained pieces %08x, zng %08x\n", whole, chained, zng_crc32(0, b, 1 << 20));
    printf("  len 0: crc32_gzip_refl(0x1234, NULL, 0) = %08x, zng_crc32(0x1234, NULL, 0) = %08x, zng_crc32(0x1234, buf, 0) = %08x\n",
           crc32_gzip_refl(0x1234, NULL, 0), zng_crc32(0x1234, NULL, 0), zng_crc32(0x1234, b, 0));
    free(b);
}

int main(int argc, char **argv)
{
    FILE *f = fopen(argc > 1 ? argv[1] : "reads.fq", "rb");
    if (!f) { perror("text"); return 1; }
    textlen = 4 << 20;
    text = malloc(textlen);
    textlen = fread(text, 1, textlen, f);
    fclose(f);
    test_streaming();
    test_stateless();
    test_deflate();
    test_crc();
    printf("\n%d check failures\n", fails);
    return 0;
}
