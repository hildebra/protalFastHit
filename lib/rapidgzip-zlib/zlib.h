// zlib.h for rapidgzip only (lib/rapidgzip.cmake): the part of zlib's API that rapidgzip's headers name,
// on protal's zlib-ng (lib/zlib-ng), which is built with its native API (zng_ names) and has no zlib.h.
// rapidgzip built with ISA-L (LIBRAPIDARCHIVE_WITH_ISAL) inflates and compresses with ISA-L; its zlib
// paths (gzip/zlib.hpp, inline) are compiled but not called. The functions are real functions, not
// macros, as zlib's are: rapidgzip has a namespace `deflate` and calls ::deflate and ::inflate.
// Only rapidgzip's include path has this directory; protal's own code includes zlib-ng.h.
#pragma once

#include <stdint.h>

#include <zlib-ng.h>

typedef zng_stream z_stream;
typedef zng_stream* z_streamp;
typedef zng_gz_header gz_header;
typedef zng_gz_header* gz_headerp;
typedef uint8_t Bytef;
typedef unsigned int uInt;
typedef unsigned long uLong;

static inline int deflateInit2(z_stream* strm, int level, int method, int windowBits, int memLevel, int strategy) {
    return zng_deflateInit2(strm, level, method, windowBits, memLevel, strategy);
}
static inline int deflateSetDictionary(z_stream* strm, const Bytef* dictionary, uInt dictLength) {
    return zng_deflateSetDictionary(strm, dictionary, dictLength);
}
static inline int deflate(z_stream* strm, int flush) { return zng_deflate(strm, flush); }
static inline int deflateEnd(z_stream* strm) { return zng_deflateEnd(strm); }

static inline int inflateInit2(z_stream* strm, int windowBits) { return zng_inflateInit2(strm, windowBits); }
static inline int inflate(z_stream* strm, int flush) { return zng_inflate(strm, flush); }
static inline int inflateEnd(z_stream* strm) { return zng_inflateEnd(strm); }
static inline int inflateSetDictionary(z_stream* strm, const Bytef* dictionary, uInt dictLength) {
    return zng_inflateSetDictionary(strm, dictionary, dictLength);
}
static inline int inflatePrime(z_stream* strm, int bits, int value) { return zng_inflatePrime(strm, bits, value); }
static inline int inflateReset2(z_stream* strm, int windowBits) { return zng_inflateReset2(strm, windowBits); }
static inline int inflateGetHeader(z_stream* strm, gz_headerp head) { return zng_inflateGetHeader(strm, head); }
