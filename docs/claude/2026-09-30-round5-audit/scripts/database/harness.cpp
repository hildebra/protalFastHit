// Probes the seek-table allocation lead: does zstd::ReadFrame / ParallelReadFrames allocate the
// seek table's claimed decompressed_size before validating it against the frame?
//   - Case A (known content size, as protal writes): getFrameContentSize should catch the forged
//     size and error out WITHOUT the big allocation.
//   - Case B (unknown content size, hand-crafted/streamed frame): the check is skipped and the code
//     resizes to the claimed size before decompressing.
// Run under a small address-space ulimit so a multi-GB resize throws std::bad_alloc, proving it.
#include "Zstd.h"
#include <cstdio>
#include <string>
using namespace protal::zstd;

static std::string CompressKnown(std::string const& in) {
    ZSTD_CCtx* c = ZSTD_createCCtx();
    ZSTD_CCtx_setPledgedSrcSize(c, in.size());  // stores content size in the header
    std::string out(ZSTD_compressBound(in.size()) + 64, '\0');
    ZSTD_inBuffer ib{in.data(), in.size(), 0};
    ZSTD_outBuffer ob{out.data(), out.size(), 0};
    size_t r; do { r = ZSTD_compressStream2(c, &ob, &ib, ZSTD_e_end); } while (r > 0 && !ZSTD_isError(r));
    out.resize(ob.pos); ZSTD_freeCCtx(c); return out;
}
static std::string CompressUnknown(std::string const& in) {
    ZSTD_CCtx* c = ZSTD_createCCtx();
    ZSTD_CCtx_setParameter(c, ZSTD_c_contentSizeFlag, 0);  // omit content size from the header
    std::string out(ZSTD_compressBound(in.size()) + 64, '\0');
    ZSTD_inBuffer ib{in.data(), in.size(), 0};
    ZSTD_outBuffer ob{out.data(), out.size(), 0};
    size_t r; do { r = ZSTD_compressStream2(c, &ob, &ib, ZSTD_e_end); } while (r > 0 && !ZSTD_isError(r));
    out.resize(ob.pos); ZSTD_freeCCtx(c); return out;
}

static void WriteForged(std::string const& path, std::string const& frame, uint64_t forged_content) {
    FrameWriter w(path);
    w.Add(frame.data(), frame.size(), forged_content);  // seek table records forged_content
    w.Finish();
    if (!w.Error().empty()) { std::fprintf(stderr, "write error: %s\n", w.Error().c_str()); }
}

static void Probe(char const* label, std::string const& frame, uint64_t forged) {
    std::string const path = std::string("/tmp/forge_") + label + ".zst";
    WriteForged(path, frame, forged);
    unsigned long long cs = ZSTD_getFrameContentSize(frame.data(), frame.size());
    std::string err;
    auto table = ReadSeekTable(path, err);
    if (!table) { std::printf("[%s] ReadSeekTable failed: %s\n", label, err.c_str()); return; }
    int fd = ::open(path.c_str(), O_RDONLY);
    std::vector<char> in, out;
    ZSTD_DCtx* d = ZSTD_createDCtx();
    ZSTD_DCtx_setParameter(d, ZSTD_d_windowLogMax, 31);
    std::printf("[%s] header content size = %s ; seek table claims %llu ...\n", label,
                cs == ZSTD_CONTENTSIZE_UNKNOWN ? "UNKNOWN" : std::to_string(cs).c_str(),
                (unsigned long long)table->frames[0].decompressed_size);
    try {
        std::string e = ReadFrame(fd, *table, 0, in, out, d);
        std::printf("[%s] ReadFrame returned: \"%s\" ; out buffer resized to %zu bytes\n",
                    label, e.c_str(), out.size());
    } catch (std::exception const& ex) {
        std::printf("[%s] ReadFrame THREW: %s  (the resize to the claimed size was attempted)\n", label, ex.what());
    }
    ZSTD_freeDCtx(d); ::close(fd);
}

int main() {
    std::string payload(6370, 'x');  // small, like reference.map
    uint64_t const forged = 0xFFFFFFF0ull;  // ~4 GB, the max a 32-bit seek-table field allows
    std::printf("== Case A: frame written WITH content size (as protal writes)\n");
    Probe("known", CompressKnown(payload), forged);
    std::printf("== Case B: frame written WITHOUT content size (streamed/hand-crafted)\n");
    Probe("unknown", CompressUnknown(payload), forged);
    return 0;
}
