// SPDX-License-Identifier: GPL-2.0-only
#include "GenomeStore.h"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <stdexcept>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zstd.h>

#include <isa-l/igzip_lib.h>

namespace fs = std::filesystem;

namespace protal::sim {

// ---- whole files -------------------------------------------------------------------------------------------------

namespace {

std::string ReadBytes(fs::path const& path) {
    int const fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
    if (fd < 0) throw std::runtime_error("cannot read " + path.string() + ": " + std::strerror(errno));
    struct stat st {};
    if (::fstat(fd, &st) != 0 || S_ISDIR(st.st_mode)) {
        std::string const why = S_ISDIR(st.st_mode) ? "it is a directory" : std::strerror(errno);
        ::close(fd);
        throw std::runtime_error("cannot read " + path.string() + ": " + why);
    }
    std::string bytes;
    bytes.resize(S_ISREG(st.st_mode) ? static_cast<std::size_t>(st.st_size) : std::size_t{1} << 20);
    std::size_t done = 0;
    for (;;) {
        if (done == bytes.size()) {
            if (S_ISREG(st.st_mode)) break;
            bytes.resize(bytes.size() * 2);  // a pipe: as long as it goes on
        }
        ssize_t const n = ::read(fd, bytes.data() + done, bytes.size() - done);
        if (n < 0 && errno == EINTR) continue;
        if (n < 0) {
            std::string const why = std::strerror(errno);
            ::close(fd);
            throw std::runtime_error("cannot read " + path.string() + ": " + why);
        }
        if (n == 0) break;
        done += static_cast<std::size_t>(n);
    }
    ::close(fd);
    bytes.resize(done);
    return bytes;
}

bool IsGzip(unsigned char const* in, std::size_t left) { return left >= 2 && in[0] == 0x1f && in[1] == 0x8b; }

// A zstd frame's magic number, or a skippable frame's (as ThreadedGzStream tells them).
bool IsZstd(unsigned char const* in, std::size_t left) {
    return left >= 4 && ((in[0] == 0x28 && in[1] == 0xb5 && in[2] == 0x2f && in[3] == 0xfd) ||
                         ((in[0] & 0xf0) == 0x50 && in[1] == 0x2a && in[2] == 0x4d && in[3] == 0x18));
}

char const* InflateError(int ret) {
    switch (ret) {
        case ISAL_INVALID_BLOCK: return "invalid block type";
        case ISAL_INVALID_SYMBOL: return "invalid code";
        case ISAL_INVALID_LOOKBACK: return "invalid distance too far back";
        case ISAL_INVALID_WRAPPER: return "incorrect header check";
        case ISAL_UNSUPPORTED_METHOD: return "unknown compression method";
        case ISAL_INCORRECT_CHECKSUM: return "incorrect data or length check";
        default: return "corrupt gzip data";
    }
}

constexpr std::size_t kChunk = std::size_t{1} << 30;  // ISA-L's counts are 32-bit

// Every gzip member of `in` (BGZF blocks are members too), each checked by its trailer; what follows a member must be
// another member or zero bytes of padding, as ThreadedGzStream reads it.
std::string Inflate(std::string const& in, fs::path const& path) {
    auto const* data = reinterpret_cast<unsigned char const*>(in.data());
    std::size_t const size = in.size();
    std::string out;
    // ISIZE of the last member: the size of a one-member file (mod 4 GB), at most 16 times the input (a damaged file's
    // last bytes can say anything); the buffer grows if it was not enough
    std::size_t const hint = size >= 4 ? (std::size_t{data[size - 4]} | std::size_t{data[size - 3]} << 8 |
                                          std::size_t{data[size - 2]} << 16 | std::size_t{data[size - 1]} << 24) : 0;
    out.resize(std::max(std::min(hint, size * 16) + 1, size * 4 + 64));
    std::size_t used = 0, pos = 0, members = 0;
    auto state = std::make_unique<inflate_state>();
    auto fail = [&](std::string const& why) {
        return std::runtime_error(path.string() + " is truncated or corrupt (" + why + ")");
    };
    while (pos < size) {
        if (!IsGzip(data + pos, size - pos)) {
            if (members > 0 && std::all_of(data + pos, data + size, [](unsigned char c) { return c == 0; })) break;
            throw fail("the data at byte " + std::to_string(pos) + ", after gzip member " + std::to_string(members) +
                       ", is no gzip member");
        }
        ++members;
        isal_inflate_init(state.get());
        state->crc_flag = ISAL_GZIP_NO_HDR_VER;
        isal_gzip_header header;
        isal_gzip_header_init(&header);  // no buffers: name, comment and extra field skipped
        state->next_in = const_cast<unsigned char*>(data + pos);
        state->avail_in = static_cast<std::uint32_t>(std::min(kChunk, size - pos));
        int ret = isal_read_gzip_header(state.get(), &header);
        if (ret != ISAL_DECOMP_OK) {
            throw fail(ret == ISAL_END_INPUT ? "the file ends inside a gzip header" : InflateError(ret));
        }
        for (;;) {
            std::size_t const consumed = reinterpret_cast<unsigned char const*>(state->next_in) - data;
            if (state->avail_in == 0 && consumed < size) {  // the next chunk of input
                state->next_in = const_cast<unsigned char*>(data + consumed);
                state->avail_in = static_cast<std::uint32_t>(std::min(kChunk, size - consumed));
            }
            if (used == out.size()) out.resize(out.size() * 2);
            state->next_out = reinterpret_cast<std::uint8_t*>(out.data() + used);
            state->avail_out = static_cast<std::uint32_t>(std::min(kChunk, out.size() - used));
            std::uint32_t const room = state->avail_out;
            std::uint32_t const had = state->avail_in;
            ret = isal_inflate(state.get());
            used += room - state->avail_out;
            if (ret != ISAL_DECOMP_OK) throw fail(std::string(InflateError(ret)) + " in gzip member " + std::to_string(members));
            if (state->block_state == ISAL_BLOCK_FINISH) break;
            if (state->avail_in == had && state->avail_out == room) {  // no progress although there is room: input ended
                throw fail("the file ends inside gzip member " + std::to_string(members));
            }
        }
        pos = reinterpret_cast<unsigned char const*>(state->next_in) - data;  // ISA-L gives back what it read past the end
    }
    out.resize(used);
    return out;
}

std::string Unzstd(std::string const& in, fs::path const& path) {
    std::unique_ptr<ZSTD_DCtx, std::size_t (*)(ZSTD_DCtx*)> dctx(ZSTD_createDCtx(), &ZSTD_freeDCtx);
    if (!dctx) throw std::runtime_error("cannot allocate a zstd decompression context");
    ZSTD_DCtx_setParameter(dctx.get(), ZSTD_d_windowLogMax, ZSTD_dParam_getBounds(ZSTD_d_windowLogMax).upperBound);
    unsigned long long const first = ZSTD_getFrameContentSize(in.data(), in.size());
    std::string out;
    bool const known = first != ZSTD_CONTENTSIZE_UNKNOWN && first != ZSTD_CONTENTSIZE_ERROR;
    out.resize(known ? std::min<std::size_t>(first, in.size() * 64) + 1 : in.size() * 4 + 64);
    ZSTD_inBuffer input{in.data(), in.size(), 0};
    std::size_t used = 0, ret = 0;
    while (input.pos < input.size) {
        if (used == out.size()) out.resize(out.size() * 2);
        ZSTD_outBuffer output{out.data() + used, out.size() - used, 0};
        ret = ZSTD_decompressStream(dctx.get(), &output, &input);
        used += output.pos;
        if (ZSTD_isError(ret)) {
            throw std::runtime_error(path.string() + " is truncated or corrupt (" + ZSTD_getErrorName(ret) + ")");
        }
    }
    while (ret != 0) {  // output still held back, or the input ended inside a frame
        if (used == out.size()) out.resize(out.size() * 2);
        ZSTD_outBuffer output{out.data() + used, out.size() - used, 0};
        ret = ZSTD_decompressStream(dctx.get(), &output, &input);
        used += output.pos;
        if (ZSTD_isError(ret)) {
            throw std::runtime_error(path.string() + " is truncated or corrupt (" + ZSTD_getErrorName(ret) + ")");
        }
        if (output.pos == 0 && ret != 0) {
            throw std::runtime_error(path.string() + " is truncated or corrupt (the file ends inside a zstd frame)");
        }
    }
    out.resize(used);
    return out;
}

}  // namespace

std::string ReadWholeFile(fs::path const& path) {
    std::string bytes = ReadBytes(path);
    auto const* data = reinterpret_cast<unsigned char const*>(bytes.data());
    if (IsGzip(data, bytes.size())) return Inflate(bytes, path);
    if (IsZstd(data, bytes.size())) return Unzstd(bytes, path);
    return bytes;  // neither: as it is
}

// ---- FASTA -------------------------------------------------------------------------------------------------------

namespace {

bool CSpace(unsigned char c) { return c == ' ' || (c >= '\t' && c <= '\r'); }  // isspace in the C locale

// A sequence line's characters, upper case (a-z), white space left out. Whole lines without white space are copied
// at once (most are); a line with some goes character by character.
void AppendSequence(std::string& seq, char const* line, std::size_t length) {
    std::size_t const old = seq.size();
    seq.resize(old + length);
    auto* out = reinterpret_cast<unsigned char*>(seq.data() + old);
    auto const* in = reinterpret_cast<unsigned char const*>(line);
    unsigned char low = 0;  // a byte of 32 or less: maybe white space
    for (std::size_t i = 0; i < length; ++i) {
        unsigned char const c = in[i];
        low |= static_cast<unsigned char>(c <= ' ');
        out[i] = static_cast<unsigned char>(c - (static_cast<unsigned char>(c - 'a') < 26 ? 32 : 0));
    }
    if (!low) return;
    seq.resize(old);
    for (std::size_t i = 0; i < length; ++i) {
        unsigned char const c = in[i];
        if (!CSpace(c)) seq.push_back(static_cast<char>(static_cast<unsigned char>(c - 'a') < 26 ? c - 32 : c));
    }
}

}  // namespace

FastaRecords ParseFasta(std::string_view text) {
    FastaRecords records;
    char const* const data = text.data();
    std::size_t const size = text.size(), none = static_cast<std::size_t>(-1);
    std::size_t pos = 0, current = none;
    while (pos < size) {
        char const* const line = data + pos;
        auto const* newline = static_cast<char const*>(std::memchr(line, '\n', size - pos));
        std::size_t const length = newline ? static_cast<std::size_t>(newline - line) : size - pos;
        pos += length + (newline ? 1 : 0);
        if (length > 0 && line[0] == '>') {
            std::size_t end = 1;
            while (end < length && line[end] != ' ' && line[end] != '\t' && line[end] != '\r') ++end;
            records.names.emplace_back(line + 1, end - 1);
            records.seqs.emplace_back();
            current = records.seqs.size() - 1;
            continue;
        }
        if (current != none) AppendSequence(records.seqs[current], line, length);
    }
    return records;
}

// ---- the genome store --------------------------------------------------------------------------------------------

// The file: a header, the FASTA's path, the contigs (Entry), their names, the runs of other characters than A, C, G and
// T (Run, by contig and start), then each contig's bases, 2 bits each (A 0, C 1, G 2, T 3; the first base in the low
// bits; a run's places 0), from an offset a multiple of 8.
namespace {

constexpr char kMagic[8] = {'p', 'r', 'o', 't', 'a', 'l', '2', 'b'};
constexpr std::uint32_t kVersion = 1;

struct Header {
    char magic[8];
    std::uint32_t version, contigs;
    std::uint64_t source_size;
    std::int64_t source_mtime;  // ns
    std::uint64_t path_bytes, names_bytes, runs, file_bytes;
};
static_assert(sizeof(Header) == 64);

std::uint64_t Align8(std::uint64_t n) { return (n + 7) & ~std::uint64_t{7}; }

// The FASTA's identity: its absolute path, size and modification time; false if it cannot be found.
bool SourceOf(fs::path const& fasta, std::string& path, std::uint64_t& size, std::int64_t& mtime) {
    std::error_code ec;
    path = fs::absolute(fasta, ec).lexically_normal().string();
    if (ec) path = fasta.string();
    struct stat st {};
    if (::stat(path.c_str(), &st) != 0) return false;
    size = static_cast<std::uint64_t>(st.st_size);
    mtime = static_cast<std::int64_t>(st.st_mtim.tv_sec) * 1'000'000'000 + st.st_mtim.tv_nsec;
    return true;
}

std::uint64_t Fnv1a(std::string_view text) {
    std::uint64_t h = 0xcbf29ce484222325ULL;
    for (unsigned char c : text) {
        h ^= c;
        h *= 0x100000001b3ULL;
    }
    return h;
}

// A, C, G, T as 0-3; any other character 4.
struct Coder {
    unsigned char code[256];
    Coder() {
        std::memset(code, 4, sizeof(code));
        code['A'] = 0;
        code['C'] = 1;
        code['G'] = 2;
        code['T'] = 3;
    }
};
Coder const kCoder;

// Four bases of a byte, as text.
struct Decoder {
    char quad[256][4];
    Decoder() {
        for (int b = 0; b < 256; ++b) {
            for (int i = 0; i < 4; ++i) quad[b][i] = "ACGT"[(b >> (2 * i)) & 3];
        }
    }
};
Decoder const kDecoder;

std::atomic<std::uint64_t> g_written{0};

}  // namespace

struct GenomeFile::Entry {
    std::uint64_t length, seq_offset, run_first, run_count, name_offset, name_length;
};

struct GenomeFile::Run {
    std::uint64_t start;
    std::uint32_t length;
    char base;
    char pad[3];
};

fs::path GenomeFile::PathOf(fs::path const& store, fs::path const& fasta) {
    std::error_code ec;
    std::string const path = fs::absolute(fasta, ec).lexically_normal().string();
    char name[32];
    std::snprintf(name, sizeof(name), "%016llx.g2b", static_cast<unsigned long long>(Fnv1a(ec ? fasta.string() : path)));
    return store / name;
}

bool GenomeFile::Write(fs::path const& store, fs::path const& fasta, FastaRecords const& records, std::string& error) {
    static_assert(sizeof(Entry) == 48 && sizeof(Run) == 16, "the file's layout");
    std::string source;
    std::uint64_t source_size = 0;
    std::int64_t source_mtime = 0;
    if (!SourceOf(fasta, source, source_size, source_mtime)) {
        error = "cannot find " + fasta.string();
        return false;
    }
    std::size_t const contigs = records.seqs.size();
    std::vector<Entry> entries(contigs);
    std::vector<Run> runs;
    std::string names;
    std::vector<std::string> packed(contigs);  // each contig's 2-bit bytes
    for (std::size_t k = 0; k < contigs; ++k) {
        std::string const& seq = records.seqs[k];
        Entry& e = entries[k];
        e.length = seq.size();
        e.name_offset = names.size();
        e.name_length = records.names[k].size();
        names += records.names[k];
        e.run_first = runs.size();
        auto const* in = reinterpret_cast<unsigned char const*>(seq.data());
        std::size_t const n = seq.size();
        std::string& bits = packed[k];
        bits.assign((n + 3) / 4, '\0');
        auto* out = reinterpret_cast<unsigned char*>(bits.data());
        auto exception = [&](std::size_t i) {  // a character other than A, C, G, T: a run of them (its place 0, 'A')
            auto const c = static_cast<char>(in[i]);
            if (runs.size() > e.run_first && runs.back().base == c && runs.back().start + runs.back().length == i &&
                runs.back().length < UINT32_MAX) {
                ++runs.back().length;
            } else {
                runs.push_back({i, 1, c, {0, 0, 0}});
            }
        };
        std::size_t i = 0;
        for (; i + 4 <= n; i += 4) {  // four bases a byte, a table lookup each
            unsigned const c0 = kCoder.code[in[i]], c1 = kCoder.code[in[i + 1]], c2 = kCoder.code[in[i + 2]],
                           c3 = kCoder.code[in[i + 3]];
            if ((c0 | c1 | c2 | c3) & 4) {
                for (std::size_t j = i; j < i + 4; ++j) {
                    unsigned const c = kCoder.code[in[j]];
                    if (c & 4) exception(j); else out[j >> 2] |= static_cast<unsigned char>(c << ((j & 3) * 2));
                }
                continue;
            }
            out[i >> 2] = static_cast<unsigned char>(c0 | c1 << 2 | c2 << 4 | c3 << 6);
        }
        for (; i < n; ++i) {
            unsigned const c = kCoder.code[in[i]];
            if (c & 4) exception(i); else out[i >> 2] |= static_cast<unsigned char>(c << ((i & 3) * 2));
        }
        e.run_count = runs.size() - e.run_first;
    }
    std::uint64_t at = sizeof(Header) + Align8(source.size());
    std::uint64_t const entries_at = at;
    at += contigs * sizeof(Entry);
    std::uint64_t const names_at = at;
    at = Align8(at + names.size());
    std::uint64_t const runs_at = at;
    at += runs.size() * sizeof(Run);
    for (auto& e : entries) {
        e.seq_offset = at;
        at = Align8(at + (e.length + 3) / 4);
    }
    std::string file(at, '\0');
    Header header{};
    std::memcpy(header.magic, kMagic, sizeof(kMagic));
    header.version = kVersion;
    header.contigs = static_cast<std::uint32_t>(contigs);
    header.source_size = source_size;
    header.source_mtime = source_mtime;
    header.path_bytes = source.size();
    header.names_bytes = names.size();
    header.runs = runs.size();
    header.file_bytes = at;
    std::memcpy(file.data(), &header, sizeof(header));
    std::memcpy(file.data() + sizeof(Header), source.data(), source.size());
    if (contigs) std::memcpy(file.data() + entries_at, entries.data(), contigs * sizeof(Entry));
    std::memcpy(file.data() + names_at, names.data(), names.size());
    if (!runs.empty()) std::memcpy(file.data() + runs_at, runs.data(), runs.size() * sizeof(Run));
    for (std::size_t k = 0; k < contigs; ++k) {
        std::memcpy(file.data() + entries[k].seq_offset, packed[k].data(), packed[k].size());
    }

    std::error_code ec;
    fs::create_directories(store, ec);
    fs::path const final_path = PathOf(store, fasta);
    fs::path const temp = final_path.string() + "." + std::to_string(::getpid()) + "." + std::to_string(g_written++) + ".tmp";
    std::FILE* out = std::fopen(temp.c_str(), "wb");
    if (!out) {
        error = "cannot write " + temp.string() + ": " + std::strerror(errno);
        return false;
    }
    bool ok = std::fwrite(file.data(), 1, file.size(), out) == file.size();
    std::string why = ok ? "" : std::strerror(errno);
    ok = (std::fclose(out) == 0) && ok;
    if (!ok) {
        error = "cannot write " + temp.string() + ": " + (why.empty() ? std::strerror(errno) : why);
        fs::remove(temp, ec);
        return false;
    }
    fs::rename(temp, final_path, ec);
    if (ec) {
        error = "cannot rename " + temp.string() + ": " + ec.message();
        fs::remove(temp, ec);
        return false;
    }
    return true;
}

std::shared_ptr<GenomeFile const> GenomeFile::Open(fs::path const& store, fs::path const& fasta) {
    std::string source;
    std::uint64_t source_size = 0;
    std::int64_t source_mtime = 0;
    if (!SourceOf(fasta, source, source_size, source_mtime)) return nullptr;
    fs::path const path = PathOf(store, fasta);
    int const fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
    if (fd < 0) return nullptr;
    struct stat st {};
    if (::fstat(fd, &st) != 0 || !S_ISREG(st.st_mode) || static_cast<std::size_t>(st.st_size) < sizeof(Header)) {
        ::close(fd);
        return nullptr;
    }
    auto const size = static_cast<std::size_t>(st.st_size);
    void* map = ::mmap(nullptr, size, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, 0);
    ::close(fd);
    if (map == MAP_FAILED) return nullptr;
    std::shared_ptr<GenomeFile> file(new GenomeFile());
    file->m_data = static_cast<char const*>(map);
    file->m_size = size;
    Header header{};
    std::memcpy(&header, file->m_data, sizeof(header));
    // made from this FASTA as it is now, and whole
    if (std::memcmp(header.magic, kMagic, sizeof(kMagic)) != 0 || header.version != kVersion ||
        header.file_bytes != size || header.source_size != source_size || header.source_mtime != source_mtime ||
        header.path_bytes != source.size() || sizeof(Header) + header.path_bytes > size ||
        std::memcmp(file->m_data + sizeof(Header), source.data(), source.size()) != 0) {
        return nullptr;
    }
    std::uint64_t const entries_at = sizeof(Header) + Align8(header.path_bytes);
    std::uint64_t const names_at = entries_at + std::uint64_t{header.contigs} * sizeof(Entry);
    std::uint64_t const runs_at = Align8(names_at + header.names_bytes);
    if (runs_at + header.runs * sizeof(Run) > size) return nullptr;
    file->m_contigs = header.contigs;
    file->m_entries = reinterpret_cast<Entry const*>(file->m_data + entries_at);
    file->m_names = file->m_data + names_at;
    file->m_runs = reinterpret_cast<Run const*>(file->m_data + runs_at);
    for (std::size_t k = 0; k < file->m_contigs; ++k) {
        Entry const& e = file->m_entries[k];
        if (e.name_offset + e.name_length > header.names_bytes || e.run_first + e.run_count > header.runs ||
            e.seq_offset < runs_at || e.seq_offset + (e.length + 3) / 4 > size) {
            return nullptr;
        }
    }
    return file;
}

GenomeFile::~GenomeFile() {
    if (m_data) ::munmap(const_cast<char*>(m_data), m_size);
}

GenomeFile::Entry const& GenomeFile::EntryOf(std::size_t contig) const { return m_entries[contig]; }

std::string_view GenomeFile::Name(std::size_t contig) const {
    Entry const& e = EntryOf(contig);
    return {m_names + e.name_offset, static_cast<std::size_t>(e.name_length)};
}

std::uint64_t GenomeFile::Length(std::size_t contig) const { return EntryOf(contig).length; }

void GenomeFile::Extract(std::size_t contig, std::uint64_t start, std::uint64_t length, std::string& out) const {
    Entry const& e = EntryOf(contig);
    if (start > e.length) throw std::out_of_range("GenomeFile::Extract: start beyond the contig");
    std::uint64_t const n = std::min(length, e.length - start);
    out.resize(n);
    auto const* bits = reinterpret_cast<unsigned char const*>(m_data + e.seq_offset);
    char* dst = out.data();
    std::uint64_t pos = start, i = 0;
    for (; i < n && (pos & 3); ++i, ++pos) dst[i] = "ACGT"[(bits[pos >> 2] >> ((pos & 3) * 2)) & 3];
    for (; i + 4 <= n; i += 4, pos += 4) std::memcpy(dst + i, kDecoder.quad[bits[pos >> 2]], 4);
    for (; i < n; ++i, ++pos) dst[i] = "ACGT"[(bits[pos >> 2] >> ((pos & 3) * 2)) & 3];
    if (e.run_count == 0 || n == 0) return;
    // the runs of other characters that reach into [start, start + n)
    Run const* first = m_runs + e.run_first;
    Run const* last = first + e.run_count;
    Run const* r = std::upper_bound(first, last, start, [](std::uint64_t at, Run const& run) { return at < run.start; });
    if (r != first) --r;
    std::uint64_t const end = start + n;
    for (; r != last && r->start < end; ++r) {
        std::uint64_t const from = std::max(r->start, start), to = std::min(r->start + r->length, end);
        if (from < to) std::memset(dst + (from - start), r->base, to - from);
    }
}

}  // namespace protal::sim
