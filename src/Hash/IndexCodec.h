// IndexCodec.h - the compressed form of the index (index.prx.zst).
//
// The raw index (Seedmap::Save) is a header, a key map and a values array. Compressed, it is a
// seekable zstd file (Zstd.h) whose frames hold the same content in columns, which compress much
// better than the raw bytes and load in parallel, one chunk per thread:
//
//   frame 0   container: "PRXSPLT1", version, the raw index header, the layout, the chunk table
//   frame i   chunk i-1: a run of control blocks and their value cells:
//               filled   one bit per block: whether it holds value cells (most blocks of a small
//                        database are empty; they cost only this bit)
//               counts   value cells per key of the filled blocks (uint16, 2 byte planes); the key
//                        map's block offsets and per-key offsets are prefix sums of these
//               flex     the flex-key cells (keys with n >= 2 cells start with (n + 2) / 3 of
//                        them, see Seedmap::FlexBlockSize), 8 byte planes
//               entries  taxid, gene id and position (byte planes of the widths this chunk needs)
//                        and the 4 flag bits (1 byte), in key order
//
// An entry cell is taxid << 40 | gene << 20 | position | flags << 60. The writer decodes every
// chunk it encodes and stores the chunk as raw cells if that does not give exactly its bytes, so
// any index round-trips. plain zstd -d does not give a raw index.prx; protal --decompress_db does.
#pragma once

#include "Zstd.h"
#include "TargetClones.h"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <optional>
#include <string>
#include <vector>

namespace protal::index_codec {
    inline constexpr char kMagic[8] = {'P', 'R', 'X', 'S', 'P', 'L', 'T', '1'};
    inline constexpr uint32_t kVersion = 1;
    inline constexpr size_t kChunkHeaderBytes = 32;
    inline constexpr uint8_t kModeRaw = 0, kModeSplit = 1;

    // An index's shape: `blocks` control blocks of ctrl_cells + keys_per_block uint16 cells (the
    // block's uint64 value offset, then one uint16 offset per key), a final control block with the
    // number of value cells, and `values` uint64 value cells.
    struct Layout {
        uint64_t blocks = 0;
        uint64_t keys_per_block = 8;
        uint64_t ctrl_cells = 4;
        uint64_t values = 0;
        uint64_t CellsPerBlock() const { return keys_per_block + ctrl_cells; }
        uint64_t KeymapCells() const { return blocks * CellsPerBlock() + ctrl_cells; }
    };

    struct Chunk {
        uint64_t first_block = 0, blocks = 0, first_value = 0, values = 0;
    };

    struct Container {
        std::string index_header;  // the raw index's header bytes (Seedmap::HeaderBytes)
        Layout layout;
        std::vector<Chunk> chunks;
    };

    namespace detail {
        inline void PutU64(std::vector<char>& out, uint64_t v) {
            char b[8];
            std::memcpy(b, &v, 8);
            out.insert(out.end(), b, b + 8);
        }

        inline uint64_t GetU64(char const* p) {
            uint64_t v;
            std::memcpy(&v, p, 8);
            return v;
        }

        inline uint64_t BlockStart(uint16_t const* keymap, Layout const& l, uint64_t block) {
            uint64_t v;
            std::memcpy(&v, keymap + block * l.CellsPerBlock(), 8);
            return v;
        }

        inline uint64_t FlexCells(uint64_t n) { return n >= 2 ? (n + 2) / 3 : 0; }

        // Bytes needed for a field value (the fields have 20 bits).
        inline int Width(uint64_t max) { return max == 0 ? 0 : max < (1u << 8) ? 1 : max < (1u << 16) ? 2 : 3; }

        // Appends the low `width` bytes of get(i), i < n, as byte planes.
        template<typename Get>
        void PutPlanes(std::vector<char>& out, size_t n, int width, Get&& get) {
            size_t const base = out.size();
            out.resize(base + n * static_cast<size_t>(width));
            for (int b = 0; b < width; b++) {
                char* plane = out.data() + base + static_cast<size_t>(b) * n;
                for (size_t i = 0; i < n; i++) plane[i] = static_cast<char>(get(i) >> (8 * b));
            }
        }

        // dst[i] |= value(first + i) << shift for i < count, with value read from `width` byte planes
        // of n bytes each starting at planes.
        inline void OrPlanes(uint64_t* dst, size_t count, unsigned char const* planes, size_t n, int width, size_t first,
                             int shift) {
            for (int b = 0; b < width; b++) {
                unsigned char const* p = planes + static_cast<size_t>(b) * n + first;
                int const s = shift + 8 * b;
                for (size_t i = 0; i < count; i++) dst[i] |= uint64_t(p[i]) << s;
            }
        }

        // n (a multiple of 4) empty control blocks of 8 keys from cells: each the value offset v and 8 zero key offsets.
        inline void FillEmptyBlocks(uint16_t* cells, uint64_t n, uint64_t v) {
            unsigned char pattern[96] = {};  // 4 blocks of 24 bytes
            for (int i = 0; i < 4; i++) std::memcpy(pattern + 24 * i, &v, 8);
            auto* out = reinterpret_cast<unsigned char*>(cells);
            for (uint64_t i = 0; i < n; i += 4) std::memcpy(out + 24 * i, pattern, sizeof(pattern));
        }

        // Decodes a chunk payload into the key map cells of its blocks (km) and its value cells (vals).
        PROTAL_CLONE_V3 inline std::string DecodeChunk(char const* data, size_t size, Layout const& l, Chunk const& c, uint16_t* km,
                                       uint64_t* vals) {
            uint64_t const cpb = l.CellsPerBlock(), kpb = l.keys_per_block;
            if (size < kChunkHeaderBytes) return "chunk shorter than its header";
            auto const* u = reinterpret_cast<unsigned char const*>(data);
            uint8_t const mode = u[0];
            if (mode == kModeRaw) {
                uint64_t const expected = kChunkHeaderBytes + 2 * c.blocks * cpb + 8 * c.values;
                if (size != expected) return "raw chunk of " + std::to_string(size) + " bytes, expected " + std::to_string(expected);
                std::memcpy(km, data + kChunkHeaderBytes, 2 * c.blocks * cpb);
                if (c.values > 0) std::memcpy(vals, data + kChunkHeaderBytes + 2 * c.blocks * cpb, 8 * c.values);
                return "";
            }
            if (mode != kModeSplit) return "unknown chunk mode " + std::to_string(mode);
            int const tw = u[1], gw = u[2], pw = u[3];
            if (tw > 3 || gw > 3 || pw > 3) return "invalid field widths";
            uint64_t const entries = GetU64(data + 8), flex = GetU64(data + 16), filled = GetU64(data + 24);
            if (entries > c.values || flex > c.values || filled > c.blocks) return "chunk counts exceed its cells";
            uint64_t const keys = filled * kpb, bitmap_bytes = (c.blocks + 7) / 8;
            uint64_t const expected = kChunkHeaderBytes + bitmap_bytes + 2 * keys + 8 * flex +
                                      static_cast<uint64_t>(tw + gw + pw + 1) * entries;
            if (size != expected) return "chunk of " + std::to_string(size) + " bytes, expected " + std::to_string(expected);
            unsigned char const* bitmap = u + kChunkHeaderBytes;
            unsigned char const* lo = bitmap + bitmap_bytes;
            unsigned char const* hi = lo + keys;
            unsigned char const* flex_planes = hi + keys;
            unsigned char const* tax = flex_planes + 8 * flex;
            unsigned char const* gene = tax + static_cast<size_t>(tw) * entries;
            unsigned char const* pos = gene + static_cast<size_t>(gw) * entries;
            unsigned char const* flags = pos + static_cast<size_t>(pw) * entries;
            auto count = [&](uint64_t k) { return uint64_t(lo[k]) | uint64_t(hi[k]) << 8; };

            // Key map: block offsets and per-key offsets are prefix sums of the counts; empty blocks
            // have all key offsets 0.
            uint64_t v = c.first_value, flex_sum = 0, entry_sum = 0, k = 0;
            for (uint64_t b = 0; b < c.blocks; b++) {
                uint16_t* cells = km + b * cpb;
                if (cpb == 12 && (b & 63) == 0 && b + 64 <= c.blocks) {
                    uint64_t word;
                    std::memcpy(&word, bitmap + (b >> 3), 8);
                    if (word == 0) {  // 64 empty blocks (most blocks of a small database): the same 24 bytes each
                        FillEmptyBlocks(cells, 64, v);
                        b += 63;
                        continue;
                    }
                }
                if (!(bitmap[b >> 3] >> (b & 7) & 1)) {
                    if (kpb == 8) {  // the only layout protal uses: one fixed-size store
                        uint16_t block[12] = {};
                        std::memcpy(block, &v, 8);
                        std::memcpy(cells, block, sizeof(block));
                    } else {
                        std::memcpy(cells, &v, 8);
                        std::fill_n(cells + l.ctrl_cells, kpb, uint16_t{0});
                    }
                    continue;
                }
                if (k == keys) return "chunk has more filled blocks than its header says";
                std::memcpy(cells, &v, 8);
                uint64_t local = 0;
                for (uint64_t j = 0; j < kpb; j++, k++) {
                    if (local > 0xffff) return "block " + std::to_string(c.first_block + b) + " holds more than 65535 cells";
                    uint64_t const n = count(k);
                    cells[l.ctrl_cells + j] = static_cast<uint16_t>(local);
                    local += n;
                    flex_sum += FlexCells(n);
                    entry_sum += n - FlexCells(n);
                }
                v += local;
            }
            if (k != keys) return "chunk has fewer filled blocks than its header says";
            if (v != c.first_value + c.values) return "chunk counts add up to " + std::to_string(v - c.first_value) +
                                                      " cells, expected " + std::to_string(c.values);
            if (flex_sum != flex || entry_sum != entries) return "chunk counts do not match its flex cells and entries";

            // Values, in key order of the filled blocks' keys: entries and flex cells are composed in
            // batches (cache-sized), then placed per key.
            constexpr uint64_t kBatch = 1 << 16;
            std::vector<uint64_t> tmp_flex, tmp_entries;
            uint64_t* dst = vals;
            uint64_t fi = 0, ei = 0;
            k = 0;
            while (k < keys) {
                uint64_t end = k, cells = 0, fcount = 0, ecount = 0;
                while (end < keys) {
                    uint64_t const n = count(end);
                    if (cells + n > kBatch && end > k) break;
                    cells += n;
                    fcount += FlexCells(n);
                    ecount += n - FlexCells(n);
                    end++;
                }
                if (cells > 0) {
                    tmp_flex.assign(fcount, 0);
                    OrPlanes(tmp_flex.data(), fcount, flex_planes, flex, 8, fi, 0);
                    tmp_entries.resize(ecount);
                    for (uint64_t i = 0; i < ecount; i++) tmp_entries[i] = uint64_t(flags[ei + i]) << 60;
                    OrPlanes(tmp_entries.data(), ecount, tax, entries, tw, ei, 40);
                    OrPlanes(tmp_entries.data(), ecount, gene, entries, gw, ei, 20);
                    OrPlanes(tmp_entries.data(), ecount, pos, entries, pw, ei, 0);
                    uint64_t fp = 0, ep = 0;
                    for (uint64_t key = k; key < end; key++) {
                        uint64_t const n = count(key), f = FlexCells(n);
                        if (n == 0) continue;
                        // Keys of one cell have no flex cells, and tmp_flex has no data pointer while it is
                        // empty. n - f >= 1, so tmp_entries is never empty here.
                        if (f) std::memcpy(dst, tmp_flex.data() + fp, f * 8);
                        std::memcpy(dst + f, tmp_entries.data() + ep, (n - f) * 8);
                        dst += n;
                        fp += f;
                        ep += n - f;
                    }
                    fi += fcount;
                    ei += ecount;
                }
                k = end;
            }
            return "";
        }

        // Encodes chunk c of the index (key map km, values vals) into out, as raw cells if the split
        // form does not decode to exactly the same bytes. Returns true if split.
        inline bool EncodeChunk(Layout const& l, Chunk const& c, uint16_t const* km, uint64_t const* vals,
                                std::vector<char>& out) {
            uint64_t const cpb = l.CellsPerBlock(), kpb = l.keys_per_block;
            std::vector<unsigned char> bitmap((c.blocks + 7) / 8, 0);
            std::vector<uint16_t> counts;
            std::vector<uint64_t> flex;
            std::vector<uint64_t> entries;
            bool canonical = true;
            uint64_t v = c.first_value;
            for (uint64_t b = c.first_block; b < c.first_block + c.blocks && canonical; b++) {
                uint16_t const* cells = km + b * cpb;
                uint64_t const start = BlockStart(km, l, b), end = BlockStart(km, l, b + 1);
                if (start != v || end < start) canonical = false;
                if (canonical && end == start) {
                    // An empty block: canonical if all its key offsets are 0.
                    for (uint64_t j = 0; j < kpb; j++) canonical = canonical && cells[l.ctrl_cells + j] == 0;
                    continue;
                }
                uint64_t const local_block = b - c.first_block;
                bitmap[local_block >> 3] |= static_cast<unsigned char>(1u << (local_block & 7));
                for (uint64_t j = 0; j < kpb && canonical; j++) {
                    uint64_t const offset = cells[l.ctrl_cells + j];
                    uint64_t const next = j + 1 < kpb ? cells[l.ctrl_cells + j + 1] : end - start;
                    if (next < offset || start + offset != v || next - offset > 0xffff) {
                        canonical = false;
                        break;
                    }
                    uint64_t const n = next - offset, f = FlexCells(n);
                    counts.push_back(static_cast<uint16_t>(n));
                    flex.insert(flex.end(), vals + v, vals + v + f);
                    entries.insert(entries.end(), vals + v + f, vals + v + n);
                    v += n;
                }
                if (canonical && v != end) canonical = false;
            }
            if (v != c.first_value + c.values) canonical = false;

            if (canonical) {
                uint64_t max_tax = 0, max_gene = 0, max_pos = 0;
                for (uint64_t e : entries) {
                    max_tax = std::max(max_tax, (e >> 40) & 0xFFFFF);
                    max_gene = std::max(max_gene, (e >> 20) & 0xFFFFF);
                    max_pos = std::max(max_pos, e & 0xFFFFF);
                }
                int const tw = Width(max_tax), gw = Width(max_gene), pw = Width(max_pos);
                out.assign(kChunkHeaderBytes, '\0');
                out[0] = static_cast<char>(kModeSplit);
                out[1] = static_cast<char>(tw);
                out[2] = static_cast<char>(gw);
                out[3] = static_cast<char>(pw);
                uint64_t const ne = entries.size(), nf = flex.size(), filled = counts.size() / kpb;
                std::memcpy(out.data() + 8, &ne, 8);
                std::memcpy(out.data() + 16, &nf, 8);
                std::memcpy(out.data() + 24, &filled, 8);
                out.insert(out.end(), bitmap.begin(), bitmap.end());
                PutPlanes(out, counts.size(), 2, [&](size_t i) { return uint64_t(counts[i]); });
                PutPlanes(out, nf, 8, [&](size_t i) { return flex[i]; });
                PutPlanes(out, ne, tw, [&](size_t i) { return (entries[i] >> 40) & 0xFFFFF; });
                PutPlanes(out, ne, gw, [&](size_t i) { return (entries[i] >> 20) & 0xFFFFF; });
                PutPlanes(out, ne, pw, [&](size_t i) { return entries[i] & 0xFFFFF; });
                PutPlanes(out, ne, 1, [&](size_t i) { return entries[i] >> 60; });

                // It must decode to exactly the chunk's bytes.
                std::vector<uint16_t> km_check(c.blocks * cpb);
                std::vector<uint64_t> vals_check(c.values);
                if (DecodeChunk(out.data(), out.size(), l, c, km_check.data(), vals_check.data()).empty() &&
                    std::memcmp(km_check.data(), km + c.first_block * cpb, km_check.size() * 2) == 0 &&
                    (vals_check.empty() ||
                     std::memcmp(vals_check.data(), vals + c.first_value, vals_check.size() * 8) == 0)) {
                    return true;
                }
            }
            out.assign(kChunkHeaderBytes, '\0');
            out[0] = static_cast<char>(kModeRaw);
            auto const* km_bytes = reinterpret_cast<char const*>(km + c.first_block * cpb);
            out.insert(out.end(), km_bytes, km_bytes + 2 * c.blocks * cpb);
            if (c.values != 0) {
                auto const* val_bytes = reinterpret_cast<char const*>(vals + c.first_value);
                out.insert(out.end(), val_bytes, val_bytes + 8 * c.values);
            }
            return false;
        }

        // Chunks: runs of blocks with about target bytes each, counting the decoded key map cells
        // (which bound the reader's and writer's buffers) and <= 10 bytes per value cell.
        inline std::optional<std::vector<Chunk>> MakeChunks(Layout const& l, uint16_t const* km, uint64_t target, std::string& error) {
            std::vector<Chunk> chunks;
            uint64_t b = 0;
            while (b < l.blocks) {
                Chunk c{b, 0, BlockStart(km, l, b), 0};
                uint64_t bytes = 0;
                while (b < l.blocks && (c.blocks == 0 || bytes < target)) {
                    uint64_t const start = BlockStart(km, l, b), end = BlockStart(km, l, b + 1);
                    if (end < start) {
                        error = "the index key map is inconsistent at block " + std::to_string(b);
                        return std::nullopt;
                    }
                    bytes += 2 * l.CellsPerBlock() + 10 * (end - start);
                    c.blocks++;
                    b++;
                }
                c.values = BlockStart(km, l, b) - c.first_value;
                chunks.push_back(c);
            }
            if (BlockStart(km, l, l.blocks) != l.values) {
                error = "the index key map ends at value " + std::to_string(BlockStart(km, l, l.blocks)) + ", not " +
                        std::to_string(l.values);
                return std::nullopt;
            }
            return chunks;
        }

        inline std::vector<char> ContainerBytes(std::string const& header, Layout const& l, std::vector<Chunk> const& chunks) {
            std::vector<char> out(kMagic, kMagic + 8);
            PutU64(out, kVersion);
            PutU64(out, header.size());
            out.insert(out.end(), header.begin(), header.end());
            for (uint64_t v : {l.blocks, l.keys_per_block, l.ctrl_cells, l.values, uint64_t(chunks.size())}) PutU64(out, v);
            for (auto const& c : chunks) {
                for (uint64_t v : {c.first_block, c.blocks, c.first_value, c.values}) PutU64(out, v);
            }
            return out;
        }
    }

    // True if the (decompressed) file starts with the container magic, e.g. a split index whose
    // seek table was cut off.
    inline bool StartsWithMagic(std::string const& path) {
        zstd::InputFile in(path);
        char magic[8] = {};
        if (!in.IsOpen()) return false;
        in.Stream().read(magic, 8);
        return in.Stream().gcount() == 8 && std::memcmp(magic, kMagic, 8) == 0;
    }

    // The container (frame 0) of a split index; nullopt with error empty if the file is not one.
    inline std::optional<Container> ReadContainer(std::string const& path, zstd::SeekTable const& table, std::string& error) {
        if (table.frames.empty()) return std::nullopt;
        int const fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
        if (fd < 0) {
            error = std::string("cannot open the file: ") + std::strerror(errno);
            return std::nullopt;
        }
        ZSTD_DCtx* dctx = ZSTD_createDCtx();
        std::vector<char> input, data;
        std::string const read_error = zstd::ReadFrame(fd, table, 0, input, data, dctx);
        ZSTD_freeDCtx(dctx);
        ::close(fd);
        if (!read_error.empty()) {
            error = read_error;
            return std::nullopt;
        }
        if (data.size() < 8 || std::memcmp(data.data(), kMagic, 8) != 0) return std::nullopt;  // a plain seekable file
        auto fail = [&](std::string const& what) {
            error = "invalid container: " + what;
            return std::nullopt;
        };
        size_t pos = 8;
        auto take = [&](uint64_t& v) {
            if (pos + 8 > data.size()) return false;
            v = detail::GetU64(data.data() + pos);
            pos += 8;
            return true;
        };
        uint64_t version = 0, header_size = 0, chunk_count = 0;
        if (!take(version) || version != kVersion) return fail("unsupported version " + std::to_string(version));
        if (!take(header_size) || pos + header_size > data.size()) return fail("header size");
        Container c;
        c.index_header.assign(data.data() + pos, header_size);
        pos += header_size;
        if (!take(c.layout.blocks) || !take(c.layout.keys_per_block) || !take(c.layout.ctrl_cells) || !take(c.layout.values) ||
            !take(chunk_count)) {
            return fail("layout");
        }
        if (c.layout.ctrl_cells != 4 || c.layout.keys_per_block == 0) return fail("layout");
        if (chunk_count != table.frames.size() - 1) {
            return fail(std::to_string(chunk_count) + " chunks, but " + std::to_string(table.frames.size() - 1) + " chunk frames");
        }
        uint64_t next_block = 0, next_value = 0;
        for (uint64_t i = 0; i < chunk_count; i++) {
            Chunk ch;
            if (!take(ch.first_block) || !take(ch.blocks) || !take(ch.first_value) || !take(ch.values)) return fail("chunk table");
            if (ch.first_block != next_block || ch.first_value != next_value || ch.blocks == 0) return fail("chunks do not tile the index");
            next_block += ch.blocks;
            next_value += ch.values;
            c.chunks.push_back(ch);
        }
        if (next_block != c.layout.blocks || next_value != c.layout.values) return fail("chunks do not cover the index");
        return c;
    }

    inline bool IsSplitIndex(std::string const& path) {
        if (!zstd::IsCompressed(path)) return false;
        std::string error;
        auto const table = zstd::ReadSeekTable(path, error);
        return table && ReadContainer(path, *table, error).has_value();
    }

    // Writes an index (its header bytes, key map and values of layout l) as a split index at path,
    // with chunks of about params.frame_size bytes. raw_chunks (if given) receives the number of
    // chunks stored raw. Returns the bytes written, or nullopt with a message in error.
    inline std::optional<uint64_t> Write(std::string const& path, std::string const& header, Layout const& l,
                                         uint16_t const* keymap, uint64_t const* values, zstd::Params const& params,
                                         std::string& error, size_t* raw_chunks = nullptr) {
        uint64_t const target = params.frame_size > 0 ? params.frame_size : uint64_t{64} << 20;
        auto const chunks = detail::MakeChunks(l, keymap, target, error);
        if (!chunks) return std::nullopt;
        std::vector<char> const container = detail::ContainerBytes(header, l, *chunks);
        std::atomic<size_t> raw{0};
        auto const written = zstd::WriteSeekable(path, chunks->size() + 1, params, [&](size_t i, std::vector<char>& out) -> std::string {
            if (i == 0) {
                out = container;
                return "";
            }
            if (!detail::EncodeChunk(l, (*chunks)[i - 1], keymap, values, out)) raw++;
            return "";
        }, error);
        if (raw_chunks) *raw_chunks = raw;
        return written;
    }

    // Decodes a split index into keymap (l.KeymapCells() cells) and values (l.values cells).
    inline std::string Decode(std::string const& path, zstd::SeekTable const& table, Container const& c, uint16_t* keymap,
                              uint64_t* values, int threads) {
        Layout const& l = c.layout;
        std::string const error = zstd::ForEachFrame(path, table, 1, threads,
                [&](size_t frame, char const* data, size_t size, size_t) -> std::string {
            Chunk const& ch = c.chunks[frame - 1];
            std::string const e = detail::DecodeChunk(data, size, l, ch, keymap + ch.first_block * l.CellsPerBlock(),
                                                      values + ch.first_value);
            return e.empty() ? e : "chunk " + std::to_string(frame) + " of " + std::to_string(c.chunks.size()) + ": " + e;
        });
        if (error.empty()) std::memcpy(keymap + l.blocks * l.CellsPerBlock(), &l.values, 8);  // final control block
        return error;
    }

    // Reads a split index back chunk by chunk and compares it with keymap and values (no second
    // copy of the index in memory). Returns an error message, empty if identical.
    inline std::string Verify(std::string const& path, std::string const& header, Layout const& l, uint16_t const* keymap,
                              uint64_t const* values, int threads) {
        std::string error;
        auto const table = zstd::ReadSeekTable(path, error);
        if (!table) return error.empty() ? "not a seekable zstd file" : error;
        auto const c = ReadContainer(path, *table, error);
        if (!c) return error.empty() ? "not a split index" : error;
        if (c->index_header != header || c->layout.blocks != l.blocks || c->layout.values != l.values ||
            c->layout.keys_per_block != l.keys_per_block) {
            return "the container does not describe this index";
        }
        std::vector<std::vector<uint16_t>> km_tmp(zstd::WorkerCount(c->chunks.size(), threads));
        std::vector<std::vector<uint64_t>> val_tmp(km_tmp.size());
        return zstd::ForEachFrame(path, *table, 1, threads, [&](size_t frame, char const* data, size_t size, size_t worker) -> std::string {
            Chunk const& ch = c->chunks[frame - 1];
            km_tmp[worker].resize(ch.blocks * l.CellsPerBlock());
            val_tmp[worker].resize(ch.values);
            std::string const e = detail::DecodeChunk(data, size, l, ch, km_tmp[worker].data(), val_tmp[worker].data());
            if (!e.empty()) return "chunk " + std::to_string(frame) + ": " + e;
            if (std::memcmp(km_tmp[worker].data(), keymap + ch.first_block * l.CellsPerBlock(), km_tmp[worker].size() * 2) != 0 ||
                (!val_tmp[worker].empty() &&
                 std::memcmp(val_tmp[worker].data(), values + ch.first_value, val_tmp[worker].size() * 8) != 0)) {
                return "chunk " + std::to_string(frame) + " differs from the index";
            }
            return "";
        });
    }
}
