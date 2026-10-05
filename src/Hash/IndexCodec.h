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

    // An index's value cells in the file's 8-byte layout, as the writer and Verify read them: an array of
    // them, or readers that give them from another layout (--build holds the values packed,
    // Seedmap::PackedLayout): a key's cells (first and n are the key's first slot and its slots), and any
    // range (for chunks kept raw, rare).
    struct Cells {
        using Reader = void (*)(void const* context, uint64_t first, uint64_t n, uint64_t* out);
        uint64_t const* array = nullptr;
        void const* context = nullptr;
        Reader read_key = nullptr, read_range = nullptr;

        Cells(uint64_t const* cells) : array(cells) {}
        Cells(void const* reader_context, Reader key_reader, Reader range_reader)
                : context(reader_context), read_key(key_reader), read_range(range_reader) {}

        // The cells of the key at slots [first, first + n): the array's own, or read into buffer.
        uint64_t const* Key(uint64_t first, uint64_t n, std::vector<uint64_t>& buffer) const {
            return Read(read_key, first, n, buffer);
        }

        // Cells [first, first + n), any range.
        uint64_t const* Range(uint64_t first, uint64_t n, std::vector<uint64_t>& buffer) const {
            return Read(read_range, first, n, buffer);
        }

    private:
        uint64_t const* Read(Reader reader, uint64_t first, uint64_t n, std::vector<uint64_t>& buffer) const {
            if (array) return array + first;
            buffer.resize(std::max<uint64_t>(n, 1));
            if (n) reader(context, first, n, buffer.data());
            return buffer.data();
        }
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

        inline std::string const kValuesDiffer = "its values differ from the index";
        inline std::string const kSinkRefused = "its values were refused";

        // Takes a chunk's value cells as DecodeChunk makes them, instead of an array of them: `key` each key's n cells
        // (its flex cells, then its entries, as the 8-byte file layout has them) from its first slot, `range` all cells
        // of a chunk stored raw. False refuses them (DecodeChunk returns kSinkRefused). Seedmap packs the cells as they
        // come, without the chunk's ~50 MB of 8-byte cells on each loading thread (1.6 GB on 32 threads at r226).
        struct ValueSink {
            using Take = bool (*)(void* context, uint64_t first_slot, uint64_t const* cells, uint64_t n);
            Take key = nullptr, range = nullptr;
            void* context = nullptr;
        };

        // Decodes a chunk payload into the key map cells of its blocks (km) and its value cells (vals). With
        // `expected` (vals unused) the value cells are compared with the index's cells instead (kValuesDiffer if
        // any differs): a check of a chunk against the index needs no copy of its values (~50 MB per chunk at r226).
        // With `sink` (vals unused) they go to it (ValueSink).
        PROTAL_CLONE_V3 inline std::string DecodeChunk(char const* data, size_t size, Layout const& l, Chunk const& c, uint16_t* km,
                                       uint64_t* vals, Cells const* expected = nullptr, ValueSink const* sink = nullptr) {
            uint64_t const cpb = l.CellsPerBlock(), kpb = l.keys_per_block;
            if (size < kChunkHeaderBytes) return "chunk shorter than its header";
            auto const* u = reinterpret_cast<unsigned char const*>(data);
            uint8_t const mode = u[0];
            if (mode == kModeRaw) {
                uint64_t const expected_size = kChunkHeaderBytes + 2 * c.blocks * cpb + 8 * c.values;
                if (size != expected_size) return "raw chunk of " + std::to_string(size) + " bytes, expected " + std::to_string(expected_size);
                std::memcpy(km, data + kChunkHeaderBytes, 2 * c.blocks * cpb);
                char const* const raw_values = data + kChunkHeaderBytes + 2 * c.blocks * cpb;
                if (c.values > 0 && expected) {
                    std::vector<uint64_t> buffer;
                    uint64_t const* const want = expected->Range(c.first_value, c.values, buffer);
                    return std::memcmp(want, raw_values, 8 * c.values) == 0 ? "" : kValuesDiffer;
                }
                if (c.values > 0 && sink) {  // rare (a chunk the split form does not give exactly): an aligned copy
                    std::vector<uint64_t> cells(c.values);
                    std::memcpy(cells.data(), raw_values, 8 * c.values);
                    return sink->range(sink->context, c.first_value, cells.data(), c.values) ? "" : kSinkRefused;
                }
                if (c.values > 0) std::memcpy(vals, raw_values, 8 * c.values);
                return "";
            }
            if (mode != kModeSplit) return "unknown chunk mode " + std::to_string(mode);
            int const tw = u[1], gw = u[2], pw = u[3];
            if (tw > 3 || gw > 3 || pw > 3) return "invalid field widths";
            uint64_t const entries = GetU64(data + 8), flex = GetU64(data + 16), filled = GetU64(data + 24);
            if (entries > c.values || flex > c.values || filled > c.blocks) return "chunk counts exceed its cells";
            uint64_t const keys = filled * kpb, bitmap_bytes = (c.blocks + 7) / 8;
            uint64_t const expected_size = kChunkHeaderBytes + bitmap_bytes + 2 * keys + 8 * flex +
                                      static_cast<uint64_t>(tw + gw + pw + 1) * entries;
            if (size != expected_size) return "chunk of " + std::to_string(size) + " bytes, expected " + std::to_string(expected_size);
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
            uint64_t want = c.first_value;  // the first slot of the next key to compare or hand to the sink
            std::vector<uint64_t> buffer;
            // A key's f flex cells and n - f entries to the values, compared with the index's cells of the key, or,
            // side by side in `buffer`, to the sink.
            auto put = [&dst, &want, &buffer, expected, sink](uint64_t const* flex_cells, uint64_t f, uint64_t const* entry_cells, uint64_t n) {
                if (sink) {
                    buffer.resize(n);
                    if (f) std::memcpy(buffer.data(), flex_cells, f * 8);
                    std::memcpy(buffer.data() + f, entry_cells, (n - f) * 8);
                    if (!sink->key(sink->context, want, buffer.data(), n)) return false;
                    want += n;
                } else if (expected) {
                    uint64_t const* const key = expected->Key(want, n, buffer);
                    if ((f && std::memcmp(key, flex_cells, f * 8) != 0) || std::memcmp(key + f, entry_cells, (n - f) * 8) != 0) return false;
                    want += n;
                } else {
                    if (f) std::memcpy(dst, flex_cells, f * 8);
                    std::memcpy(dst + f, entry_cells, (n - f) * 8);
                    dst += n;
                }
                return true;
            };
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
                        if (!put(f ? tmp_flex.data() + fp : nullptr, f, tmp_entries.data() + ep, n)) return sink ? kSinkRefused : kValuesDiffer;
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

        // Encodes chunk c of the index (key map km, value cells vals) into out, as raw cells if the split
        // form does not decode to exactly the same bytes. Returns true if split. The planes are written
        // from the index's cells a key at a time and the check compares the decoded values with them, so
        // a chunk costs its encoded bytes and its key map cells, not two more copies of its values (64
        // threads held ~6 GB of such copies writing an r226 index).
        inline bool EncodeChunk(Layout const& l, Chunk const& c, uint16_t const* km, Cells const& vals,
                                std::vector<char>& out) {
            uint64_t const cpb = l.CellsPerBlock(), kpb = l.keys_per_block;
            std::vector<unsigned char> bitmap((c.blocks + 7) / 8, 0);
            std::vector<uint16_t> counts;
            std::vector<uint64_t> buffer;  // a key's cells, unless vals is an array
            uint64_t nf = 0, ne = 0, max_tax = 0, max_gene = 0, max_pos = 0;
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
                    nf += f;
                    ne += n - f;
                    uint64_t const* const key_cells = vals.Key(v, n, buffer);
                    for (uint64_t i = f; i < n; i++) {
                        max_tax = std::max(max_tax, (key_cells[i] >> 40) & 0xFFFFF);
                        max_gene = std::max(max_gene, (key_cells[i] >> 20) & 0xFFFFF);
                        max_pos = std::max(max_pos, key_cells[i] & 0xFFFFF);
                    }
                    v += n;
                }
                if (canonical && v != end) canonical = false;
            }
            if (v != c.first_value + c.values) canonical = false;

            if (canonical) {
                int const tw = Width(max_tax), gw = Width(max_gene), pw = Width(max_pos);
                uint64_t const keys = counts.size(), filled = keys / kpb;
                // Header, bitmap, then the planes: counts (2), flex cells (8), taxids, genes, positions, flags (1).
                size_t const counts_at = kChunkHeaderBytes + bitmap.size(), flex_at = counts_at + 2 * keys;
                size_t const tax_at = flex_at + 8 * nf, gene_at = tax_at + static_cast<size_t>(tw) * ne;
                size_t const pos_at = gene_at + static_cast<size_t>(gw) * ne, flags_at = pos_at + static_cast<size_t>(pw) * ne;
                out.assign(flags_at + ne, '\0');
                out[0] = static_cast<char>(kModeSplit);
                out[1] = static_cast<char>(tw);
                out[2] = static_cast<char>(gw);
                out[3] = static_cast<char>(pw);
                std::memcpy(out.data() + 8, &ne, 8);
                std::memcpy(out.data() + 16, &nf, 8);
                std::memcpy(out.data() + 24, &filled, 8);
                std::copy(bitmap.begin(), bitmap.end(), out.begin() + kChunkHeaderBytes);
                char* const lo = out.data() + counts_at;
                for (size_t k = 0; k < keys; k++) {
                    lo[k] = static_cast<char>(counts[k]);
                    lo[keys + k] = static_cast<char>(counts[k] >> 8);
                }
                char* const flex = out.data() + flex_at;
                char* const tax = out.data() + tax_at;
                char* const gene = out.data() + gene_at;
                char* const pos = out.data() + pos_at;
                char* const flags = out.data() + flags_at;
                size_t fi = 0, ei = 0;
                v = c.first_value;
                for (uint16_t const n : counts) {
                    uint64_t const f = FlexCells(n);
                    uint64_t const* const key_cells = vals.Key(v, n, buffer);
                    for (uint64_t i = 0; i < f; i++, fi++) {
                        for (int b = 0; b < 8; b++) flex[b * nf + fi] = static_cast<char>(key_cells[i] >> (8 * b));
                    }
                    for (uint64_t i = f; i < n; i++, ei++) {
                        uint64_t const e = key_cells[i];
                        for (int b = 0; b < tw; b++) tax[b * ne + ei] = static_cast<char>(((e >> 40) & 0xFFFFF) >> (8 * b));
                        for (int b = 0; b < gw; b++) gene[b * ne + ei] = static_cast<char>(((e >> 20) & 0xFFFFF) >> (8 * b));
                        for (int b = 0; b < pw; b++) pos[b * ne + ei] = static_cast<char>((e & 0xFFFFF) >> (8 * b));
                        flags[ei] = static_cast<char>(e >> 60);
                    }
                    v += n;
                }

                // It must decode to exactly the chunk's bytes.
                std::vector<uint16_t> km_check(c.blocks * cpb);
                if (DecodeChunk(out.data(), out.size(), l, c, km_check.data(), nullptr, &vals).empty() &&
                    std::memcmp(km_check.data(), km + c.first_block * cpb, km_check.size() * 2) == 0) {
                    return true;
                }
            }
            out.assign(kChunkHeaderBytes, '\0');
            out[0] = static_cast<char>(kModeRaw);
            auto const* km_bytes = reinterpret_cast<char const*>(km + c.first_block * cpb);
            out.insert(out.end(), km_bytes, km_bytes + 2 * c.blocks * cpb);
            if (c.values != 0) {
                auto const* val_bytes = reinterpret_cast<char const*>(vals.Range(c.first_value, c.values, buffer));
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
                                         uint16_t const* keymap, Cells const& values, zstd::Params const& params,
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

    // What each chunk of a split index writes when it is decoded: its key map cells (2 bytes each) and its values,
    // value_bits_32 / 32 bits each (64 in the file's layout, Seedmap::PackedLayout::SlotBits32 packed), for
    // ForEachFrame's budget (zstd::LoadBudget).
    inline zstd::LoadBudget ChunkOutput(Container const& c, uint64_t value_bits_32) {
        zstd::LoadBudget budget;
        budget.output.reserve(c.chunks.size());
        for (Chunk const& ch : c.chunks) {
            budget.output.push_back(ch.blocks * c.layout.CellsPerBlock() * sizeof(uint16_t) + ((ch.values * value_bits_32 >> 5) + 7) / 8);
        }
        return budget;
    }

    // Decodes a split index into keymap (l.KeymapCells() cells) and values (l.values cells), both memory that becomes
    // resident as it is written (calloc'd): near the end fewer threads decode (zstd::LoadBudget). chunks_on_fewer, if
    // given: the chunks decoded after a thread had stopped for memory.
    inline std::string Decode(std::string const& path, zstd::SeekTable const& table, Container const& c, uint16_t* keymap,
                              uint64_t* values, int threads, size_t* chunks_on_fewer = nullptr) {
        Layout const& l = c.layout;
        zstd::LoadBudget budget = ChunkOutput(c, 64 * 32);
        std::string const error = zstd::ForEachFrame(path, table, 1, threads,
                [&](size_t frame, char const* data, size_t size, size_t) -> std::string {
            Chunk const& ch = c.chunks[frame - 1];
            std::string const e = detail::DecodeChunk(data, size, l, ch, keymap + ch.first_block * l.CellsPerBlock(),
                                                      values + ch.first_value);
            return e.empty() ? e : "chunk " + std::to_string(frame) + " of " + std::to_string(c.chunks.size()) + ": " + e;
        }, &budget);
        if (error.empty()) std::memcpy(keymap + l.blocks * l.CellsPerBlock(), &l.values, 8);  // final control block
        if (chunks_on_fewer) *chunks_on_fewer = budget.frames_on_fewer;
        return error;
    }

    // Reads a split index back chunk by chunk and compares it with keymap and values (no second
    // copy of the index in memory, nor of a chunk's values: they are compared as they are decoded).
    // Returns an error message, empty if identical.
    inline std::string Verify(std::string const& path, std::string const& header, Layout const& l, uint16_t const* keymap,
                              Cells const& values, int threads) {
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
        return zstd::ForEachFrame(path, *table, 1, threads, [&](size_t frame, char const* data, size_t size, size_t worker) -> std::string {
            Chunk const& ch = c->chunks[frame - 1];
            km_tmp[worker].resize(ch.blocks * l.CellsPerBlock());
            std::string const e = detail::DecodeChunk(data, size, l, ch, km_tmp[worker].data(), nullptr, &values);
            if (e == detail::kValuesDiffer) return "chunk " + std::to_string(frame) + " differs from the index";
            if (!e.empty()) return "chunk " + std::to_string(frame) + ": " + e;
            if (std::memcmp(km_tmp[worker].data(), keymap + ch.first_block * l.CellsPerBlock(), km_tmp[worker].size() * 2) != 0) {
                return "chunk " + std::to_string(frame) + " differs from the index";
            }
            return "";
        });
    }
}
