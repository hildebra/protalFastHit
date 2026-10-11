// BinaryTable.h - a table's arrays as one binary member of the single-file database, read straight into the table's
// arrays: no text parsed and no copy of the whole, decompressed frame by frame on several threads (column_weights.bin,
// strain_alleles.bin). At GTDB r226 size the column weights' text took 12 s of every run on one thread and ~7.6 GB
// beside the table while it was read (docs/claude/2026-10-11-performance-audit).
//
// A binary table is derived from its text table (column_weights.tsv, strain_alleles.tsv), which the database keeps
// beside it for the tools that read the text and for --unpack_db: the packer writes it (protal --build, --compress_db,
// --add_tables), the header records the text table's size, and a run uses the binary table only where the database's
// text table still has that size, as gene_table.bin is used (GeneTableFile.h). A database folder reads the text.
//
// Layout (little-endian, as protal's index): magic (8 bytes), version, the text table's size, the number of arrays,
// each array's size in bytes (uint64 each), then the arrays one after another, each from an offset that is a multiple
// of 8 (zero bytes between).
#pragma once

#include "Database.h"

#include <cstdint>
#include <cstring>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

namespace protal::binary_table {
    // An array to write: its bytes.
    struct ConstArray {
        char const* data = nullptr;
        uint64_t size = 0;
    };

    // An array to read: the memory it goes to, of the size the header gives.
    struct Array {
        char* data = nullptr;
        uint64_t size = 0;
    };

    template<typename T>
    ConstArray Of(std::vector<T> const& v) {
        return { reinterpret_cast<char const*>(v.data()), v.size() * sizeof(T) };
    }

    template<typename T>
    Array Into(std::vector<T>& v) {
        return { reinterpret_cast<char*>(v.data()), v.size() * sizeof(T) };
    }

    inline uint64_t Padded(uint64_t n) {
        return (n + 7) / 8 * 8;
    }

    namespace detail {
        inline void PutU64(std::ostream& os, uint64_t v) {
            char bytes[8];
            for (int i = 0; i < 8; i++) bytes[i] = static_cast<char>((v >> (8 * i)) & 0xff);
            os.write(bytes, 8);
        }

        inline uint64_t GetU64(char const* p) {
            uint64_t v = 0;
            for (int i = 0; i < 8; i++) v |= uint64_t(static_cast<unsigned char>(p[i])) << (8 * i);
            return v;
        }
    }

    // Writes a binary table; false on a write error.
    inline bool Write(std::ostream& os, char const (&magic)[8], uint64_t version, uint64_t text_size, std::vector<ConstArray> const& arrays) {
        os.write(magic, 8);
        detail::PutU64(os, version);
        detail::PutU64(os, text_size);
        detail::PutU64(os, arrays.size());
        for (auto const& a : arrays) detail::PutU64(os, a.size);
        uint64_t const header = 32 + 8 * arrays.size();
        static constexpr char kZeros[8] = {};
        os.write(kZeros, static_cast<std::streamsize>(Padded(header) - header));
        for (auto const& a : arrays) {
            if (a.size > 0) os.write(a.data, static_cast<std::streamsize>(a.size));
            os.write(kZeros, static_cast<std::streamsize>(Padded(a.size) - a.size));
        }
        return static_cast<bool>(os);
    }

    // A binary table's header: its version, the text table's size it was made from, its arrays' sizes and where they are.
    struct Header {
        uint64_t version = 0;
        uint64_t text_size = 0;
        std::vector<uint64_t> sizes;
        std::vector<uint64_t> offsets;  // each array's offset in the content
    };

    // The header of the binary table `file`; nullopt with error set if it is none (another magic, too short, or a size
    // that does not add up to the file's).
    inline std::optional<Header> ReadHeader(db::DbFile const& file, char const (&magic)[8], std::string& error) {
        static constexpr uint64_t kMaxArrays = 64;
        auto const start = file.ReadPrefix(32, error);
        if (!start) return std::nullopt;
        if (start->size() < 32 || std::memcmp(start->data(), magic, 8) != 0) {
            error = "not a binary table of this kind (wrong magic)";
            return std::nullopt;
        }
        Header h;
        h.version = detail::GetU64(start->data() + 8);
        h.text_size = detail::GetU64(start->data() + 16);
        uint64_t const count = detail::GetU64(start->data() + 24);
        if (count > kMaxArrays) {
            error = "invalid header (" + std::to_string(count) + " arrays)";
            return std::nullopt;
        }
        uint64_t const header = 32 + 8 * count;
        auto const all = file.ReadPrefix(static_cast<size_t>(header), error);
        if (!all) return std::nullopt;
        if (all->size() < header) {
            error = "the header is cut off (truncated file?)";
            return std::nullopt;
        }
        uint64_t offset = Padded(header);
        for (uint64_t i = 0; i < count; i++) {
            uint64_t const size = detail::GetU64(all->data() + 32 + 8 * i);
            h.sizes.push_back(size);
            h.offsets.push_back(offset);
            offset += Padded(size);
        }
        auto const total = file.KnownSize();
        if (!total || *total != offset) {
            error = "its arrays add up to " + std::to_string(offset) + " bytes, the file holds " +
                    (total ? std::to_string(*total) : std::string("an unknown number"));
            return std::nullopt;
        }
        return h;
    }

    // Reads the arrays of `header` into `arrays` (one per array, sized as the header says) with up to `threads` threads;
    // "" or what failed.
    inline std::string ReadArrays(db::DbFile const& file, Header const& header, std::vector<Array> const& arrays, int threads) {
        if (arrays.size() != header.sizes.size()) return "the table holds " + std::to_string(header.sizes.size()) + " arrays, not " +
                                                         std::to_string(arrays.size());
        std::vector<db::Section> sections;
        for (size_t i = 0; i < arrays.size(); i++) {
            if (arrays[i].size != header.sizes[i]) return "array " + std::to_string(i) + " was sized wrongly";
            if (arrays[i].size > 0) sections.push_back({ header.offsets[i], arrays[i].size, arrays[i].data });
        }
        return file.ReadSections(std::move(sections), threads);
    }
}
