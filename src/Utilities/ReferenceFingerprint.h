#pragma once

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include "xxhash64.h"
#include "Zstd.h"
#include "Database.h"

namespace protal {
    // Identifies the reference an index was built against: the xxhash64 of reference.map (every gene's
    // id and byte range in reference.fna) and the size of reference.fna. An index queried with other
    // reference files would align reads against the wrong sequences, so the two are compared at load.
    // The size is that of the uncompressed content, so the fingerprint stays the same when
    // reference.fna is compressed or packed into a single-file database.
    struct ReferenceFingerprint {
        uint64_t map_hash = 0;
        uint64_t fna_size = 0;

        bool operator==(ReferenceFingerprint const&) const = default;

        static ReferenceFingerprint Of(db::DbFile const& map, db::DbFile const& fna) {
            ReferenceFingerprint fingerprint;
            XXHash64 hash(0);
            auto input = map.Open();
            std::istream& is = input->Stream();
            std::vector<char> buffer(1 << 20);
            try {
                while (is) {
                    is.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
                    if (is.gcount() > 0) hash.add(buffer.data(), static_cast<uint64_t>(is.gcount()));
                }
            } catch (std::exception const&) {
                is.setstate(std::ios_base::badbit);
            }
            // An unreadable map cannot match any index.
            fingerprint.map_hash = is.bad() || !map.Exists() ? 0 : hash.hash();
            fingerprint.fna_size = fna.Size().value_or(static_cast<uint64_t>(-1));
            return fingerprint;
        }

        // fna_path may be compressed or have a .zst sibling.
        static ReferenceFingerprint Of(std::string const& map_path, std::string const& fna_path) {
            return Of(db::DbFile::OnDisk(map_path), db::DbFile::OnDisk(zstd::Resolve(fna_path)));
        }
    };
}
