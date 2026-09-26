#pragma once

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include "xxhash64.h"

namespace protal {
    // Identifies the reference an index was built against: the xxhash64 of reference.map (every gene's
    // id and byte range in reference.fna) and the size of reference.fna. An index queried with other
    // reference files would align reads against the wrong sequences, so the two are compared at load.
    struct ReferenceFingerprint {
        uint64_t map_hash = 0;
        uint64_t fna_size = 0;

        bool operator==(ReferenceFingerprint const&) const = default;

        static ReferenceFingerprint Of(std::string const& map_path, std::string const& fna_path) {
            ReferenceFingerprint fingerprint;
            XXHash64 hash(0);
            std::ifstream is(map_path, std::ios::binary);
            std::vector<char> buffer(1 << 20);
            while (is) {
                is.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
                if (is.gcount() > 0) hash.add(buffer.data(), static_cast<uint64_t>(is.gcount()));
            }
            fingerprint.map_hash = hash.hash();
            std::error_code ec;
            fingerprint.fna_size = std::filesystem::file_size(fna_path, ec);
            return fingerprint;
        }
    };
}
