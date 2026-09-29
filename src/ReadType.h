// ReadType.h - the kinds of reads protal profiles: short reads, paired-end or single-end, and PacBio
// long reads. Each kind is aligned in its own way and has its own presence model (the model's
// features are counted per read), and a SAM names the kind of its reads in a header comment.
// Tokens and model files are those of the database build (branch zstd-compression's
// db::kReadTypeModels).
#pragma once

#include <array>
#include <optional>
#include <string>
#include <vector>

namespace protal {
    enum class ReadType { Paired, Single, PacBio };

    struct ReadTypeInfo {
        ReadType type;
        std::string token;         // --read_type, a map's READ_TYPE, the SAM header comment
        std::string name;          // in messages
        std::string model_file;    // the database's model of these reads
        std::string model_option;  // the option naming another model for them
    };

    inline constexpr size_t kReadTypeCount = 3;
    inline const std::array<ReadTypeInfo, kReadTypeCount> kReadTypes = {{
            { ReadType::Paired, "pe", "paired-end", "model_pe.xml", "--model" },
            { ReadType::Single, "se", "single-end", "model_se.xml", "--model_se" },
            { ReadType::PacBio, "pb", "PacBio", "model_PB.xml", "--model_pb" } }};

    // Databases from before read types hold one model, for paired-end reads.
    inline const std::vector<std::string> kLegacyModelFiles = { "model.xml", "random_forest.xml" };

    inline ReadTypeInfo const& Info(ReadType type) {
        return kReadTypes[static_cast<size_t>(type)];
    }

    inline std::string const& ReadTypeName(ReadType type) {
        return Info(type).name;
    }

    inline std::optional<ReadType> ReadTypeFromToken(std::string const& token) {
        for (auto const& info : kReadTypes) {
            if (info.token == token) return info.type;
        }
        return std::nullopt;
    }

    // The SAM header line naming the reads of the SAM: kSamReadTypeComment + the type's token.
    inline const std::string kSamReadTypeComment = "@CO\tprotal read type: ";
}
