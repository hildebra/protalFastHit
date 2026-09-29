// ReadType.h - the kinds of reads protal profiles: short reads, paired-end or single-end, and PacBio
// long reads. Each kind is aligned in its own way and has its own presence model (the model's
// features are counted per read), and a SAM names the kind of its reads in a header comment.
#pragma once

#include <array>
#include <optional>
#include <string>

namespace protal {
    enum class ReadType { Paired, Single, PacBio };
    inline constexpr std::array<ReadType, 3> kReadTypes = { ReadType::Paired, ReadType::Single, ReadType::PacBio };

    // The name in messages and in the SAM header comment.
    inline std::string ReadTypeName(ReadType type) {
        switch (type) {
            case ReadType::Paired: return "paired-end";
            case ReadType::Single: return "single-end";
            case ReadType::PacBio: return "PacBio";
        }
        return "";
    }

    inline std::optional<ReadType> ReadTypeFromName(std::string const& name) {
        for (auto type : kReadTypes) {
            if (ReadTypeName(type) == name) return type;
        }
        return std::nullopt;
    }

    // The SAM header line naming the reads of the SAM: kSamReadsComment + ReadTypeName.
    inline const std::string kSamReadsComment = "@CO\tprotal reads: ";
}
