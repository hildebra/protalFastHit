// ReadType.h - the kinds of reads protal profiles: short reads, paired-end or single-end, and PacBio
// and ONT long reads. Each kind is aligned in its own way and has its own presence model (the
// model's features are counted per read), which a database holds as model_<...>.xml (--build packs
// them, --add_model stores one), and a SAM names the kind of its reads in a header comment.
#pragma once

#include <array>
#include <optional>
#include <string>
#include <vector>

namespace protal {
    enum class ReadType { Paired, Single, PacBio, ONT };

    struct ReadTypeInfo {
        ReadType type;
        std::string token;         // --read_type, a map's READ_TYPE, the SAM header comment
        std::string name;          // in messages
        std::string model_file;    // the database's model of these reads
        std::string model_option;  // the option naming another model for them
        int fasta_quality;         // the Phred quality of every base of a read without qualities (FASTA)
        // Defaults of these reads unless -a / --snp_min_af are given (none: the options' defaults).
        // ONT errors are mostly indels, which the alignment identity counts twice (see
        // AlignmentInfo::GetProxyANI), and they put a low-frequency second allele at many positions.
        std::optional<double> max_score_ani;
        std::optional<double> snp_min_af;
    };

    inline constexpr size_t kReadTypeCount = 4;
    inline const std::array<ReadTypeInfo, kReadTypeCount> kReadTypes = {{
            { ReadType::Paired, "pe", "paired-end", "model_pe.xml", "--model", 30, std::nullopt, std::nullopt },
            { ReadType::Single, "se", "single-end", "model_se.xml", "--model_se", 30, std::nullopt, std::nullopt },
            { ReadType::PacBio, "pb", "PacBio", "model_PB.xml", "--model_pb", 30, std::nullopt, std::nullopt },
            { ReadType::ONT, "ont", "ONT", "model_ONT.xml", "--model_ont", 18, 0.85, 0.2 } }};

    // Long reads, aligned per gene (LongReads.h).
    inline bool IsLongReadType(ReadType type) {
        return type == ReadType::PacBio || type == ReadType::ONT;
    }

    // Databases from before read types hold one model, for paired-end reads.
    inline const std::vector<std::string> kLegacyModelFiles = { "model.xml", "random_forest.xml" };

    inline ReadTypeInfo const& Info(ReadType type) {
        return kReadTypes[static_cast<size_t>(type)];
    }

    inline std::string const& ReadTypeName(ReadType type) {
        return Info(type).name;
    }

    // The quality character (Phred+33) given to every base of a read of `type` without qualities.
    inline char FastaQualityChar(ReadType type) {
        return static_cast<char>(33 + Info(type).fasta_quality);
    }

    inline std::optional<ReadType> ReadTypeFromToken(std::string const& token) {
        for (auto const& info : kReadTypes) {
            if (info.token == token) return info.type;
        }
        return std::nullopt;
    }

    // "pe, se, pb, ont", for messages.
    inline std::string ReadTypeTokens() {
        std::string tokens;
        for (auto const& info : kReadTypes) tokens += (tokens.empty() ? "" : ", ") + info.token;
        return tokens;
    }

    // The model files a database may hold for reads of `type`, in order of precedence: the type's
    // (Info(type).model_file), and for paired-end reads those of databases from before read types.
    inline std::vector<std::string> ModelCandidates(ReadType type) {
        std::vector<std::string> names = { Info(type).model_file };
        if (type == ReadType::Paired) names.insert(names.end(), kLegacyModelFiles.begin(), kLegacyModelFiles.end());
        return names;
    }

    // Every model file a database may hold.
    inline std::vector<std::string> AllModelFiles() {
        std::vector<std::string> names;
        for (auto const& info : kReadTypes) names.push_back(info.model_file);
        names.insert(names.end(), kLegacyModelFiles.begin(), kLegacyModelFiles.end());
        return names;
    }

    // The SAM header line naming the reads of the SAM: kSamReadTypeComment + the type's token.
    inline const std::string kSamReadTypeComment = "@CO\tprotal read type: ";
}
