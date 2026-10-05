// ReadTypeDetection.h - the kind of reads in a single read file, from its first reads, for a sample whose kind
// is not given (--read_type, a map's READ_TYPE): Options::ResolveReadTypes takes it and says so.
//
// - Short reads (at most a quarter of the first reads longer than the longest short read, so that a stray long
//   read leaves short ones to the reader, which stops at it) are single-end reads.
// - Long reads are PacBio or ONT reads by their names where the instruments' software names them its own way:
//   MinKNOW and dorado give each read a UUID (MinKNOW adds runid=...), PacBio names a read by its movie and ZMW
//   (m64011_190830_220126/123/ccs).
// - Else by their base qualities: PacBio HiFi reads are Q30 and better, ONT reads around Q20 (Q12-25), so
//   kPacBioMinQuality (Q25) tells them apart. A read's quality is that of its bases' mean error probability, as
//   the instruments' software reports it (the mean of Phred values would let a few Q93 bases outweigh many bad
//   ones); the reads' median decides.
// - Long reads without qualities (FASTA, or Q0 at every base, which means unknown: pbsim3's PacBio reads, some
//   converters) are PacBio reads unless their names say ONT.
#pragma once

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>
#include "ReadType.h"
#include "SequenceUtils/SeqReader.h"
#include "ThreadedGzStream.h"

namespace protal {
    // The reads looked at, and what they show.
    struct ReadSample {
        size_t reads = 0;
        size_t longest = 0;
        size_t median_length = 0;
        std::vector<size_t> lengths;           // of the reads, in file order
        std::optional<double> median_quality;  // Phred; none without qualities (FASTA, Q0 throughout)
        size_t ont_names = 0;                  // reads named as MinKNOW or dorado names them
        size_t pacbio_names = 0;               // reads named as PacBio names them
    };

    struct ReadTypeGuess {
        ReadType type = ReadType::Single;
        std::string evidence;  // for the message: what the reads showed
    };

    inline constexpr double kPacBioMinQuality = 25.0;

    namespace read_type_detection {
        inline bool IsHex(char c) { return std::isxdigit(static_cast<unsigned char>(c)) != 0; }

        // A UUID, 8-4-4-4-12 hex digits: how MinKNOW and dorado name ONT reads.
        inline bool IsUuid(std::string const& id) {
            if (id.size() != 36) return false;
            for (size_t i = 0; i < id.size(); i++) {
                bool const dash = i == 8 || i == 13 || i == 18 || i == 23;
                if (dash ? id[i] != '-' : !IsHex(id[i])) return false;
            }
            return true;
        }

        inline bool IsOntName(FastxRecord const& record) {
            return IsUuid(record.id) || record.header.find(" runid=") != std::string::npos ||
                   record.header.find("\trunid=") != std::string::npos;
        }

        // A PacBio movie and ZMW: m, digits, an underscore in the movie name, then /digits, and the end or a slash
        // (m64011_190830_220126/123/ccs, m84011_220902_175841_s1/4567/ccs/fwd).
        inline bool IsPacBioName(std::string const& id) {
            if (id.size() < 5 || id[0] != 'm' || !std::isdigit(static_cast<unsigned char>(id[1]))) return false;
            size_t const slash = id.find('/');
            if (slash == std::string::npos || id.find('_') > slash) return false;
            size_t end = slash + 1;
            while (end < id.size() && std::isdigit(static_cast<unsigned char>(id[end]))) end++;
            return end > slash + 1 && (end == id.size() || id[end] == '/');
        }

        // A read's quality: -10 log10 of its bases' mean error probability (Phred+33).
        inline double ReadQuality(std::string const& quality) {
            static std::array<double, 94> const error = [] {
                std::array<double, 94> e{};
                for (size_t q = 0; q < e.size(); q++) e[q] = std::pow(10.0, -static_cast<double>(q) / 10.0);
                return e;
            }();
            if (quality.empty()) return 0;
            double sum = 0;
            for (char c : quality) sum += error[static_cast<size_t>(std::clamp(static_cast<int>(c) - 33, 0, 93))];
            return std::max(0.0, -10.0 * std::log10(std::max(sum / static_cast<double>(quality.size()), 1e-10)));
        }

        template<typename T>
        T Median(std::vector<T> values) {
            std::nth_element(values.begin(), values.begin() + static_cast<long>(values.size() / 2), values.end());
            return values[values.size() / 2];
        }
    }

    // The first max_reads reads of a FASTQ or FASTA file (plain, gzip or zstd), or fewer if they hold max_bases
    // bases. None read (reads 0) if the file is not a regular file (a pipe would lose them) or holds no reads.
    inline ReadSample SampleReads(std::string const& path, size_t max_reads = 200, size_t max_bases = 5'000'000) {
        using namespace read_type_detection;
        ReadSample sample;
        std::error_code ec;
        if (!std::filesystem::is_regular_file(path, ec)) return sample;
        ThreadedGzIstream is(path.c_str());
        SeqReaderSE reader(is);
        FastxRecord record;
        auto& lengths = sample.lengths;
        std::vector<double> qualities;
        size_t bases = 0;
        bool any_quality = false;  // a base above Q0: all Q0 ('!') is no quality, as pbsim3 and some converters write
        while (sample.reads < max_reads && bases < max_bases && reader(record)) {
            sample.reads++;
            bases += record.sequence.size();
            lengths.push_back(record.sequence.size());
            sample.longest = std::max(sample.longest, record.sequence.size());
            if (record.quality.size() == record.sequence.size() && !record.quality.empty()) {
                qualities.push_back(ReadQuality(record.quality));
                any_quality = any_quality || std::any_of(record.quality.begin(), record.quality.end(), [](char c) { return c > '!'; });
            }
            sample.ont_names += IsOntName(record);
            sample.pacbio_names += IsPacBioName(record.id);
        }
        if (!lengths.empty()) sample.median_length = Median(lengths);
        if (any_quality && qualities.size() * 2 > sample.reads) sample.median_quality = Median(qualities);
        return sample;
    }

    // The kind of reads `sample` shows, short reads being those of up to short_read_length bases; none if no reads
    // were read.
    inline std::optional<ReadTypeGuess> GuessReadType(ReadSample const& sample, size_t short_read_length) {
        if (sample.reads == 0) return std::nullopt;
        auto kb = [](size_t bases) {
            char text[32];
            std::snprintf(text, sizeof text, "%.1f kb", static_cast<double>(bases) / 1000.0);
            return std::string(text);
        };
        std::string const lengths = "the first " + std::to_string(sample.reads) + " reads up to " + kb(sample.longest) +
                                    " long, median " + kb(sample.median_length);
        size_t const long_reads = static_cast<size_t>(std::count_if(sample.lengths.begin(), sample.lengths.end(),
                                                                    [short_read_length](size_t n) { return n > short_read_length; }));
        if (4 * long_reads <= sample.reads) return ReadTypeGuess{ ReadType::Single, lengths };
        if (sample.ont_names * 2 > sample.reads) return ReadTypeGuess{ ReadType::ONT, lengths + ", named as MinKNOW and dorado name reads" };
        if (sample.pacbio_names * 2 > sample.reads) return ReadTypeGuess{ ReadType::PacBio, lengths + ", named as PacBio names reads" };
        if (sample.median_quality) {
            char text[96];
            std::snprintf(text, sizeof text, ", median read quality Q%.1f (PacBio from Q%.0f, ONT below)", *sample.median_quality,
                          kPacBioMinQuality);
            return ReadTypeGuess{ *sample.median_quality >= kPacBioMinQuality ? ReadType::PacBio : ReadType::ONT, lengths + text };
        }
        return ReadTypeGuess{ ReadType::PacBio, lengths + ", without base qualities (FASTA, or Q0 throughout) to tell PacBio from ONT reads" };
    }
}
