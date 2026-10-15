// ReadTypeDetection.h - the kind of reads in a single read file, from its first reads, for a sample whose kind
// is not given (--read_type, a map's READ_TYPE): Options::ResolveReadTypes takes it and says so.
//
// - Short reads (at most a quarter of the first reads longer than the longest short read, so that a stray long
//   read leaves short ones to the reader, which stops at it) are single-end reads.
// - Long reads are PacBio or ONT reads by their names where the instruments' software names them its own way:
//   MinKNOW and dorado give each read a UUID (MinKNOW adds runid=...; a dorado duplex read is "uuid;uuid"), PacBio
//   names a read by its movie and ZMW (m64011_190830_220126/123/ccs). Anywhere in the header: SRA and ENA put their
//   own id first (SRR123.1) and keep the original name after it, if at all.
// - Else by their base qualities: PacBio HiFi reads are Q30 and better, ONT reads around Q20 (Q12-25), so
//   kPacBioMinQuality (Q25) tells them apart. A read's quality is that of its bases' mean error probability, as
//   the instruments' software reports it (the mean of Phred values would let a few Q93 bases outweigh many bad
//   ones); the reads' median decides. This can be wrong (ONT reads of R10.4.1 sup and duplex basecalling reach
//   Q25, SRA Lite gives every base Q30), so Options::ResolveReadTypes warns.
// - Long reads without qualities (FASTA, or Q0 at every base, which means unknown: pbsim3's PacBio reads, some
//   converters) are PacBio reads unless their names say ONT, with a warning too.
#pragma once

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
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
        std::optional<double> median_quality;  // Phred, of the reads with qualities where they are most; none
                                               // without (FASTA, Q0 throughout)
        size_t ont_names = 0;                  // reads named as MinKNOW or dorado names them
        size_t pacbio_names = 0;               // reads named as PacBio names them
    };

    // What decided the kind of reads: their lengths (single-end reads), their names, their median base quality, or the
    // lack of qualities (PacBio reads). The last two can be wrong, so DescribeReadType makes them a warning.
    enum class ReadTypeBasis { Length, Names, Quality, NoQualities };

    struct ReadTypeGuess {
        ReadType type = ReadType::Single;
        std::string evidence;  // for the message: what the reads showed
        ReadTypeBasis basis = ReadTypeBasis::Length;
        // Single-end reads that look like long ones all the same: a median length beyond Illumina's
        // (kLongLikeMedianLength) and a median quality below PacBio's, as short ONT reads have.
        bool long_like = false;
    };

    inline constexpr double kPacBioMinQuality = 25.0;
    // Illumina's reads are at most 300 bp (MiSeq), Ultima's mostly shorter: single-end reads of a longer median length
    // and with ONT's qualities look like long reads (that are mostly up to 1 kb, so their lengths made them short ones).
    inline constexpr size_t kLongLikeMedianLength = 400;

    namespace read_type_detection {
        inline bool IsHex(char c) { return std::isxdigit(static_cast<unsigned char>(c)) != 0; }

        // A UUID, 8-4-4-4-12 hex digits: how MinKNOW and dorado name ONT reads.
        inline bool IsUuid(std::string_view id) {
            if (id.size() != 36) return false;
            for (size_t i = 0; i < id.size(); i++) {
                bool const dash = i == 8 || i == 13 || i == 18 || i == 23;
                if (dash ? id[i] != '-' : !IsHex(id[i])) return false;
            }
            return true;
        }

        // Whether `text` holds a UUID (IsUuid) between characters that are no letters or digits: a MinKNOW or dorado
        // read id, also after an archive's own (SRA, ENA: "SRR123.1 <uuid> ..."), in a dorado duplex read's
        // "<uuid>;<uuid>" and before a suffix ("<uuid>_Basecall_1D_template").
        inline bool HoldsUuid(std::string_view text) {
            auto word = [](char c) { return std::isalnum(static_cast<unsigned char>(c)) != 0; };
            for (size_t i = 0; i + 36 <= text.size(); i++) {
                if ((i > 0 && word(text[i - 1])) || (i + 36 < text.size() && word(text[i + 36]))) continue;
                if (IsUuid(text.substr(i, 36))) return true;
            }
            return false;
        }

        // A read's header line without its '@' or '>', or its id where it has none (a record not read from a file).
        inline std::string_view HeaderOf(FastxRecord const& record) {
            std::string_view header = record.header.empty() ? std::string_view(record.id) : std::string_view(record.header);
            if (!record.header.empty() && (header.front() == '@' || header.front() == '>')) header.remove_prefix(1);
            return header;
        }

        // A UUID or MinKNOW's runid= anywhere in the header.
        inline bool IsOntName(FastxRecord const& record) {
            std::string_view const header = HeaderOf(record);
            return HoldsUuid(header) || header.find("runid=") != std::string_view::npos;
        }

        // A PacBio movie and ZMW: m, digits, an underscore in the movie name, then /digits, and the end or a slash
        // (m64011_190830_220126/123/ccs, m84011_220902_175841_s1/4567/ccs/fwd).
        inline bool IsPacBioName(std::string_view id) {
            if (id.size() < 5 || id[0] != 'm' || !std::isdigit(static_cast<unsigned char>(id[1]))) return false;
            size_t const slash = id.find('/');
            if (slash == std::string_view::npos || id.find('_') > slash) return false;
            size_t end = slash + 1;
            while (end < id.size() && std::isdigit(static_cast<unsigned char>(id[end]))) end++;
            return end > slash + 1 && (end == id.size() || id[end] == '/');
        }

        // A PacBio name (IsPacBioName) as any word of the header: the id, or one after an archive's id
        // ("SRR123.1 m64011_190830_220126/123/ccs length=...").
        inline bool HoldsPacBioName(FastxRecord const& record) {
            std::string_view const header = HeaderOf(record);
            for (size_t start = 0; start < header.size();) {
                size_t const end = std::min(header.find_first_of(" \t", start), header.size());
                if (end > start && IsPacBioName(header.substr(start, end - start))) return true;
                start = end + 1;
            }
            return false;
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
        while (sample.reads < max_reads && bases < max_bases && reader(record)) {
            sample.reads++;
            bases += record.sequence.size();
            lengths.push_back(record.sequence.size());
            sample.longest = std::max(sample.longest, record.sequence.size());
            // Q0 at every base is no quality (IsPlaceholderQuality), as the alignment's reader takes it.
            if (record.quality.size() == record.sequence.size() && !record.quality.empty() &&
                !IsPlaceholderQuality(record.quality)) {
                qualities.push_back(ReadQuality(record.quality));
            }
            sample.ont_names += IsOntName(record);
            sample.pacbio_names += HoldsPacBioName(record);
        }
        if (!lengths.empty()) sample.median_length = Median(lengths);
        if (qualities.size() * 2 > sample.reads) sample.median_quality = Median(qualities);
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
        if (4 * long_reads <= sample.reads) {
            ReadTypeGuess guess{ ReadType::Single, lengths, ReadTypeBasis::Length };
            guess.long_like = sample.median_length > kLongLikeMedianLength && sample.median_quality &&
                              *sample.median_quality < kPacBioMinQuality;
            if (guess.long_like) {
                char text[96];
                std::snprintf(text, sizeof text, ", at most a quarter over %s, median read quality Q%.1f",
                              kb(short_read_length).c_str(), *sample.median_quality);
                guess.evidence += text;
            }
            return guess;
        }
        if (sample.ont_names * 2 > sample.reads) {
            return ReadTypeGuess{ ReadType::ONT, lengths + ", named as MinKNOW and dorado name reads", ReadTypeBasis::Names };
        }
        if (sample.pacbio_names * 2 > sample.reads) {
            return ReadTypeGuess{ ReadType::PacBio, lengths + ", named as PacBio names reads", ReadTypeBasis::Names };
        }
        if (sample.median_quality) {
            char text[96];
            std::snprintf(text, sizeof text, ", median read quality Q%.1f (PacBio from Q%.0f, ONT below)", *sample.median_quality,
                          kPacBioMinQuality);
            return ReadTypeGuess{ *sample.median_quality >= kPacBioMinQuality ? ReadType::PacBio : ReadType::ONT, lengths + text,
                                  ReadTypeBasis::Quality };
        }
        return ReadTypeGuess{ ReadType::PacBio, lengths + ", without base qualities (FASTA, or Q0 throughout) to tell PacBio from ONT reads",
                              ReadTypeBasis::NoQualities };
    }

    struct ReadTypeMessage {
        std::string text;
        bool warning = false;  // else a note
    };

    // What Options::ResolveReadTypes says of `guess` for the reads of sample `sample` in `file`. A note where the reads'
    // names decided, or where single-end reads look like long ones (long_like). A warning where their base quality, or
    // the lack of it, decided a long-read type alone: ONT reads downloaded from SRA or ENA have lost their names, and
    // ONT reads reach PacBio's qualities (R10.4.1 sup and duplex basecalling; SRA Lite gives every base Q30). None for
    // other single-end reads.
    inline std::optional<ReadTypeMessage> DescribeReadType(ReadTypeGuess const& guess, std::string const& sample,
                                                           std::string const& file) {
        std::string const found = "Sample " + sample + ": " + ReadTypeName(guess.type) + " reads (" + guess.evidence + " in " +
                                  file + "), aligned and profiled as such";
        std::string const names = "their names show no UUID (MinKNOW, dorado) and no movie and ZMW (PacBio)";
        switch (guess.basis) {
            case ReadTypeBasis::Length:
                if (!guess.long_like) return std::nullopt;
                return ReadTypeMessage{ found + ", though their lengths and qualities are those of long reads rather than "
                                        "Illumina's; --read_type ont or pb (or a map's READ_TYPE) for ONT or PacBio reads" };
            case ReadTypeBasis::Names:
                return ReadTypeMessage{ found + "; --read_type (or a map's READ_TYPE) sets the kind of reads" };
            case ReadTypeBasis::Quality:
                return ReadTypeMessage{ found + " by their base quality alone: " + names + ", as after a download from SRA or ENA" +
                                        (guess.type == ReadType::PacBio ? ", and ONT reads reach PacBio's Q25 with recent "
                                                                          "basecalling (R10.4.1 sup, duplex) or as SRA Lite "
                                                                          "(Q30 throughout)" : "") +
                                        "; --read_type ont or pb (or a map's READ_TYPE) sets the kind of reads", true };
            case ReadTypeBasis::NoQualities:
                return ReadTypeMessage{ found + " for want of base qualities: " + names + ", and long reads without "
                                        "qualities are taken for PacBio's; --read_type ont or pb (or a map's READ_TYPE) sets "
                                        "the kind of reads", true };
        }
        return std::nullopt;
    }
}
