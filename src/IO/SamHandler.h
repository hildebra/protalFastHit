//
// Created by fritsche on 26/08/22.
//

#ifndef PROTAL_SAMHANDLER_H
#define PROTAL_SAMHANDLER_H


#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include "LineSplitter.h"
#include "ReadType.h"
#include <iostream>
#include <omp.h>

namespace protal {

    using QNAME_t = std::string;
    using FLAG_t = uint16_t;
    using RNAME_t = std::string;
    using POS_t = uint32_t;
    using MAPQ_t = uint8_t;
    using CIGAR_t = std::string;
    using TLEN_t = int64_t;
    using SEQ_t = std::string;
    using QUAL_t = std::string;

    class Flag {
    public:
        static void SetPairedEnd(FLAG_t &flag, bool ispaired, bool bothalign, bool is_read1=false, bool is_read2=false) {
            flag |= ispaired;
//            flag |= (bothalign << 1);
            flag |= (is_read1 << 6);
            flag |= (is_read2 << 7);
        }

        static void SetPairBothAlign(FLAG_t &flag, bool both_align) {
            flag |= (both_align << 1);
        }

        static void SetMateUnmapped(FLAG_t &flag, bool unmapped) {
            flag |= (unmapped << 3);
        }

        static void SetReadUnmapped(FLAG_t &flag, bool unmapped) {
            flag |= (unmapped << 2);
        }

        static void SetReadReverseComplement(FLAG_t &flag, bool rc) {
            flag |= (rc << 4);
        }

        static void SetMateReverseComplement(FLAG_t &flag, bool rc) {
            flag |= (rc << 5);
        }

        static void SetNotPrimaryAlignment(FLAG_t &flag, bool not_primary) {
            flag |= (not_primary << 8);
        }

        static void SetAlignmentFailsQuality(FLAG_t &flag, bool fails) {
            flag |= (fails << 9);
        }

        static void SetDuplicate(FLAG_t &flag, bool duplicate) {
            flag |= (duplicate << 10);
        }

        static void SetSupplementaryAlignment(FLAG_t &flag, bool is_supplementary) {
            flag |= (is_supplementary << 11);
        }


        static bool IsPaired(FLAG_t flag) {
            return flag & 1;
        }

        static bool IsPairBothAlign(FLAG_t flag) {
            return flag & (1 << 1);
        }

        static bool IsUnmapped(FLAG_t flag) {
            return flag & (1 << 2);
        }

        static bool IsMateUnmapped(FLAG_t flag) {
            return flag & (1 << 3);
        }

        static bool IsReverseComplement(FLAG_t flag) {
            return flag & (1 << 4);
        }

        static bool IsMateReverseComplement(FLAG_t flag) {
            return flag & (1 << 5);
        }

        static bool IsRead1(FLAG_t flag) {
            return flag & (1 << 6);
        }

        static bool IsRead2(FLAG_t flag) {
            return flag & (1 << 7);
        }

        static bool IsNotPrimaryAlignment(FLAG_t flag) {
            return flag & (1 << 8);
        }

        static bool IsAlignmentFailsQuality(FLAG_t flag) {
            return flag & (1 << 9);
        }

        static bool IsDuplicate(FLAG_t &flag) {
            return flag & (1 << 10);
        }

        static bool IsSupplementaryAlignment(FLAG_t &flag) {
            return flag & (1 << 11);
        }
    };

    struct ReducedSam {
        size_t m_qid;
        size_t m_rid;
        FLAG_t m_flag;
        POS_t m_pos;
        MAPQ_t m_mapq;
        RNAME_t m_rnext;
        POS_t m_pnext;
        TLEN_t m_tlen;
        QUAL_t m_qual;
        SEQ_t m_seq;
        CIGAR_t m_cigar;

        [[nodiscard]] std::string ToString() const {
            return  std::to_string(m_qid) + '\t'
                    + std::to_string(m_flag) + '\t'
                    + std::to_string(m_rid) + '\t'
                    + std::to_string(m_pos) + '\t'
                    + std::to_string(m_mapq) + '\t'
                    + m_cigar + '\t'
                    + m_rnext + '\t'
                    + std::to_string(m_pnext) + '\t'
                    + std::to_string(m_tlen) + '\t'
                    + m_seq + '\t'
                    + m_qual;
        }

        // Strand of THIS record (0x10), for read1 and read2 alike; 0x20 is the mate's strand.
        bool IsReversed() const {
            return Flag::IsReverseComplement(m_flag);
        }

    };

    struct SamEntry {
        QNAME_t m_qname;
        FLAG_t m_flag;
        RNAME_t m_rname;
        POS_t m_pos;
        MAPQ_t m_mapq;
        RNAME_t m_rnext;
        POS_t m_pnext;
        TLEN_t m_tlen;
        uint16_t m_uniques;
        uint16_t m_uniques_two;
        QUAL_t m_qual;
        SEQ_t m_seq;
        CIGAR_t m_cigar;
        // The hard clips at the CIGAR's start and end, which NormalizeCigar drops from m_cigar: where on its read a long
        // read's record lies (UnusableRecord sets them).
        uint32_t m_hard_clip_start = 0;
        uint32_t m_hard_clip_end = 0;
        // ZA tag of a primary or supplementary record: the read's alternative alignments to other taxa,
        // "<taxid>:<edits more than this one>" comma-separated (AlternativesTag), or "*" for none. Empty:
        // not written (secondary records, SAM files of older protal versions).
        std::string m_alternatives;

        [[nodiscard]] std::string ToString() const {
            return  m_qname + '\t'
                    + std::to_string(m_flag) + '\t'
                    + m_rname + '\t'
                    + std::to_string(m_pos) + '\t'
                    + std::to_string(m_mapq) + '\t'
                    + m_cigar + '\t'
                    + m_rnext + '\t'
                    + std::to_string(m_pnext) + '\t'
                    + std::to_string(m_tlen) + '\t'
                    + m_seq + '\t'
                    + m_qual + '\t'
                    + "ZU:i:" + std::to_string(m_uniques) + '\t'
                    + "ZT:i:" + std::to_string(m_uniques_two)
                    + (m_alternatives.empty() ? std::string() : "\tZA:Z:" + m_alternatives);
        }

        // Strand of THIS record (0x10), for read1 and read2 alike; 0x20 is the mate's strand.
        bool IsReversed() const {
            return Flag::IsReverseComplement(m_flag);
        }
    };


    // A SAM line that cannot be parsed: too few fields or a non-numeric FLAG, POS, MAPQ, PNEXT or TLEN.
    struct SamFormatError : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    namespace sam_detail {
        inline bool ParseUnsigned(std::string_view s, uint64_t& value) {
            auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), value);
            return !s.empty() && ec == std::errc() && ptr == s.data() + s.size();
        }

        inline bool ParseSigned(std::string_view s, int64_t& value) {
            auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), value);
            return !s.empty() && ec == std::errc() && ptr == s.data() + s.size();
        }

        // Value of an integer tag ("ZU:i:12") among the optional fields, or nullopt if it is absent
        // or not a non-negative integer.
        inline std::optional<uint64_t> IntTag(std::vector<std::string_view> const& tokens, std::string_view name) {
            for (size_t i = 11; i < tokens.size(); i++) {
                auto const t = tokens[i];
                if (t.size() > 5 && t.substr(0, 2) == name && t.substr(2, 3) == ":i:") {
                    uint64_t value;
                    if (ParseUnsigned(t.substr(5), value)) return value;
                    return std::nullopt;
                }
            }
            return std::nullopt;
        }

        // Value of a string tag ("ZA:Z:12:0,40:1") among the optional fields, or nullopt if it is absent.
        inline std::optional<std::string_view> StringTag(std::vector<std::string_view> const& tokens, std::string_view name) {
            for (size_t i = 11; i < tokens.size(); i++) {
                auto const t = tokens[i];
                if (t.size() >= 5 && t.substr(0, 2) == name && t.substr(2, 3) == ":Z:") return t.substr(5);
            }
            return std::nullopt;
        }

        // The tab-separated fields of a line, as LineSplitter::Split gives them (empty ones kept, none for an
        // empty line), as views into it.
        inline void SplitFields(std::string_view line, std::vector<std::string_view>& fields) {
            fields.clear();
            if (line.empty()) return;
            size_t start = 0;
            for (size_t tab; (tab = line.find('\t', start)) != std::string_view::npos; start = tab + 1) {
                fields.push_back(line.substr(start, tab - start));
            }
            fields.push_back(line.substr(start));
        }

        // protal references are genes named "<taxid>_<gene id>".
        inline bool IsProtalGene(std::string const& rname) {
            auto underscore = rname.find('_');
            if (underscore == std::string::npos) return false;
            auto end = rname.find('_', underscore + 1);
            if (end == std::string::npos) end = rname.size();
            auto digits = [&](size_t from, size_t to) {
                return to > from && std::all_of(rname.begin() + from, rname.begin() + to,
                                                [](char c) { return std::isdigit(static_cast<unsigned char>(c)); });
            };
            return digits(0, underscore) && digits(underscore + 1, end);
        }
    }

    // Whether NormalizeCigar leaves a CIGAR as it is: ops M, X, I, D and S only, no two in a row the same,
    // each count without leading zeros (of up to 9 digits here; longer ones are left to NormalizeCigar), and the
    // bases of SEQ.
    inline bool IsNormalCigar(std::string const& cigar, size_t seq_length) {
        if (cigar.empty()) return false;
        size_t query = 0;
        char last = 0;
        for (size_t i = 0; i < cigar.size();) {
            size_t j = i;
            size_t count = 0;
            while (j < cigar.size() && cigar[j] >= '0' && cigar[j] <= '9') count = count * 10 + static_cast<size_t>(cigar[j++] - '0');
            if (j == i || j == cigar.size() || j - i > 9 || cigar[i] == '0') return false;
            char const op = cigar[j];
            if (!(op == 'M' || op == 'X' || op == 'I' || op == 'D' || op == 'S') || op == last) return false;
            if (op != 'D') query += count;
            last = op;
            i = j + 1;
        }
        return query == seq_length;
    }

    // Rewrites a CIGAR into the ops the profiler walks, which are those protal writes: M for exact
    // matches, X, I, D and S. '=' (sequence match) becomes M, and hard clips are dropped because they
    // consume neither SEQ nor the reference (clip_start and clip_end, if given, get those at the start and
    // at the end). Returns false for a CIGAR that cannot be walked: malformed, with N or P ops, or not
    // covering exactly the bases in SEQ.
    inline bool NormalizeCigar(std::string& cigar, size_t seq_length, uint32_t* clip_start = nullptr, uint32_t* clip_end = nullptr) {
        // Most CIGARs are left as they are (protal writes them so): only checked, not rebuilt.
        if (IsNormalCigar(cigar, seq_length)) {
            if (clip_start) *clip_start = 0;
            if (clip_end) *clip_end = 0;
            return true;
        }
        uint64_t hard_start = 0, hard_end = 0;  // hard_end: those since the last other op
        bool other = false;
        std::string out;
        size_t query = 0;
        char last_op = 0;
        uint64_t last_count = 0;
        auto flush = [&]() {
            if (last_count > 0) out += std::to_string(last_count) + last_op;
        };
        for (size_t i = 0; i < cigar.size();) {
            size_t j = i;
            while (j < cigar.size() && std::isdigit(static_cast<unsigned char>(cigar[j]))) j++;
            uint64_t count;
            if (j == cigar.size() || !sam_detail::ParseUnsigned(std::string_view(cigar).substr(i, j - i), count) || count == 0) return false;
            char op = cigar[j];
            i = j + 1;
            switch (op) {
                case '=': op = 'M'; break;
                case 'M': case 'X': case 'I': case 'D': case 'S': break;
                case 'H':
                    (other ? hard_end : hard_start) += count;
                    continue;
                default: return false;
            }
            other = true;
            hard_end = 0;
            if (op != 'D') query += count;
            if (op == last_op) {
                last_count += count;
            } else {
                flush();
                last_op = op;
                last_count = count;
            }
        }
        flush();
        if (out.empty() || query != seq_length) return false;
        cigar = std::move(out);
        if (clip_start) *clip_start = static_cast<uint32_t>(std::min<uint64_t>(hard_start, UINT32_MAX));
        if (clip_end) *clip_end = static_cast<uint32_t>(std::min<uint64_t>(hard_end, UINT32_MAX));
        return true;
    }

    // Fills sam from the fields of one SAM line; throws SamFormatError if the line cannot be parsed.
    // ZU and ZT (protal's unique k-mer counts) are looked up by name and are 0 when absent; ZA (the
    // read's alternatives in other taxa) is empty when absent.
    inline void SamFromTokens(std::vector<std::string_view> const& tokens, SamEntry &sam) {
        if (tokens.size() < 11) {
            throw SamFormatError("expected at least 11 tab-separated fields, found " + std::to_string(tokens.size()));
        }
        uint64_t flag, pos, mapq, pnext;
        int64_t tlen;
        auto text = [](std::string_view s) { return std::string(s); };
        if (!sam_detail::ParseUnsigned(tokens[1], flag) || flag > UINT16_MAX) throw SamFormatError("FLAG is not a 16-bit integer: " + text(tokens[1]));
        if (!sam_detail::ParseUnsigned(tokens[3], pos) || pos > UINT32_MAX) throw SamFormatError("POS is not a position: " + text(tokens[3]));
        if (!sam_detail::ParseUnsigned(tokens[4], mapq) || mapq > 255) throw SamFormatError("MAPQ is not between 0 and 255: " + text(tokens[4]));
        if (!sam_detail::ParseUnsigned(tokens[7], pnext) || pnext > UINT32_MAX) throw SamFormatError("PNEXT is not a position: " + text(tokens[7]));
        if (!sam_detail::ParseSigned(tokens[8], tlen)) throw SamFormatError("TLEN is not an integer: " + text(tokens[8]));

        sam.m_qname = tokens[0];
        sam.m_flag = static_cast<FLAG_t>(flag);
        sam.m_rname = tokens[2];
        sam.m_pos = static_cast<POS_t>(pos);
        sam.m_mapq = static_cast<MAPQ_t>(mapq);
        sam.m_cigar = tokens[5];
        sam.m_rnext = tokens[6];
        sam.m_pnext = static_cast<POS_t>(pnext);
        sam.m_tlen = tlen;
        sam.m_seq = tokens[9];
        sam.m_qual = tokens[10];
        sam.m_uniques = static_cast<uint16_t>(std::min<uint64_t>(sam_detail::IntTag(tokens, "ZU").value_or(0), UINT16_MAX));
        sam.m_uniques_two = static_cast<uint16_t>(std::min<uint64_t>(sam_detail::IntTag(tokens, "ZT").value_or(0), UINT16_MAX));
        sam.m_alternatives = sam_detail::StringTag(tokens, "ZA").value_or(std::string_view());
    }

    // Why the profiler cannot use a parsed record, or nullptr if it can. Normalizes the CIGAR.
    inline char const* UnusableRecord(SamEntry& sam) {
        if (Flag::IsUnmapped(sam.m_flag) || sam.m_rname == "*") return "unmapped";
        if (!sam_detail::IsProtalGene(sam.m_rname)) return "reference is not a protal gene (<taxid>_<gene id>)";
        if (sam.m_pos == 0) return "no position";
        if (sam.m_seq == "*") return "no sequence";
        if (sam.m_qual.size() != sam.m_seq.size()) return "no base qualities, or QUAL and SEQ lengths differ";
        if (!NormalizeCigar(sam.m_cigar, sam.m_seq.size(), &sam.m_hard_clip_start, &sam.m_hard_clip_end)) return "CIGAR is missing, malformed, has N/P ops or does not match SEQ";
        return nullptr;
    }

    // Reads the alignment records of a SAM stream. Header lines and blank lines are skipped, and so
    // are records the profiler cannot use (counted by reason in Skipped()). Next() returns a read1
    // record together with its mate when the mate's record follows, or a single record: a read2
    // without its read1 in sam2, any other one (an orphan read1, a single-end read) in sam1. A line
    // that cannot be parsed throws SamFormatError naming the line. The lines come from a stream or from a
    // text in memory (the profiler's chunks of a SAM); the fields are views into the line.
    class SamReader {
        std::istream* m_is = nullptr;
        std::string_view m_text;
        size_t m_text_pos = 0;
        std::string m_line;
        std::vector<std::string_view> m_tokens;
        SamEntry m_next;
        bool m_has_next = false;  // m_next holds a record read ahead
        size_t m_line_no = 0;
        size_t m_records = 0;
        size_t m_paired_records = 0;
        size_t m_without_tags = 0;
        size_t m_primary_records = 0;  // not secondary (0x100)
        size_t m_primary_without_alternatives = 0;  // of those, without a ZA tag
        std::map<std::string, size_t> m_skipped;
        std::function<void(std::string const&)> m_on_header;  // sees every header line

        // The next line, without its newline (as getline reads it); false at the end.
        bool NextLine(std::string_view& line) {
            if (m_is) {
                if (!std::getline(*m_is, m_line)) return false;
                line = m_line;
                return true;
            }
            if (m_text_pos >= m_text.size()) return false;
            size_t end = m_text.find('\n', m_text_pos);
            if (end == std::string_view::npos) end = m_text.size();
            line = m_text.substr(m_text_pos, end - m_text_pos);
            m_text_pos = end + 1;
            return true;
        }

        // Reads up to and including the next usable record.
        bool Advance(SamEntry& sam) {
            std::string_view line;
            while (NextLine(line)) {
                m_line_no++;
                if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
                if (line.empty()) continue;
                if (line[0] == '@') {
                    if (m_on_header) m_on_header(std::string(line));
                    continue;
                }
                sam_detail::SplitFields(line, m_tokens);
                try {
                    SamFromTokens(m_tokens, sam);
                } catch (SamFormatError const& e) {
                    throw SamFormatError("line " + std::to_string(m_line_no) + ": " + e.what());
                }
                if (auto reason = UnusableRecord(sam)) {
                    m_skipped[reason]++;
                    continue;
                }
                m_records++;
                m_paired_records += Flag::IsPaired(sam.m_flag);
                m_without_tags += !sam_detail::IntTag(m_tokens, "ZU").has_value();
                if (!Flag::IsNotPrimaryAlignment(sam.m_flag)) {
                    m_primary_records++;
                    m_primary_without_alternatives += sam.m_alternatives.empty();
                }
                return true;
            }
            if (m_is && m_is->bad()) throw SamFormatError("read error after line " + std::to_string(m_line_no));
            return false;
        }

    public:
        // `on_header`, if given, is called with each header line and may throw SamFormatError to stop.
        explicit SamReader(std::istream& is, std::function<void(std::string const&)> on_header = {}) :
                m_is(&is), m_on_header(std::move(on_header)) {}

        // The lines of `text`, which must outlive the reader; `first_line`: the lines before it, for the line
        // numbers of errors when the text is a part of a file (the profiler's chunks, SamChunks.h).
        SamReader(std::string_view text, std::function<void(std::string const&)> on_header, size_t first_line) :
                m_text(text), m_line_no(first_line), m_on_header(std::move(on_header)) {}

        bool Next(SamEntry &sam1, SamEntry &sam2, bool &has_sam1, bool &has_sam2) {
            has_sam1 = false;
            has_sam2 = false;
            if (!m_has_next && !Advance(m_next)) return false;
            m_has_next = false;

            if (Flag::IsPaired(m_next.m_flag) && !Flag::IsRead1(m_next.m_flag)) {
                std::swap(sam2, m_next);
                has_sam2 = true;
                return true;
            }
            std::swap(sam1, m_next);
            has_sam1 = true;

            // The mate follows only when it aligned too. SAM files from older protal versions do not
            // set 0x8 on an orphan read1, so the next record is also checked to be this read's read2.
            if (!Flag::IsPaired(sam1.m_flag) || Flag::IsMateUnmapped(sam1.m_flag)) return true;
            if (!Advance(m_next)) return true;
            if (Flag::IsRead2(m_next.m_flag) && m_next.m_qname == sam1.m_qname) {
                std::swap(sam2, m_next);
                has_sam2 = true;
            } else {
                m_has_next = true;
            }
            return true;
        }

        size_t Records() const { return m_records; }
        // Records of paired reads (0x1); single-end reads have none.
        size_t PairedRecords() const { return m_paired_records; }
        size_t Lines() const { return m_line_no; }
        size_t RecordsWithoutTags() const { return m_without_tags; }
        size_t PrimaryRecords() const { return m_primary_records; }
        // Primary and supplementary records without ZA: all of them in a SAM file of an older protal.
        size_t PrimaryRecordsWithoutAlternatives() const { return m_primary_without_alternatives; }
        std::map<std::string, size_t> const& Skipped() const { return m_skipped; }
    };

    // The reads a SAM stream holds: the kind its header names (kSamReadTypeComment, as protal writes it),
    // and whether its first usable record is of paired reads (0x1). protal writes a sample's reads
    // of one kind only. Throws SamFormatError as SamReader.
    struct SamReads {
        std::optional<ReadType> declared;
        std::optional<bool> paired;  // nullopt: no usable record
    };

    inline SamReads ReadsOfSam(std::istream& is) {
        SamReads reads;
        SamReader reader(is, [&reads](std::string const& line) {
            if (line.compare(0, kSamReadTypeComment.size(), kSamReadTypeComment) == 0) {
                reads.declared = ReadTypeFromToken(line.substr(kSamReadTypeComment.size()));
            }
        });
        SamEntry sam1, sam2;
        bool has_sam1 = false, has_sam2 = false;
        if (reader.Next(sam1, sam2, has_sam1, has_sam2)) reads.paired = reader.PairedRecords() > 0;
        return reads;
    }

    // Whether a SAM stream holds paired reads, judged by its first usable record (see ReadsOfSam).
    inline std::optional<bool> HoldsPairedReads(std::istream& is) {
        return ReadsOfSam(is).paired;
    }


    class SamHandler {

    };
}


#endif //PROTAL_SAMHANDLER_H
