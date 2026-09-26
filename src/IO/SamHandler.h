//
// Created by fritsche on 26/08/22.
//

#ifndef PROTAL_SAMHANDLER_H
#define PROTAL_SAMHANDLER_H


#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstdint>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include "LineSplitter.h"
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
                    + "ZT:i:" + std::to_string(m_uniques_two);
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
        inline bool ParseUnsigned(std::string const& s, uint64_t& value) {
            auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), value);
            return !s.empty() && ec == std::errc() && ptr == s.data() + s.size();
        }

        inline bool ParseSigned(std::string const& s, int64_t& value) {
            auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), value);
            return !s.empty() && ec == std::errc() && ptr == s.data() + s.size();
        }

        // Value of an integer tag ("ZU:i:12") among the optional fields, or nullopt if it is absent
        // or not a non-negative integer.
        inline std::optional<uint64_t> IntTag(std::vector<std::string> const& tokens, std::string_view name) {
            for (size_t i = 11; i < tokens.size(); i++) {
                auto const& t = tokens[i];
                if (t.size() > 5 && t.compare(0, 2, name) == 0 && t.compare(2, 3, ":i:") == 0) {
                    uint64_t value;
                    if (ParseUnsigned(t.substr(5), value)) return value;
                    return std::nullopt;
                }
            }
            return std::nullopt;
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

    // Rewrites a CIGAR into the ops the profiler walks, which are those protal writes: M for exact
    // matches, X, I, D and S. '=' (sequence match) becomes M, and hard clips are dropped because they
    // consume neither SEQ nor the reference. Returns false for a CIGAR that cannot be walked:
    // malformed, with N or P ops, or not covering exactly the bases in SEQ.
    inline bool NormalizeCigar(std::string& cigar, size_t seq_length) {
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
            if (j == cigar.size() || !sam_detail::ParseUnsigned(cigar.substr(i, j - i), count) || count == 0) return false;
            char op = cigar[j];
            i = j + 1;
            switch (op) {
                case '=': op = 'M'; break;
                case 'M': case 'X': case 'I': case 'D': case 'S': break;
                case 'H': continue;
                default: return false;
            }
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
        return true;
    }

    // Fills sam from the fields of one SAM line; throws SamFormatError if the line cannot be parsed.
    // ZU and ZT (protal's unique k-mer counts) are looked up by name and are 0 when absent.
    inline void SamFromTokens(std::vector<std::string> const& tokens, SamEntry &sam) {
        if (tokens.size() < 11) {
            throw SamFormatError("expected at least 11 tab-separated fields, found " + std::to_string(tokens.size()));
        }
        uint64_t flag, pos, mapq, pnext;
        int64_t tlen;
        if (!sam_detail::ParseUnsigned(tokens[1], flag) || flag > UINT16_MAX) throw SamFormatError("FLAG is not a 16-bit integer: " + tokens[1]);
        if (!sam_detail::ParseUnsigned(tokens[3], pos) || pos > UINT32_MAX) throw SamFormatError("POS is not a position: " + tokens[3]);
        if (!sam_detail::ParseUnsigned(tokens[4], mapq) || mapq > 255) throw SamFormatError("MAPQ is not between 0 and 255: " + tokens[4]);
        if (!sam_detail::ParseUnsigned(tokens[7], pnext) || pnext > UINT32_MAX) throw SamFormatError("PNEXT is not a position: " + tokens[7]);
        if (!sam_detail::ParseSigned(tokens[8], tlen)) throw SamFormatError("TLEN is not an integer: " + tokens[8]);

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
    }

    // Why the profiler cannot use a parsed record, or nullptr if it can. Normalizes the CIGAR.
    inline char const* UnusableRecord(SamEntry& sam) {
        if (Flag::IsUnmapped(sam.m_flag) || sam.m_rname == "*") return "unmapped";
        if (!sam_detail::IsProtalGene(sam.m_rname)) return "reference is not a protal gene (<taxid>_<gene id>)";
        if (sam.m_pos == 0) return "no position";
        if (sam.m_seq == "*") return "no sequence";
        if (sam.m_qual.size() != sam.m_seq.size()) return "no base qualities, or QUAL and SEQ lengths differ";
        if (!NormalizeCigar(sam.m_cigar, sam.m_seq.size())) return "CIGAR is missing, malformed, has N/P ops or does not match SEQ";
        return nullptr;
    }

    // Reads the alignment records of a SAM stream. Header lines and blank lines are skipped, and so
    // are records the profiler cannot use (counted by reason in Skipped()). Next() returns a read1
    // record together with its mate when the mate's record follows, or a single record. A line that
    // cannot be parsed throws SamFormatError naming the line.
    class SamReader {
        std::istream& m_is;
        std::string m_line;
        std::vector<std::string> m_tokens;
        SamEntry m_next;
        bool m_has_next = false;  // m_next holds a record read ahead
        size_t m_line_no = 0;
        size_t m_records = 0;
        size_t m_without_tags = 0;
        std::map<std::string, size_t> m_skipped;

        // Reads up to and including the next usable record.
        bool Advance(SamEntry& sam) {
            static const std::string delim = "\t";
            while (std::getline(m_is, m_line)) {
                m_line_no++;
                if (!m_line.empty() && m_line.back() == '\r') m_line.pop_back();
                if (m_line.empty() || m_line[0] == '@') continue;
                LineSplitter::Split(m_line, delim, m_tokens);
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
                m_without_tags += !sam_detail::IntTag(m_tokens, "ZU").has_value();
                return true;
            }
            if (m_is.bad()) throw SamFormatError("read error after line " + std::to_string(m_line_no));
            return false;
        }

    public:
        explicit SamReader(std::istream& is) : m_is(is) {}

        bool Next(SamEntry &sam1, SamEntry &sam2, bool &has_sam1, bool &has_sam2) {
            has_sam1 = false;
            has_sam2 = false;
            if (!m_has_next && !Advance(m_next)) return false;
            m_has_next = false;

            if (!Flag::IsRead1(m_next.m_flag)) {
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
        size_t Lines() const { return m_line_no; }
        size_t RecordsWithoutTags() const { return m_without_tags; }
        std::map<std::string, size_t> const& Skipped() const { return m_skipped; }
    };


    class SamHandler {

    };
}


#endif //PROTAL_SAMHANDLER_H
