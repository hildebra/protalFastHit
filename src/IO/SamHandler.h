//
// Created by fritsche on 26/08/22.
//

#ifndef PROTAL_SAMHANDLER_H
#define PROTAL_SAMHANDLER_H


#include <cstdint>
#include <string>
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


    static void SamFromTokens(std::vector<std::string>& tokens, SamEntry &sam) {
        if (tokens.size() < 11) {
            std::cerr << omp_get_thread_num() << " Broken SAM file "<< std::endl;
            exit(120);
        }
        sam.m_qname = tokens[0];
        sam.m_flag = std::stoul(tokens[1]);
        sam.m_rname = tokens[2];
        sam.m_pos = std::stoul(tokens[3]);
        sam.m_mapq = std::stoul(tokens[4]);
        sam.m_cigar = tokens[5];
        sam.m_rnext = tokens[6];
        sam.m_pnext = std::stoul(tokens[7]);
        sam.m_tlen = std::stol(tokens[8]);
        sam.m_seq = tokens[9];
        sam.m_qual = tokens[10];

        if (tokens.size() < 13) {
            // std::cerr << omp_get_thread_num() << " Sam file lacks unique information" << std::endl;
            sam.m_uniques = 0;
            sam.m_uniques_two = 0;
            return;
        }

        const auto s = tokens[11];
        sam.m_uniques = s.size() > 4 && s[4] == ':' ? std::stoul(s.substr(5)) : 0;
        const auto s2 = tokens[12];
        sam.m_uniques_two = s2.size() > 4 && s2[4] == ':' ? std::stoul(s2.substr(5)) : 0;

    }



    static bool GetSam(std::istream &file, std::string &line, std::vector<std::string> &tokens, SamEntry &sam) {
        static std::string delim = "\t";
        if (!std::getline(file, line)) return false;
        LineSplitter::Split(line, delim, tokens);
        SamFromTokens(tokens, sam);
        return true;
    }

    // Reads the next alignment: a read1 record together with its mate when the mate record follows,
    // or a single record. `line_loaded` is the read-ahead state shared across calls: true on entry
    // means `line` already holds the next record; true on return means a record was read ahead
    // (a paired read1's next line that turned out not to be its mate) and is kept for the next call.
    static bool GetSamPair(std::istream &file, std::string &line, std::vector<std::string> &tokens, SamEntry &sam1, SamEntry &sam2, bool &has_sam1, bool &has_sam2, bool &line_loaded) {
        static std::string delim = "\t";
        has_sam1 = false;
        has_sam2 = false;
        if (!line_loaded && !std::getline(file, line)) return false;
        line_loaded = false;

        LineSplitter::Split(line, delim, tokens);
        if (tokens.size() < 2) return false;

        if (!Flag::IsRead1(stoul(tokens[1]))) {
            SamFromTokens(tokens, sam2);
            has_sam2 = true;
            return true;
        }

        SamFromTokens(tokens, sam1);
        has_sam1 = true;

        // The mate follows only when it aligned too. SAM files from older protal versions do not
        // set 0x8 on an orphan read1, so the next line is also checked to be this read's read2.
        if (!Flag::IsPaired(sam1.m_flag) || Flag::IsMateUnmapped(sam1.m_flag)) return true;
        if (!std::getline(file, line)) return true;

        LineSplitter::Split(line, delim, tokens);
        if (tokens.size() >= 2 && Flag::IsRead2(stoul(tokens[1])) && tokens[0] == sam1.m_qname) {
            SamFromTokens(tokens, sam2);
            has_sam2 = true;
        } else {
            line_loaded = true;
        }
        return true;
    }



    class SamHandler {

    };
}


#endif //PROTAL_SAMHANDLER_H
