//
// Created by fritsche on 17/08/22.
//

#pragma once


#include "Constants.h"
#include <cstdio>
#include <string>
#include "wfa2-lib/bindings/cpp/WFAligner.hpp"


using namespace::wfa;

namespace protal {
    // Gap-affine WFA2-lib (v2.3.6) aligner as protal uses it: the reference is the pattern and the
    // read the text, so in the returned operations 'I' is a read base missing from the reference and
    // 'D' a reference base missing from the read (as in SAM). WFA2's default wf-adaptive heuristic
    // stays on, as it was with the previously vendored v2.3.
    class WFA2Wrapper2 {
        WFAlignerGapAffine m_aligner;

        int m_mismatch;
        int m_gap_opening;
        int m_gap_extension;
        size_t m_x_drop;

        WFAligner::AlignmentStatus m_status = WFAligner::StatusMaxStepsReached;
    public:
        WFA2Wrapper2(int mismatch, int gap_opening, int gap_extension, size_t x_drop) :
                m_aligner(mismatch, gap_opening, gap_extension, WFAligner::Alignment, WFAligner::MemoryHigh),
                m_mismatch(mismatch),
                m_gap_opening(gap_opening),
                m_gap_extension(gap_extension),
                m_x_drop(x_drop) {
        }

        WFA2Wrapper2(const WFA2Wrapper2 &other) :
                m_aligner(other.m_mismatch, other.m_gap_opening, other.m_gap_extension, WFAligner::Alignment,
                          WFAligner::MemoryHigh),
                m_mismatch(other.m_mismatch),
                m_gap_opening(other.m_gap_opening),
                m_gap_extension(other.m_gap_extension),
                m_x_drop(other.m_x_drop) {
        };

        // Per-base operations (M, X, I, D) of the last successful alignment.
        inline std::string Cigar() {
            return m_aligner.getAlignment();
        }

        void Reset() {
            m_status = WFAligner::StatusMaxStepsReached;
        }

        bool Success() const {
            return m_status == WFAligner::StatusAlgCompleted;
        }

        // Ends-free alignment of the read `query` (text) against `ref` (pattern). An alignment whose
        // score reaches max_score is abandoned: WFA2 counts this limit in steps, which for the
        // gap-affine model is the same cut-off as the max score of WFA2-lib v2.3.
        void Alignment(std::string const& query, std::string const& ref, int q_begin_free, int q_end_free,
                       int r_begin_free, int r_end_free, int max_score=INT32_MAX) {
            m_aligner.setMaxAlignmentSteps(max_score);
            m_status = m_aligner.alignEndsFree(ref, r_begin_free, r_end_free, query, q_begin_free, q_end_free);
        }

        int GetAlignmentScore() {
            return m_aligner.getAlignmentScore();
        }

        WFAlignerGapAffine &GetAligner() {
            return m_aligner;
        }

        // Debugging: the last alignment of `query` against `ref`, drawn by WFA2.
        void PrintAlignment(std::string const& query, std::string const& ref, FILE* stream = stderr) {
            std::fprintf(stream, "%s (score %d)\n", m_aligner.getAlignment().c_str(), m_aligner.getAlignmentScore());
            m_aligner.printPretty(stream, ref.c_str(), static_cast<int>(ref.length()),
                                  query.c_str(), static_cast<int>(query.length()));
        }
    };
}
