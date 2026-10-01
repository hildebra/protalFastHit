//
// Created by fritsche on 17/08/22.
//

#pragma once


#include "Constants.h"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <string>
#include "wfa2-lib/bindings/cpp/WFAligner.hpp"


using namespace::wfa;

namespace protal {
    // Gap-affine WFA2-lib (v2.3.6) aligner as protal uses it: the reference is the pattern and the
    // read the text, so in the returned operations 'I' is a read base missing from the reference and
    // 'D' a reference base missing from the read (as in SAM). WFA2's default wf-adaptive heuristic
    // stays on, as it was with the previously vendored v2.3.
    //
    // x_drop (--x_drop) adds WFA2's X-drop to wf-adaptive; 0 leaves it off. The two share WFA2's
    // step counter, so X-drop prunes only at steps where wf-adaptive did not (short wavefronts).
    // With protal's scores (a match costs 0), X-drop alone stops good long alignments (a 2 kb read
    // at 2% divergence with 50); added to wf-adaptive, 50-1000 changed no 150 bp alignment in tests,
    // and 200 or less gave an 8 kb one a worse score. The default, 1000, is for short reads: over the
    // gene-long windows of long reads it lost the own species' alignment of genes an ONT read ends
    // in, so RunProtal aligns long reads without X-drop.
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
            ApplyXDrop();
        }

        WFA2Wrapper2(const WFA2Wrapper2 &other) :
                m_aligner(other.m_mismatch, other.m_gap_opening, other.m_gap_extension, WFAligner::Alignment,
                          WFAligner::MemoryHigh),
                m_mismatch(other.m_mismatch),
                m_gap_opening(other.m_gap_opening),
                m_gap_extension(other.m_gap_extension),
                m_x_drop(other.m_x_drop) {
            ApplyXDrop();
        };

        size_t XDrop() const {
            return m_x_drop;
        }

        // Per-base operations (M, X, I, D) of the last successful alignment.
        inline std::string Cigar() {
            return m_aligner.getAlignment();
        }

        // The same operations written into `out` (replaced, or appended to with append), from WFA2's own buffer:
        // no string of its own per alignment.
        void CigarInto(std::string& out, bool append = false) {
            char* operations = nullptr;
            int length = 0;
            m_aligner.getAlignment(&operations, &length);
            if (append) out.append(operations, static_cast<size_t>(length));
            else out.assign(operations, static_cast<size_t>(length));
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
            // X-drop can report an alignment complete whose operations do not fit the read (151 read
            // bases for 150, seen with --x_drop 50): count it as dropped, as WFA2 reports a drop otherwise.
            if (m_x_drop > 0 && m_status == WFAligner::StatusAlgCompleted && !CoversText(query.size())) {
                m_status = WFAligner::StatusAlgPartial;
            }
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

    private:
        // Every score step, as wf-adaptive (1 step between cut-offs: the steps are shared).
        void ApplyXDrop() {
            if (m_x_drop > 0) m_aligner.setHeuristicXDrop(static_cast<int>(std::min<size_t>(m_x_drop, INT32_MAX)), 1);
        }

        // Whether the last alignment consumed exactly text_length read bases (every operation but 'D').
        bool CoversText(size_t text_length) {
            char* operations = nullptr;
            int length = 0;
            m_aligner.getAlignment(&operations, &length);
            size_t consumed = 0;
            for (int i = 0; i < length; i++) consumed += operations[i] != 'D';
            return consumed == text_length;
        }
    };
}
