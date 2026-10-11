//
// Created by fritsche on 17/08/22.
//

#pragma once


#include "Constants.h"
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <numeric>
#include <string>
#include "WfaTraceback.h"
#include "wfa2-lib/bindings/cpp/WFAligner.hpp"


using namespace::wfa;

namespace protal {
    // Gap-affine WFA2-lib aligner as protal uses it (lib/wfa2-lib.cmake: v2.3.6 with the Singletrack
    // traceback): the reference is the pattern and the read the text, so in the returned operations 'I'
    // is a read base missing from the reference and 'D' a reference base missing from the read (as in
    // SAM). WFA2's default wf-adaptive heuristic stays on, as it was with the previously vendored v2.3.
    //
    // traceback (--wfa_traceback, WfaTraceback.h): WFA2's high-memory mode (classic), or Singletrack,
    // which keeps only the M wavefronts. The scores are the same; among equally good alignments the
    // operations may differ.
    //
    // x_drop (--x_drop) adds WFA2's X-drop to wf-adaptive; 0 leaves it off. The two share WFA2's
    // step counter, so X-drop prunes only at steps where wf-adaptive did not (short wavefronts).
    // With protal's scores (a match costs 0), X-drop alone stops good long alignments (a 2 kb read
    // at 2% divergence with 50); added to wf-adaptive, 50-1000 changed no 150 bp alignment in tests,
    // and 200 or less gave an 8 kb one a worse score. The default, 1000, is for short reads: over the
    // gene-long windows of long reads it lost the own species' alignment of genes an ONT read ends
    // in, so RunProtal aligns long reads without X-drop.
    //
    // Singletrack keeps only the M wavefronts, so its backtrace can fail where wf-adaptive pruned the
    // wavefronts a path ran through (lib/wfa2-lib.cmake: protal's patch makes that rare and reports it
    // instead of exiting). Such an alignment is aligned again with the classic traceback, which gives
    // the same score; SingletrackFallbacks() counts them.
    //
    // Two aligners, the same alignments (docs/claude/2026-10-11-wfa2-singletrack-speed):
    // - The fast one (most calls) aligns with the penalties divided by their greatest common divisor
    //   (4/6/2 -> 2/3/1) and without X-drop. With 4/6/2 every score is even, so WFA2 spent every odd
    //   score step on empty wavefronts; with the score scaled back and the budget ceil(max/g), the
    //   alignments are the same (wf-adaptive counts bases and steps that reach a wavefront, not scores).
    //   X-drop is left out where it cannot prune: with a match score of 0, WFA2's X-drop scores a state
    //   (-(v+h)-score)/2 (a "match bonus" of -1), so the drop it measures is at most (v+h+score)/2, below
    //   x_drop whenever the read part, the reference part and the budget add up to less than 2*x_drop
    //   (every short read at the default 1000). It does not change when wf-adaptive runs either.
    // - The X-drop one, with the penalties as given, for the rest: X-drop's scores do not scale.
    class WFA2Wrapper2 {
        int m_mismatch;
        int m_gap_opening;
        int m_gap_extension;
        size_t m_x_drop;
        WfaTraceback m_traceback;
        int m_scale;  // the penalties' greatest common divisor: the fast aligner's penalties are divided by it

        WFAlignerGapAffine m_aligner;                         // the fast aligner
        std::unique_ptr<WFAlignerGapAffine> m_xdrop_aligner;  // the X-drop aligner, made when first needed
        std::unique_ptr<WFAlignerGapAffine> m_classic;        // Singletrack's fallbacks for either, made when needed
        std::unique_ptr<WFAlignerGapAffine> m_xdrop_classic;
        WFAlignerGapAffine* m_last = &m_aligner;              // the aligner of the last alignment
        int m_last_scale = 1;                                 // and the factor of its scores

        WFAligner::AlignmentStatus m_status = WFAligner::StatusMaxStepsReached;

        static inline std::atomic<size_t> s_singletrack_fallbacks{ 0 };

        static WFAligner::MemoryModel Memory(WfaTraceback traceback) {
            return traceback == WfaTraceback::Singletrack ? WFAligner::MemorySingletrack : WFAligner::MemoryHigh;
        }

        static int Scale(int mismatch, int gap_opening, int gap_extension) {
            int const scale = std::gcd(std::gcd(mismatch, gap_opening), gap_extension);
            return scale > 0 ? scale : 1;
        }
    public:
        WFA2Wrapper2(int mismatch, int gap_opening, int gap_extension, size_t x_drop,
                     WfaTraceback traceback = WfaTraceback::Classic) :
                m_mismatch(mismatch),
                m_gap_opening(gap_opening),
                m_gap_extension(gap_extension),
                m_x_drop(x_drop),
                m_traceback(traceback),
                m_scale(Scale(mismatch, gap_opening, gap_extension)),
                m_aligner(mismatch / m_scale, gap_opening / m_scale, gap_extension / m_scale, WFAligner::Alignment,
                          Memory(traceback)) {
        }

        WFA2Wrapper2(const WFA2Wrapper2 &other) :
                WFA2Wrapper2(other.m_mismatch, other.m_gap_opening, other.m_gap_extension, other.m_x_drop, other.m_traceback) {
        };

        // Alignments of all aligners whose Singletrack backtrace failed and that were aligned again classic.
        static size_t SingletrackFallbacks() {
            return s_singletrack_fallbacks.load(std::memory_order_relaxed);
        }

        size_t XDrop() const {
            return m_x_drop;
        }

        WfaTraceback Traceback() const {
            return m_traceback;
        }

        // The gap-affine penalties (a match costs 0).
        int Mismatch() const { return m_mismatch; }
        int GapOpening() const { return m_gap_opening; }
        int GapExtension() const { return m_gap_extension; }

        // Per-base operations (M, X, I, D) of the last successful alignment.
        inline std::string Cigar() {
            return m_last->getAlignment();
        }

        // The same operations written into `out` (replaced, or appended to with append), from WFA2's own buffer:
        // no string of its own per alignment.
        void CigarInto(std::string& out, bool append = false) {
            char* operations = nullptr;
            int length = 0;
            m_last->getAlignment(&operations, &length);
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
            bool const xdrop = XDropMayPrune(query.size(), ref.size(), max_score);
            WFAlignerGapAffine& aligner = Select(xdrop);
            int const steps = Steps(max_score);
            aligner.setMaxAlignmentSteps(steps);
            m_status = aligner.alignEndsFree(ref, r_begin_free, r_end_free, query, q_begin_free, q_end_free);
            if (BacktraceFailed()) {
                WFAlignerGapAffine& classic = Classic(xdrop);
                classic.setMaxAlignmentSteps(steps);
                m_status = classic.alignEndsFree(ref, r_begin_free, r_end_free, query, q_begin_free, q_end_free);
            }
            // X-drop can report an alignment complete whose operations do not fit the read (151 read
            // bases for 150, seen with --x_drop 50): count it as dropped, as WFA2 reports a drop otherwise.
            if (xdrop && m_status == WFAligner::StatusAlgCompleted && !CoversText(query.size())) {
                m_status = WFAligner::StatusAlgPartial;
            }
        }

        // End-to-end alignment of the read `query` (text) against `ref` (pattern), both whole, with WFA2's end-to-end
        // kernels (no ends-free bookkeeping); abandoned at max_score as Alignment is.
        void EndToEnd(std::string const& query, std::string const& ref, int max_score=INT32_MAX) {
            bool const xdrop = XDropMayPrune(query.size(), ref.size(), max_score);
            WFAlignerGapAffine& aligner = Select(xdrop);
            int const steps = Steps(max_score);
            aligner.setMaxAlignmentSteps(steps);
            m_status = aligner.alignEnd2End(ref, query);
            if (BacktraceFailed()) {
                WFAlignerGapAffine& classic = Classic(xdrop);
                classic.setMaxAlignmentSteps(steps);
                m_status = classic.alignEnd2End(ref, query);
            }
            if (xdrop && m_status == WFAligner::StatusAlgCompleted && !CoversText(query.size())) {
                m_status = WFAligner::StatusAlgPartial;
            }
        }

        // The score of the last successful alignment (negative: WFA2's penalty), on the scale of the penalties as given.
        int GetAlignmentScore() {
            return m_last->getAlignmentScore() * m_last_scale;
        }

        // The fast aligner (penalties divided by their greatest common divisor, no X-drop), e.g. to change its heuristic.
        WFAlignerGapAffine &GetAligner() {
            return m_aligner;
        }

        // Whether X-drop could prune an alignment of `text_length` read and `pattern_length` reference bases within
        // max_score (see the class comment): then the X-drop aligner aligns it, else the fast one. WFA2 prunes at the
        // scores it extends, below max_score and 0 (always extended).
        bool XDropMayPrune(size_t text_length, size_t pattern_length, int max_score) const {
            if (m_x_drop == 0) return false;
            return static_cast<int64_t>(text_length) + static_cast<int64_t>(pattern_length) + std::max(max_score, 0) >=
                   2 * static_cast<int64_t>(std::min<size_t>(m_x_drop, INT32_MAX));
        }

        // Debugging: the last alignment of `query` against `ref`, drawn by WFA2.
        void PrintAlignment(std::string const& query, std::string const& ref, FILE* stream = stderr) {
            std::fprintf(stream, "%s (score %d)\n", m_last->getAlignment().c_str(), GetAlignmentScore());
            m_last->printPretty(stream, ref.c_str(), static_cast<int>(ref.length()),
                                query.c_str(), static_cast<int>(query.length()));
        }

    private:
        // WFA2's step limit for a budget of max_score on the given penalties' scale. WFA2 stops when the score reaches
        // the limit (score >= steps), so an alignment of score s (a multiple of m_scale) succeeds iff s < max_score,
        // i.e. iff s/m_scale <= (max_score-1)/m_scale. A budget of 0 or less stops every alignment but one of score 0
        // at either scale. Only the fast aligner is scaled; the X-drop one takes max_score as it is.
        int Steps(int max_score) const {
            if (m_last_scale == 1 || max_score <= 0) return max_score;
            return (max_score - 1) / m_last_scale + 1;
        }

        // The aligner of the next alignment (the X-drop one when X-drop could prune, else the fast one); the last one
        // from now.
        WFAlignerGapAffine& Select(bool xdrop) {
            if (!xdrop) {
                m_last = &m_aligner;
                m_last_scale = m_scale;
                return m_aligner;
            }
            if (!m_xdrop_aligner) m_xdrop_aligner = MakeXDropAligner(Memory(m_traceback));
            m_last = m_xdrop_aligner.get();
            m_last_scale = 1;
            return *m_xdrop_aligner;
        }

        // An aligner with the penalties as given and X-drop at every score step, as wf-adaptive (1 step between
        // cut-offs: the steps are shared).
        std::unique_ptr<WFAlignerGapAffine> MakeXDropAligner(WFAligner::MemoryModel memory) const {
            auto aligner = std::make_unique<WFAlignerGapAffine>(m_mismatch, m_gap_opening, m_gap_extension, WFAligner::Alignment,
                                                                memory);
            aligner->setHeuristicXDrop(static_cast<int>(std::min<size_t>(m_x_drop, INT32_MAX)), 1);
            return aligner;
        }

        // Whether the Singletrack backtrace of the last alignment found no path (the patched WFA2 reports it as
        // unattainable); then the classic aligner is the last one, and the failure is counted.
        bool BacktraceFailed() {
            if (m_traceback != WfaTraceback::Singletrack || static_cast<int>(m_status) != WF_STATUS_UNATTAINABLE) return false;
            s_singletrack_fallbacks.fetch_add(1, std::memory_order_relaxed);
            return true;
        }

        // The classic aligner of Singletrack's fallback, configured as the aligner that failed; the last aligner from now
        // (the score scale stays).
        WFAlignerGapAffine& Classic(bool xdrop) {
            auto& classic = xdrop ? m_xdrop_classic : m_classic;
            if (!classic) {
                classic = xdrop ? MakeXDropAligner(WFAligner::MemoryHigh)
                                : std::make_unique<WFAlignerGapAffine>(m_mismatch / m_scale, m_gap_opening / m_scale,
                                                                       m_gap_extension / m_scale, WFAligner::Alignment,
                                                                       WFAligner::MemoryHigh);
            }
            m_last = classic.get();
            return *classic;
        }

        // Whether the last alignment consumed exactly text_length read bases (every operation but 'D').
        bool CoversText(size_t text_length) {
            char* operations = nullptr;
            int length = 0;
            m_last->getAlignment(&operations, &length);
            size_t consumed = 0;
            for (int i = 0; i < length; i++) consumed += operations[i] != 'D';
            return consumed == text_length;
        }
    };
}
