#pragma once

#include <algorithm>
#include <cstddef>
#include <string>
#include "ChainingStrategy.h"
#include "WFA2Wrapper2.h"

namespace protal {

    // The window a read is aligned into, as SimpleAlignmentHandler::AlignAnchor sets it up around
    // its anchor: gene positions [ref_start, ref_end), how many reference bases are free at each
    // end (the dovetails), how many read bases are free at each end (where the read runs past the
    // gene), and the penalty the alignment must stay below.
    struct AlignmentWindow {
        size_t ref_start = 0;
        size_t ref_end = 0;
        int ref_begin_free = 0;
        int ref_end_free = 0;
        int read_begin_free = 0;
        int read_end_free = 0;
        int max_score = 0;
    };

    // Aligns a read into its window from the exact matches of its anchor chain, instead of aligning
    // the whole read into the whole window. The links are taken as they are; the read left of the
    // first link is aligned to the window left of it (reversed, so that the alignment is anchored at
    // the link and its free end is the window's start), the read right of the last link likewise,
    // and between two links the bases in between. Most of a read's anchors are one link of 20-50
    // bases on a relative's gene: anchored at it, WFA gives up on such a flank as soon as it exceeds
    // the window's budget, where aligning the whole read started from every free diagonal at once.
    //
    // The result is the operations the whole-window alignment would give when its best path runs
    // through the links (M, X, I, D over the window, free ends included), so the caller
    // post-processes it the same way. The budget is shared: the pieces' penalties add up to less than
    // the window's max_score, as for a whole alignment. Chains this does not handle are left to the
    // whole-window alignment (NotApplicable): links on different diagonals (the seeds imply an
    // indel), out of order, outside the window, or not exact matches apart from Ns.
    class AnchoredAligner {
    public:
        enum class Status { Aligned, Failed, NotApplicable };

        Status Align(std::string const& read, std::string_view const gene, ChainList const& chain, AlignmentWindow const& w,
                     WFA2Wrapper2& aligner, std::string& ops) {
            ops.clear();
            if (chain.empty() || w.ref_end > gene.size() || w.ref_start > w.ref_end) return Status::NotApplicable;
            long const diagonal = static_cast<long>(chain.front().genepos) - static_cast<long>(chain.front().readpos);
            for (size_t i = 0; i < chain.size(); i++) {
                auto const& link = chain[i];
                if (link.length == 0 || static_cast<long>(link.genepos) - static_cast<long>(link.readpos) != diagonal) return Status::NotApplicable;
                if (link.ReadEnd() > read.size() || link.GeneStart() < w.ref_start || link.GeneEnd() > w.ref_end) return Status::NotApplicable;
                if (i > 0 && link.ReadStart() < chain[i - 1].ReadEnd()) return Status::NotApplicable;
            }

            int used = 0;
            // Left of the first link, reversed: anchored at the link, free at the window's start.
            m_left.clear();
            auto const& first = chain.front();
            if (Flank(read, gene, 0, first.ReadStart(), w.ref_start, first.GeneStart(), true, w.ref_begin_free, w.read_begin_free,
                      w.max_score, aligner, used, m_left) != Status::Aligned) {
                return m_status;
            }
            ops.append(m_left.rbegin(), m_left.rend());

            for (size_t i = 0; i < chain.size(); i++) {
                auto const& link = chain[i];
                // The link: exact but for Ns, which count as mismatches as in any alignment.
                for (uint32_t k = 0; k < link.length; k++) {
                    char const q = read[link.readpos + k], r = gene[link.genepos + k];
                    if (q == r) { ops += 'M'; continue; }
                    if (q != 'N' && r != 'N') return Status::NotApplicable;
                    ops += 'X';
                    if ((used += kMismatch) >= w.max_score) return Status::Failed;
                }
                if (i + 1 < chain.size()) {
                    // Up to the next link on the same diagonal: as many read as gene bases.
                    size_t const from = link.ReadEnd(), to = chain[i + 1].ReadStart();
                    if (Between(read, gene, from, to, link.GeneEnd(), w.max_score, aligner, used, ops) != Status::Aligned) return m_status;
                }
            }

            // Right of the last link: anchored at the link, free at the window's end.
            auto const& last = chain.back();
            m_right.clear();
            if (Flank(read, gene, last.ReadEnd(), read.size(), last.GeneEnd(), w.ref_end, false, w.ref_end_free, w.read_end_free,
                      w.max_score, aligner, used, m_right) != Status::Aligned) {
                return m_status;
            }
            ops += m_right;

            // The operations cover the read and the window exactly, as a whole-window alignment's do.
            size_t read_bases = 0, ref_bases = 0;
            for (char c : ops) {
                read_bases += c != 'D';
                ref_bases += c != 'I';
            }
            if (read_bases != read.size() || ref_bases != w.ref_end - w.ref_start) return Status::NotApplicable;
            return Status::Aligned;
        }

    private:
        static constexpr int kMismatch = 4;  // as SimpleAlignmentHandler sets up WFA2Wrapper2 (4, 6, 2)

        std::string m_query, m_ref, m_left, m_right, m_between;
        Status m_status = Status::Aligned;

        // A flank: read [read_from, read_to) against gene [ref_from, ref_to), anchored at the link
        // (the flank's end for the left flank, its start for the right one) and free at the other
        // end by up to ref_free reference and read_free read bases. The left flank is aligned
        // reversed; its operations come back reversed too, for the caller to turn around.
        Status Flank(std::string const& read, std::string_view const gene, size_t read_from, size_t read_to, size_t ref_from,
                     size_t ref_to, bool reversed, int ref_free, int read_free, int max_score, WFA2Wrapper2& aligner,
                     int& used, std::string& out) {
            size_t const r = read_to - read_from, g = ref_to - ref_from;
            if (r == 0 && g == 0) return m_status = Status::Aligned;
            if (r == 0 || g == 0) {
                // Only reference or only read bases: free if the window allows that many at this end.
                if (r == 0 && g <= static_cast<size_t>(std::max(ref_free, 0))) { out.append(g, 'D'); return m_status = Status::Aligned; }
                if (g == 0 && r <= static_cast<size_t>(std::max(read_free, 0))) { out.append(r, 'I'); return m_status = Status::Aligned; }
                return m_status = Status::NotApplicable;
            }
            if (reversed) {
                m_query.assign(read.rbegin() + static_cast<long>(read.size() - read_to), read.rbegin() + static_cast<long>(read.size() - read_from));
                m_ref.assign(gene.rbegin() + static_cast<long>(gene.size() - ref_to), gene.rbegin() + static_cast<long>(gene.size() - ref_from));
            } else {
                m_query.assign(read, read_from, r);
                m_ref.assign(gene, ref_from, g);
            }
            aligner.Reset();
            aligner.Alignment(m_query, m_ref, 0, std::min<int>(std::max(read_free, 0), static_cast<int>(r)),
                              0, std::min<int>(std::max(ref_free, 0), static_cast<int>(g)), max_score - used);
            if (!aligner.Success()) return m_status = Status::Failed;
            used += -aligner.GetAlignmentScore();
            out = aligner.Cigar();
            return m_status = Status::Aligned;
        }

        // The bases between two links on one diagonal: read [from, to) against as many gene bases
        // from gene_from, both ends fixed. With up to 3 mismatches the ungapped alignment is the best
        // (12 against at least 16 for the insertion and deletion any gapped one needs); else WFA.
        Status Between(std::string const& read, std::string_view const gene, size_t from, size_t to, size_t gene_from, int max_score,
                       WFA2Wrapper2& aligner, int& used, std::string& ops) {
            size_t const n = to - from;
            if (n == 0) return m_status = Status::Aligned;
            int mismatches = 0;
            for (size_t k = 0; k < n; k++) mismatches += read[from + k] != gene[gene_from + k];
            if (mismatches <= 3) {
                for (size_t k = 0; k < n; k++) ops += read[from + k] == gene[gene_from + k] ? 'M' : 'X';
                used += kMismatch * mismatches;
                return m_status = used < max_score ? Status::Aligned : Status::Failed;
            }
            m_query.assign(read, from, n);
            m_ref.assign(gene, gene_from, n);
            aligner.Reset();
            aligner.Alignment(m_query, m_ref, 0, 0, 0, 0, max_score - used);
            if (!aligner.Success()) return m_status = Status::Failed;
            used += -aligner.GetAlignmentScore();
            m_between = aligner.Cigar();
            ops += m_between;
            return m_status = Status::Aligned;
        }
    };
}
