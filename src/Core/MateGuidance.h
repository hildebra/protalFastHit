#pragma once

// A fragment whose mates did not align together: the mate that is sure of its alignment guides the other one to its
// taxon, and on its gene, where the fragment can reach. A pair is one organism; the other mate otherwise counts for a
// relative (it has no candidate of the right taxon among the few that are aligned), or not at all (no anchor: too
// divergent, or partly beyond the gene's end), and its bases are lost to the strain's SNPs.

#include <algorithm>
#include <cstdint>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include "AcrossGenes.h"
#include "AlignmentOutputHandler.h"
#include "ChainingStrategy.h"
#include "GenomeLoader.h"

namespace protal {
    inline constexpr int kGuideMinMapq = 20;              // a mate this sure of its best alignment guides the other
    inline constexpr int64_t kRescueMaxFragment = 1000;   // the farthest a fragment reaches from its guiding mate
    inline constexpr int64_t kRescueMinBases = 30;        // the rest of the gene the fragment can reach, at least
    inline constexpr size_t kRescueK = 12;                // k-mers that place the other mate on the gene
    inline constexpr size_t kRescueMinHits = 2;           // of them on one diagonal, at least

    // What mate guidance did, over the fragments of a run.
    struct MateGuidanceCounts {
        size_t guided = 0;    // fragments with a guiding mate and another without a candidate of its taxon
        size_t from_anchor = 0;  // other mates aligned from their own anchor of that taxon
        size_t looked_for = 0;   // other mates looked for on the guiding mate's gene
        size_t rescued = 0;      // of these, found and aligned there
        size_t rescued_on_neighbour = 0;  // of those, on a gene next to the guiding mate's (the database's gene neighbours)
        size_t paired_across_genes = 0;   // fragments whose best candidate has its mates on two neighbouring genes

        void Join(MateGuidanceCounts const& other) {
            guided += other.guided;
            from_anchor += other.from_anchor;
            looked_for += other.looked_for;
            rescued += other.rescued;
            rescued_on_neighbour += other.rescued_on_neighbour;
            paired_across_genes += other.paired_across_genes;
        }
    };

    namespace mate_guidance {
        inline int Base(char c) {
            switch (c) {
                case 'A': case 'a': return 0;
                case 'C': case 'c': return 1;
                case 'G': case 'g': return 2;
                case 'T': case 't': return 3;
                default: return -1;
            }
        }

        // (code, position) of each k-mer of `seq` without an ambiguous base.
        inline void Kmers(std::string_view seq, size_t k, std::vector<std::pair<uint32_t, uint32_t>>& out) {
            out.clear();
            uint32_t const mask = static_cast<uint32_t>((uint64_t{1} << (2 * k)) - 1);
            uint32_t code = 0;
            size_t valid = 0;
            for (size_t i = 0; i < seq.size(); i++) {
                int const b = Base(seq[i]);
                if (b < 0) {
                    valid = 0;
                    code = 0;
                    continue;
                }
                code = ((code << 2) | static_cast<uint32_t>(b)) & mask;
                if (++valid >= k) out.emplace_back(code, static_cast<uint32_t>(i + 1 - k));
            }
        }

        // The diagonal (gene position minus read position) on which most k-mers of `read` match `window` (the gene's
        // bases from gene position `window_start`): its offset, the matches on it, and the read position of the first.
        struct Diagonal {
            int64_t offset = 0;
            size_t hits = 0;
            uint32_t readpos = 0;
        };

        inline Diagonal BestDiagonal(std::string_view read, std::string_view window, int64_t window_start, size_t k) {
            static thread_local std::vector<std::pair<uint32_t, uint32_t>> window_kmers, read_kmers;
            static thread_local std::vector<std::pair<int64_t, uint32_t>> diagonals;  // offset, read position
            Kmers(window, k, window_kmers);
            std::sort(window_kmers.begin(), window_kmers.end());
            Kmers(read, k, read_kmers);
            diagonals.clear();
            auto by_code = [](std::pair<uint32_t, uint32_t> const& a, std::pair<uint32_t, uint32_t> const& b) { return a.first < b.first; };
            for (auto const& [code, readpos] : read_kmers) {
                auto range = std::equal_range(window_kmers.begin(), window_kmers.end(), std::pair<uint32_t, uint32_t>{ code, 0 }, by_code);
                for (auto it = range.first; it != range.second; ++it) {
                    diagonals.emplace_back(window_start + static_cast<int64_t>(it->second) - static_cast<int64_t>(readpos), readpos);
                }
            }
            std::sort(diagonals.begin(), diagonals.end());
            Diagonal best;
            for (size_t i = 0; i < diagonals.size();) {
                size_t j = i;
                while (j < diagonals.size() && diagonals[j].first == diagonals[i].first) j++;
                if (j - i > best.hits) best = { diagonals[i].first, j - i, diagonals[i].second };
                i = j;
            }
            return best;
        }
    }

    // After the candidates of a fragment are joined and sorted (best first): if the best has only one mate, that mate
    // is sure of it (MAPQ kGuideMinMapq or more against its own other alignments) and the other mate has no candidate
    // of its taxon, the other mate is aligned to that taxon:
    //  1. from its own longest anchor of the taxon (only a read's longest anchors are aligned; a mate on another gene
    //     of the taxon, e.g. a neighbour in an operon, often has one further down);
    //  2. else on the guiding mate's gene, where the fragment can reach (up to kRescueMaxFragment bases from the
    //     guiding mate, on its strand's side, and at least kRescueMinBases of them): the diagonal on which most of its
    //     k-mers match that part of the gene (kRescueMinHits or more) anchors an alignment as any other; a mate that
    //     runs past the gene's end aligns in part, clipped, as other reads at gene ends do. Where the fragment reaches
    //     past the gene's end, the genes the taxon's clade has there (the database's gene neighbours) are looked on
    //     too, over the stretch the fragment reaches across the gap between them.
    // A new alignment on the guiding mate's gene in the opposite orientation, or on a neighbouring gene across the
    // facing ends (across_genes::PairAcrossNeighbours), completes the best candidate; one on another gene is added as
    // a candidate of its own, so that the mates are written on their consensus taxon
    // (ProtalPairedOutputHandler::WriteSplitMates). Returns whether the candidates changed (sort them again).
    template<typename AlignmentHandler>
    bool GuideMate(PairedAlignmentResultList& pairs, AlignmentAnchorList const& anchors1, AlignmentAnchorList const& anchors2,
                   FastxRecord& record1, FastxRecord& record2, std::string const& reverse1, std::string const& reverse2,
                   AlignmentHandler& handler, GenomeLoader& genomes, MateGuidanceCounts& counts) {
        if (pairs.empty()) return false;
        if (pairs.front().first.IsSet() == pairs.front().second.IsSet()) return false;  // both mates aligned (or none)
        bool const guide_is_first = pairs.front().first.IsSet();
        AlignmentResult const guide = guide_is_first ? pairs.front().first : pairs.front().second;

        // The guiding mate's MAPQ against its other alignments; whether the other mate has one of its taxon.
        auto const& guide_info = guide.GetAlignmentInfo();
        int const own = guide_info.Score(2, 3, 1, 2);
        int second = 0;
        bool other_has_taxon = false;
        for (auto const& [a1, a2] : pairs) {
            auto const& g = guide_is_first ? a1 : a2;
            auto const& o = guide_is_first ? a2 : a1;
            if (g.IsSet() && !(g.Taxid() == guide.Taxid() && g.GeneId() == guide.GeneId() && g.Forward() == guide.Forward() &&
                               g.GetAlignmentInfo().gene_alignment_start == guide_info.gene_alignment_start)) {
                second = std::max(second, g.GetAlignmentInfo().Score(2, 3, 1, 2));
            }
            if (o.IsSet() && o.Taxid() == guide.Taxid()) other_has_taxon = true;
        }
        if (own <= 0 || other_has_taxon || MAPQv2(own, std::max(0, second)) < kGuideMinMapq) return false;
        counts.guided++;

        auto& record = guide_is_first ? record2 : record1;
        auto const& reverse = guide_is_first ? reverse2 : reverse1;
        auto const& anchors = guide_is_first ? anchors2 : anchors1;
        static thread_local AlignmentAnchorList one;
        static thread_local AlignmentResultList found;

        auto attach = [&](AlignmentResult&& alignment) {
            bool const orientation = guide_is_first ? CorrectOrientation(guide, alignment) : CorrectOrientation(alignment, guide);
            bool const across = guide_is_first ? across_genes::PairAcrossNeighbours(guide, alignment, genomes, kRescueMaxFragment)
                                               : across_genes::PairAcrossNeighbours(alignment, guide, genomes, kRescueMaxFragment);
            if (alignment.Taxid() == guide.Taxid() && ((alignment.GeneId() == guide.GeneId() && orientation) || across)) {
                (guide_is_first ? pairs.front().second : pairs.front().first) = std::move(alignment);
            } else if (guide_is_first) {
                pairs.emplace_back(AlignmentResult(), std::move(alignment));
            } else {
                pairs.emplace_back(std::move(alignment), AlignmentResult());
            }
        };

        // 1. The other mate's own anchor of the taxon (anchors come longest first).
        auto anchor = std::find_if(anchors.begin(), anchors.end(), [&guide](CAlignmentAnchor const& a) { return a.taxid == guide.Taxid(); });
        if (anchor != anchors.end()) {
            one.assign(1, *anchor);
            found.clear();
            handler(one, found, record.sequence, reverse, 1, record.id);
            if (!found.empty()) {
                counts.from_anchor++;
                attach(std::move(found.front()));
                return true;
            }
        }

        // 2. On the guiding mate's gene, where the fragment can reach: downstream of a forward mate, upstream of a
        // reverse one; the other mate aligns in the other orientation. Where the reach runs past the gene's end,
        // also on the genes the taxon's clade has at that end (the database's gene neighbours, AcrossGenes.h), in
        // the orientation their strands give; the place whose diagonal has the most k-mers is aligned.
        auto& gene = genomes.GetGenome(guide.Taxid()).GetGeneOMP(guide.GeneId());
        int64_t const gene_length = static_cast<int64_t>(gene.GetLength());
        int64_t const start = guide_info.gene_alignment_start;
        int64_t const end = start + static_cast<int64_t>(guide_info.alignment_length);
        int64_t const from = guide.Forward() ? std::max<int64_t>(0, start) : std::max<int64_t>(0, end - kRescueMaxFragment);
        int64_t const to = guide.Forward() ? std::min(gene_length, start + kRescueMaxFragment) : std::min(gene_length, end);
        bool const forward = !guide.Forward();
        static thread_local std::vector<across_genes::Stretch> places, neighbours;
        places.clear();
        if (to - from >= kRescueMinBases) places.push_back({ static_cast<uint32_t>(guide.GeneId()), forward, from, to });
        int64_t const past = guide.Forward() ? start + kRescueMaxFragment - gene_length : kRescueMaxFragment - end;
        if (past >= kRescueMinBases) {
            across_genes::NeighbourStretches(genomes, guide.Taxid(), guide.GeneId(), gene_neighbours::EndAhead(guide.Forward()), forward,
                                             0, past, kRescueMinBases, 2, neighbours);
            places.insert(places.end(), neighbours.begin(), neighbours.end());
        }
        if (places.empty()) return false;
        counts.looked_for++;
        mate_guidance::Diagonal diagonal;
        across_genes::Stretch const* place = nullptr;
        for (auto const& candidate : places) {
            auto& on = genomes.GetGenome(guide.Taxid()).GetGeneOMP(candidate.gene);
            std::string const& query = candidate.forward ? record.sequence : reverse;
            auto const window = on.Window(static_cast<size_t>(candidate.first), static_cast<size_t>(candidate.last));
            auto const d = mate_guidance::BestDiagonal(query, window.View().substr(static_cast<size_t>(candidate.first),
                                                       static_cast<size_t>(candidate.last - candidate.first)), candidate.first, kRescueK);
            if (d.hits > diagonal.hits) {
                diagonal = d;
                place = &candidate;
            }
        }
        if (!place || diagonal.hits < kRescueMinHits) return false;
        CAlignmentAnchor rescue(guide.Taxid(), place->gene, place->forward);
        rescue.chain.emplace_back(ChainLink(static_cast<uint32_t>(diagonal.offset + diagonal.readpos), static_cast<uint16_t>(diagonal.readpos),
                                            static_cast<uint16_t>(kRescueK)));
        rescue.total_length = static_cast<uint16_t>(kRescueK);
        one.assign(1, std::move(rescue));
        found.clear();
        handler(one, found, record.sequence, reverse, 1, record.id);
        if (found.empty()) return false;
        counts.rescued++;
        counts.rescued_on_neighbour += place->gene != guide.GeneId();
        attach(std::move(found.front()));
        return true;
    }
}
