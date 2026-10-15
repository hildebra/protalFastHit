#pragma once

// A fragment whose mates did not align together as its best candidate (--mate_guidance, MateGuidanceMode.h). A pair is
// one organism; a mate otherwise counts for a relative (it has no candidate of the right taxon among the few that are
// aligned), or not at all (no anchor: too divergent, or partly beyond the gene's end), and its bases are lost to the
// strain's SNPs.
//  - consensus (VoteMates): both mates vote. Every taxon that either mate aligned to is scored by the sum of the two
//    mates' best scores on it; a mate without an alignment there is looked for as guidance looks for it (RescueMate),
//    and where it is not found scores as an alignment at the identity floor would (at most its best elsewhere). The
//    taxon of the best sum takes the fragment, with the MAPQ of its sum against the next taxon's, so that mates that
//    disagree get a low MAPQ.
//  - guide (GuideMate): the mate that is sure of its alignment guides the other one to its taxon, and on its gene,
//    where the fragment can reach; the other mate's own alignments elsewhere are not weighed.

#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include "AcrossGenes.h"
#include "AlignmentOutputHandler.h"
#include "ChainingStrategy.h"
#include "GenomeLoader.h"
#include "ReadConsensus.h"

namespace protal {
    inline constexpr int kGuideMinMapq = 20;              // a mate this sure of its best alignment guides the other
    inline constexpr int64_t kRescueMaxFragment = 1000;   // the farthest a fragment reaches from its guiding mate
    inline constexpr int64_t kRescueMinBases = 30;        // the rest of the gene the fragment can reach, at least
    inline constexpr size_t kRescueK = 12;                // k-mers that place the other mate on the gene
    inline constexpr size_t kRescueMinHits = 2;           // of them on one diagonal, at least

    // The consensus vote (VoteMates).
    inline constexpr int kVoteUnaligned = 0;      // a mate's score where it aligned nowhere, as in a pair's Bitscore
    inline constexpr size_t kVoteMaxRescues = 4;  // mates looked for per fragment, on the taxa of the highest scores first
    // A mate's own alignment of another taxon is clearly better than its alignment of the vote's taxon at this MAPQ
    // between the two, the one at which guidance calls a mate sure: ~12 mismatches more on 150 bases, or ~60 bases
    // clipped at a gene's end. Such a mate (a chimera, a mate past its gene's end that another taxon's neighbouring gene
    // has) keeps its own alignment, with MAPQ 0, instead of being moved to the vote's taxon; fewer differences are the
    // strain's, and the pair, one organism, outweighs them.
    inline constexpr int kVoteClearlyBetterMapq = kGuideMinMapq;
    inline constexpr int kVoteLowMapq = 4;  // the profiler's MAPQ cutoff (Profiler.h m_min_mapq), for the log

    // What mate guidance did, over the fragments of a run.
    struct MateGuidanceCounts {
        size_t guided = 0;    // fragments with a guiding mate and another without a candidate of its taxon
        size_t from_anchor = 0;  // other mates aligned from their own anchor of that taxon
        size_t looked_for = 0;   // other mates looked for on the guiding mate's gene
        size_t rescued = 0;      // of these, found and aligned there
        size_t rescued_on_neighbour = 0;  // of those, on a gene next to the guiding mate's (the database's gene neighbours)
        size_t paired_across_genes = 0;   // fragments whose best candidate has its mates on two neighbouring genes
        // The consensus vote (VoteMates).
        size_t voted = 0;          // fragments voted on: their best candidate had one mate only
        size_t vote_disagreed = 0;  // of them, with each mate's own best alignment on another taxon
        size_t vote_moved = 0;      // of them, given another taxon than their best candidate's
        size_t vote_looked_for = 0;  // mates looked for on a taxon they had no alignment to
        size_t vote_found = 0;      // of them, found
        size_t vote_kept = 0;       // mates that kept a clearly better alignment of their own (MAPQ 0)
        size_t vote_uncertain = 0;  // fragments whose vote has a MAPQ under kVoteLowMapq

        void Join(MateGuidanceCounts const& other) {
            guided += other.guided;
            from_anchor += other.from_anchor;
            looked_for += other.looked_for;
            rescued += other.rescued;
            rescued_on_neighbour += other.rescued_on_neighbour;
            paired_across_genes += other.paired_across_genes;
            voted += other.voted;
            vote_disagreed += other.vote_disagreed;
            vote_moved += other.vote_moved;
            vote_looked_for += other.vote_looked_for;
            vote_found += other.vote_found;
            vote_kept += other.vote_kept;
            vote_uncertain += other.vote_uncertain;
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

        // Every pair of a k-mer of the read and an equal one of the window gives a diagonal; the k-mers of the shorter
        // of the two are sorted (and their codes marked in a small filter) and those of the other looked up in them.
        // Which side is sorted changes nothing: the diagonals are the same, sorted before they are counted.
        inline Diagonal BestDiagonal(std::string_view read, std::string_view window, int64_t window_start, size_t k) {
            static thread_local std::vector<std::pair<uint32_t, uint32_t>> sorted_kmers, other_kmers;
            static thread_local std::vector<std::pair<int64_t, uint32_t>> diagonals;  // offset, read position
            bool const read_sorted = read.size() <= window.size();
            Kmers(read_sorted ? read : window, k, sorted_kmers);
            std::sort(sorted_kmers.begin(), sorted_kmers.end());
            constexpr uint32_t kFilterBits = 4096;
            auto const slot = [](uint32_t code) { return (code ^ (code >> 11) ^ (code >> 22)) & (kFilterBits - 1); };
            std::array<uint64_t, kFilterBits / 64> filter{};
            for (auto const& kmer : sorted_kmers) filter[slot(kmer.first) >> 6] |= uint64_t{1} << (slot(kmer.first) & 63);
            Kmers(read_sorted ? window : read, k, other_kmers);
            diagonals.clear();
            auto by_code = [](std::pair<uint32_t, uint32_t> const& a, std::pair<uint32_t, uint32_t> const& b) { return a.first < b.first; };
            for (auto const& [code, pos] : other_kmers) {
                uint32_t const s = slot(code);
                if (!((filter[s >> 6] >> (s & 63)) & 1)) continue;
                auto range = std::equal_range(sorted_kmers.begin(), sorted_kmers.end(), std::pair<uint32_t, uint32_t>{ code, 0 }, by_code);
                for (auto it = range.first; it != range.second; ++it) {
                    uint32_t const readpos = read_sorted ? it->second : pos;
                    uint32_t const windowpos = read_sorted ? pos : it->second;
                    diagonals.emplace_back(window_start + static_cast<int64_t>(windowpos) - static_cast<int64_t>(readpos), readpos);
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

    // The other mate's alignment to the taxon of `guide` (a mate's alignment), looked for as mate guidance does:
    //  1. from its own longest anchor of the taxon (only a read's longest anchors are aligned; a mate on another gene
    //     of the taxon, e.g. a neighbour in an operon, often has one further down);
    //  2. else on the guiding mate's gene, where the fragment can reach (up to kRescueMaxFragment bases from the
    //     guiding mate, on its strand's side, and at least kRescueMinBases of them): the diagonal on which most of its
    //     k-mers match that part of the gene (kRescueMinHits or more) anchors an alignment as any other; a mate that
    //     runs past the gene's end aligns in part, clipped, as other reads at gene ends do. Where the fragment reaches
    //     past the gene's end, the genes the taxon's clade has there (the database's gene neighbours) are looked on
    //     too, over the stretch the fragment reaches across the gap between them.
    // Nothing if neither aligns (within the alignment handler's identity floor). `anchors`, `record` and `reverse` are
    // the other mate's; `guide` must stay valid while it runs.
    template<typename AlignmentHandler>
    std::optional<AlignmentResult> RescueMate(AlignmentResult const& guide, AlignmentAnchorList const& anchors, FastxRecord& record,
                                              std::string const& reverse, AlignmentHandler& handler, GenomeLoader& genomes,
                                              MateGuidanceCounts& counts) {
        static thread_local AlignmentAnchorList one;
        static thread_local AlignmentResultList found;
        auto const& guide_info = guide.GetAlignmentInfo();

        // 1. The other mate's own anchor of the taxon (anchors come longest first).
        auto anchor = std::find_if(anchors.begin(), anchors.end(), [&guide](CAlignmentAnchor const& a) { return a.taxid == guide.Taxid(); });
        if (anchor != anchors.end()) {
            one.assign(1, *anchor);
            found.clear();
            handler(one, found, record.sequence, reverse, 1, record.id);
            if (!found.empty()) {
                counts.from_anchor++;
                return std::move(found.front());
            }
        }

        // 2. On the guiding mate's gene, where the fragment can reach: downstream of a forward mate, upstream of a
        // reverse one; the other mate aligns in the other orientation. Where the reach runs past the gene's end,
        // also on the genes the taxon's clade has at that end (the database's gene neighbours, AcrossGenes.h), in
        // the orientation their strands give; the place whose diagonal has the most k-mers is aligned.
        auto& gene = genomes.GetGeneOMP(guide.Taxid(), guide.GeneId());
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
        if (places.empty()) return std::nullopt;
        counts.looked_for++;
        mate_guidance::Diagonal diagonal;
        across_genes::Stretch const* place = nullptr;
        for (auto const& candidate : places) {
            auto& on = genomes.GetGeneOMP(guide.Taxid(), candidate.gene);
            std::string const& query = candidate.forward ? record.sequence : reverse;
            auto const window = on.Window(static_cast<size_t>(candidate.first), static_cast<size_t>(candidate.last));
            auto const d = mate_guidance::BestDiagonal(query, window.View().substr(static_cast<size_t>(candidate.first),
                                                       static_cast<size_t>(candidate.last - candidate.first)), candidate.first, kRescueK);
            if (d.hits > diagonal.hits) {
                diagonal = d;
                place = &candidate;
            }
        }
        if (!place || diagonal.hits < kRescueMinHits) return std::nullopt;
        CAlignmentAnchor rescue(guide.Taxid(), place->gene, place->forward);
        rescue.chain.emplace_back(ChainLink(static_cast<uint32_t>(diagonal.offset + diagonal.readpos), static_cast<uint16_t>(diagonal.readpos),
                                            static_cast<uint16_t>(kRescueK)));
        rescue.total_length = static_cast<uint16_t>(kRescueK);
        one.assign(1, std::move(rescue));
        found.clear();
        handler(one, found, record.sequence, reverse, 1, record.id);
        if (found.empty()) return std::nullopt;
        counts.rescued++;
        counts.rescued_on_neighbour += place->gene != guide.GeneId();
        return std::move(found.front());
    }

    // Whether a mate 1 and a mate 2 alignment are one fragment, as classify::JoinAlignmentPairs pairs them: on one
    // gene of a taxon in opposite orientations, or on two genes of it across their facing ends.
    inline bool MatesAreOneFragment(AlignmentResult const& first, AlignmentResult const& second, GenomeLoader& genomes) {
        if (first.Taxid() != second.Taxid()) return false;
        if (first.GeneId() == second.GeneId()) return CorrectOrientation(first, second);
        return across_genes::PairAcrossNeighbours(first, second, genomes, kRescueMaxFragment);
    }

    // --mate_guidance guide. After the candidates of a fragment are joined and sorted (best first): if the best has
    // only one mate, that mate is sure of it (MAPQ kGuideMinMapq or more against its own other alignments) and the
    // other mate has no candidate of its taxon, the other mate is aligned to that taxon (RescueMate). A new alignment
    // on the guiding mate's gene in the opposite orientation, or on a neighbouring gene across the facing ends
    // (across_genes::PairAcrossNeighbours), completes the best candidate; one on another gene is added as a candidate
    // of its own, so that the mates are written on their consensus taxon (ProtalPairedOutputHandler::WriteSplitMates).
    // Returns whether the candidates changed (sort them again).
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

        auto alignment = guide_is_first ? RescueMate(guide, anchors2, record2, reverse2, handler, genomes, counts)
                                        : RescueMate(guide, anchors1, record1, reverse1, handler, genomes, counts);
        if (!alignment) return false;
        bool const orientation = guide_is_first ? CorrectOrientation(guide, *alignment) : CorrectOrientation(*alignment, guide);
        bool const across = guide_is_first ? across_genes::PairAcrossNeighbours(guide, *alignment, genomes, kRescueMaxFragment)
                                           : across_genes::PairAcrossNeighbours(*alignment, guide, genomes, kRescueMaxFragment);
        if (alignment->Taxid() == guide.Taxid() && ((alignment->GeneId() == guide.GeneId() && orientation) || across)) {
            (guide_is_first ? pairs.front().second : pairs.front().first) = std::move(*alignment);
        } else if (guide_is_first) {
            pairs.emplace_back(AlignmentResult(), std::move(*alignment));
        } else {
            pairs.emplace_back(std::move(*alignment), AlignmentResult());
        }
        return true;
    }

    // --mate_guidance consensus. After the candidates of a fragment are joined and sorted (best first), if the best has
    // one mate only (the mates did not align together as their best: one did not align, they prefer two taxa, or they
    // lie on two genes that are no neighbours), the two mates vote:
    //  - Every taxon of the candidates is scored by the sum of the two mates' best scores on it (Score(2, 3, 1, 2), as a
    //    pair's Bitscore). A mate without an alignment of the taxon is looked for there (RescueMate, from the other
    //    mate's best alignment of it): on the taxa of the highest known scores first, at most kVoteMaxRescues of them,
    //    and only where finding it could change the best or the next sum (its read aligned perfectly).
    //  - A mate not found on a taxon, but aligned to another, scores what an alignment of its read at the identity
    //    floor (-a) would, 1.5 a base at 0.9, or its best alignment's score if that is lower: the database may lack
    //    the part of the taxon's genome it is from (a gene, or a gene's neighbour, missing from the reference) or the
    //    strain may differ there, so its absence counts as the poorest alignment protal accepts, never as more than the
    //    mate achieved anywhere. Scored 0 instead (as in a pair's Bitscore), a congener whose two mates both align, if
    //    poorly, beat the right taxon of a mate aligned perfectly whose other mate the reference lacks (2026-10-15
    //    evaluation). A mate aligned nowhere scores kVoteUnaligned (0) everywhere.
    //  - The taxon of the highest sum wins (of equal sums, the one of the better candidate). The vote's MAPQ is that of
    //    the winner's sum against the next taxon's (MAPQv2; against 0 without another taxon), over the mates aligned to
    //    either of them (a mate aligned to neither scores the same on both): mates that disagree get a low MAPQ, and a
    //    mate found on the winner does not double the pair's score against its competitors. As for a pair on one gene,
    //    a mate that fits both taxa equally well is part of the sums, and lowers the MAPQ of the other mate's margin.
    //  - A mate's MAPQ is also at most its margin within the winner: its best alignment there against its next one at
    //    another place of the taxon (a paralog, another gene), as MAPQv1 bounded it over the candidates (2026-10-15).
    //  - Each mate is written at its best alignment of the winner, its own or one found, with the vote's MAPQ; unless
    //    its own best alignment, of another taxon, is clearly better (MAPQ kVoteClearlyBetterMapq or more between the
    //    two): then it keeps that one with MAPQ 0, as does a mate without an alignment of the winner, so that it counts
    //    for no taxon. A mate found on the winner completes the other mate's candidate where the two are one fragment
    //    (MatesAreOneFragment), else it is a candidate of its own.
    // The candidate to write first (a mate's on the winner) is moved to the front, the others keep their order, and
    // `vote` says how to write them: one candidate (a pair, or a mate alone), or two mates on two candidates
    // (ProtalPairedOutputHandler::WriteSplitMates). Returns whether the vote ran. Thread- and order-independent: it
    // depends on the fragment's candidates alone.
    template<typename AlignmentHandler>
    bool VoteMates(PairedAlignmentResultList& pairs, AlignmentAnchorList const& anchors1, AlignmentAnchorList const& anchors2,
                   FastxRecord& record1, FastxRecord& record2, std::string const& reverse1, std::string const& reverse2,
                   AlignmentHandler& handler, GenomeLoader& genomes, MateGuidanceCounts& counts, PairVote& vote,
                   double min_identity = 0.9) {
        vote = PairVote{};
        if (pairs.empty() || (pairs.front().first.IsSet() && pairs.front().second.IsSet())) return false;
        constexpr int kNone = INT_MIN;
        struct Taxon {
            uint32_t taxid = 0;
            std::array<int, 2> score{ kNone, kNone };            // each mate's best score on the taxon
            std::array<size_t, 2> index{ SIZE_MAX, SIZE_MAX };   // the candidate holding it (SIZE_MAX: found by the vote)
            std::array<int, 2> next{ kNone, kNone };             // each mate's next score on the taxon, at another place
            std::array<std::pair<uint32_t, int32_t>, 2> place{};  // the gene and position of its best
            // The sum of the mates' scores, a mate without an alignment of the taxon scoring `missing`.
            int Sum(std::array<int, 2> const& missing = { kVoteUnaligned, kVoteUnaligned }) const {
                return (score[0] == kNone ? missing[0] : score[0]) + (score[1] == kNone ? missing[1] : score[1]);
            }
            int Of(size_t m, std::array<int, 2> const& missing) const { return score[m] == kNone ? missing[m] : score[m]; }
        };
        static thread_local std::vector<Taxon> taxa;
        static thread_local std::vector<size_t> order;
        taxa.clear();
        auto mate = [&pairs](size_t i, size_t m) -> AlignmentResult& { return m == 0 ? pairs[i].first : pairs[i].second; };
        auto vote_score = [](AlignmentResult const& ar) { return ar.GetAlignmentInfo().Score(2, 3, 1, 2); };

        // The taxa in the order of their best candidates, and each mate's best alignment of all.
        std::array<int, 2> own_score{ kNone, kNone };
        std::array<size_t, 2> own_index{ SIZE_MAX, SIZE_MAX };
        for (size_t i = 0; i < pairs.size(); i++) {
            for (size_t m = 0; m < 2; m++) {
                auto const& ar = mate(i, m);
                if (!ar.IsSet()) continue;
                auto const taxid = static_cast<uint32_t>(ar.Taxid());
                auto it = std::find_if(taxa.begin(), taxa.end(), [taxid](Taxon const& t) { return t.taxid == taxid; });
                if (it == taxa.end()) it = taxa.insert(taxa.end(), Taxon{ taxid });
                int const score = vote_score(ar);
                std::pair<uint32_t, int32_t> const place{ ar.GeneId(), ar.GenePos() };
                if (score > it->score[m]) {
                    if (it->score[m] != kNone && it->place[m] != place) it->next[m] = it->score[m];
                    it->score[m] = score;
                    it->index[m] = i;
                    it->place[m] = place;
                } else if (place != it->place[m] && score > it->next[m]) {
                    it->next[m] = score;
                }
                if (score > own_score[m]) {
                    own_score[m] = score;
                    own_index[m] = i;
                }
            }
        }
        counts.voted++;
        counts.vote_disagreed += own_index[0] != SIZE_MAX && own_index[1] != SIZE_MAX &&
                                 mate(own_index[0], 0).Taxid() != mate(own_index[1], 1).Taxid();

        // The mates looked for: on taxa that have one mate's alignment, the highest first, while finding the other
        // mate (at most a perfect alignment of its read) could still raise the taxon above the second-highest sum.
        order.clear();
        for (size_t t = 0; t < taxa.size(); t++) {
            if ((taxa[t].score[0] == kNone) != (taxa[t].score[1] == kNone)) order.push_back(t);
        }
        auto known = [&](size_t t) { return std::max(taxa[t].score[0], taxa[t].score[1]); };
        std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return known(a) > known(b); });
        auto second_sum = [&]() {
            int first = INT_MIN, second = INT_MIN;
            for (auto const& t : taxa) {
                int const s = t.Sum();
                if (s > first) {
                    second = first;
                    first = s;
                } else if (s > second) {
                    second = s;
                }
            }
            return second;
        };
        std::optional<AlignmentResult> found;  // the winner's, if a mate was found on it
        size_t found_taxon = SIZE_MAX;
        static thread_local std::vector<std::optional<AlignmentResult>> found_on;
        found_on.assign(taxa.size(), std::nullopt);
        size_t tried = 0;
        for (size_t t : order) {
            if (tried == kVoteMaxRescues) break;
            size_t const m = taxa[t].score[0] == kNone ? 0 : 1;  // the mate without an alignment of the taxon
            size_t const k = 1 - m;
            auto& record = m == 0 ? record1 : record2;
            int const perfect = 2 * static_cast<int>(record.sequence.size());
            if (static_cast<long>(taxa[t].score[k]) + perfect <= second_sum()) continue;
            tried++;
            counts.vote_looked_for++;
            auto alignment = RescueMate(mate(taxa[t].index[k], k), m == 0 ? anchors1 : anchors2, record, m == 0 ? reverse1 : reverse2,
                                        handler, genomes, counts);
            if (!alignment || alignment->Taxid() != taxa[t].taxid) continue;
            counts.vote_found++;
            taxa[t].score[m] = vote_score(*alignment);
            found_on[t] = std::move(alignment);
        }

        // What a mate scores on a taxon it has no alignment to: an alignment of its read at the identity floor (2 a
        // match, -3 a mismatch: 5 * floor - 3 a base), at most its best alignment's score; 0 if it aligned nowhere.
        std::array<int, 2> missing{ kVoteUnaligned, kVoteUnaligned };
        for (size_t m = 0; m < 2; m++) {
            int best_any = kNone;
            for (auto const& t : taxa) best_any = std::max(best_any, t.score[m]);
            if (best_any == kNone) continue;
            auto const& record = m == 0 ? record1 : record2;
            int const floor = static_cast<int>((5.0 * min_identity - 3.0) * static_cast<double>(record.sequence.size()));
            missing[m] = std::min(floor, best_any);
        }
        // The winner, and its MAPQ against the next taxon, over the mates aligned to either: a mate aligned to neither
        // scores the same on both and is no evidence between them.
        size_t win = 0;
        for (size_t t = 1; t < taxa.size(); t++) {
            if (taxa[t].Sum(missing) > taxa[win].Sum(missing)) win = t;
        }
        size_t runner = SIZE_MAX;
        for (size_t t = 0; t < taxa.size(); t++) {
            if (t != win && (runner == SIZE_MAX || taxa[t].Sum(missing) > taxa[runner].Sum(missing))) runner = t;
        }
        int best = 0, next = 0;
        for (size_t m = 0; m < 2; m++) {
            if (taxa[win].score[m] == kNone && (runner == SIZE_MAX || taxa[runner].score[m] == kNone)) continue;
            best += taxa[win].Of(m, missing);
            if (runner != SIZE_MAX) next += taxa[runner].Of(m, missing);
        }
        int const mapq = best > 0 ? MAPQv2(best, std::clamp(next, 0, best)) : 0;
        Taxon const winner = taxa[win];
        if (found_on[win]) {
            found = std::move(found_on[win]);
            found_taxon = win;
        }
        auto const& front = pairs.front().first.IsSet() ? pairs.front().first : pairs.front().second;
        counts.vote_moved += front.Taxid() != winner.taxid;
        counts.vote_uncertain += mapq < kVoteLowMapq;

        // Each mate's record: its best of the winner, unless its own of another taxon is clearly better.
        std::array<size_t, 2> index{ SIZE_MAX, SIZE_MAX };
        std::array<int, 2> mate_mapq{ 0, 0 };
        std::array<bool, 2> on_winner{ false, false };
        for (size_t m = 0; m < 2; m++) {
            int const on = winner.score[m];
            bool const elsewhere = own_index[m] != SIZE_MAX && mate(own_index[m], m).Taxid() != winner.taxid;
            bool const clearly = on != kNone && elsewhere && own_score[m] > 0 && own_score[m] > on &&
                                 MAPQv2(own_score[m], std::max(0, on)) >= kVoteClearlyBetterMapq;
            if (on != kNone && !clearly) {
                on_winner[m] = true;
                // At most its margin within the winner: a read that fits two places of the taxon equally well (a
                // paralog, another gene) is no surer of its place than MAPQv1 made it before the vote.
                mate_mapq[m] = winner.next[m] == kNone || on <= 0 ? mapq
                                                                 : std::min(mapq, MAPQv2(on, std::clamp(winner.next[m], 0, on)));
                index[m] = winner.index[m];  // SIZE_MAX: the mate found (below)
            } else if (own_index[m] != SIZE_MAX) {
                index[m] = own_index[m];
                counts.vote_kept++;
            }
        }
        // A mate found on the winner: into the other mate's candidate of the winner if they are one fragment.
        for (size_t m = 0; m < 2; m++) {
            if (!on_winner[m] || index[m] != SIZE_MAX) continue;
            if (!found || found_taxon != win) {  // cannot happen: a winner's score without a candidate is a found mate's
                on_winner[m] = false;
                mate_mapq[m] = 0;
                index[m] = own_index[m];
                continue;
            }
            size_t const k = 1 - m;
            bool const joins = on_winner[k] && index[k] != SIZE_MAX && !mate(index[k], m).IsSet() &&
                               (m == 0 ? MatesAreOneFragment(*found, mate(index[k], k), genomes)
                                       : MatesAreOneFragment(mate(index[k], k), *found, genomes));
            if (joins) {
                mate(index[k], m) = std::move(*found);
                index[m] = index[k];
            } else {
                if (m == 0) {
                    pairs.emplace_back(std::move(*found), AlignmentResult());
                } else {
                    pairs.emplace_back(AlignmentResult(), std::move(*found));
                }
                index[m] = pairs.size() - 1;
            }
        }
        if (index[0] == SIZE_MAX && index[1] == SIZE_MAX) return true;  // nothing to write (cannot happen): as without

        // The candidate written first: a mate's of the winner (mate 1's first), moved to the front.
        size_t const first = on_winner[0] && index[0] != SIZE_MAX ? index[0] : on_winner[1] && index[1] != SIZE_MAX ? index[1]
                           : index[0] != SIZE_MAX ? index[0] : index[1];
        std::rotate(pairs.begin(), pairs.begin() + static_cast<std::ptrdiff_t>(first), pairs.begin() + static_cast<std::ptrdiff_t>(first) + 1);
        auto moved = [first](size_t i) { return i == SIZE_MAX ? i : i == first ? size_t{0} : i < first ? i + 1 : i; };
        vote.decided = true;
        vote.taxid = winner.taxid;
        vote.mapq = mapq;
        vote.index = { moved(index[0]), moved(index[1]) };
        vote.mate_mapq = mate_mapq;
        vote.split = vote.index[0] != SIZE_MAX && vote.index[1] != SIZE_MAX && vote.index[0] != vote.index[1];
        return true;
    }
}
