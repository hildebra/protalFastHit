// AcrossGenes.h - following a fragment or a long read from one marker gene into the next, with the database's
// gene neighbours (GeneNeighbours.h; none: nothing here does anything): which genes a species' clade has at a
// gene's end and where they lie past it, and whether mates on two genes of one taxon are one fragment across
// their facing ends. Used by mate guidance (MateGuidance.h), the joining of mates into pairs
// (classify::JoinAlignmentPairs) and long reads (LongReads.h).
#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>
#include "AlignmentUtils.h"
#include "GeneNeighbours.h"
#include "GenomeLoader.h"

namespace protal::across_genes {
    using gene_neighbours::End;

    // The expected neighbours of species taxid at end `end` of `gene` (Table::Partners: their share smoothed over
    // its clades), the most common first, at most `max` of them, among the genes the species has in the database:
    // the rules of the nearest clades that saw them.
    inline void Neighbours(GenomeLoader const& genomes, uint32_t taxid, uint32_t gene, End end, size_t max,
                           std::vector<gene_neighbours::Rule const*>& out) {
        static thread_local std::vector<gene_neighbours::Partner> partners;
        out.clear();
        auto const& table = genomes.GetGeneNeighbours();
        if (table.Empty()) return;
        table.Partners(taxid, gene, end, partners);
        std::stable_sort(partners.begin(), partners.end(), [](auto const& a, auto const& b) { return a.share > b.share; });
        for (auto const& p : partners) {
            if (p.rule->HasPartner() && p.verdict == gene_neighbours::Verdict::Expected && genomes.HasGene(taxid, p.rule->partner)) {
                out.push_back(p.rule);
            }
        }
        if (out.size() > max) out.resize(max);
    }

    // A stretch [first, last) of a gene to look for a read on, in orientation `forward` there.
    struct Stretch {
        uint32_t gene = 0;
        bool forward = true;
        int64_t first = 0;
        int64_t last = 0;
    };

    // The stretches of the neighbours at end `end` of `gene` (Neighbours, at most `max`) that a read in orientation
    // `forward` on the gene, reaching `from` to `to` bases past that end, covers, with the clade's gap range and
    // kGapSlack either side; those shorter than min_bases are left out.
    inline void NeighbourStretches(GenomeLoader& genomes, uint32_t taxid, uint32_t gene, End end, bool forward,
                                   int64_t from, int64_t to, int64_t min_bases, size_t max, std::vector<Stretch>& out) {
        static thread_local std::vector<gene_neighbours::Rule const*> rules;
        out.clear();
        Neighbours(genomes, taxid, gene, end, max, rules);
        for (auto const* rule : rules) {
            int64_t const length = static_cast<int64_t>(genomes.GeneLength(taxid, rule->partner));
            auto const [first, last] = gene_neighbours::PartnerStretch(from, to, rule->gap_min - gene_neighbours::kGapSlack,
                                                                       rule->gap_max + gene_neighbours::kGapSlack, rule->partner_end, length);
            if (last - first < min_bases) continue;
            out.push_back({ rule->partner, gene_neighbours::OrientationOnPartner(forward, end, rule->partner_end), first, last });
        }
    }

    // The bases of the gene an alignment covers (its CIGAR without insertions and clips).
    inline int64_t ReferenceSpan(AlignmentInfo const& info) {
        return static_cast<int64_t>(info.alignment_length) - static_cast<int64_t>(info.insertions);
    }

    // Whether mates aligned to two genes of one taxon (a1 of mate 1, a2 of mate 2) are one fragment across the genes'
    // facing ends: each mate runs towards the end of its gene that faces the other gene, those ends are expected
    // neighbours in the taxon's clades (Table::Assess), and the fragment, from the start of one mate over the
    // shortest gap of the nearest clade that saw them (less kGapSlack) to the start of the other, is at most
    // max_fragment long.
    inline bool PairAcrossNeighbours(AlignmentResult const& a1, AlignmentResult const& a2, GenomeLoader& genomes, int64_t max_fragment) {
        if (a1.Taxid() != a2.Taxid() || a1.GeneId() == a2.GeneId()) return false;
        auto const& table = genomes.GetGeneNeighbours();
        if (table.Empty()) return false;
        End const end1 = gene_neighbours::EndAhead(a1.Forward());
        End const end2 = gene_neighbours::EndAhead(a2.Forward());
        auto const assessment = table.Assess(a1.Taxid(), a1.GeneId(), end1, a2.GeneId(), end2);
        if (assessment.verdict != gene_neighbours::Verdict::Expected) return false;
        auto const& i1 = a1.GetAlignmentInfo();
        auto const& i2 = a2.GetAlignmentInfo();
        int64_t const length1 = static_cast<int64_t>(genomes.GeneLength(a1.Taxid(), a1.GeneId()));
        int64_t const length2 = static_cast<int64_t>(genomes.GeneLength(a2.Taxid(), a2.GeneId()));
        int64_t const reach = gene_neighbours::ReachToEnd(i1.gene_alignment_start, ReferenceSpan(i1), length1, end1) +
                              gene_neighbours::ReachToEnd(i2.gene_alignment_start, ReferenceSpan(i2), length2, end2);
        return reach + assessment.rule->gap_min - gene_neighbours::kGapSlack <= max_fragment;
    }
}
