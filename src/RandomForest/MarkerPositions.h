// SPDX-License-Identifier: GPL-2.0-only
#pragma once

// Where the marker genes are in the genomes simulated, so that simulate_metagenomes counts the reads that overlap one
// (--gene_positions): a species whose marker genes got no read cannot be seen by a marker-gene profiler, and an
// evaluation can leave it out of the truth (collect_training_data.py meta_marker_reads, machine_learning_cmdline.py).
//
// The positions are a protal database's gene_positions.tsv (scripts/mini_db/gene_neighbours.py: each marker gene
// placed in each genome of the genome table; columns accession, contig, start, end among others, 1-based, inclusive;
// comment lines start with #). An in-silico strain (insilico_strains.py) is a mutated copy of its representative, base
// for base: its genome is named insilico_<accession with _ for .>, its contigs <genome>_<the representative's contig>,
// and it has its representative's positions.

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace protal::sim {

class MarkerPositions {
public:
    // The rows of the given accessions (every row if none are given), from a plain-text gene_positions.tsv. A read
    // overlaps a marker gene if it shares min_overlap bases with one (all its bases, if it is shorter). Throws if the
    // file cannot be read or lacks a column.
    static std::shared_ptr<MarkerPositions const> Load(std::filesystem::path const& file,
                                                       std::unordered_set<std::string> const& accessions,
                                                       std::uint32_t min_overlap);

    // The accession a genome's positions are under: an in-silico strain's representative's
    // (insilico_GCA_000018565_1: GCA_000018565.1), else the genome's name.
    static std::string AccessionOf(std::string const& genome);

    // Whether the genome has positions (its reads can be counted; a genome without, say one whose genes could not all
    // be placed, cannot).
    bool Has(std::string const& genome) const;

    // Whether bases [start, end) of a contig of the genome overlap a marker gene by min_overlap bases (or by all of
    // them, if fewer).
    bool Overlaps(std::string const& genome, std::string const& contig, std::uint64_t start, std::uint64_t end) const;

    std::uint32_t MinOverlap() const { return m_min_overlap; }
    std::size_t Genes() const { return m_genes; }

private:
    // A contig's genes, by start: [starts[i], ends[i]), and the largest end of genes 0..i (max_ends).
    struct Contig {
        std::vector<std::uint64_t> starts, ends, max_ends;
    };
    std::unordered_map<std::string, Contig> m_contigs;  // accession \t contig -> its genes
    std::unordered_set<std::string> m_accessions;
    std::uint32_t m_min_overlap = 1;
    std::size_t m_genes = 0;
};

}  // namespace protal::sim
