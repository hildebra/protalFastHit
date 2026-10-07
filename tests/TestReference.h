#pragma once

// A reference for the unit tests: taxa with genes 1..n each, written as protal's reference.fna and reference.map into a
// scratch directory of its own, and loaded.

#include <cstdint>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>
#include "SequenceUtils/GenomeLoader.h"
#include "TestUtil.h"

namespace protal::test {

    class LoadedReference {
    public:
        using Taxa = std::map<uint32_t, std::vector<std::string>>;  // taxid -> gene i at [i - 1]

        ScratchDir dir;
        Taxa genes;
        std::unique_ptr<GenomeLoader> loader;

        explicit LoadedReference(Taxa taxa, std::string const& name = "reference") : dir(name), genes(std::move(taxa)) {
            std::string fna, map;
            for (auto const& [taxid, seqs] : genes) {
                for (size_t i = 0; i < seqs.size(); i++) {
                    fna += ">" + std::to_string(taxid) + "_" + std::to_string(i + 1) + "\n";
                    map += std::to_string(taxid) + '\t' + std::to_string(i + 1) + '\t' + std::to_string(fna.size()) + '\t' +
                           std::to_string(fna.size() + seqs[i].size()) + '\n';
                    fna += seqs[i] + '\n';
                }
            }
            loader = std::make_unique<GenomeLoader>(dir.Write("reference.fna", fna), dir.Write("reference.map", map));
            loader->LoadAllGenomes();
        }

        // A gene's sequence.
        std::string const& Gene(uint32_t taxid, size_t gene) const { return genes.at(taxid).at(gene - 1); }

        // A SAM header of the genes (@HD and an @SQ line per gene).
        std::string Header() const {
            std::string header = "@HD\tVN:1.6\n";
            for (auto const& [taxid, seqs] : genes) {
                for (size_t i = 0; i < seqs.size(); i++) {
                    header += "@SQ\tSN:" + std::to_string(taxid) + "_" + std::to_string(i + 1) + "\tLN:" + std::to_string(seqs[i].size()) + "\n";
                }
            }
            return header;
        }

        // Writes a file into the reference's directory; returns its path.
        std::string Write(std::string const& name, std::string const& content) const { return dir.Write(name, content); }
    };

    // Taxon 1 with two genes of 50 bases (genes 1 and 2), loaded.
    struct TinyReference : LoadedReference {
        static constexpr char kGene1[] = "ACGTTGCAAGGCTTACCGATGACTGAAACCGGTTTACGATCGGTAGCATG";
        static constexpr char kGene2[] = "TTGACCAGTCAGGATCCATTGCAGGTACTTGACCGTAAGCTGCATTGACA";
        std::string const gene = kGene1, gene2 = kGene2;

        TinyReference() : LoadedReference({ { 1, { kGene1, kGene2 } } }, "tiny reference") {}
    };

    // Random genes of the given lengths, one draw of rng per base, gene after gene.
    inline std::vector<std::string> RandomGenes(std::vector<size_t> const& lengths, std::mt19937& rng) {
        std::vector<std::string> genes;
        for (auto length : lengths) genes.push_back(RandomSequence(length, rng));
        return genes;
    }
}
