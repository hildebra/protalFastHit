// SpikeCheck.h - a sample's spike check (SpikeGenes.h) against the database: its records per gene narrowed to those the
// profile counts, every taxon's spike genes by the lengths of its genes, and the sample check. The profiling makes it
// before it profiles a SAM (Profiler::ProfileSam); the alignment makes it once a SAM it wrote is complete, from the
// records counted as they were written (spike_genes::Tally, AlignedChecks), so that the profiling need not read the SAM
// once more for it. Both give the same check from the same records.
#pragma once

#include <cstdint>
#include <map>
#include <utility>
#include <vector>

#include "GenomeLoader.h"
#include "SpikeGenes.h"

namespace protal::spike_genes {
    // The records of `counts` that the profile can count (CheckSam): those on a gene of the database with unique k-mers.
    // Records on a gene without them are skipped as no evidence of the taxon (in a real D6331 sample 12,683 of an absent
    // congener's, on one gene), and on a gene the database lacks rejected, so neither the spike rule nor the sample check
    // sees them.
    inline Counts CountedRecords(GenomeLoader& genomes, Counts counts) {
        std::erase_if(counts, [&](auto const& entry) {
            auto const taxid = static_cast<int>(entry.first >> 32);
            auto const geneid = static_cast<int>(static_cast<uint32_t>(entry.first));
            return !genomes.HasGene(taxid, geneid) || !genomes.GetGenome(taxid).IsGeneHittable(geneid);
        });
        return counts;
    }

    // The spike genes of a sample whose records per gene are `counts` (CountedRecords): each taxon's (SpikesOf) by the
    // lengths of its hittable genes in the database; taxa and genes the database lacks are left alone (their records are
    // rejected). By taxid, then gene.
    inline std::vector<Spike> FindSpikeGenes(GenomeLoader& genomes, Counts const& counts) {
        std::map<uint32_t, std::vector<std::pair<uint32_t, uint64_t>>> by_taxon;
        for (auto const& [key, n] : counts) by_taxon[static_cast<uint32_t>(key >> 32)].emplace_back(static_cast<uint32_t>(key), n);
        std::vector<Spike> spikes;
        auto const& genome_map = genomes.GetGenomeMap();
        for (auto& [taxid, records] : by_taxon) {
            if (records.size() <= kMinRestGenes + 1) continue;
            auto const it = genome_map.find(taxid);
            if (it == genome_map.end()) continue;
            auto& genome = genomes.GetGenome(taxid);
            std::vector<std::pair<uint32_t, uint64_t>> lengths;
            for (auto const gene : genome.GetHittableGenes()) {
                if (!genome.HasGene(gene)) continue;
                lengths.emplace_back(gene, genome.GetGene(gene).GetLength());
            }
            auto found = SpikesOf(taxid, std::move(records), lengths);
            spikes.insert(spikes.end(), found.begin(), found.end());
        }
        return spikes;
    }

    // The spike check of a sample whose records per gene (of its whole SAM: CountRecords) are `counts`.
    inline Found Check(GenomeLoader& genomes, Counts counts) {
        Found found;
        auto const counted = CountedRecords(genomes, std::move(counts));
        found.spikes = FindSpikeGenes(genomes, counted);
        found.check = Judge(counted, found.spikes);
        return found;
    }
}
