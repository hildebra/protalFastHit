// Composition.h - what the called species explain of a sample. A species' depth (Taxon::VerticalCoverage: its fragment
// bases per base of its marker genes) times its genome size (species_priors.tsv) is the bases its whole genome put into
// the sample; their sum over the called species, against the bases the aligner read (ScannedReads.h, from the SAM
// header; a pair's overlapping bases counted once, as in the depth), is the share of the sample they explain. The rest
// came from organisms the profile does not name: species the database lacks or the model did not call, plasmids,
// viruses, eukaryotes and host. Their genome equivalents, at the called species' average genome size, give the
// profile's unknown share ("?"), so that the abundances are shares of all genomes sequenced; divided by the depth of a
// typical called species (the median) and of the least abundant one, they give how many species they would be.
//
// The average genome size is the mean of the called species' sizes weighted by their depth, which is proportional to
// their cells: the genome of the average cell (as MicrobeCensus estimates it from single-copy genes, here from the
// species). Its spread is given as the depth-weighted standard deviation and quantiles.
#pragma once

#include "ScannedReads.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <optional>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace protal::composition {
    inline constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    // The depth-weighted quantiles of the genome sizes Result gives.
    inline constexpr double kQuantiles[5] = { 0.1, 0.25, 0.5, 0.75, 0.9 };

    // A called species: its depth and its genome size in bp (0 or less: unknown).
    struct Species {
        uint32_t taxid = 0;
        double depth = 0;
        double genome_size = -1;
    };

    struct Result {
        std::optional<ScannedReads> scanned;  // what the SAM header says
        double fragment_base_share = 1;       // fragment bases per read base: below 1 by the bases a pair's mates share
        size_t species = 0;                   // called
        size_t sized = 0;                     // of them, with a genome size
        double scanned_fragment_bases = kNaN; // the scanned bases times fragment_base_share
        double attributed_bases = kNaN;       // the called species' depth times genome size, summed
        double attributed_share = kNaN;       // of the scanned fragment bases; above 1 where depths or sizes are overestimated
        double attributed_fragments = kNaN;   // the scanned fragments times attributed_share
        double average_genome_size = kNaN;    // bp, the called species' sizes weighted by depth
        double genome_size_sd = kNaN;         // depth-weighted
        double genome_size_quantiles[5] = { kNaN, kNaN, kNaN, kNaN, kNaN };  // at kQuantiles, depth-weighted
        double genome_equivalents = 0;        // the called species' depths, summed: genomes sequenced at 1x
        double unknown_genome_equivalents = kNaN;  // of the bases not explained, at the average genome size
        double unknown_share = kNaN;          // of all genome equivalents: the profile's "?"
        double median_depth = kNaN;           // of the called species
        double lowest_depth = kNaN;
        double missing_species_at_median_depth = kNaN;  // the unknown genome equivalents over the median depth
        double missing_species_at_lowest_depth = kNaN;  // over the lowest depth

        // Whether the profile gets its "?" line: the reads scanned and the genome sizes are known.
        bool HasUnknownShare() const { return !std::isnan(unknown_share); }

        // The fragments the depth of a taxon implies over its genome (the average genome size for one without a size);
        // NaN without the scanned reads or a size.
        double FragmentsOf(double depth, double genome_size) const {
            double const size = genome_size > 0 ? genome_size : average_genome_size;
            if (!scanned || !(scanned_fragment_bases > 0) || std::isnan(size)) return kNaN;
            return depth * size / scanned_fragment_bases * static_cast<double>(scanned->fragments);
        }
    };

    // The smallest size whose cumulative weight reaches q of the total, of (size, weight) pairs sorted by size.
    inline double WeightedQuantile(std::vector<std::pair<double, double>> const& sorted, double total, double q) {
        if (sorted.empty() || !(total > 0)) return kNaN;
        double cumulative = 0;
        for (auto const& [size, weight] : sorted) {
            cumulative += weight;
            if (cumulative >= q * total * (1 - 1e-12)) return size;
        }
        return sorted.back().first;
    }

    // The composition of a sample of the called `species`, from its `scanned` reads (nullopt: not known) and the
    // fragment bases per read base of its aligned reads. `database_has_sizes`: the database gives genome sizes at all
    // (without, a sample without called species gets no unknown share either, as the others).
    inline Result Compute(std::vector<Species> const& species, std::optional<ScannedReads> const& scanned,
                          double fragment_base_share, bool database_has_sizes) {
        Result r;
        r.scanned = scanned;
        r.fragment_base_share = fragment_base_share > 0 && fragment_base_share <= 1 ? fragment_base_share : 1;
        r.species = species.size();

        std::vector<std::pair<double, double>> sized;  // (size, depth)
        std::vector<double> depths;
        double sized_depth = 0, sized_bases = 0;
        for (auto const& s : species) {
            if (!(s.depth > 0)) continue;
            r.genome_equivalents += s.depth;
            depths.push_back(s.depth);
            if (s.genome_size > 0) {
                sized.emplace_back(s.genome_size, s.depth);
                sized_depth += s.depth;
                sized_bases += s.depth * s.genome_size;
            }
        }
        for (auto const& s : species) r.sized += s.genome_size > 0;
        if (!depths.empty()) {
            std::sort(depths.begin(), depths.end());
            size_t const mid = depths.size() / 2;
            r.median_depth = depths.size() % 2 ? depths[mid] : (depths[mid - 1] + depths[mid]) / 2;
            r.lowest_depth = depths.front();
        }
        if (sized_depth > 0) {
            std::sort(sized.begin(), sized.end());
            r.average_genome_size = sized_bases / sized_depth;
            double variance = 0;
            for (auto const& [size, depth] : sized) variance += depth * (size - r.average_genome_size) * (size - r.average_genome_size);
            r.genome_size_sd = std::sqrt(variance / sized_depth);
            for (size_t i = 0; i < 5; i++) r.genome_size_quantiles[i] = WeightedQuantile(sized, sized_depth, kQuantiles[i]);
            // A species without a size counts at the average one.
            r.attributed_bases = sized_bases + (r.genome_equivalents - sized_depth) * r.average_genome_size;
        } else if (depths.empty()) {
            r.attributed_bases = 0;
        }

        if (!scanned || scanned->bases == 0) return r;
        r.scanned_fragment_bases = static_cast<double>(scanned->bases) * r.fragment_base_share;
        if (std::isnan(r.attributed_bases)) return r;  // called species, none with a size
        r.attributed_share = r.attributed_bases / r.scanned_fragment_bases;
        r.attributed_fragments = r.attributed_share * static_cast<double>(scanned->fragments);
        if (!database_has_sizes) return r;
        if (depths.empty()) {
            r.unknown_share = 1;  // nothing called: every genome sequenced is unknown
            return r;
        }
        double const unexplained = std::max(0.0, r.scanned_fragment_bases - r.attributed_bases);
        r.unknown_genome_equivalents = unexplained / r.average_genome_size;
        r.unknown_share = r.unknown_genome_equivalents / (r.genome_equivalents + r.unknown_genome_equivalents);
        r.missing_species_at_median_depth = r.unknown_genome_equivalents / r.median_depth;
        r.missing_species_at_lowest_depth = r.unknown_genome_equivalents / r.lowest_depth;
        return r;
    }

    // The columns of <profile>.composition, one line per sample.
    inline void WriteHeader(std::ostream& os) {
        os << "Sample\tScannedFragments\tScannedReads\tScannedBases\tFragmentBaseShare\tCalledSpecies\tSpeciesWithGenomeSize\t"
              "AttributedShare\tAttributedFragments\tAverageGenomeSize\tGenomeSizeSD\tGenomeSizeP10\tGenomeSizeP25\tGenomeSizeMedian\t"
              "GenomeSizeP75\tGenomeSizeP90\tGenomeEquivalents\tUnknownGenomeEquivalents\tUnknownShare\tMedianDepth\tLowestDepth\t"
              "MissingSpeciesAtMedianDepth\tMissingSpeciesAtLowestDepth\n";
    }

    // A number of the table: NA if it is not known.
    inline std::string Value(double v, int precision = 6) {
        if (std::isnan(v) || std::isinf(v)) return "NA";
        std::ostringstream os;
        os << std::setprecision(precision) << v;
        return os.str();
    }

    inline void WriteRow(std::ostream& os, std::string const& sample, Result const& r) {
        auto const count = [&r](uint64_t ScannedReads::*field) { return r.scanned ? std::to_string((*r.scanned).*field) : std::string("NA"); };
        os << sample << '\t' << count(&ScannedReads::fragments) << '\t' << count(&ScannedReads::reads) << '\t'
           << count(&ScannedReads::bases) << '\t' << Value(r.fragment_base_share) << '\t' << r.species << '\t' << r.sized << '\t'
           << Value(r.attributed_share) << '\t' << Value(std::round(r.attributed_fragments), 15) << '\t'
           << Value(std::round(r.average_genome_size), 15) << '\t' << Value(std::round(r.genome_size_sd), 15);
        for (double q : r.genome_size_quantiles) os << '\t' << Value(std::round(q), 15);
        os << '\t' << Value(r.genome_equivalents) << '\t' << Value(r.unknown_genome_equivalents) << '\t' << Value(r.unknown_share)
           << '\t' << Value(r.median_depth) << '\t' << Value(r.lowest_depth) << '\t' << Value(r.missing_species_at_median_depth)
           << '\t' << Value(r.missing_species_at_lowest_depth) << '\n';
    }

    // One line for the run's log, after "Sample <name>: ".
    inline std::string Summary(Result const& r) {
        std::ostringstream os;
        os << std::fixed;
        if (!r.scanned) {
            os << "its SAM does not say how many reads were scanned (written by protal before 2026-10-08, or with "
                  "--full_sam_header): no composition";
            return os.str();
        }
        if (std::isnan(r.attributed_share)) {
            os << r.scanned->fragments << " fragments scanned; " << r.species
               << " species called, none with a genome size in the database: no composition";
            return os.str();
        }
        os << std::setprecision(1) << 100 * r.attributed_share << "% of " << r.scanned->fragments << " fragments explained by the "
           << r.species << " species called";
        if (!std::isnan(r.average_genome_size)) {
            os << "; average genome size " << std::setprecision(2) << r.average_genome_size / 1e6 << " Mb (10-90%: "
               << r.genome_size_quantiles[0] / 1e6 << "-" << r.genome_size_quantiles[4] / 1e6 << " Mb)";
        }
        if (r.HasUnknownShare()) {
            os << "; unknown " << std::setprecision(1) << 100 * r.unknown_share << "% of genomes";
            if (!std::isnan(r.missing_species_at_median_depth)) {
                os << ", about " << std::setprecision(0) << r.missing_species_at_median_depth << " species at the median depth ("
                   << r.missing_species_at_lowest_depth << " at the lowest)";
            }
        }
        return os.str();
    }
}
