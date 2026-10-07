// SampleContext.h - what a taxon's call takes from the other taxa of its sample.
//
// Amplicon denoisers (AmpliconNoise, DADA2, UNOISE) judge a candidate sequence against the more abundant sequences
// that could have produced it, at an error rate that falls with the distance between them: their bar rises with the
// abundance of the possible parent, not with the depth of the sample (docs/claude/2026-10-02-amplicon-denoising). In
// protal a taxon's reads may be those of a more abundant congener that align a little better to it (spill-over),
// and false calls grew with depth. What this file gives the profiler (MicrobialProfile::ApplySampleContext):
//   - how abundant a taxon's relatives in the sample are, by rank and by the distance between their references
//     (CongenerDistances: Mash distances of sketched marker genes, the build's congener comparison at k = 12), and
//     how many of its reads they would spill onto it at SpillRate(distance) per read;
//   - how much of its reads an abundance-weighted assignment leaves to it: an EM over the reads' alternatives (the ZA
//     tag), each read's taxon weighted by the taxa's abundance, as AmpliconNoise weighs its clusters
//     (AbundanceWeightedShares);
//   - calibrated calls that keep a sample's expected share of false calls at a target (FalseCallKnob), the prior
//     adjusted to the sample (SampleAdjusted): a deep sample has many more absent candidates, so each of its taxa is
//     less likely present at the same score, as DADA2 corrects its test for the number of sequences it tests.
#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "GeneConservation.h"
#include "GenomeLoader.h"

namespace protal::profiler::context {

    // ---- distances between congeners' references ---------------------------------------------------------

    inline constexpr size_t kSketchSize = 64;      // the smallest hashes of a gene's k-mers kept (Mash's bottom sketch)
    inline constexpr size_t kMinSharedGenes = 10;  // genes two references must share for a distance (as the build's)
    inline constexpr size_t kMaxSketches = 5000;   // sketches cached in a run before the cache starts over (~25 kB each)
    inline constexpr size_t kMaxPairs = 200000;    // and compared pairs (~1 kB each); both are recomputed alike
    inline constexpr double kFarDistance = 1;      // two references without kMinSharedGenes shared sketched genes

    // The hash of a k-mer, a gene's bottom sketch of kSketchSize hashes and the Mash distance of two sketches: the
    // build's (GeneConservation.h), which GeneIncongruence.h uses too.
    using gene_conservation::MixKmer;
    using gene_conservation::SketchDistance;

    inline std::vector<uint32_t> GeneSketch(std::string_view seq) {
        return gene_conservation::BottomSketch(seq, kSketchSize);
    }

    // A reference's sketches by gene id, ascending.
    using TaxonSketch = std::vector<std::pair<uint32_t, std::vector<uint32_t>>>;

    // Two references compared gene by gene: the SketchDistance of each gene both have (by gene id), and their distance,
    // the median of those, kFarDistance if they share fewer than kMinSharedGenes (the build's typical distance of two
    // congeners, gene_congeners.tsv).
    struct PairDistances {
        double distance = kFarDistance;
        std::vector<std::pair<uint32_t, float>> genes;
    };

    inline PairDistances SketchedTaxonDistances(TaxonSketch const& a, TaxonSketch const& b) {
        PairDistances result;
        auto i = a.begin();
        auto j = b.begin();
        while (i != a.end() && j != b.end()) {
            if (i->first < j->first) {
                ++i;
            } else if (j->first < i->first) {
                ++j;
            } else {
                result.genes.emplace_back(i->first, static_cast<float>(SketchDistance(i->second, j->second)));
                ++i;
                ++j;
            }
        }
        if (result.genes.size() >= kMinSharedGenes) {
            std::vector<double> distances;
            distances.reserve(result.genes.size());
            for (auto const& [_, d] : result.genes) distances.push_back(d);
            result.distance = gene_conservation::Median(std::move(distances));
        }
        return result;
    }

    // Distances between the database's references (SketchedTaxonDistances of their hittable genes), computed when
    // first asked for and kept for the run; any thread may ask. A distance depends only on the two references, so the
    // features computed from them are the same on any number of threads.
    class CongenerDistances {
    public:
        explicit CongenerDistances(GenomeLoader& loader) : m_loader(&loader) {}

        // The distance of the two references (PairDistances::distance).
        double Between(uint32_t a, uint32_t b) {
            if (a == b) return 0;
            return Compared(a, b)->distance;
        }

        // Makes the sketch of taxid now if the run has none yet (SketchOf): for one pass over a sample's taxa on
        // all its threads before their pairs are compared, each sketch once.
        void Sketch(uint32_t taxid) {
            SketchOf(taxid);
        }

        // The two references gene by gene (PairDistances).
        std::shared_ptr<PairDistances const> Compared(uint32_t a, uint32_t b) {
            uint64_t const key = (static_cast<uint64_t>(std::min(a, b)) << 32) | std::max(a, b);
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                auto const found = m_pairs.find(key);
                if (found != m_pairs.end()) return found->second;
            }
            auto const sa = SketchOf(a);
            auto const sb = SketchOf(b);
            auto compared = std::make_shared<PairDistances const>(SketchedTaxonDistances(*sa, *sb));
            std::lock_guard<std::mutex> lock(m_mutex);
            if (m_pairs.size() >= kMaxPairs) m_pairs.clear();
            return m_pairs.emplace(key, std::move(compared)).first->second;
        }

    private:
        std::shared_ptr<TaxonSketch const> SketchOf(uint32_t taxid) {
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                auto const found = m_sketches.find(taxid);
                if (found != m_sketches.end()) return found->second;
            }
            auto sketch = std::make_shared<TaxonSketch const>(MakeSketch(taxid));
            std::lock_guard<std::mutex> lock(m_mutex);
            if (m_sketches.size() >= kMaxSketches) m_sketches.clear();
            return m_sketches.emplace(taxid, std::move(sketch)).first->second;
        }

        TaxonSketch MakeSketch(uint32_t taxid) {
            TaxonSketch sketch;
            if (!m_loader->GetGenomeMap().contains(taxid)) return sketch;
            auto& genome = m_loader->GetGenome(taxid);
            for (uint32_t const id : genome.GetHittableGenes()) {
                if (!genome.HasGene(id)) continue;
                auto const sequence = genome.GetGeneOMP(id).Sequence();  // keeps the bases its view points to
                auto hashes = GeneSketch(sequence.View());
                if (!hashes.empty()) sketch.emplace_back(id, std::move(hashes));
            }
            std::sort(sketch.begin(), sketch.end(), [](auto const& x, auto const& y) { return x.first < y.first; });
            return sketch;
        }

        GenomeLoader* m_loader;
        std::mutex m_mutex;
        std::unordered_map<uint64_t, std::shared_ptr<PairDistances const>> m_pairs;
        std::unordered_map<uint32_t, std::shared_ptr<TaxonSketch const>> m_sketches;
    };

    // ---- spill-over from relatives ------------------------------------------------------------------------

    // The share of a relative's fragments expected on a taxon whose reference is `distance` from the relative's (Mash,
    // about the share of bases that differ): kSpillAtZero for identical references, ten times less per kSpillDecade.
    // A fixed prior, as UNOISE's skew 1/2^(2d+1) per difference is (one decade per 5% here, as a shotgun read of a
    // congener 5% away lands on the wrong reference far less often than an amplicon one error away); the model learns
    // from relative_skew and relative_distance themselves what they mean.
    inline constexpr double kSpillAtZero = 1e-2;
    inline constexpr double kSpillDecade = 0.05;

    inline double SpillRate(double distance) {
        return kSpillAtZero * std::pow(10.0, -std::max(distance, 0.0) / kSpillDecade);
    }

    // genus_spill's rates by rank alone (no distances): a congener's fragments, and those of a species of another genus
    // of the family.
    inline constexpr double kGenusSpill = 1e-3;
    inline constexpr double kFamilySpill = 1e-4;

    // The congeners a taxon is compared with by distance: those with more fragments than it, the most first, at most
    // this many.
    inline constexpr size_t kMaxRelatives = 4;

    // The singleton rule's default (--singleton_congener): a taxon of a single fragment beside a congener of this many
    // fragments or more is not reported if its read looks like the congener's: the abundance-weighted assignment gives
    // it mostly to others (em_own_share below kSingletonOwnShare), or it is further from the taxon's reference than a
    // strain of the species would be (identity below kSingletonIdentity). Of the r226 v3 build's paired-end calls of
    // such taxa, 40 of 40 (species held out) and 19 of 19 (test set) were false; but on the benchmark world's samples
    // with congener groups the rule without the read's condition also removed 2-11 true calls per test set, of minor
    // congeners whose read fits their own reference (identity ~0.975, EM share ~0.98), which the condition keeps
    // (docs/claude/2026-10-03-denoising-implementation).
    inline constexpr size_t kSingletonCongener = 0;  // off: on r226 training data the model called none of the vetoed taxa
    inline constexpr double kSingletonOwnShare = 0.5;
    inline constexpr double kSingletonIdentity = 0.95;

    // ---- the sample's complexity ---------------------------------------------------------------------------

    // What all the taxa of a sample say together, before any is called (SampleEvidence::sample_log_taxa,
    // sample_low_identity, sample_identity): a soil sample has thousands of taxa and much of its reads on relatives the
    // database lacks, a gut sample hundreds, mostly on species the database holds. Without them a model can tell the
    // kind of community only by the sample's depth, which the simulated scenarios draw from a narrow band around each
    // preset: at GTDB r226 (v14), with whole samples held out, pe's shallowest shallow-soil sample, below every other soil
    // sample's depth, scored F1 0.758 against 0.920 with species held out; with these three as features the six shallow
    // samples scored 0.935 against 0.908 (docs/claude/2026-10-06-r226-v14).
    inline constexpr size_t kSampleIdentityMinFragments = 10;  // the taxa sample_identity is the median of, if any

    // A taxon as SampleComplexityOf sees it: its fragments, identity (BaseIdentity) and share of low-identity bases
    // (LowIdentityShare).
    struct TaxonSummary {
        size_t fragments = 0;
        double identity = 0;
        double low_identity_share = 0;
    };

    struct SampleComplexity {
        double log_taxa = 0;      // log10 of the taxa with fragments, at least 0
        double low_identity = 0;  // the fragments' share on low-identity bases: the taxa's low_identity_share weighted by
                                  // their fragments; 0 without fragments
        double identity = 0;      // the fragment-weighted median identity of the taxa with at least
                                  // kSampleIdentityMinFragments fragments (of every taxon with fragments if none has
                                  // that many); 0 without fragments
    };

    // The sample's complexity from its taxa, given in a fixed order (taxids ascending): taxa of equal identity are
    // taken in that order for the median, so the values do not depend on the order of the reads or the threads. The
    // median is the identity of the first taxon, by identity, at which the cumulative fragments reach half of the
    // taxa's (as numpy's searchsorted over the cumulative sum finds it).
    inline SampleComplexity SampleComplexityOf(std::vector<TaxonSummary> const& taxa) {
        SampleComplexity c;
        size_t with_fragments = 0, deep = 0;
        double fragments = 0, low = 0;
        for (auto const& t : taxa) {
            if (t.fragments == 0) continue;
            with_fragments++;
            deep += t.fragments >= kSampleIdentityMinFragments;
            fragments += static_cast<double>(t.fragments);
            low += static_cast<double>(t.fragments) * t.low_identity_share;
        }
        c.log_taxa = std::log10(static_cast<double>(std::max<size_t>(with_fragments, 1)));
        if (fragments == 0) return c;
        c.low_identity = low / fragments;
        size_t const min_fragments = deep > 0 ? kSampleIdentityMinFragments : 1;
        std::vector<std::pair<double, size_t>> by_identity;  // (identity, fragments)
        double total = 0;
        for (auto const& t : taxa) {
            if (t.fragments < min_fragments) continue;
            by_identity.emplace_back(t.identity, t.fragments);
            total += static_cast<double>(t.fragments);
        }
        std::stable_sort(by_identity.begin(), by_identity.end(), [](auto const& a, auto const& b) { return a.first < b.first; });
        double cumulative = 0;
        for (auto const& [identity, n] : by_identity) {
            cumulative += static_cast<double>(n);
            if (cumulative >= total / 2) {
                c.identity = identity;
                break;
            }
        }
        return c;
    }

    // ---- abundance-weighted assignment of ambiguous reads -------------------------------------------------

    // One kind of best record with alternatives (RecordEvidenceCollector::NoteRecord): its taxon, whether the MAPQ and
    // length filters keep it (1) or not (0), then each alternative taxon (ZA) and the edits more it takes, by taxon.
    using AmbiguityKey = std::vector<uint32_t>;
    using AmbiguityClasses = std::map<AmbiguityKey, uint64_t>;

    // The likelihood of a read under a taxon whose reference it fits one edit worse, against one it fits best: the
    // chance that one of that taxon's read bases differs from its reference (its read errors and its strain's
    // divergence), over the chance that it does not, θ / (1 - θ), θ the taxon's own differences per aligned base in the
    // sample (RecordCounts), within kMinDivergence..kMaxDivergence; kEditRatio for a taxon without aligned bases. A
    // strain of an abundant species far from its reference differs from it at many bases of its own reads, and so
    // explains a read that fits a congener's reference a few edits better more readily than a strain close to it does.
    inline constexpr double kEditRatio = 0.05;
    inline constexpr double kMinDivergence = 0.002;
    inline constexpr double kMaxDivergence = 0.2;
    // The sweeps stop when no taxon's own share (its reads left to it over its records, the feature the EM is for)
    // changed by kEmTolerance or more in the last sweep, or after kEmIterations. The first rule was a relative change
    // of any weight below 1e-10, which a sample that converges reached at about twice the sweeps (a 5M-pair sample:
    // 197 against 105) for the same shares to 1e-6; a taxon whose reads a congener explains about as well (its weight
    // times its edit factor close to the taxon's reads) loses them at a rate within 1e-4 per sweep of standing still
    // and runs to the cap under either rule (docs/claude/2026-10-03-performance-review). The cap was 200 sweeps to
    // 0.7.6; at GTDB r226 a paired-end sample's EM ran to it every time (319k classes, 1.1 s of the profiling), so it
    // is 100 since 0.7.7: about what a converging sample needs, and half the time for those that do not converge,
    // whose slowest taxa (reads a congener explains about as well) are left with a somewhat larger share.
    inline constexpr size_t kEmIterations = 100;
    inline constexpr double kEmTolerance = 1e-6;

    // A taxon's best records (all, and those the filters keep) and their differences and aligned bases, for
    // AbundanceWeightedShares.
    struct RecordCounts {
        uint64_t records = 0;
        uint64_t kept = 0;
        uint64_t differences = 0;
        uint64_t aligned = 0;
    };

    // The edit ratio (kEditRatio) of a taxon of these records.
    inline double EditRatio(RecordCounts const& c) {
        if (c.aligned == 0) return kEditRatio;
        double const theta = std::clamp(static_cast<double>(c.differences) / static_cast<double>(c.aligned),
                                        kMinDivergence, kMaxDivergence);
        return theta / (1 - theta);
    }

    // Of a taxon's best records (all, and those kept), the share an abundance-weighted assignment leaves to it.
    struct OwnShares {
        double all = 1;
        double kept = 1;
    };

    // Each read is given to the taxa it fits in proportion to the taxon's reads times its EditRatio^(edits more): an EM
    // over the reads' alternatives (ZA lists those within 5 edits), as a mixture model's clusters share their reads by
    // their weights (AmpliconNoise). A read that fits a congener a hundred times as abundant two edits worse is mostly
    // the congener's (100 x 0.05^2 = 0.25 against 1 at 5% divergence of the congener's reads), one it fits as well
    // nearly all. The taxa start at their own record counts; reads without alternatives stay with their taxon. Classes
    // are visited in key order and the result depends only on the records, not on the threads that read them.
    //
    // The sweeps run on dense vectors: every taxon of the counts or of a class gets an index once, and a class is
    // resolved once into its taxon's index, its alternatives' indices and their constant factors
    // EditRatio^(edits more) (the ratio of a taxon does not change between sweeps). A sweep is then a multiply per
    // alternative and the additions, in the order the first version did them with hash maps and a pow per visit,
    // so the shares are the same to the last bit (test ReferenceShares) at a small fraction of the instructions
    // (docs/claude/2026-10-03-performance-review). A taxon without records (an alternative the counts lack) has
    // weight 0 throughout, as it had no entry before.
    //
    // What the last AbundanceWeightedShares did, for the run's log (`stats`).
    struct EmStats {
        size_t classes = 0, alternatives = 0, taxa = 0, named = 0, sweeps = 0;  // named: the taxa classes name
        double setup = 0, sweeping = 0, sharing = 0;  // seconds: the indices and classes, the sweeps, the shares
    };

    // `tolerance` and `max_sweeps` are kEmTolerance and kEmIterations (tests pass others: 0 runs every sweep).
    inline std::unordered_map<uint32_t, OwnShares> AbundanceWeightedShares(AmbiguityClasses const& classes,
                                                                           std::map<uint32_t, RecordCounts> const& counts,
                                                                           double tolerance = kEmTolerance,
                                                                           size_t max_sweeps = kEmIterations,
                                                                           EmStats* stats = nullptr) {
        auto clock = std::chrono::steady_clock::now();
        auto lap = [&clock]() {
            auto const now = std::chrono::steady_clock::now();
            double const seconds = std::chrono::duration<double>(now - clock).count();
            clock = now;
            return seconds;
        };
        size_t sweeps = 0;
        // Indices: the taxa of the counts first, then any other taxon a class names, in class order.
        std::unordered_map<uint32_t, uint32_t> index;
        index.reserve(counts.size() * 2);
        std::vector<double> plain, weight, ratio, records;
        auto index_of = [&](uint32_t taxid) {
            auto const [it, added] = index.emplace(taxid, static_cast<uint32_t>(plain.size()));
            if (added) {
                plain.push_back(0);
                weight.push_back(0);
                ratio.push_back(0);
                records.push_back(0);
            }
            return it->second;
        };
        for (auto const& [taxid, c] : counts) {
            auto const i = index_of(taxid);
            plain[i] = static_cast<double>(c.records);
            weight[i] = static_cast<double>(c.records);
            ratio[i] = EditRatio(c);
            records[i] = static_cast<double>(c.records);
        }

        // The classes, resolved: own index, reads, kept flag, and the alternatives' (index, factor).
        struct Alternative {
            uint32_t index;
            double factor;  // ratio[index] ^ edits more, 0 for a taxon without records (its weight is 0 anyway)
        };
        struct Class {
            uint32_t own;
            bool kept;
            uint64_t reads;
            double n;        // reads, as a double
            uint32_t first;  // its alternatives in `alternatives`
            uint32_t count;
        };
        std::vector<Class> resolved;
        std::vector<Alternative> alternatives;
        resolved.reserve(classes.size());
        for (auto const& [key, n] : classes) {
            Class c{ index_of(key[0]), key[1] != 0, n, static_cast<double>(n), static_cast<uint32_t>(alternatives.size()), 0 };
            for (size_t i = 2; i + 1 < key.size(); i += 2) {
                auto const a = index_of(key[i]);
                bool const has_records = counts.contains(key[i]);
                alternatives.push_back({ a, has_records ? std::pow(ratio[a], static_cast<double>(key[i + 1])) : 0.0 });
                c.count++;
            }
            resolved.push_back(c);
        }
        for (auto const& c : resolved) plain[c.own] -= c.n;
        double const setup_seconds = lap();

        // A class's shares (own first, then its alternatives'), in a buffer on the stack (ZA lists a few alternatives; a
        // longer list goes to a vector).
        constexpr uint32_t kStackShares = 16;
        double stack_shares[1 + kStackShares];
        std::vector<double> long_shares;
        double* shares = stack_shares;
        auto posterior = [&](Class const& c) {
            shares = c.count <= kStackShares ? stack_shares : (long_shares.resize(1 + c.count), long_shares.data());
            shares[0] = weight[c.own];
            double total = shares[0];
            for (uint32_t k = 0; k < c.count; k++) {
                auto const& alt = alternatives[c.first + k];
                shares[1 + k] = weight[alt.index] * alt.factor;
                total += shares[1 + k];
            }
            if (total <= 0) {
                std::fill(shares, shares + 1 + c.count, 0.0);
                shares[0] = 1;
                return;
            }
            for (uint32_t k = 0; k <= c.count; k++) shares[k] /= total;
        };
        // Only the taxa a class names change in a sweep: the others keep their records as weight and their share. At r226 a
        // sample's taxa are all 143,614 species (each has failed candidates), its classes name a fraction of them; the
        // sweeps reset and check those alone, which gives the same weights and changes as resetting and checking all.
        std::vector<char> named(plain.size(), 0), owns(plain.size(), 0);
        for (auto const& c : resolved) {
            owns[c.own] = 1;
            named[c.own] = 1;
            for (uint32_t k = 0; k < c.count; k++) named[alternatives[c.first + k].index] = 1;
        }
        std::vector<uint32_t> named_taxa, owners;
        for (uint32_t i = 0; i < plain.size(); i++) {
            if (named[i]) named_taxa.push_back(i);
            if (owns[i]) owners.push_back(i);
        }
        std::vector<double> next = plain, own(plain.size(), 0.0), previous_share(plain.size(), 0.0);
        for (size_t iteration = 0; iteration < max_sweeps && !resolved.empty(); iteration++) {
            for (auto const i : named_taxa) next[i] = plain[i];
            for (auto const i : owners) own[i] = 0.0;
            for (auto const& c : resolved) {
                posterior(c);
                double const to_own = c.n * shares[0];
                next[c.own] += to_own;
                own[c.own] += to_own;
                for (uint32_t k = 0; k < c.count; k++) {
                    if (shares[1 + k] > 0) next[alternatives[c.first + k].index] += c.n * shares[1 + k];
                }
            }
            // The sweep's largest change of any taxon's own share: converged when below the tolerance. A taxon that owns
            // no class keeps its share (`own` stays 0), so only the owners can change.
            double change = 0;
            for (auto const i : owners) {
                double const share = (plain[i] + own[i]) / std::max(1.0, records[i]);
                if (iteration > 0) change = std::max(change, std::abs(share - previous_share[i]));
                previous_share[i] = share;
            }
            weight.swap(next);
            sweeps = iteration + 1;
            if (iteration > 0 && change < tolerance) break;
        }

        double const sweep_seconds = lap();
        std::vector<double> own_all(plain.size(), 0.0), own_kept(plain.size(), 0.0);
        std::vector<uint64_t> ambiguous_kept(plain.size(), 0);
        for (auto const& c : resolved) {
            posterior(c);
            own_all[c.own] += c.n * shares[0];
            if (c.kept) {
                own_kept[c.own] += c.n * shares[0];
                ambiguous_kept[c.own] += c.reads;
            }
        }
        std::unordered_map<uint32_t, OwnShares> result;
        result.reserve(counts.size());
        for (auto const& [taxid, c] : counts) {
            auto const i = index.at(taxid);
            OwnShares s;
            if (c.records > 0) s.all = std::clamp((plain[i] + own_all[i]) / static_cast<double>(c.records), 0.0, 1.0);
            if (c.kept > 0) {
                double const plain_kept = static_cast<double>(c.kept) - static_cast<double>(ambiguous_kept[i]);
                s.kept = std::clamp((plain_kept + own_kept[i]) / static_cast<double>(c.kept), 0.0, 1.0);
            }
            result[taxid] = s;
        }
        if (stats) *stats = { resolved.size(), alternatives.size(), plain.size(), named_taxa.size(), sweeps, setup_seconds, sweep_seconds, lap() };
        return result;
    }

    // ---- calibrated calls with a target share of false calls ----------------------------------------------

    // The model's score -> the probability that a taxon is present, piecewise linear, x increasing, y not decreasing
    // (machine_learning_cmdline.py --fdr-calls: an isotonic fit on species held out in training).
    using CalibrationCurve = std::vector<std::pair<double, double>>;

    // The curve at `score`: linear between its points, its end points' beyond them.
    inline double Calibrated(CalibrationCurve const& curve, double score) {
        if (curve.empty()) return score;
        if (score <= curve.front().first) return curve.front().second;
        if (score >= curve.back().first) return curve.back().second;
        auto const next = std::upper_bound(curve.begin(), curve.end(), score,
                                           [](double v, auto const& point) { return v < point.first; });
        auto const& [x1, y1] = *next;
        auto const& [x0, y0] = *std::prev(next);
        return x1 > x0 ? y0 + (y1 - y0) * (score - x0) / (x1 - x0) : y1;
    }

    inline constexpr double kMinProbability = 1e-6;
    inline constexpr size_t kPriorIterations = 1000;
    inline constexpr double kPriorTolerance = 1e-10;
    // Candidates at the training rate added to every sample's in SampleAdjusted: a sample of a few candidates, all
    // likely, would otherwise reach a rate of 1, and every probability 1 with it.
    inline constexpr double kPriorPseudoCount = 20;

    // The calibrated probabilities `q` of a sample's taxa, made for candidates that were present at rate `prior` (the
    // training rows'), adjusted to the sample's own rate: the fixed point of the EM of Saerens, Latinne and Decaestecker
    // (2002), the rate the mean of the adjusted probabilities (with kPriorPseudoCount candidates at `prior`), each q's
    // odds times the rates' odds ratio. A deep sample's many absent candidates lower its rate and with it every
    // probability. The rate goes to `sample_prior`.
    inline std::vector<double> SampleAdjusted(std::vector<double> q, double prior, double* sample_prior = nullptr) {
        prior = std::clamp(prior, kMinProbability, 1 - kMinProbability);
        for (auto& v : q) v = std::clamp(v, kMinProbability, 1 - kMinProbability);
        double rate = prior;
        std::vector<double> adjusted(q.size());
        for (size_t iteration = 0; iteration < kPriorIterations && !q.empty(); iteration++) {
            double const a = rate / prior, b = (1 - rate) / (1 - prior);
            double sum = 0;
            for (size_t i = 0; i < q.size(); i++) {
                adjusted[i] = a * q[i] / (a * q[i] + b * (1 - q[i]));
                sum += adjusted[i];
            }
            double const next = std::clamp((sum + kPriorPseudoCount * prior) / (static_cast<double>(q.size()) + kPriorPseudoCount),
                                           kMinProbability, 1 - kMinProbability);
            bool const done = std::abs(next - rate) < kPriorTolerance;
            rate = next;
            if (done) break;
        }
        double const a = rate / prior, b = (1 - rate) / (1 - prior);
        for (size_t i = 0; i < q.size(); i++) adjusted[i] = a * q[i] / (a * q[i] + b * (1 - q[i]));
        if (sample_prior) *sample_prior = rate;
        return adjusted;
    }

    // What FalseCallKnob decided for a sample.
    struct FalseCalls {
        double knob = 0.5;          // the lowest score called; above every score when none is
        size_t called = 0;          // taxa called
        double sample_prior = 0;    // the sample's rate of present candidates (SampleAdjusted)
        double expected_false = 0;  // the expected false calls among them: their summed 1 - probability
    };

    // The threshold for a sample's model scores (its taxa that may be called: not vetoed) at which the expected share
    // of false calls stays at or below `fdr`: the taxa by score, the highest first, each with its probability (the
    // calibrated score adjusted to the sample, SampleAdjusted), called as long as the mean of their (1 - probability)
    // is at most fdr. As the probabilities fall with the score, the calls are the highest-scoring taxa.
    inline FalseCalls FalseCallKnob(std::vector<double> scores, CalibrationCurve const& curve, double prior, double fdr) {
        FalseCalls result;
        std::sort(scores.begin(), scores.end(), std::greater<>());
        std::vector<double> q(scores.size());
        for (size_t i = 0; i < scores.size(); i++) q[i] = Calibrated(curve, scores[i]);
        auto const p = SampleAdjusted(q, prior, &result.sample_prior);
        double false_sum = 0;
        for (size_t i = 0; i < p.size(); i++) {
            double const sum = false_sum + (1 - p[i]);
            if (sum > fdr * static_cast<double>(i + 1)) break;
            false_sum = sum;
            result.called = i + 1;
        }
        result.expected_false = false_sum;
        double const top = scores.empty() ? 1.0 : scores.front();
        result.knob = result.called > 0 ? scores[result.called - 1]
                                        : std::nextafter(std::max(top, 1.0), std::numeric_limits<double>::infinity());
        return result;
    }
}
