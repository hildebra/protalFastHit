#pragma once

// What Run.cpp, Profile.cpp and Strains.cpp share of profiling (RunStages.h says what each unit holds).

#include <atomic>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include "Profiler.h"
#include "RunStages.h"

namespace protal {

    using Profile = profiler::MicrobialProfile;
    using Profiles = std::vector<Profile>;

    // The models the samples are profiled with (Options::SampleModelDbFile: of their read types, trained at the gates
    // nearest theirs), by ModelKey, each loaded once per read type.
    using ReadTypeModels = std::map<std::string, profiler::TaxonFilterObj>;

    // A model's key in ReadTypeModels: the read type it profiles and its file (one file given for several read types,
    // --model, is loaded and named for each).
    inline std::string ModelKey(ReadType type, db::DbFile const& file) {
        return Info(type).token + '\t' + file.Name();
    }

    // A taxon that a sample's profile leaves out (score below --knob) although its own reads are strong
    // evidence that it is present (profiler::StrongOwnEvidence).
    struct UnreportedTaxon {
        std::string species;
        std::string row;  // of unreported_species.tsv
    };

    // What profiling the samples of a run shares: the taxonomy, every taxid's genus and family (which of a read's
    // alternatives (ZA) are congeners of its taxon, and what the other taxa of a sample say of each:
    // MicrobialProfile::ApplySampleContext), the distances of the references (made as they are asked for, kept for
    // the run), the models, and a slot per sample for its profile and for the taxa it leaves out although their own
    // reads are strong evidence. Made before the alignment when samples are profiled while others are aligned
    // (ProfilingAhead), else by the profiling stage.
    struct ProfilingContext {
        Options& options;
        GenomeLoader& genomes;
        taxonomy::IntTaxonomy* taxonomy = nullptr;
        std::shared_ptr<std::vector<uint32_t> const> genera, families;
        std::shared_ptr<profiler::context::CongenerDistances> distances;
        ReadTypeModels const& models;
        std::vector<size_t> range;  // the samples; a sample's slot is its position here
        // MicrobialProfile holds a reference member so it is not assignable; a slot is filled in place.
        std::vector<std::optional<profiler::MicrobialProfile>> profile_slots;
        std::vector<std::vector<UnreportedTaxon>> unreported_slots;
        DeferredTables* deferred = nullptr;  // the tables still loading beside the alignment, if any (ProfileSample waits)
        // The samples whose profiling has not reached its end yet, and whether one is returning freed memory to the
        // system (malloc_trim): ProfileSample skips it after the run's last profile, and while another's runs.
        std::atomic<size_t> profiles_left;
        std::atomic<bool> trimming{ false };

        ProfilingContext(Options& options, ProtalDB& db, ReadTypeModels const& models) :
                options(options), genomes(db.GetGenomes()), models(models), range(options.GetRange()),
                profile_slots(range.size()), unreported_slots(range.size()), profiles_left(range.size()) {
            if (!db.IsTaxonomyLoaded()) db.LoadTaxonomy(options.TaxonomyDbFile());
            taxonomy = &db.GetTaxonomy();
            genera = profiler::GeneraOf(*taxonomy);
            families = profiler::FamiliesOf(*taxonomy);
            distances = std::make_shared<profiler::context::CongenerDistances>(genomes);
        }
    };

    // The parameters of the strain MSA's items of a sample with reads of `read_type` and a taxon whose own reads are
    // at least `min_identity` identical (StrainLevelContainer::MSAItem): the SNP filters of the options.
    inline profiler::MSAItemParameters MSAItemFor(Options const& options, ReadType read_type, double min_identity) {
        profiler::MSAItemParameters p;
        p.min_identity = min_identity;
        p.min_cov = options.GetSNPMinCov();
        p.min_af = options.GetSNPMinAF(read_type);
        p.min_mean_qual = options.GetSNPMinMeanQual();
        p.min_qual_sum = options.GetSNPMinPhredSum();
        p.require_strand = options.GetSNPRequireStrand();
        return p;
    }

    // What a sample's taxa that enter the MSAs keep of their genes (Taxon::ReleaseReadData): the items at the identity
    // of each taxon's own reads, and the SNP counts the strain stage reads (SelectGenesForTaxon: (2, 60) and (0, 0);
    // GetMSAForTaxon: the SNP filters').
    inline profiler::MSAReleaseParameters MSAReleaseFor(Options const& options, ReadType read_type) {
        profiler::MSAReleaseParameters p;
        p.identity_margin = options.GetMSAIdentityMargin();
        p.item = MSAItemFor(options, read_type, 0);
        p.count_keys = { { 2, 60 }, { 0, 0 }, { static_cast<uint32_t>(options.GetSNPMinCov()), options.GetSNPMinPhredSum() } };
        return p;
    }

    // Profile.cpp. ProfilingAhead (--profile_ahead) profiles samples while others are still being aligned; Run holds it
    // through these, so that only Profile.cpp compiles the profiling itself.
    class ProfilingAhead;
    struct ProfilingAheadDeleter {
        void operator()(ProfilingAhead* ahead) const;
    };
    using ProfilingAheadPtr = std::unique_ptr<ProfilingAhead, ProfilingAheadDeleter>;
    ProfilingAheadPtr MakeProfilingAhead(ProfilingContext& ctx, size_t threads);
    // The sample of slot idx has its SAM: the worker profiles it when it is free.
    void SubmitProfilingAhead(ProfilingAhead& ahead, size_t idx);
    // Profiles the samples of ctx that are not yet profiled, each with the model of its kind of reads, and returns
    // every sample's profile (in the order of the samples). With `ahead`, the samples it profiled or is profiling
    // while the reads were aligned are left to it, and it is waited for.
    Profiles ProfileWrapper(ProfilingContext& ctx, ProfilingAhead* ahead = nullptr);

    // Strains.cpp.
    using TaxidList = std::vector<uint32_t>;
    // Taxa in at least `min_samples` profiles (entering the MSAs by `filter`, if given), the most frequent first.
    TaxidList ExtractTaxa(Profiles const& profiles, std::optional<profiler::TaxonFilterObj> filter = {}, size_t min_samples = 2);
    // The taxids of --msa_species (s__Genus_species names of the taxonomy's species; exits 2 if one is not); none without it.
    std::vector<uint32_t> ResolveMSASpecies(Options& options, taxonomy::IntTaxonomy& taxonomy);
    // The strain-level MSAs of the taxa in `msa_taxids` (if empty, those ExtractTaxa finds) over the profiles.
    void StrainWrapper2(Options& options, Profiles& profiles, GenomeLoader& loader, taxonomy::IntTaxonomy& taxonomy,
                        std::vector<uint32_t> msa_taxids = {}, std::optional<profiler::TaxonFilterObj> filter = {});
}
