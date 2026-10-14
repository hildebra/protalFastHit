    // Prototype split (compile-time study): protal::Run's stages, each defined in a translation unit of its own.
    namespace split {
        void Align(Options& options, ProtalDB& db, std::function<void(size_t)> const& sample_done);
        Profiles Profile(ProfilingContext& ctx, ProfilingAhead* ahead);
        std::unique_ptr<ProfilingAhead> MakeAhead(ProfilingContext& ctx, size_t threads);
        void Strains(Options& options, Profiles& profiles, GenomeLoader& loader, taxonomy::IntTaxonomy& taxonomy,
                     std::vector<uint32_t> msa_taxids, std::optional<profiler::TaxonFilterObj> filter);
        int BuildMaint(Options& options);
    }
    namespace split {
        Options ParseOptions(int argc, char* argv[]);
        void Help(Options& options, bool dev);
        std::string OptionsText(Options const& options);
    }
