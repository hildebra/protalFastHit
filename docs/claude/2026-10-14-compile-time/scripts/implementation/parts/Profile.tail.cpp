
    void ProfilingAheadDeleter::operator()(ProfilingAhead* ahead) const {
        delete ahead;
    }

    ProfilingAheadPtr MakeProfilingAhead(ProfilingContext& ctx, size_t threads) {
        return ProfilingAheadPtr(new ProfilingAhead(ctx, threads));
    }

    void SubmitProfilingAhead(ProfilingAhead& ahead, size_t idx) {
        ahead.Submit(idx);
    }

    std::shared_ptr<std::vector<uint32_t> const> TaxonomyGenera(taxonomy::IntTaxonomy const& taxonomy) {
        return profiler::GeneraOf(taxonomy);
    }
}
