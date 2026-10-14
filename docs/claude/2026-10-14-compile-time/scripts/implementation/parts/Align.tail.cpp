
    void RunAlignment(Options& options, ProtalDB& db, std::function<void(size_t)> const& sample_done) {
        using AlignmentBenchmark = CoreBenchmark;
        // Untangle Template options that need to be written out specifically.
        if (options.BenchmarkAlignment()) {
            AlignmentBenchmark alignment_benchmark{};
            if (!options.GetBenchmarkAlignmentOutputFile().empty()) {
                alignment_benchmark.SetOutput(options.GetBenchmarkAlignmentOutputFile());
            }
            RunWrapper(options, db, alignment_benchmark, sample_done);
            if (!options.GetBenchmarkAlignmentOutputFile().empty()) {
                alignment_benchmark.DestroyOutput();
            }
        } else {
            RunWrapper(options, db, NoBenchmark{}, sample_done);
        }
    }
}
