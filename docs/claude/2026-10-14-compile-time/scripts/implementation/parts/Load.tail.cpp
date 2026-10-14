
    void StartIndexLoading(Options const& options, ProtalDB& db) {
        // The index's field widths from the genes, before the preload changes them.
        auto const maxima = db.GetGenomes().IndexFieldMaxima();
        auto const fingerprint = db.GetGenomes().MapFingerprint();
        db.SetIndexLoading(std::async(std::launch::async, [&options, maxima, fingerprint]() {
            LoadedIndex loaded;
            loaded.timer.Start();
            std::ostringstream log;
            loaded.map = LoadIndex(options, maxima, fingerprint, log);
            loaded.timer.Stop();
            loaded.log = log.str();
            return loaded;
        }));
    }
}
