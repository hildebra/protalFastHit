            {   // probe: per-sweep convergence (experiment build only)
                std::vector<double> oa(plain.size(), 0.0);
                for (auto const& c : resolved) { posterior(c); oa[c.own] += c.n * shares[0]; }
                static thread_local std::vector<double> prev;
                if (iteration == 0) prev.assign(oa.size(), -1.0);
                double ds = 0, dr = 0;
                for (size_t i = 0; i < oa.size(); i++) {
                    double const r = plain[i] + oa[i];
                    if (prev[i] >= -0.5) { ds = std::max(ds, std::abs(r - prev[i])); dr = std::max(dr, std::abs(r - prev[i]) / std::max(1.0, r)); }
                    prev[i] = r;
                }
                std::cerr << "[sweep] " << iteration + 1 << " change " << change << " max_own_reads_change " << ds << " max_own_relative_change " << dr << std::endl;
            }
