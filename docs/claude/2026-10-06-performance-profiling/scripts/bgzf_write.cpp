// bgzf_write.cpp: how fast protal writes BGZF (bgzf::CompressFile, Bgzf.h), and how large: compresses SRC to DST
// with `threads` threads, --repeats times, and prints the median wall and CPU seconds, MB/s of input and the size.
// Built against a protal source tree (libdeflate before, ISA-L after; build_read_bench.sh).
//
// usage: bgzf_write SRC DST [threads] [--repeats 5]
#include <sys/resource.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "IO/Bgzf.h"

namespace {
    double Cpu() {
        rusage r {};
        getrusage(RUSAGE_SELF, &r);
        return static_cast<double>(r.ru_utime.tv_sec + r.ru_stime.tv_sec) +
               static_cast<double>(r.ru_utime.tv_usec + r.ru_stime.tv_usec) / 1e6;
    }

    double Median(std::vector<double> v) {
        std::sort(v.begin(), v.end());
        return v[v.size() / 2];
    }
}

int main(int argc, char** argv) {
    std::vector<std::string> args;
    int repeats = 5;
    for (int i = 1; i < argc; i++) {
        std::string const arg = argv[i];
        if (arg == "--repeats" && i + 1 < argc) repeats = std::stoi(argv[++i]);
        else args.push_back(arg);
    }
    if (args.size() < 2) {
        std::fprintf(stderr, "usage: bgzf_write SRC DST [threads] [--repeats 5]\n");
        return 2;
    }
    int const threads = args.size() > 2 ? std::stoi(args[2]) : 1;
    std::vector<double> walls, cpus;
    for (int r = 0; r < repeats; r++) {
        auto const start = std::chrono::steady_clock::now();
        double const cpu = Cpu();
        std::string const error = protal::bgzf::CompressFile(args[0], args[1], threads);
        if (!error.empty()) {
            std::fprintf(stderr, "%s\n", error.c_str());
            return 1;
        }
        walls.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
        cpus.push_back(Cpu() - cpu);
    }
    double const in = static_cast<double>(std::filesystem::file_size(args[0]));
    std::printf("threads\twall_s\tMB_per_s\tcpu_s\tbytes_in\tbytes_out\n%d\t%.3f\t%.0f\t%.3f\t%.0f\t%ju\n", threads, Median(walls),
                in / Median(walls) / 1e6, Median(cpus), in, static_cast<uintmax_t>(std::filesystem::file_size(args[1])));
    return 0;
}
