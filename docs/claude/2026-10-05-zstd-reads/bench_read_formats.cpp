// How fast protal's input stream (ThreadedGzStream) hands over the reads of one file in each format: the reader
// takes FASTQ records in batches with TakeLines, as BufferedFastxReader::LoadBatch does, while the stream's own
// thread decompresses. Prints the median of --repeats runs per file: wall seconds, MB/s of FASTQ, and the
// process's CPU seconds (the decompressing thread's and the reader's).
//
// usage: bench_read_formats FILE... [--repeats 5]
#include <sys/resource.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <string>
#include <vector>

#include "IO/ThreadedGzStream.h"

namespace {
    double Cpu() {
        rusage r {};
        getrusage(RUSAGE_SELF, &r);
        return static_cast<double>(r.ru_utime.tv_sec + r.ru_stime.tv_sec) +
               static_cast<double>(r.ru_utime.tv_usec + r.ru_stime.tv_usec) / 1e6;
    }
}

int main(int argc, char** argv) {
    std::vector<std::string> files;
    int repeats = 5;
    for (int i = 1; i < argc; i++) {
        std::string const arg = argv[i];
        if (arg == "--repeats" && i + 1 < argc) repeats = std::stoi(argv[++i]);
        else files.push_back(arg);
    }
    std::printf("file\twall_s\tMB_per_s\tcpu_s\tbytes\terror\n");
    for (auto const& file : files) {
        std::vector<double> walls, cpus;
        size_t bytes = 0;
        std::string error;
        for (int r = 0; r < repeats; r++) {
            auto const start = std::chrono::steady_clock::now();
            double const cpu = Cpu();
            protal::ThreadedGzIstream is(file.c_str());
            std::string batch;
            bytes = 0;
            while (true) {
                batch.clear();
                if (is.rdbuf()->TakeLines(4 * 32, batch) == 0) break;  // 32 records, as SeqReaderSE takes them
                bytes += batch.size();
            }
            error = is.rdbuf()->read_failed() ? is.rdbuf()->read_error_message() : "";
            is.close();
            walls.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
            cpus.push_back(Cpu() - cpu);
        }
        std::sort(walls.begin(), walls.end());
        std::sort(cpus.begin(), cpus.end());
        double const wall = walls[walls.size() / 2];
        std::printf("%s\t%.3f\t%.0f\t%.3f\t%zu\t%s\n", file.c_str(), wall, static_cast<double>(bytes) / 1e6 / wall,
                    cpus[cpus.size() / 2], bytes, error.c_str());
    }
}
