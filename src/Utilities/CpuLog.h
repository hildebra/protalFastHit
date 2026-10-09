// CpuLog.h - misc/cpu.tsv of a run: for each stage, its wall-clock seconds and the CPU seconds the whole process used
// meanwhile (all threads, user and system time), so that cpu / wall tells how many cores a stage kept busy. The stages
// follow each other without gaps: the start-up (from the process's start: the database's index and tables loaded), each
// sample's alignment (up to its SAM written), the profiling stage, and the rest of the run (the strains, the outputs).
// With --profile_ahead the profiling worker's time falls into the samples' rows. A file of its own, not on the console:
// build_gtdb_database.py's collector copies it into its protal runs' CPU table (docs/claude/2026-10-09-build-parallelism),
// where a stage that leaves cores idle on a 120-core node shows.
#pragma once

#include <sys/resource.h>

#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>

namespace protal {
    // The process's CPU seconds so far: user and system time of all its threads.
    inline double ProcessCpuSeconds() {
        rusage usage{};
        if (getrusage(RUSAGE_SELF, &usage) != 0) return 0;
        auto seconds = [](timeval const& t) { return static_cast<double>(t.tv_sec) + static_cast<double>(t.tv_usec) / 1e6; };
        return seconds(usage.ru_utime) + seconds(usage.ru_stime);
    }

    class CpuLog {
    public:
        // The process's start, as near as a static's initialisation gets to it.
        static std::chrono::steady_clock::time_point ProcessStart() {
            static auto const start = std::chrono::steady_clock::now();
            return start;
        }

        static CpuLog& Get() {
            static CpuLog log;
            return log;
        }

        // Starts the table at `path` (its folder made if missing; an earlier table replaced). Rows are written by Mark;
        // nothing is written if the file cannot be opened.
        void Open(std::filesystem::path const& path, size_t threads) {
            std::error_code ec;
            std::filesystem::create_directories(path.parent_path(), ec);
            m_os.open(path, std::ios::out | std::ios::trunc);
            if (!m_os) return;
            m_threads = threads;
            m_os << "stage\twall_seconds\tcpu_seconds\tthreads\n" << std::fixed << std::setprecision(3);
            m_os.flush();
        }

        // A row for the stage that ends now: the time and CPU since the last row (the first: since the process started).
        void Mark(std::string const& stage) {
            auto const now = std::chrono::steady_clock::now();
            double const cpu = ProcessCpuSeconds();
            if (m_os) {
                double const wall = std::chrono::duration<double>(now - m_last).count();
                m_os << stage << '\t' << wall << '\t' << (cpu - m_last_cpu) << '\t' << m_threads << '\n';
                m_os.flush();
            }
            m_last = now;
            m_last_cpu = cpu;
        }

        // Marks `stage` when it goes out of scope (a sample's alignment, whichever way its loop iteration ends).
        class Stage {
        public:
            explicit Stage(std::string name) : m_name(std::move(name)) {}
            ~Stage() { CpuLog::Get().Mark(m_name); }
            Stage(Stage const&) = delete;
            Stage& operator=(Stage const&) = delete;

        private:
            std::string m_name;
        };

    private:
        CpuLog() : m_last(ProcessStart()) {}

        std::ofstream m_os;
        size_t m_threads = 0;
        std::chrono::steady_clock::time_point m_last;
        double m_last_cpu = 0;
    };

    // Takes the process's start when the program's statics are initialised (before main), not at the first Get().
    inline auto const kCpuLogProcessStart = CpuLog::ProcessStart();
}
