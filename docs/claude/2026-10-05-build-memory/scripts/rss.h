// rss.h - resident memory of this process for the harnesses: VmRSS and VmHWM from /proc/self/status,
// and a sampler thread that keeps the highest VmRSS seen while a stage runs (VmHWM cannot be reset
// without privileges on WSL, so a stage's own peak is taken from samples every 5 ms).
#pragma once

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <thread>

inline long StatusKb(char const* field) {
    std::ifstream is("/proc/self/status");
    std::string line;
    size_t const n = std::strlen(field);
    while (std::getline(is, line)) {
        if (line.compare(0, n, field) == 0) return std::stol(line.substr(n + 1));
    }
    return -1;
}

inline double RssMb() { return StatusKb("VmRSS") / 1024.0; }
inline double HwmMb() { return StatusKb("VmHWM") / 1024.0; }

// The highest VmRSS while it lives.
class PeakSampler {
public:
    PeakSampler() : m_peak(RssMb()), m_thread([this] {
        while (!m_stop) {
            double const now = RssMb();
            double seen = m_peak.load();
            while (now > seen && !m_peak.compare_exchange_weak(seen, now)) {}
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
        }
    }) {}
    ~PeakSampler() { Stop(); }
    double Stop() {
        if (!m_stop.exchange(true)) m_thread.join();
        return m_peak;
    }

private:
    std::atomic<bool> m_stop{false};
    std::atomic<double> m_peak;
    std::thread m_thread;
};
