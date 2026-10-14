#pragma once

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#ifdef __GLIBC__
#include <malloc.h>
#endif

// The process's memory, as the build and a run report it (from Build.h, so that the units of a run that do not build
// need not include it).
namespace protal::build {
    // The build's resident memory now and its peak so far (VmRSS and VmHWM in /proc/self/status), after
    // a phase: which phase sets a GTDB build's peak (docs/claude/2026-10-05-build-memory). Nothing
    // where /proc is not there.
    inline void PrintMemory(std::string const& after) {
        std::ifstream status("/proc/self/status");
        std::string line;
        double rss = -1, hwm = -1;
        while (std::getline(status, line)) {
            if (line.rfind("VmRSS:", 0) == 0) rss = std::stod(line.substr(6));
            else if (line.rfind("VmHWM:", 0) == 0) hwm = std::stod(line.substr(6));
        }
        if (rss < 0 || hwm < 0) return;
        char text[96];
        std::snprintf(text, sizeof(text), "%.2f GB resident, peak %.2f GB", rss / (1024.0 * 1024.0), hwm / (1024.0 * 1024.0));
        std::cout << "Memory after " << after << ": " << text << std::endl;
    }

    // Gives the memory freed so far back to the system. glibc keeps freed small blocks in its
    // per-thread arenas: the suspect-copy scan's millions of sketches left ~12 GB resident at r226
    // scale (64 threads) that no later phase reused, under the index.
    inline void ReleaseFreeMemory() {
#ifdef __GLIBC__
        malloc_trim(0);
#endif
    }

    inline std::string HumanBytes(uint64_t bytes) {
        char buffer[32];
        if (bytes >= (uint64_t{1} << 30)) std::snprintf(buffer, sizeof(buffer), "%.2f GB", bytes / double(uint64_t{1} << 30));
        else if (bytes >= (uint64_t{1} << 20)) std::snprintf(buffer, sizeof(buffer), "%.2f MB", bytes / double(uint64_t{1} << 20));
        else std::snprintf(buffer, sizeof(buffer), "%llu bytes", static_cast<unsigned long long>(bytes));
        return buffer;
    }
}
