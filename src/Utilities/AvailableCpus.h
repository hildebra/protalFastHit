// AvailableCpus.h - the CPUs this process may run on: its affinity (a SLURM job's cpuset, taskset's), else the online
// CPUs. For a run's threads against them (RunProtal's input line, ThreadedGzStreambuf's reading threads).
#pragma once

#include <sched.h>

#include <cstddef>
#include <thread>

namespace protal {

    inline size_t AvailableCpus() {
        cpu_set_t set;
        CPU_ZERO(&set);
        if (sched_getaffinity(0, sizeof(set), &set) == 0 && CPU_COUNT(&set) > 0) return static_cast<size_t>(CPU_COUNT(&set));
        unsigned const n = std::thread::hardware_concurrency();
        return n > 0 ? n : 1;
    }

}  // namespace protal
