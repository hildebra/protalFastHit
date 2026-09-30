//
// Created by joachim on 05/06/2020.
//

#pragma once

#include <string>
#include <chrono>
#include <iostream>
#include <cstdint>

using namespace std::chrono;
using namespace std;

namespace protal {
    enum Time {
        days, hours, minutes, seconds, milliseconds, microseconds
    };

    // Time spent in a stage: the sum of its intervals (Start to Stop) and, after Join, of other
    // threads' sums. Intervals are summed in nanoseconds, as most per-read stages take well under a
    // microsecond; flooring each one to whole microseconds lost most of their time.
    class Benchmark {
        string name;
        int64_t time_sum = 0;  // nanoseconds
        size_t samplings = 0;
        size_t threads_sampled = 0;
        time_point<std::chrono::steady_clock> start_time;

    public:
        Benchmark(string name, size_t threads_sampled=0) : name(name), threads_sampled(threads_sampled) {
            //Start();
        }

        std::string GetName() const {
            return name;
        }

        void AddObservation() {
            samplings++;
        }

        void Start(bool new_sample=true) {
            samplings += new_sample;
            start_time = std::chrono::steady_clock::now();
        }

        void Stop() {
            auto const stop_time = std::chrono::steady_clock::now();
            time_sum += duration_cast<std::chrono::nanoseconds>(stop_time - start_time).count();
        }

        // The summed time in whole units of format.
        uint64_t GetDuration(Time format) const {
            uint64_t const us = static_cast<uint64_t>(time_sum) / 1'000;
            switch (format) {
                case Time::microseconds:
                    return us;
                case Time::milliseconds:
                    return us / 1'000;
                case Time::seconds:
                    return us / 1'000'000;
                case Time::minutes:
                    return us / 60'000'000;
                case Time::hours:
                    return us / 3'600'000'000;
                case Time::days:
                    return us / 86'400'000'000;
            }
            return 0L;
        }

        // Threads whose time is in the sum: those joined, and this one if it timed anything itself.
        size_t Threads() const {
            return threads_sampled + (samplings > 0);
        }

        // The summed time in seconds, and its mean per thread (what PrintResults shows).
        double Seconds() const {
            return static_cast<double>(time_sum) / 1e9;
        }

        double MeanSeconds() const {
            return Threads() > 1 ? Seconds() / static_cast<double>(Threads()) : Seconds();
        }

        void Join(Benchmark const& other, bool add_thread=true) {
            if (add_thread && other.samplings > 0) threads_sampled++;
            time_sum += other.time_sum;
        }

        // Prints the time (the mean per thread if several threads were joined). A benchmark started
        // but never stopped is stopped now. Printing does not change the sum, so it can be printed
        // and written (Seconds) in any order.
        void PrintResults() {
            if (time_sum == 0 && samplings > 0) {
                Stop();
            }
            size_t const threads = Threads();
            uint64_t time_sum_local = static_cast<uint64_t>(time_sum) / 1'000 / (threads > 1 ? threads : 1);  // µs

            uint64_t hours = time_sum_local / 3600000000;
            time_sum_local -= hours * 3600000000;
            uint64_t mins = time_sum_local / 60000000;
            time_sum_local -= mins * 60000000;
            uint64_t secs = time_sum_local / 1000000;
            time_sum_local -= secs * 1000000;
            uint64_t msecs = time_sum_local / 1000;

            std::cout << name << " took ";
            if (hours) std::cout << hours << "h ";
            if (mins) std::cout << mins << "m ";
            if (secs) std::cout << secs << "s ";
            if (msecs) std::cout << msecs << "ms";
            if (!hours && !mins && !secs && !msecs) std::cout << "less than 1ms";
            if (samplings > 1) std::cout << " (" << std::to_string(samplings) << ")";
            if (threads > 1) std::cout << " mean over " << threads << " threads";
            std::cout << std::endl;
        }
    };
}
