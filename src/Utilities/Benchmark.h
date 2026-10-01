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
    //
    // A stage that runs several times per read (kPerRead) times only every kPerRead-th call: reading
    // the clock costs ~16 ns, and a dozen such stages cost a read pair 10-15% of its alignment time.
    // Its time is estimated as the mean of the timed intervals times the number of calls. The period
    // is odd so that a stage called once per mate, or once per anchor, is not always timed on the
    // same mate; the first call is always timed, so short runs have a value.
    class Benchmark {
        string name;
        int64_t time_sum = 0;  // nanoseconds: whole intervals, and the estimates of joined benchmarks
        size_t samplings = 0;
        size_t threads_sampled = 0;
        time_point<std::chrono::steady_clock> start_time;
        bool running = false;  // between a Start that took the time and its Stop

        uint32_t sample_every = 1;  // 1: every call is timed
        uint32_t until_sample = 0;  // calls to skip before the next timed one
        uint64_t calls = 0;         // calls to Start, when sampling
        uint64_t timed_calls = 0;   // intervals that were timed (all of them, unless sampling)
        int64_t sampled_ns = 0;

        // The time in nanoseconds: the whole intervals and, when sampling, the estimate for the calls.
        int64_t Estimate() const {
            if (sample_every == 1 || timed_calls == 0) return time_sum;
            return time_sum + static_cast<int64_t>(static_cast<double>(sampled_ns) / static_cast<double>(timed_calls) *
                                                   static_cast<double>(calls));
        }

    public:
        // The period for the stages that run per read or per anchor.
        static constexpr uint32_t kPerRead = 61;

        Benchmark(string name, size_t threads_sampled=0, uint32_t sample_every=1) :
                name(name), threads_sampled(threads_sampled), sample_every(sample_every < 1 ? 1 : sample_every) {
            //Start();
        }

        std::string GetName() const {
            return name;
        }

        // Whether the last Start took the time (it is the one the next Stop adds), and how many calls were
        // timed, when sampling: for tests, which can then check the sampling without a clock.
        bool Timing() const {
            return running;
        }

        uint64_t TimedCalls() const {
            return timed_calls;
        }

        void AddObservation() {
            samplings++;
        }

        void Start(bool new_sample=true) {
            samplings += new_sample;
            if (sample_every > 1) {
                calls++;
                if (until_sample > 0) {
                    until_sample--;
                    running = false;
                    return;
                }
                until_sample = sample_every - 1;
            }
            running = true;
            start_time = std::chrono::steady_clock::now();
        }

        void Stop() {
            if (!running) return;
            running = false;
            auto const stop_time = std::chrono::steady_clock::now();
            int64_t const ns = duration_cast<std::chrono::nanoseconds>(stop_time - start_time).count();
            timed_calls++;
            if (sample_every > 1) {
                sampled_ns += ns;
            } else {
                time_sum += ns;
            }
        }

        // The summed time in whole units of format.
        uint64_t GetDuration(Time format) const {
            uint64_t const us = static_cast<uint64_t>(Estimate()) / 1'000;
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
            return static_cast<double>(Estimate()) / 1e9;
        }

        double MeanSeconds() const {
            return Threads() > 1 ? Seconds() / static_cast<double>(Threads()) : Seconds();
        }

        void Join(Benchmark const& other, bool add_thread=true) {
            if (add_thread && other.samplings > 0) threads_sampled++;
            time_sum += other.Estimate();
        }

        // Prints the time (the mean per thread if several threads were joined). A benchmark started
        // but never stopped is stopped now. Printing does not change the sum, so it can be printed
        // and written (Seconds) in any order.
        void PrintResults() {
            if (running) {
                Stop();
            }
            size_t const threads = Threads();
            uint64_t time_sum_local = static_cast<uint64_t>(Estimate()) / 1'000 / (threads > 1 ? threads : 1);  // µs

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
