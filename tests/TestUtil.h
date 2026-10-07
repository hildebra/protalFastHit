#pragma once

// Helpers the unit tests share: a scratch directory, whole files read and written, random and mutated sequences.
// TestReference.h adds a reference of genes written as protal's reference.fna and reference.map and loaded.

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <random>
#include <string>
#include <system_error>
#include <unistd.h>

namespace protal::test {

    // A directory of its own in the temporary directory, removed with its content when the object goes. Every object
    // has its own, also two at once in one process. The name holds a space: protal never splits a path.
    class ScratchDir {
    public:
        std::filesystem::path path;

        explicit ScratchDir(std::string const& name = "test") {
            static std::atomic<unsigned> count{ 0 };
            path = std::filesystem::temp_directory_path() /
                   ("protal " + name + " " + std::to_string(::getpid()) + "_" + std::to_string(count++));
            std::error_code ec;
            std::filesystem::remove_all(path, ec);
            std::filesystem::create_directories(path);
        }
        ~ScratchDir() {
            std::error_code ec;
            std::filesystem::remove_all(path, ec);
        }
        ScratchDir(ScratchDir const&) = delete;
        ScratchDir& operator=(ScratchDir const&) = delete;

        // The path of a file in it.
        std::string operator/(std::string const& name) const { return (path / name).string(); }

        // Writes a file in it; returns its path.
        std::string Write(std::string const& name, std::string const& content) const {
            auto const file = path / name;
            std::ofstream(file, std::ios::binary) << content;
            return file.string();
        }
    };

    inline std::string Slurp(std::filesystem::path const& path) {
        std::ifstream is(path, std::ios::binary);
        return { std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>() };
    }

    // Writes a file; returns its path.
    inline std::string Spit(std::filesystem::path const& path, std::string const& content) {
        std::ofstream(path, std::ios::binary) << content;
        return path.string();
    }

    // `length` random bases, one draw of rng each.
    inline std::string RandomSequence(size_t length, std::mt19937& rng) {
        static constexpr char kBases[] = "ACGT";
        std::string seq(length, 'A');
        for (auto& c : seq) c = kBases[rng() % 4];
        return seq;
    }

    // seq with a share `rate` of its bases changed to another base.
    inline std::string Mutated(std::string seq, double rate, std::mt19937& rng) {
        static constexpr char kBases[] = "ACGT";
        std::bernoulli_distribution change(rate);
        for (auto& c : seq) {
            if (!change(rng)) continue;
            char other = c;
            while (other == c) other = kBases[rng() % 4];
            c = other;
        }
        return seq;
    }

    // `size` bases that compress like a genome's: random, with repeats of up to 300 bases of what came before.
    inline std::string TestData(size_t size, unsigned seed = 1) {
        std::mt19937 rng(seed);
        std::string s;
        s.reserve(size);
        while (s.size() < size) {
            if (s.size() > 1000 && rng() % 4 == 0) {
                size_t const from = rng() % (s.size() - 500);
                s.append(s, from, std::min<size_t>(300, size - s.size()));
            } else {
                s.push_back("ACGT"[rng() % 4]);
            }
        }
        return s;
    }
}
