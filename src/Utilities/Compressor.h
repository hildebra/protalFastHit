#pragma once

#include <filesystem>
#include <stdexcept>
#include <string>
#include <system_error>

#include "Bgzf.h"

class Compressor {
public:
    // Compresses a file in place: `inputFile` becomes `inputFile`.gz, as BGZF (libdeflate, with
    // `threads` threads; Bgzf.h), and the original is removed. The output holds no file name or
    // time stamp and does not depend on the thread count, so the same input always gives
    // byte-identical output. On failure nothing is removed and the partial .gz is deleted.
    static void compressInPlace(const std::filesystem::path& inputFile, int threads)
    {
        namespace fs = std::filesystem;
        if (threads <= 0)
            throw std::invalid_argument("Thread count must be greater than zero.");
        if (inputFile.empty())
            throw std::invalid_argument("Input file path is empty.");
        if (!fs::exists(inputFile))
            throw std::invalid_argument("Input file does not exist: " + inputFile.string());

        fs::path gzFile = inputFile;
        gzFile += ".gz";
        std::string const error = protal::bgzf::CompressFile(inputFile.string(), gzFile.string(), threads);
        std::error_code ec;
        if (!error.empty()) {
            fs::remove(gzFile, ec);
            throw std::runtime_error("Compressing " + inputFile.string() + " failed: " + error);
        }
        fs::remove(inputFile, ec);
        if (ec) throw std::runtime_error("Cannot remove " + inputFile.string() + " after compressing it: " + ec.message());
    }
};
